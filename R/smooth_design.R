# Penalized smooth terms (s(), te(), ti(), t2() of mgcv) in model formulas.
#
# A formula with smooth terms is split by mgcv::interpret.gam() into its
# parametric part, which goes through model.frame()/model.matrix() as
# before, and its smooth specifications. Each smooth is set up by
# mgcv::smoothCon() with its identifiability constraint absorbed and its
# penalty scaled, and its basis columns are appended to the design matrix
# of its component (mean, a distribution parameter, the zero model, ...).
# The smooth objects are kept in the blueprint, so mgcv::PredictMat()
# rebuilds the same bases (same knots) for new data.
#
# Each smooth entry records:
#   label  the term label, e.g. "s(lnaadt)"
#   cols   its columns in the component's design matrix
#   S      its penalty matrices (one per smoothing parameter)
#   sm     the mgcv smooth object (without its training basis)

.fcr_smooth_fns <- c("s", "te", "ti", "t2")

# Does a formula contain smooth terms? Only calls count: a variable named
# `s` (as in ~ x + s) is not a smooth.
.fcr_has_smooth <- function(fml) {
  if (is.null(fml)) return(FALSE)
  walk <- function(e) {
    if (!is.call(e)) return(FALSE)
    fn <- e[[1]]
    if (is.name(fn) && as.character(fn) %in% .fcr_smooth_fns) return(TRUE)
    if (.fcr_is_tv(e) &&
        identical(.fcr_tv_basis_of(.fcr_tv_match(e)), "smooth")) {
      return(TRUE)
    }
    any(vapply(as.list(e)[-1], walk, logical(1)))
  }
  walk(stats::as.formula(fml))
}

# Split a formula into its parametric formula and smooth specifications.
.fcr_smooth_split <- function(fml) {
  rlang::check_installed("mgcv", reason = "for smooth terms (s(), te()).")
  fml <- .fcr_tv_smooth_rewrite(fml)
  g <- mgcv::interpret.gam(fml)
  pf <- g$pf
  environment(pf) <- environment(stats::as.formula(fml))
  list(pf = pf, specs = g$smooth.spec)
}

# Build the bases of smooth specifications on `data` (the estimation rows).
# Returns list(B = basis matrix, smooths = list of entries with `cols`
# relative to the basis matrix).
.fcr_smooth_build <- function(specs, data) {
  blocks <- list()
  entries <- list()
  pos <- 0L
  for (spec in specs) {
    sms <- mgcv::smoothCon(spec, data = data, knots = NULL,
                           absorb.cons = TRUE, scale.penalty = TRUE)
    for (sm in sms) {
      k <- ncol(sm$X)
      B <- sm$X
      colnames(B) <- paste0(sm$label, ".", seq_len(k))
      blocks[[length(blocks) + 1L]] <- B
      S <- sm$S
      sm$X <- NULL
      entries[[length(entries) + 1L]] <- list(
        label = sm$label, cols = pos + seq_len(k), S = S, sm = sm,
        vars = sm$term, by = if (!identical(sm$by, "NA")) sm$by)
      pos <- pos + k
    }
  }
  list(B = do.call(cbind, blocks), smooths = entries)
}

# Basis matrix of stored smooth entries for new data; rows with missing
# values in the smooth's variables are NA rows.
.fcr_smooth_predict <- function(entries, data) {
  if (!length(entries)) return(NULL)
  n <- nrow(data)
  blocks <- lapply(entries, function(e) {
    v <- unique(c(e$vars, e$by))
    miss <- setdiff(v, names(data))
    if (length(miss)) {
      stop("`newdata` is missing required column(s): ", toString(miss),
           call. = FALSE)
    }
    ok <- stats::complete.cases(data[v])
    out <- matrix(NA_real_, n, length(e$cols))
    if (any(ok)) {
      out[ok, ] <- mgcv::PredictMat(e$sm, data[ok, , drop = FALSE])
    }
    colnames(out) <- paste0(e$label, ".", seq_along(e$cols))
    out
  })
  do.call(cbind, blocks)
}

# Shift the column indices of smooth entries by `offset`.
.fcr_smooth_shift <- function(entries, offset) {
  lapply(entries, function(e) {
    e$cols <- e$cols + offset
    e
  })
}

# Data columns a formula needs: for a formula with smooth terms, those of
# its parametric part and the smooths' variables and `by` variables (not
# objects in their arguments, such as the neighborhood list in
# s(region, bs = "mrf", xt = list(nb = nb))).
.fcr_formula_vars <- function(fml) {
  if (is.null(fml)) return(character(0))
  if (!.fcr_has_smooth(fml)) return(all.vars(stats::as.formula(fml)))
  sp <- .fcr_smooth_split(fml)
  sv <- unlist(lapply(sp$specs, function(s) c(s$term,
    if (!identical(s$by, "NA")) s$by)))
  unique(c(all.vars(sp$pf), sv))
}
