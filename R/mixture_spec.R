# Finite mixtures and latent classes: control settings, family checks and the
# free-parameter map.
#
# Naming: internal helpers of the finite-mixture engine use the prefix
# `.fcr_fm_`. (The prefix `.fcr_mix_` belongs to the numerical integration of
# the mixed-Poisson families, such as PLN, in mixture_quad.R: those families
# mix over a continuous random effect *within* one component.)

#' Control settings for finite-mixture and latent-class count models
#'
#' Settings of the estimation algorithm used by [countreg_mix()] and
#' [countreg_lc()]. `lc_control()` is the same function under the name used
#' with the latent-class interface; both return the same object.
#'
#' @param algorithm `"em"` (default): each start is improved by the
#'   expectation-maximization (EM) algorithm, and the best solutions are then
#'   polished by maximizing the observed-data log-likelihood directly (with
#'   the optimizer chosen by `method`). `"direct"`: every start is maximized
#'   directly, without EM iterations.
#' @param em_maxit maximum number of EM iterations per start.
#' @param em_tol relative convergence tolerance of the EM iterations. EM stops
#'   when the relative change in the observed log-likelihood is below
#'   `em_tol` *and* the largest change in a parameter is below
#'   `sqrt(em_tol)`.
#' @param mstep_maxit maximum number of quasi-Newton iterations in each
#'   numerical M-step (a generalized EM step: it need not reach the
#'   conditional maximum, but it must increase the likelihood).
#' @param n_polish number of distinct EM solutions (the best ones) that are
#'   polished by direct maximization.
#' @param boot_starts number of starts used in each refit of a bootstrap
#'   (`stderr = "boot"`), in addition to the original estimates.
#' @param small_class posterior class share below which a component is
#'   flagged as small in the output. This is a reporting threshold only: it
#'   never changes the fit.
#' @param start_spread spread of the response-based starting partitions (the
#'   standard deviation, in units of the standardized score, of the noise
#'   added to the score that orders the units into starting classes).
#' @param ... arguments passed from `lc_control()` to `mixture_control()`.
#'
#' @returns An object of class `fcr_mixture_control`.
#' @seealso [countreg_mix()], [countreg_lc()]
#' @examples
#' mixture_control(em_maxit = 200)
#' lc_control(algorithm = "direct")
#' @export
mixture_control <- function(algorithm = c("em", "direct"), em_maxit = 200L,
                            em_tol = 1e-6, mstep_maxit = 25L, n_polish = 3L,
                            boot_starts = 5L, small_class = 0.05,
                            start_spread = 0.5) {
  algorithm <- match.arg(algorithm)
  whole <- function(v, nm, lo) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < lo ||
        v != round(v)) {
      stop("`", nm, "` must be a single whole number of at least ", lo, ".",
           call. = FALSE)
    }
    as.integer(v)
  }
  pos <- function(v, nm, hi = Inf) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0 ||
        v >= hi) {
      stop("`", nm, "` must be a single positive number",
           if (is.finite(hi)) paste0(" below ", hi), ".", call. = FALSE)
    }
    as.numeric(v)
  }
  structure(list(
    algorithm = algorithm,
    em_maxit = whole(em_maxit, "em_maxit", 0),
    em_tol = pos(em_tol, "em_tol", 1),
    mstep_maxit = whole(mstep_maxit, "mstep_maxit", 1),
    n_polish = whole(n_polish, "n_polish", 1),
    boot_starts = whole(boot_starts, "boot_starts", 1),
    small_class = pos(small_class, "small_class", 1),
    start_spread = if (identical(start_spread, 0)) 0 else
      pos(start_spread, "start_spread")
  ), class = "fcr_mixture_control")
}

#' @rdname mixture_control
#' @export
lc_control <- function(...) mixture_control(...)

#' @export
print.fcr_mixture_control <- function(x, ...) {
  cat("Finite-mixture control settings\n")
  for (nm in names(x)) cat(" ", nm, ":", format(x[[nm]]), "\n")
  invisible(x)
}

# ---- Families ----------------------------------------------------------------

# Status of each registry family as a mixture component. "validated":
# analytic scores, independent likelihood checks and simulation recovery in
# the test suite. "experimental": the component density, mean, variance,
# distribution function and random generator come from the same registry
# entries as countreg() (one-component agreement is tested), but recovery of
# mixtures of this family has not been studied.
.fcr_fm_validated <- c("POISSON", "NB1", "NB2")

.fcr_fm_family_status <- function(families) {
  ifelse(families %in% .fcr_fm_validated, "validated", "experimental")
}

# Normalize the family argument and the number of components. Returns the
# family vector (length K) and the display labels.
.fcr_fm_families <- function(family, components, components_missing,
                             label_prefix = "Component") {
  if (!is.character(family) || !length(family) || anyNA(family)) {
    stop("`family` must be a character vector of family codes.",
         call. = FALSE)
  }
  labels <- names(family)
  fam <- toupper(gsub("[^[:alnum:]]", "", family))
  bad <- setdiff(fam, .fcr_families)
  if (length(bad)) {
    stop("Unknown `family` ", toString(paste0('"', bad, '"')),
         ". Must be one of: ", toString(.fcr_families), ".", call. = FALSE)
  }
  if (!is.numeric(components) || length(components) != 1L ||
      !is.finite(components) || components < 1 ||
      components != round(components)) {
    stop("The number of components (classes) must be a single whole ",
         "number of at least 1.", call. = FALSE)
  }
  K <- as.integer(components)
  if (length(fam) > 1L) {
    if (components_missing) {
      K <- length(fam)
    } else if (length(fam) != K) {
      stop("`family` has ", length(fam), " entries but ", K,
           " components (classes) were requested; give one family, or one ",
           "per component.", call. = FALSE)
    }
  } else {
    fam <- rep(fam, K)
    labels <- NULL
  }
  if (is.null(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) {
    if (!is.null(labels)) {
      stop("Names of `family` (component labels) must be non-empty and ",
           "unique.", call. = FALSE)
    }
    labels <- paste0(label_prefix, seq_len(K))
  }
  list(families = unname(fam), labels = unname(labels), K = K)
}

# ---- Free-parameter map --------------------------------------------------------
#
# All parameters are stored in one named free vector `theta`:
#   Component<k>:Mean:<column>           component-specific mean coefficient
#   Shared:Mean:<column>                 mean coefficient shared by all
#   Component<k>:<label>[:<column>]      distribution parameter (label as in
#                                        countreg(), e.g. ln(alpha))
#   Shared:<label>[:<column>]            shared distribution parameter block
#   Mixing:Component<k>:<column>         multinomial-logit coefficient of
#                                        component k (k < K; component K is
#                                        the reference with coefficients 0)
# Coefficient names always use the immutable component number, never the
# display label. The map gives, for each component, the positions of its
# mean, dis1 and dis2 coefficients in theta (shared positions repeat across
# components) and the positions of the mixing coefficients.

# Logical vector over the columns of X: which belong to shared terms.
.fcr_fm_shared_cols <- function(shared_terms, X, tt) {
  p <- ncol(X)
  out <- rep(FALSE, p)
  if (!length(shared_terms)) return(out)
  if (!is.character(shared_terms) || anyNA(shared_terms)) {
    stop("`shared$mean` must be a character vector of formula terms.",
         call. = FALSE)
  }
  labels <- attr(tt, "term.labels")
  asg <- attr(X, "assign")
  for (s in unique(shared_terms)) {
    if (s == "(Intercept)") {
      if (!"(Intercept)" %in% colnames(X)) {
        stop("`shared$mean` names \"(Intercept)\" but the model has no ",
             "intercept.", call. = FALSE)
      }
      out[colnames(X) == "(Intercept)"] <- TRUE
    } else {
      j <- match(s, labels)
      if (is.na(j)) {
        stop("`shared$mean` term \"", s, "\" is not a term of the formula. ",
             "Terms: ", toString(c("(Intercept)", labels)), ". Name whole ",
             "terms (for a factor, the factor, not one of its dummy ",
             "columns).", call. = FALSE)
      }
      out[asg == j] <- TRUE
    }
  }
  out
}

.fcr_fm_check_shared <- function(shared) {
  if (is.null(shared)) shared <- list()
  if (!is.list(shared)) {
    stop("`shared` must be a list with elements `mean`, `dis1` and `dis2`.",
         call. = FALSE)
  }
  bad <- setdiff(names(shared), c("mean", "dis1", "dis2"))
  if (length(bad) || (length(shared) && is.null(names(shared)))) {
    stop("Unknown element(s) of `shared`: ", toString(bad), ". Use `mean` ",
         "(term names), `dis1` and `dis2` (TRUE or FALSE).", call. = FALSE)
  }
  flag <- function(v, nm) {
    if (is.null(v)) return(FALSE)
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop("`shared$", nm, "` must be TRUE or FALSE.", call. = FALSE)
    }
    v
  }
  list(mean = as.character(shared$mean %||% character(0)),
       dis1 = flag(shared$dis1, "dis1"), dis2 = flag(shared$dis2, "dis2"))
}

# Build the map. D1, D2: design matrices of the distribution parameters (a
# column of ones when the parameter has no formula); d1_formula, d2_formula:
# whether a formula was given (for the coefficient names). Z: mixing design
# at the unit level (NULL for fixed mixing or K = 1).
.fcr_fm_map <- function(families, X, tt, D1, D2, d1_formula, d2_formula, Z,
                        shared) {
  K <- length(families)
  xn <- colnames(X)
  p <- length(xn)
  params <- lapply(families, get_params)
  has1 <- vapply(params, function(pr) !is.null(pr[[1]]), logical(1))
  has2 <- vapply(params, function(pr) !is.null(pr[[2]]), logical(1))
  sh_cols <- if (K > 1L) .fcr_fm_shared_cols(shared$mean, X, tt) else
    rep(FALSE, p)
  for (s in c("dis1", "dis2")) {
    if (!isTRUE(shared[[s]]) || K == 1L) next
    has <- if (s == "dis1") has1 else has2
    if (!all(has)) {
      stop("`shared$", s, "` = TRUE, but not every component has that ",
           "distribution parameter.", call. = FALSE)
    }
    if (length(unique(families)) > 1L) {
      stop("`shared$", s, "` = TRUE needs the same family in every ",
           "component: distribution parameters of different families do ",
           "not have a common meaning, even when their labels agree.",
           call. = FALSE)
    }
  }
  if (K > 1L && length(unique(families)) == 1L && all(sh_cols) &&
      (!has1[1] || isTRUE(shared$dis1)) && (!has2[1] || isTRUE(shared$dis2))) {
    stop("Every parameter is shared, so all components are identical and ",
         "the mixing probabilities are not identified. Leave some ",
         "coefficients component-specific.", call. = FALSE)
  }

  nm <- character(0)
  info <- list()
  add <- function(name, block, comp, col, label, fam) {
    nm <<- c(nm, name)
    info[[length(info) + 1L]] <<- data.frame(
      name = name, block = block, component = comp, column = col,
      label = label, family = fam, stringsAsFactors = FALSE)
    length(nm)
  }
  dis_name <- function(prefix, label, D, j, has_formula) {
    if (has_formula) paste0(prefix, label, ":", colnames(D)[j]) else
      paste0(prefix, label)
  }
  mean_idx <- lapply(seq_len(K), function(k) integer(p))
  d1_idx <- lapply(seq_len(K), function(k) integer(0))
  d2_idx <- d1_idx
  comp_specific <- d1_idx          # component-specific positions, in order
  for (k in seq_len(K)) {
    pre <- paste0("Component", k, ":")
    for (j in which(!sh_cols)) {
      mean_idx[[k]][j] <- add(paste0(pre, "Mean:", xn[j]), "mean", k, xn[j],
                              "Mean", families[k])
    }
    if (has1[k] && !(isTRUE(shared$dis1) && K > 1L)) {
      d1_idx[[k]] <- vapply(seq_len(ncol(D1)), function(j)
        add(dis_name(pre, params[[k]][[1]], D1, j, d1_formula), "dis1", k,
            colnames(D1)[j], params[[k]][[1]], families[k]), integer(1))
    }
    if (has2[k] && !(isTRUE(shared$dis2) && K > 1L)) {
      d2_idx[[k]] <- vapply(seq_len(ncol(D2)), function(j)
        add(dis_name(pre, params[[k]][[2]], D2, j, d2_formula), "dis2", k,
            colnames(D2)[j], params[[k]][[2]], families[k]), integer(1))
    }
    comp_specific[[k]] <- c(mean_idx[[k]][!sh_cols], d1_idx[[k]],
                            d2_idx[[k]])
  }
  for (j in which(sh_cols)) {
    i <- add(paste0("Shared:Mean:", xn[j]), "mean", NA_integer_, xn[j],
             "Mean", families[1])
    for (k in seq_len(K)) mean_idx[[k]][j] <- i
  }
  if (K > 1L && isTRUE(shared$dis1)) {
    i <- vapply(seq_len(ncol(D1)), function(j)
      add(dis_name("Shared:", params[[1]][[1]], D1, j, d1_formula), "dis1",
          NA_integer_, colnames(D1)[j], params[[1]][[1]], families[1]),
      integer(1))
    for (k in seq_len(K)) d1_idx[[k]] <- i
  }
  if (K > 1L && isTRUE(shared$dis2)) {
    i <- vapply(seq_len(ncol(D2)), function(j)
      add(dis_name("Shared:", params[[1]][[2]], D2, j, d2_formula), "dis2",
          NA_integer_, colnames(D2)[j], params[[1]][[2]], families[1]),
      integer(1))
    for (k in seq_len(K)) d2_idx[[k]] <- i
  }
  mix_idx <- NULL
  if (!is.null(Z) && K > 1L) {
    mix_idx <- matrix(0L, K - 1L, ncol(Z))
    for (k in seq_len(K - 1L)) {
      for (j in seq_len(ncol(Z))) {
        mix_idx[k, j] <- add(paste0("Mixing:Component", k, ":",
                                    colnames(Z)[j]), "mixing", k,
                             colnames(Z)[j], "Mixing", NA_character_)
      }
    }
  }
  info <- do.call(rbind, info)
  rownames(info) <- NULL
  comp_par <- which(info$block != "mixing")
  list(K = K, families = families, names = nm, n_par = length(nm),
       info = info, mean_idx = mean_idx, d1_idx = d1_idx, d2_idx = d2_idx,
       mix_idx = mix_idx, comp_specific = comp_specific,
       comp_par = comp_par, shared_mean = sh_cols,
       shared = shared, has1 = has1, has2 = has2,
       any_shared = any(sh_cols) || (K > 1L &&
                                       (isTRUE(shared$dis1) ||
                                          isTRUE(shared$dis2))))
}

# Split theta into component and mixing parameters.
.fcr_fm_unpack <- function(theta, map) {
  K <- map$K
  gamma <- NULL
  if (!is.null(map$mix_idx)) {
    gamma <- cbind(t(matrix(theta[map$mix_idx], nrow = K - 1L)), 0)
  }
  list(beta = lapply(map$mean_idx, function(i) theta[i]),
       d1 = lapply(map$d1_idx, function(i) if (length(i)) theta[i]),
       d2 = lapply(map$d2_idx, function(i) if (length(i)) theta[i]),
       gamma = gamma)
}

# Reorder components: new component k is old component perm[k]. Only valid
# between components with the same family (exchangeable blocks of the same
# structure). Mixing coefficients are re-expressed against the new
# reference component (the last one); this is a linear map, so the same
# function transforms estimates, and the covariance follows from its
# matrix (.fcr_fm_perm_matrix()).
.fcr_fm_permute <- function(theta, map, perm) {
  K <- map$K
  if (identical(as.integer(perm), seq_len(K))) return(theta)
  if (!identical(map$families[perm], map$families)) {
    stop("Components can only be relabeled among components of the same ",
         "family.", call. = FALSE)
  }
  out <- theta
  for (k in seq_len(K)) {
    out[map$comp_specific[[k]]] <- theta[map$comp_specific[[perm[k]]]]
  }
  if (!is.null(map$mix_idx)) {
    G <- .fcr_fm_unpack(theta, map)$gamma[, perm, drop = FALSE]
    G <- G - G[, K]
    out[map$mix_idx] <- t(G[, seq_len(K - 1L), drop = FALSE])
  }
  out
}

# Matrix A with permute(theta) = A %*% theta (for transforming covariances)
.fcr_fm_perm_matrix <- function(map, perm) {
  n <- map$n_par
  A <- vapply(seq_len(n), function(j) {
    e <- numeric(n)
    e[j] <- 1
    .fcr_fm_permute(e, map, perm)
  }, numeric(n))
  dimnames(A) <- list(map$names, map$names)
  A
}

# Permutations that keep every component's family (all of them when the
# components are exchangeable). Fixed mixing probabilities must also agree.
.fcr_fm_allowed_perms <- function(families, fixed = NULL) {
  K <- length(families)
  perms <- .fcr_fm_all_perms(K)
  ok <- apply(perms, 1, function(p) {
    identical(families[p], families) &&
      (is.null(fixed) || isTRUE(all.equal(fixed[p], fixed)))
  })
  perms[ok, , drop = FALSE]
}

.fcr_fm_all_perms <- function(K) {
  if (K == 1L) return(matrix(1L, 1, 1))
  if (K > 7L) {
    stop("Relabeling is limited to at most 7 components.", call. = FALSE)
  }
  sub <- .fcr_fm_all_perms(K - 1L)
  do.call(rbind, lapply(seq_len(K), function(i) {
    cbind(i, matrix(setdiff(seq_len(K), i)[sub], ncol = K - 1L))
  }))
}
