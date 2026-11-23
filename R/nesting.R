# Likelihood-ratio tests only for valid nesting.
#
# A model "spec" describes what a likelihood-ratio comparison needs: the
# estimator, count family, normalized coefficient names, and the model
# components (zero model, underreporting, random parameters) whose absence
# would put a parameter at a nonregular boundary. Two specs are compared with
# .fcr_nesting(); .fcr_lr_pvalue() then gives the p-value for the kind of
# restriction found.

# Count families nested in other families. "interior": the restriction is in
# the interior of the parameter space (ordinary chi-square). "boundary": one
# parameter is on the boundary of its space (a 50:50 mixture of chi-square
# distributions); valid only when that parameter is a scalar in the full
# model (with a dispersion formula the restriction is not regular).
# PL, PLG and PLL do not nest Poisson (their variance is at least
# mu + mu^2 / 2), NBP and SI have a parameter that is not identified at the
# Poisson limit, and PW, PGE and GW reach Poisson only as a parameter goes to
# infinity, so none of these comparisons is listed.
# DPO, HYPERPOIS and GAMMACOUNT (all mean-parameterized) are Poisson at
# ln(phi) = 0, ln(gamma) = 0 and ln(shape) = 0, interior points at which the
# models are identified. DELAPORTE, HERMITE, PA and NBGE reach Poisson (or
# PGE) only as a parameter goes to the edge of its space and are not
# listed.
.fcr_family_nesting <- data.frame(
  restricted = c(rep("POISSON", 9), "NB1", "NB2", "PL", "PL",
                 rep("POISSON", 3)),
  full = c("NB1", "NB2", "PLN", "PIG1", "PIG2", "PIG", "GP1", "GP2", "COM",
           "NBP", "NBP", "PLG", "PLL", "DPO", "HYPERPOIS", "GAMMACOUNT"),
  type = c(rep("boundary", 6), rep("interior", 5), "boundary", "boundary",
           rep("interior", 3)),
  param = c("ln(alpha)", "ln(alpha)", "ln(sigma)", "ln(eta)", "ln(eta)",
            "ln(eta)", "ln(phi+1)", "ln(alpha+1)", "ln(nu)", "ln(p)", "ln(p)",
            "ln(alpha)", "ln(sigma)", "ln(phi)", "ln(gamma)", "ln(shape)"),
  stringsAsFactors = FALSE)

# Coefficient names with the "(Intercept)" of a distribution-parameter formula
# folded into the scalar name, so "ln(alpha):(Intercept)" matches "ln(alpha)".
.fcr_norm_coef_names <- function(nm, family) {
  labels <- if (!is.null(family) && family %in% .fcr_families) {
    unlist(get_params(family))
  }
  for (lab in labels) {
    nm[nm == paste0(lab, ":(Intercept)")] <- lab
  }
  nm
}

.fcr_model_spec <- function(object) {
  m <- object$model
  family <- m$family
  if (is.character(family)) family <- toupper(family)
  rp <- m$rpar %||% if (inherits(object, "countreg_rp"))
    sub(":Mean$", "", grep(":Mean$", names(stats::coef(object)), value = TRUE))
  rp <- as.character(rp)
  rd <- m$rpardists
  list(class = class(object)[1],
       family = family,
       coefs = .fcr_norm_coef_names(names(stats::coef(object)), family),
       # Mean coefficients: their nesting is checked on the design matrices
       mean = (m$blueprint %||% list())$x_names,
       zero = if (!is.null(m$zero_formula)) m$zero_type else "none",
       underreport = !is.null(m$underreport_formula),
       rpar = sort(rp),
       rpardists = if (length(rp) && !is.null(rd)) unname(rd[sort(rp)]),
       correlated = isTRUE(m$correlated),
       converged = m$converged)
}

# Spec of a countreg() base model of `family` with the given mean
# coefficients and scalar distribution parameters.
.fcr_base_spec <- function(family, mean_coefs) {
  family <- toupper(family)
  list(class = "countreg", family = family,
       coefs = c(mean_coefs, unlist(get_params(family))),
       mean = mean_coefs,
       zero = "none", underreport = FALSE, rpar = character(0),
       rpardists = NULL, correlated = FALSE,
       converged = TRUE)
}

# Is `r` (restricted) nested in `f` (full)? Returns list(nested, type, df,
# reason).
.fcr_nesting <- function(r, f) {
  no <- function(reason) list(nested = FALSE, type = NA_character_,
                              df = NA_integer_, reason = reason)
  if (!identical(r$class, f$class)) {
    return(no("the models use different estimators"))
  }
  if (!identical(r$zero, f$zero) || !identical(r$underreport, f$underreport)) {
    return(no(paste("the models differ in their zero or underreporting",
                    "components, whose restriction is not regular")))
  }
  if (!identical(r$rpar, f$rpar)) {
    return(no(paste("the models have different random parameters (a",
                    "variance of zero is on the boundary)")))
  }
  if (!identical(r$rpardists, f$rpardists) ||
      !identical(r$correlated, f$correlated)) {
    return(no(paste("the random parameters have different distributions",
                    "or correlation structures")))
  }
  # Mean coefficients are compared through their design matrices (see
  # .fcr_design_nested()), so equivalent bases are recognized; the other
  # parameters by name
  r_other <- setdiff(r$coefs, r$mean)
  f_other <- setdiff(f$coefs, f$mean)
  df <- length(f$coefs) - length(r$coefs)
  if (!all(r_other %in% f_other) || df < 1L ||
      (is.null(r$mean) != is.null(f$mean))) {
    return(no(paste("the parameters of the smaller model are not a subset",
                    "of those of the larger model")))
  }
  type <- "interior"
  if (!identical(r$family, f$family)) {
    k <- which(.fcr_family_nesting$restricted == r$family &
                 .fcr_family_nesting$full == f$family)
    if (!length(k)) {
      return(no(paste0("the ", r$family, " family is not nested in the ",
                       f$family, " family with a regular restriction")))
    }
    type <- .fcr_family_nesting$type[k]
    if (type == "boundary") {
      par <- .fcr_family_nesting$param[k]
      if (any(startsWith(f$coefs, paste0(par, ":")))) {
        return(no(paste0("the restriction sets ", par, " to its boundary, ",
                         "which is not regular when ", par, " has a formula")))
      }
    }
  }
  list(nested = TRUE, type = type, df = df, reason = NA_character_)
}

# LR statistic and p-value for a nesting result; NA (with a reason) when not
# nested, when a log-likelihood is not finite, or when the larger model has a
# materially lower log-likelihood than the restricted one (LR below
# -max(0.01, 2e-6 |logLik|)): the larger fit then did not reach its maximum,
# and the test is unavailable. Smaller negative values, which arise when an
# optimizer stops on a flat boundary (e.g. NB2 with ln(alpha) -> -Inf), are
# rounded to 0; they change no inference (p > 0.9).
.fcr_lr_pvalue <- function(ll_full, ll_restricted, nest) {
  na <- function(reason) list(LR = NA_real_, p = NA_real_, reason = reason)
  if (!isTRUE(nest$nested)) return(na(nest$reason %||% "not nested"))
  lr <- 2 * (ll_full - ll_restricted)
  if (!is.finite(lr)) return(na("a log-likelihood is not finite"))
  if (lr < -max(0.01, 2e-6 * abs(ll_full))) {
    return(na(paste("the larger model has a lower log-likelihood than the",
                    "restricted model, so it did not reach its maximum")))
  }
  lr <- max(lr, 0)
  df <- nest$df
  p <- if (nest$type == "boundary") {
    # One parameter on the boundary: 0.5 chi2(df - 1) + 0.5 chi2(df)
    0.5 * stats::pchisq(lr, df - 1, lower.tail = FALSE) +
      0.5 * stats::pchisq(lr, df, lower.tail = FALSE)
  } else {
    stats::pchisq(lr, df, lower.tail = FALSE)
  }
  list(LR = lr, p = p, reason = NA_character_)
}

# ---- Sample and design evidence for LR tests ---------------------------------

# What an LR comparison needs to know about the analysis sample of a fit:
# the rows used (of nrow(data)), the response, weights, offset and panel of
# each row, and the raw values of every variable the model uses on those
# rows. Stored in the fitted model at estimation.
.fcr_lr_sample <- function(prep, data, group_var = NULL) {
  rows <- prep$rows_used
  bp <- prep$blueprint
  forms <- c(list(prep$formula), unname(prep$aux_formulas))
  v <- unique(c(unlist(lapply(forms, all.vars)), bp$offset_cols,
                if (is.character(group_var)) group_var))
  v <- intersect(v, names(data))
  vars <- as.data.frame(data)[rows, v, drop = FALSE]
  rownames(vars) <- NULL
  list(version = 1L, n_data = nrow(data), rows_used = rows,
       y = as.numeric(prep$y), weights = prep$weights, offset = prep$offset,
       group = if (!is.null(prep$group)) as.integer(prep$group),
       vars = vars)
}

# "" when two fits used the same observations with the same data, otherwise
# the reason they did not. Row numbers alone are not enough (they repeat
# across unrelated data sets), so the response, weights, offsets, panels and
# the values of every variable both models use must agree too. Variables
# only one model uses are not compared.
.fcr_same_sample <- function(sa, sb) {
  # A sample stored without its data (keep_data = FALSE) holds fingerprints
  # of its values (.fcr_lr_sample_strip()); the other sample's values are
  # then compared through their fingerprints
  fp <- isTRUE(sa$fingerprinted) || isTRUE(sb$fingerprinted)
  val <- function(s, nm) {
    if (fp && !isTRUE(s$fingerprinted)) .fcr_fingerprint(s[[nm]]) else s[[nm]]
  }
  var <- function(s, v) {
    if (fp && !isTRUE(s$fingerprinted)) .fcr_fingerprint(s$vars[[v]]) else
      s$vars[[v]]
  }
  eq <- function(a, b) isTRUE(all.equal(a, b, check.attributes = FALSE))
  if (length(sa$rows_used) != length(sb$rows_used) ||
      !identical(sa$rows_used, sb$rows_used))
    return("the models were fit to different rows")
  if (!eq(val(sa, "y"), val(sb, "y"))) return("the models have different outcomes")
  if (!eq(val(sa, "weights"), val(sb, "weights")))
    return("the models have different weights")
  if (!eq(val(sa, "offset"), val(sb, "offset")))
    return("the models have different offsets")
  if (!is.null(sa$group) || !is.null(sb$group)) {
    ga <- sa$group %||% seq_along(sa$rows_used)
    gb <- sb$group %||% seq_along(sb$rows_used)
    pairs <- unique(cbind(ga, gb))
    if (nrow(pairs) != length(unique(ga)) || nrow(pairs) != length(unique(gb)))
      return("the models have different panels")
  }
  shared <- intersect(names(sa$vars), names(sb$vars))
  for (v in shared) {
    if (!eq(var(sa, v), var(sb, v))) {
      return(paste0("the models use different values of `", v, "`"))
    }
  }
  ""
}

# Summary of a vector that identifies its values without storing them:
# numeric vectors keep their length, number of missing values and three
# weighted sums (compared with all.equal()'s tolerance, like the values);
# other vectors keep a hash of their values.
.fcr_fingerprint <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.numeric(x) || is.logical(x)) {
    x <- as.numeric(x)
    ok <- !is.na(x)
    t <- seq_along(x) / max(1, length(x))
    return(c(length(x), sum(!ok), sum(x[ok]), sum(x[ok]^2),
             sum(t[ok] * x[ok])))
  }
  rlang::hash(as.character(x))
}

# Sample evidence without the data (for keep_data = FALSE): row numbers and
# panels are kept, values are replaced by fingerprints.
.fcr_lr_sample_strip <- function(s) {
  if (is.null(s) || isTRUE(s$fingerprinted)) return(s)
  for (nm in c("y", "weights", "offset")) s[[nm]] <- .fcr_fingerprint(s[[nm]])
  s$vars <- lapply(s$vars, .fcr_fingerprint)
  s$fingerprinted <- TRUE
  s
}

# Is every column of A in the column space of B (relative tolerance)?
.fcr_in_span <- function(A, B, tol = 1e-8) {
  if (is.null(A) || !ncol(A)) return(TRUE)
  if (is.null(B) || !ncol(B)) return(FALSE)
  q <- qr(B, tol = 1e-10)
  res <- qr.resid(q, A)
  all(abs(res) <= tol * pmax(1, abs(A)) * max(1, sqrt(nrow(A))))
}

# Design matrices of a fit on its analysis rows: the mean design and each
# auxiliary design (a column of ones for a scalar distribution parameter).
.fcr_fit_designs <- function(object) {
  m <- object$model
  bp <- m$blueprint
  if (isTRUE(m$lr_sample$fingerprinted)) {
    if (is.null(object$data)) {
      stop("A model was fit with keep_data = FALSE.", call. = FALSE)
    }
    vars <- as.data.frame(object$data)[m$lr_sample$rows_used, , drop = FALSE]
  } else {
    vars <- m$lr_sample$vars
  }
  fd <- forge_new_data(bp, vars)
  out <- list(mean = fd$X)
  n <- nrow(fd$X)
  ones <- matrix(1, n, 1)
  fam <- m$family
  par_labels <- if (is.character(fam) && toupper(fam) %in% .fcr_families)
    get_params(toupper(fam)) else list(NULL, NULL)
  for (k in 1:2) {
    nm <- paste0("dis", k)
    if (!is.null(fd$aux[[nm]])) out[[nm]] <- fd$aux[[nm]]
    else if (!is.null(par_labels[[k]])) out[[nm]] <- ones
  }
  for (nm in setdiff(names(fd$aux), c("dis1", "dis2"))) {
    out[[nm]] <- fd$aux[[nm]]
  }
  out
}

# "" when every design block of the restricted fit is in the column space of
# the matching block of the full fit (so the restriction is linear in the
# parameters), otherwise the reason. A block only the full fit has is the
# extra parameters.
.fcr_design_nested <- function(obj_r, obj_f) {
  dr <- tryCatch(.fcr_fit_designs(obj_r), error = function(e) NULL)
  df <- tryCatch(.fcr_fit_designs(obj_f), error = function(e) NULL)
  if (is.null(dr) || is.null(df)) {
    if (is.null(obj_r$data) || is.null(obj_f$data)) {
      return(paste("the design matrices of the models could not be rebuilt",
                   "(a model was fit with keep_data = FALSE)"))
    }
    return("the design matrices of the models could not be rebuilt")
  }
  for (nm in names(dr)) {
    if (is.null(df[[nm]])) {
      return(paste0("the smaller model has a `", nm, "` component that the ",
                    "larger model does not have"))
    }
    if (!.fcr_in_span(dr[[nm]], df[[nm]])) {
      return(paste0("the `", nm, "` design of the smaller model is not in ",
                    "the column space of the larger model's"))
    }
  }
  ""
}

# Full eligibility of an LR test between two fitted models: stored sample
# evidence, convergence, finite log-likelihoods, the same sample, nested
# specifications and nested designs. Returns a nesting result as
# .fcr_nesting() does, plus `sample_ok` (FALSE when the samples differ).
.fcr_lr_eligibility <- function(obj_r, obj_f) {
  no <- function(reason, sample_ok = TRUE)
    list(nested = FALSE, type = NA_character_, df = NA_integer_,
         reason = reason, sample_ok = sample_ok)
  sr <- obj_r$model$lr_sample
  sf <- obj_f$model$lr_sample
  if (is.null(sr) || is.null(sf)) {
    return(no(paste("a model was fit with an older version of flexCountReg",
                    "and lacks the sample information an LR test needs;",
                    "refit it")))
  }
  why <- .fcr_same_sample(sr, sf)
  if (nzchar(why)) return(no(why, sample_ok = FALSE))
  if (identical(obj_r$model$weight_type, "sampling") ||
      identical(obj_f$model$weight_type, "sampling")) {
    return(no(paste("a model uses sampling weights, so its log-likelihood is",
                    "a pseudo-likelihood")))
  }
  if (!isTRUE(obj_r$model$converged) || !isTRUE(obj_f$model$converged)) {
    return(no(paste("a model did not report convergence (or its",
                    "convergence is unknown)")))
  }
  ll <- c(as.numeric(stats::logLik(obj_r)), as.numeric(stats::logLik(obj_f)))
  if (!all(is.finite(ll))) return(no("a log-likelihood is not finite"))
  nest <- .fcr_nesting(.fcr_model_spec(obj_r), .fcr_model_spec(obj_f))
  if (!isTRUE(nest$nested)) return(c(nest, sample_ok = TRUE))
  why <- .fcr_design_nested(obj_r, obj_f)
  if (nzchar(why)) return(no(why))
  c(nest, sample_ok = TRUE)
}
