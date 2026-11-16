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
.fcr_family_nesting <- data.frame(
  restricted = c(rep("POISSON", 9), "NB1", "NB2", "PL", "PL"),
  full = c("NB1", "NB2", "PLN", "PIG1", "PIG2", "PIG", "GP1", "GP2", "COM",
           "NBP", "NBP", "PLG", "PLL"),
  type = c(rep("boundary", 6), rep("interior", 5), "boundary", "boundary"),
  param = c("ln(alpha)", "ln(alpha)", "ln(sigma)", "ln(eta)", "ln(eta)",
            "ln(eta)", "ln(phi+1)", "ln(alpha+1)", "ln(nu)", "ln(p)", "ln(p)",
            "ln(alpha)", "ln(sigma)"),
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
  list(class = class(object)[1],
       family = family,
       coefs = .fcr_norm_coef_names(names(stats::coef(object)), family),
       zero = if (!is.null(m$zero_formula)) m$zero_type else "none",
       underreport = !is.null(m$underreport_formula),
       rpar = sort(as.character(rp)),
       converged = m$converged)
}

# Spec of a countreg() base model of `family` with the given mean
# coefficients and scalar distribution parameters.
.fcr_base_spec <- function(family, mean_coefs) {
  family <- toupper(family)
  list(class = "countreg", family = family,
       coefs = c(mean_coefs, unlist(get_params(family))),
       zero = "none", underreport = FALSE, rpar = character(0),
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
  extra <- setdiff(f$coefs, r$coefs)
  if (!all(r$coefs %in% f$coefs) || !length(extra)) {
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
  list(nested = TRUE, type = type, df = length(extra), reason = NA_character_)
}

# LR statistic and p-value for a nesting result; NA when not nested.
.fcr_lr_pvalue <- function(ll_full, ll_restricted, nest) {
  if (!isTRUE(nest$nested)) return(list(LR = NA_real_, p = NA_real_))
  lr <- 2 * (ll_full - ll_restricted)
  if (!is.finite(lr)) return(list(LR = NA_real_, p = NA_real_))
  lr <- max(lr, 0)
  df <- nest$df
  p <- if (nest$type == "boundary") {
    # One parameter on the boundary: 0.5 chi2(df - 1) + 0.5 chi2(df)
    0.5 * stats::pchisq(lr, df - 1, lower.tail = FALSE) +
      0.5 * stats::pchisq(lr, df, lower.tail = FALSE)
  } else {
    stats::pchisq(lr, df, lower.tail = FALSE)
  }
  list(LR = lr, p = p)
}
