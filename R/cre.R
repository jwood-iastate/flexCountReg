# Correlated random effects (Mundlak 1978; Chamberlain 1982): the group
# means of time-varying regressors are added to the mean model of a
# random-effects model. The random effect is then allowed to depend on the
# regressors through their means, and a Wald test that the coefficients of
# the means are zero is a robust Hausman-type test of random against fixed
# effects (Wooldridge 2019).
#
# The means are computed on the estimation rows only and stored with the
# fit, so predictions for a group in the estimation data use its training
# means: held-out rows never change them. A group not in the estimation data
# gets the means of its own rows in `newdata` (its covariates only).

# Add the group-mean columns to `data` and to `formula`. `formulas` are the
# model's other formulas (their variables define the complete rows).
.fcr_cre_prepare <- function(cre, formula, data, group_var,
                             formulas = list(), extra_vars = character(0)) {
  if (!inherits(cre, "formula") || length(cre) != 2L) {
    stop("`cre` must be a one-sided formula of time-varying regressors, ",
         "e.g. cre = ~ x1 + x2.", call. = FALSE)
  }
  if (is.null(group_var)) {
    stop("`cre` needs `group_var`: the group means are taken within ",
         "groups.", call. = FALSE)
  }
  vars <- attr(stats::terms(cre), "term.labels")
  miss <- setdiff(vars, names(data))
  if (!length(vars) || length(miss)) {
    stop("`cre` must list columns of `data`",
         if (length(miss)) paste0("; not found: ", toString(miss)), ".",
         call. = FALSE)
  }
  num <- vapply(data[vars], is.numeric, logical(1))
  if (!all(num)) {
    stop("`cre` variables must be numeric (make indicator columns for ",
         "factors): ", toString(vars[!num]), ".", call. = FALSE)
  }
  key <- .fcr_cre_key(data, group_var)
  used <- unique(c(all.vars(formula),
                   unlist(lapply(Filter(Negate(is.null), formulas),
                                 all.vars)),
                   vars, group_var, extra_vars))
  used <- intersect(used, names(data))
  ok <- stats::complete.cases(data[used])
  cols <- character(0)
  means <- data.frame(.key = unique(key[ok]), stringsAsFactors = FALSE)
  for (v in vars) {
    col <- .fcr_unused_name(c(names(data), cols), make.names(paste0("mean_",
                                                                     v)))
    m <- tapply(data[[v]][ok], key[ok], mean)
    data[[col]] <- as.vector(m[key])
    means[[col]] <- as.vector(m[means$.key])
    cols <- c(cols, col)
  }
  formula <- stats::update(
    formula, stats::as.formula(paste(". ~ . +", paste(cols, collapse = " + ")),
                               env = environment(formula)))
  list(data = data, formula = formula,
       info = list(vars = vars, cols = stats::setNames(cols, vars),
                   group_var = group_var, means = means))
}

# Group key of each row (one or more grouping columns).
.fcr_cre_key <- function(data, group_var) {
  if (length(group_var) == 1L) return(as.character(data[[group_var]]))
  do.call(paste, c(lapply(data[group_var], as.character), sep = "\r"))
}

# The group-mean columns for `newdata` of a fit with `cre`: training means
# for known groups, the newdata rows' own means for new groups.
.fcr_cre_newdata <- function(cre, data) {
  miss <- setdiff(c(cre$group_var, cre$vars), names(data))
  if (length(miss)) {
    stop("This model has correlated random effects (`cre`); `newdata` ",
         "needs the column(s) ", toString(miss), ".", call. = FALSE)
  }
  key <- .fcr_cre_key(data, cre$group_var)
  idx <- match(key, cre$means$.key)
  for (v in cre$vars) {
    col <- cre$cols[[v]]
    val <- cre$means[[col]][idx]
    new <- is.na(idx)
    if (any(new)) {
      own <- tapply(data[[v]][new], key[new], mean, na.rm = TRUE)
      val[new] <- as.vector(own[key[new]])
    }
    data[[col]] <- val
  }
  data
}

#' Test of correlated random effects (random against fixed effects)
#'
#' Wald test that the coefficients of the group means added by
#' \code{cre} are all zero. Under the null hypothesis the random effects
#' are uncorrelated with the time-varying regressors, so the random-effects
#' estimates are consistent; a rejection favors fixed effects (or keeping
#' the means in the model). With robust or cluster-robust standard errors
#' the test is a robust version of the Hausman test (Wooldridge 2019).
#'
#' @param object a model fitted with \code{cre} (\code{\link{countreg_rp}()},
#'   \code{\link{renb}()} or \code{\link{repois}()}).
#' @param vcov optional covariance matrix of the estimates (default: the
#'   fit's, \code{vcov(object)}).
#'
#' @returns An object of class \code{"htest"}.
#' @references
#' Mundlak, Y. (1978). On the pooling of time series and cross section
#' data. \emph{Econometrica}, 46(1), 69-85.
#'
#' Wooldridge, J.M. (2019). Correlated random effects models with
#' unbalanced panels. \emph{Journal of Econometrics}, 211(1), 137-150.
#' @examples
#' \donttest{
#' data("washington_roads")
#' m <- repois(Total_crashes ~ lnaadt + offset(lnlength), group_var = "ID",
#'             data = washington_roads, cre = ~ lnaadt, stderr = "robust")
#' cre_test(m)
#' }
#' @export
cre_test <- function(object, vcov = NULL) {
  cre <- object$model$cre
  if (is.null(cre)) {
    stop("The model was not fitted with `cre`.", call. = FALSE)
  }
  cols <- unname(cre$cols)
  b <- stats::coef(object)[cols]
  V <- if (is.null(vcov)) stats::vcov(object) else vcov
  V <- V[cols, cols, drop = FALSE]
  W <- as.numeric(crossprod(b, solve(V, b)))
  df <- length(cols)
  structure(list(
    statistic = c("Wald chi-squared" = W), parameter = c(df = df),
    p.value = stats::pchisq(W, df, lower.tail = FALSE),
    method = "Wald test of correlated random effects (group means)",
    data.name = paste(deparse(object$call), collapse = " "),
    estimate = b), class = "htest")
}
