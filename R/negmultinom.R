#' Negative Multinomial Family Panel Regression (NM1, NM2, NM-P)
#'
#' @name negmultinom
#' @aliases negmultinom
#'
#' @description
#' Estimate a negative multinomial (NM) regression model for clustered or panel
#' count data. Within each cluster (e.g., a site observed over several years),
#' the counts are conditionally independent Poisson variables that share a
#' common gamma-distributed multiplicative frailty with mean 1. Integrating
#' over the frailty gives the negative multinomial likelihood in closed form.
#'
#' Three members of the family are available, mirroring the NB-1, NB-2, and
#' NB-P forms of the negative binomial in \code{\link{countreg}}:
#' \itemize{
#'   \item \code{"NM2"} (default): the conventional negative multinomial. Every
#'     observation, and every within-cluster aggregate, has an NB-2 marginal
#'     with the same \eqn{\alpha}. This is the model estimated by
#'     \code{\link{repois}}.
#'   \item \code{"NM1"}: the frailty variance is calibrated so that the cluster
#'     total follows an NB-1 distribution.
#'   \item \code{"NMP"}: the frailty variance is calibrated so that the cluster
#'     total follows an NB-P distribution with power parameter \eqn{p}.
#' }
#'
#' Both \eqn{\ln(\alpha)} and \eqn{p} can be specified as functions of
#' cluster-level variables.
#'
#' @param formula an R formula for the mean, e.g. \code{y ~ x1 + x2}. An
#'   \code{offset()} term may be included in the formula.
#' @param group_var character vector naming the variable(s) that identify the
#'   clusters (e.g., a site or panel ID). If more than one variable is given,
#'   clusters are defined by their interaction.
#' @param data a data frame containing all variables used in the model.
#' @param family the member of the negative multinomial family to estimate:
#'   \code{"NM2"} (default), \code{"NM1"}, or \code{"NMP"}. Case and
#'   non-alphanumeric characters are ignored, so \code{"NM-P"} and
#'   \code{"nm2"} are accepted.
#' @param offset optional character vector naming offset variable(s) in
#'   \code{data}. Multiple offsets are summed. Do not also use \code{offset()}
#'   in \code{formula} for the same variable.
#' @param dis_param_formula_1 optional one-sided formula for
#'   \eqn{\ln(\alpha)}, e.g. \code{~ x3}. If \code{NULL}, a single constant
#'   \eqn{\ln(\alpha)} is estimated. See Details on within-cluster variation.
#' @param dis_param_formula_2 optional one-sided formula for \eqn{p} (only
#'   used when \code{family = "NMP"}). The identity link is used, so
#'   \eqn{p = W\delta}. If \code{NULL}, a single constant \eqn{p} is
#'   estimated.
#' @param method optimization method passed to \code{\link[maxLik]{maxLik}}.
#'   One of \code{"BFGS"} (default), \code{"NR"}, \code{"BHHH"}, \code{"NM"},
#'   \code{"CG"}, \code{"BFGSR"}, or \code{"SANN"}. An analytic gradient is
#'   supplied to the gradient-based methods.
#' @param max.iters maximum number of iterations for the optimizer.
#' @param start.vals optional numeric vector of starting values, ordered as
#'   the mean coefficients, then the \eqn{\ln(\alpha)} coefficient(s), then
#'   (for \code{"NMP"}) the \eqn{p} coefficient(s).
#' @param stderr type of standard errors: \code{"normal"} (default; inverse
#'   of the negative Hessian), \code{"robust"} (sandwich estimator using
#'   cluster-level scores, i.e., cluster-robust), or \code{"boot"} (cluster
#'   bootstrap; requires \code{bootstraps}).
#' @param bootstraps integer number of cluster bootstrap replications used when
#'   \code{stderr = "boot"}.
#' @param verbose logical. If \code{TRUE}, optimizer progress is printed.
#'
#' @details
#' ## Model structure
#'
#' Let \eqn{i} index clusters and \eqn{t} index the observations within a
#' cluster. The mean is modeled as
#' \deqn{\mu_{it} = \exp(X_{it}\beta + \text{offset}_{it}).}
#' Conditional on a cluster frailty \eqn{U_i},
#' \deqn{Y_{it} \mid U_i \sim \text{Poisson}(\mu_{it} U_i),}
#' independently over \eqn{t}, and
#' \eqn{U_i \sim \text{Gamma}(r_i, 1/r_i)} so that \eqn{E(U_i) = 1} and
#' \eqn{\text{Var}(U_i) = 1/r_i}.
#'
#' Define the cluster totals over the reference collection \eqn{R} (all of a
#' cluster's observations used in estimation):
#' \deqn{y_{iR} = \sum_t y_{it}, \qquad \mu_{iR} = \sum_t \mu_{it}.}
#'
#' The three family members differ only in the gamma shape parameter
#' \eqn{r_i}:
#' \deqn{\text{NM2: } r_i = \frac{1}{\alpha_i}, \qquad
#'       \text{NM1: } r_i = \frac{\mu_{iR}}{\alpha_i}, \qquad
#'       \text{NM-P: } r_i = \frac{\mu_{iR}^{2-p_i}}{\alpha_i}.}
#' NM2 and NM1 are the special cases \eqn{p = 2} and \eqn{p = 1} of NM-P.
#'
#' ## Probability mass function
#'
#' Integrating over \eqn{U_i} gives the negative multinomial PMF:
#' \deqn{P(\mathbf{y}_i \mid \boldsymbol{\mu}_i) =
#'   \frac{\Gamma(r_i + y_{iR})}{\Gamma(r_i)\prod_t y_{it}!}
#'   \left(\frac{r_i}{r_i + \mu_{iR}}\right)^{r_i}
#'   \prod_t \left(\frac{\mu_{it}}{r_i + \mu_{iR}}\right)^{y_{it}}.}
#'
#' The log-likelihood for cluster \eqn{i} is:
#' \deqn{LL_i = \ln\Gamma(r_i + y_{iR}) - \ln\Gamma(r_i) + r_i\ln(r_i)
#'   - (r_i + y_{iR})\ln(r_i + \mu_{iR})
#'   + \sum_t y_{it}\ln(\mu_{it}) - \sum_t \ln(y_{it}!).}
#'
#' ## Mean and variance
#'
#' For every family member, \eqn{E(Y_{it}) = \mu_{it}}. The variance and
#' covariance are:
#' \deqn{\text{Var}(Y_{it}) = \mu_{it} + \frac{\mu_{it}^2}{r_i}, \qquad
#'       \text{Cov}(Y_{is}, Y_{it}) = \frac{\mu_{is}\mu_{it}}{r_i}.}
#' The cluster total therefore satisfies
#' \deqn{\text{Var}(Y_{iR}) = \mu_{iR} + \alpha_i\mu_{iR}^{p_i},}
#' i.e., it is NB-2, NB-1, or NB-P for \code{"NM2"}, \code{"NM1"}, and
#' \code{"NMP"}, respectively.
#'
#' ## Interpretation of NM1 and NM-P
#'
#' For NM2, \eqn{r_i = 1/\alpha_i} does not depend on the mean, so each
#' individual observation (and any within-cluster subtotal) is NB-2 with the
#' same \eqn{\alpha}. The model is coherent across temporal aggregation.
#'
#' For \eqn{p \neq 2}, the frailty distribution depends on \eqn{\mu_{iR}}.
#' Only the cluster total over the reference collection has the NB-P (or NB-1)
#' variance function; individual observations have
#' \eqn{\text{Var}(Y_{it}) = \mu_{it} + \alpha_i\mu_{iR}^{p_i - 2}\mu_{it}^2},
#' which is not NB-P with the same \eqn{\alpha} and \eqn{p}. The reference
#' collection used here is all rows of a cluster that remain after removing
#' missing values. Changing the set of observations in a cluster (e.g.,
#' adding or dropping years) changes the implied frailty distribution, so
#' NM1 and NM-P should be interpreted as reference-window-specific joint
#' models.
#'
#' ## Parameters as functions of variables
#'
#' When \code{dis_param_formula_1} or \code{dis_param_formula_2} is supplied,
#' \deqn{\ln(\alpha_i) = \bar{Z}_i\gamma, \qquad p_i = \bar{W}_i\delta,}
#' where \eqn{\bar{Z}_i} and \eqn{\bar{W}_i} are the within-cluster means of
#' the corresponding design matrices. Because \eqn{\alpha} and \eqn{p}
#' describe the shared cluster frailty, they must be constant within a
#' cluster. Cluster-level variables (constant within a cluster) are
#' therefore recommended; if a variable varies within a cluster, a warning
#' is issued and its cluster mean is used.
#'
#' The \eqn{\ln(\alpha)} parameters use a log link so that \eqn{\alpha > 0}.
#' The \eqn{p} parameters use the identity link, so \eqn{p} is unrestricted.
#'
#' ## Empirical Bayes estimates
#'
#' Given the cluster history, the posterior of the frailty is
#' \eqn{U_i \mid \mathbf{y}_i \sim
#'   \text{Gamma}(r_i + y_{iR}, 1/(r_i + \mu_{iR}))}. The returned model
#' includes the empirical Bayes (EB) expected frequency for each observation:
#' \deqn{E(\Lambda_{it} \mid \mathbf{y}_i) =
#'   \mu_{it}\frac{r_i + y_{iR}}{r_i + \mu_{iR}}.}
#' For NM2 this is
#' \eqn{\mu_{it}(1 + \alpha y_{iR})/(1 + \alpha\mu_{iR})}. EB estimates for
#' new data can be obtained with \code{predict(model, newdata, method = "EB")}
#' when \code{newdata} contains the outcome and grouping variables.
#'
#' ## Standard errors
#'
#' The log-likelihood is evaluated cluster by cluster, so the \code{"robust"}
#' sandwich estimator and the \code{"BHHH"} optimizer both use cluster-level
#' scores. The bootstrap resamples whole clusters.
#'
#' @returns
#' An object of class \code{flexCountReg}, a list with components:
#' \itemize{
#'   \item \code{model}: the fitted model. Along with the \code{maxLik}
#'     output, this includes \code{estimate}, \code{se}, \code{vcov},
#'     \code{family}, \code{alpha} and \code{p} (cluster-level values),
#'     \code{fitted.values} (marginal means \eqn{\mu_{it}}),
#'     \code{eb.fitted} (EB expected frequencies),
#'     \code{frailty.posterior} (cluster posterior mean of \eqn{U_i}),
#'     \code{n_groups}, \code{nobs}, \code{aic}, and \code{bic}.
#'   \item \code{data}: the data frame used to fit the model.
#'   \item \code{call}: the matched call.
#'   \item \code{formula}: the formula used to fit the model.
#' }
#'
#' @examples
#' \donttest{
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT > 10000, 1, 0)
#'
#' # Conventional negative multinomial (NM2)
#' nm2 <- negmultinom(Animal ~ lnaadt + speed50 + ShouldWidth04 + AADTover10k,
#'                    group_var = "ID", data = washington_roads,
#'                    offset = "lnlength")
#' summary(nm2)
#'
#' # NM1: cluster totals follow an NB-1 distribution
#' nm1 <- negmultinom(Animal ~ lnaadt + speed50 + ShouldWidth04 + AADTover10k,
#'                    group_var = "ID", data = washington_roads,
#'                    family = "NM1", offset = "lnlength")
#' summary(nm1)
#'
#' # NM-P with ln(alpha) as a function of a site-level variable and
#' # cluster-robust standard errors
#' nmp <- negmultinom(Animal ~ lnaadt + speed50 + ShouldWidth04 + AADTover10k,
#'                    group_var = "ID", data = washington_roads,
#'                    family = "NMP", offset = "lnlength",
#'                    dis_param_formula_1 = ~ speed50,
#'                    stderr = "robust")
#' summary(nmp)
#'
#' # Marginal means and empirical Bayes expected frequencies
#' head(predict(nm2))
#' head(predict(nm2, method = "EB"))
#' }
#'
#' @references
#' Cameron, A. C., & Trivedi, P. K. (2013). Regression Analysis of Count Data
#' (2nd ed.). Cambridge University Press.
#'
#' Greene, W. (2008). Functional forms for the negative binomial model for count
#' data. Economics Letters, 99(3), 585-590.
#'
#' Hausman, J., Hall, B. H., & Griliches, Z. (1984). Econometric models for
#' count data with an application to the patents-R&D relationship.
#' Econometrica, 52(4), 909-938.
#'
#' @seealso \code{\link{repois}}, \code{\link{countreg}},
#'   \code{\link{predict.flexCountReg}}
#'
#' @importFrom maxLik maxLik maxControl
#' @importFrom stats model.frame model.matrix model.response model.offset
#'   terms delete.response complete.cases glm.fit poisson as.formula vcov
#' @export
negmultinom <- function(formula,
                        group_var,
                        data,
                        family = "NM2",
                        offset = NULL,
                        dis_param_formula_1 = NULL,
                        dis_param_formula_2 = NULL,
                        method = "BFGS",
                        max.iters = 1000,
                        start.vals = NULL,
                        stderr = "normal",
                        bootstraps = NULL,
                        verbose = FALSE) {

  call <- match.call()

  family <- .nm_clean_family(family)
  method <- .nm_clean_method(method)
  stderr <- match.arg(stderr, c("normal", "robust", "boot"))
  if (stderr == "boot" &&
      (!is.numeric(bootstraps) || length(bootstraps) != 1L || bootstraps < 2)) {
    stop("`bootstraps` must be an integer >= 2 when `stderr = \"boot\"`.")
  }
  if (family != "NMP" && !is.null(dis_param_formula_2)) {
    warning("`dis_param_formula_2` is only used when family = \"NMP\"; ",
            "it will be ignored.")
    dis_param_formula_2 <- NULL
  }

  prep <- .nm_prepare(formula = formula, group_var = group_var, data = data,
                      offset = offset,
                      dis_param_formula_1 = dis_param_formula_1,
                      dis_param_formula_2 = dis_param_formula_2,
                      family = family)

  start <- if (is.null(start.vals)) {
    .nm_start_values(prep, family)
  } else {
    if (length(start.vals) != length(prep$par_names)) {
      stop("`start.vals` must have length ", length(prep$par_names),
           " (", paste(prep$par_names, collapse = ", "), ").")
    }
    as.numeric(start.vals)
  }
  names(start) <- prep$par_names

  fit <- .nm_fit_engine(prep = prep, start = start, family = family,
                        method = method, max.iters = max.iters,
                        print.level = if (isTRUE(verbose)) 2L else 0L)

  # ---- Standard errors -----------------------------------------------------
  k <- length(fit$estimate)
  vc_hess <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  if (is.null(vc_hess) || any(!is.finite(vc_hess))) {
    vc_hess <- matrix(NA_real_, k, k)
    warning("The Hessian could not be inverted; standard errors are NA.")
  }
  dimnames(vc_hess) <- list(prep$par_names, prep$par_names)

  vc <- vc_hess
  if (stderr == "robust") {
    G <- .nm_grad(fit$estimate, prep = prep, family = family)
    meat <- crossprod(G)
    vc <- vc_hess %*% meat %*% vc_hess
    dimnames(vc) <- list(prep$par_names, prep$par_names)
  } else if (stderr == "boot") {
    vc <- .nm_bootstrap_vcov(prep = prep, start = fit$estimate,
                             family = family, method = method,
                             max.iters = max.iters, B = bootstraps)
  }
  se <- sqrt(pmax(diag(vc), 0))
  names(se) <- prep$par_names

  # ---- Fitted quantities ---------------------------------------------------
  comp <- .nm_components(fit$estimate, prep = prep, family = family)
  post_mean <- (comp$r + comp$Y) / (comp$r + comp$M)

  fit$estimate <- stats::setNames(as.numeric(fit$estimate), prep$par_names)
  fit$coefficients <- fit$estimate
  fit$beta_pred <- fit$estimate[seq_len(prep$n_beta)]
  fit$se <- se
  fit$vcov <- vc
  fit$stderr <- stderr
  fit$bootstraps <- bootstraps
  fit$bootstrapped_se <- NULL
  fit$alpha <- stats::setNames(comp$alpha, prep$group_levels)
  fit$p <- stats::setNames(comp$p, prep$group_levels)
  fit$r <- stats::setNames(comp$r, prep$group_levels)
  fit$frailty.posterior <- stats::setNames(post_mean, prep$group_levels)
  fit$fitted.values <- unname(comp$mu)
  fit$predictions <- comp$mu
  fit$eb.fitted <- unname(comp$mu * post_mean[prep$g])
  fit$linear.predictors <- comp$eta
  fit$observed <- prep$y
  fit$residuals <- prep$y - comp$mu
  fit$LL <- fit$maximum
  fit$logLik <- fit$maximum
  fit$nobs <- length(prep$y)
  fit$n_groups <- prep$G
  fit$aic <- -2 * fit$maximum + 2 * k
  fit$bic <- -2 * fit$maximum + log(length(prep$y)) * k
  fit$converged <- fit$code %in% c(0L, 1L, 2L, 8L)
  fit$family <- family
  fit$model <- family
  fit$modelType <- "NegMultinom"
  fit$method <- method
  fit$formula <- formula
  fit$dis_param_formula_1 <- dis_param_formula_1
  fit$dis_param_formula_2 <- dis_param_formula_2
  fit$offset <- offset
  fit$group_var <- group_var
  fit$n_beta <- prep$n_beta
  fit$n_alpha <- prep$n_alpha
  fit$n_p <- prep$n_p
  fit$rows_used <- prep$rows_used
  fit$call <- call

  .createFlexCountReg(model = fit, data = data, call = call, formula = formula)
}


# ============================================================================
# Internal helpers
# ============================================================================

.nm_clean_family <- function(family) {
  if (!is.character(family) || length(family) != 1L) {
    stop("`family` must be a single string: \"NM2\", \"NM1\", or \"NMP\".")
  }
  fam <- toupper(gsub("[^[:alnum:]]", "", family))
  if (!fam %in% c("NM1", "NM2", "NMP")) {
    stop("`family` must be one of \"NM2\", \"NM1\", or \"NMP\".")
  }
  fam
}

.nm_clean_method <- function(method) {
  m <- toupper(gsub("[^[:alnum:]]", "", method))
  if (m == "NELDERMEAD") m <- "NM"
  if (!m %in% c("NM", "BFGS", "BFGSR", "NR", "BHHH", "CG", "SANN")) {
    stop("Unsupported optimization method: ", method)
  }
  m
}

# Build a design matrix for a one-sided formula on the rows used, and
# collapse it to the cluster level (within-cluster means).
.nm_cluster_design <- function(fml, data_used, g, G, label) {
  fml <- stats::as.formula(fml)
  if (length(fml) == 3L) {
    stop("`", label, "` must be a one-sided formula, e.g. ~ x1 + x2.")
  }
  Z <- stats::model.matrix(fml, stats::model.frame(fml, data_used,
                                                   na.action = stats::na.pass))
  n_g <- tabulate(g, nbins = G)
  Zc <- rowsum(Z, g, reorder = TRUE) / n_g
  dev <- abs(Z - Zc[g, , drop = FALSE])
  if (any(dev > 1e-8 * pmax(1, abs(Z)))) {
    bad <- colnames(Z)[apply(dev > 1e-8 * pmax(1, abs(Z)), 2, any)]
    warning("Variable(s) in `", label, "` vary within clusters (",
            paste(bad, collapse = ", "), "). The within-cluster mean of the ",
            "linear predictor is used because the parameter applies to the ",
            "shared cluster frailty.", call. = FALSE)
  }
  colnames(Zc) <- colnames(Z)
  Zc
}

.nm_prepare <- function(formula, group_var, data, offset,
                        dis_param_formula_1, dis_param_formula_2, family) {
  if (!inherits(data, "data.frame")) stop("`data` must be a data frame.")
  if (missing(group_var) || length(group_var) < 1L) {
    stop("`group_var` must be supplied.")
  }
  group_var <- as.character(group_var)
  miss <- setdiff(group_var, names(data))
  if (length(miss) > 0L) {
    stop("Grouping variable(s) not found in `data`: ",
         paste(miss, collapse = ", "))
  }
  if (!is.null(offset)) {
    if (!is.character(offset)) {
      stop("`offset` must be a character vector naming variable(s) in `data`.")
    }
    miss <- setdiff(offset, names(data))
    if (length(miss) > 0L) {
      stop("Offset variable(s) not found in `data`: ",
           paste(miss, collapse = ", "))
    }
  }

  # Complete cases across every variable used anywhere in the model so all
  # design matrices refer to the same rows.
  vars <- unique(c(all.vars(formula),
                   if (!is.null(dis_param_formula_1))
                     all.vars(stats::as.formula(dis_param_formula_1)),
                   if (!is.null(dis_param_formula_2))
                     all.vars(stats::as.formula(dis_param_formula_2)),
                   group_var, offset))
  vars <- intersect(vars, names(data))
  ok <- stats::complete.cases(data[, vars, drop = FALSE])
  rows_used <- which(ok)
  data_used <- data[ok, , drop = FALSE]
  if (nrow(data_used) == 0L) stop("No complete observations remain.")

  mf <- stats::model.frame(formula, data_used, na.action = stats::na.fail)
  y <- as.numeric(stats::model.response(mf))
  if (any(y < 0)) stop("The response variable must be nonnegative.")
  if (any(abs(y - round(y)) > 1e-8)) {
    warning("The response variable is not integer-valued.")
  }
  X <- stats::model.matrix(stats::terms(mf), mf)

  off <- stats::model.offset(mf)
  if (is.null(off)) off <- rep(0, length(y))
  if (!is.null(offset)) {
    off <- off + rowSums(as.matrix(data_used[, offset, drop = FALSE]))
  }

  gvals <- lapply(group_var, function(v) data_used[[v]])
  gfac <- if (length(gvals) == 1L) {
    factor(gvals[[1L]])
  } else {
    interaction(gvals, drop = TRUE, sep = ":")
  }
  gfac <- droplevels(gfac)
  g <- as.integer(gfac)
  G <- nlevels(gfac)

  n_g <- tabulate(g, nbins = G)
  one_col <- matrix(1, nrow = G, ncol = 1L, dimnames = list(NULL, "(Intercept)"))

  if (is.null(dis_param_formula_1)) {
    Zc <- one_col
    alpha_names <- "ln(alpha)"
  } else {
    Zc <- .nm_cluster_design(dis_param_formula_1, data_used, g, G,
                             "dis_param_formula_1")
    alpha_names <- paste0("ln(alpha):", colnames(Zc))
  }

  if (family == "NMP") {
    if (is.null(dis_param_formula_2)) {
      Wc <- one_col
      p_names <- "p"
    } else {
      Wc <- .nm_cluster_design(dis_param_formula_2, data_used, g, G,
                               "dis_param_formula_2")
      p_names <- paste0("p:", colnames(Wc))
    }
  } else {
    Wc <- NULL
    p_names <- character(0)
  }

  list(
    y = y, X = X, offset = off, g = g, G = G, n_g = n_g,
    group_levels = levels(gfac),
    Zc = Zc, Wc = Wc,
    n_beta = ncol(X), n_alpha = ncol(Zc),
    n_p = if (is.null(Wc)) 0L else ncol(Wc),
    par_names = c(colnames(X), alpha_names, p_names),
    Y = as.numeric(rowsum(y, g, reorder = TRUE)),
    lgy = as.numeric(rowsum(lgamma(y + 1), g, reorder = TRUE)),
    rows_used = rows_used
  )
}

# Compute all per-cluster quantities at a parameter vector.
.nm_components <- function(par, prep, family) {
  kb <- prep$n_beta
  ka <- prep$n_alpha
  beta <- par[seq_len(kb)]
  gam <- par[kb + seq_len(ka)]

  eta <- drop(prep$X %*% beta) + prep$offset
  mu <- exp(eta)
  M <- as.numeric(rowsum(mu, prep$g, reorder = TRUE))

  lnalpha <- drop(prep$Zc %*% gam)
  p <- switch(family,
              NM2 = rep(2, prep$G),
              NM1 = rep(1, prep$G),
              NMP = drop(prep$Wc %*% par[kb + ka + seq_len(prep$n_p)]))

  logM <- log(M)
  logr <- (2 - p) * logM - lnalpha
  r <- exp(logr)

  list(eta = eta, mu = mu, M = M, Y = prep$Y, lnalpha = lnalpha,
       alpha = exp(lnalpha), p = p, logM = logM, logr = logr, r = r)
}

# Cluster-level log-likelihood contributions (vector of length G).
.nm_loglik <- function(par, prep, family) {
  cp <- .nm_components(par, prep, family)
  if (any(!is.finite(cp$eta)) || any(cp$eta > 700) ||
      any(!is.finite(cp$logM)) || any(!is.finite(cp$r)) || any(cp$r <= 0)) {
    return(NA_real_)
  }
  syeta <- as.numeric(rowsum(prep$y * cp$eta, prep$g, reorder = TRUE))
  ll <- lgamma(cp$r + cp$Y) - lgamma(cp$r) + cp$r * cp$logr -
    (cp$r + cp$Y) * log(cp$r + cp$M) + syeta - prep$lgy
  if (any(!is.finite(ll))) return(NA_real_)
  ll
}

# Analytic gradient: G x k matrix of cluster-level scores.
.nm_grad <- function(par, prep, family) {
  cp <- .nm_components(par, prep, family)
  r <- cp$r
  rpM <- r + cp$M
  # d LL / d r and d LL / d M (cluster level)
  dr <- digamma(r + cp$Y) - digamma(r) + cp$logr + 1 - log(rpM) -
    (r + cp$Y) / rpM
  dM <- -(r + cp$Y) / rpM
  # Total derivative w.r.t. M, including r's dependence on M:
  # d r / d M = (2 - p) r / M
  dMtot <- dM + dr * (2 - cp$p) * r / cp$M
  deta <- prep$y + cp$mu * dMtot[prep$g]

  g_beta <- rowsum(prep$X * deta, prep$g, reorder = TRUE)
  g_alpha <- (-dr * r) * prep$Zc           # d r / d ln(alpha) = -r
  out <- cbind(g_beta, g_alpha)
  if (family == "NMP") {
    g_p <- (-dr * r * cp$logM) * prep$Wc   # d r / d p = -r ln(M)
    out <- cbind(out, g_p)
  }
  out <- unname(out)
  colnames(out) <- prep$par_names
  out
}

.nm_start_values <- function(prep, family) {
  fit0 <- tryCatch(
    stats::glm.fit(prep$X, prep$y, family = stats::poisson(),
                   offset = prep$offset),
    error = function(e) NULL)
  beta0 <- if (is.null(fit0)) rep(0, prep$n_beta) else fit0$coefficients
  beta0[!is.finite(beta0)] <- 0

  mu0 <- exp(drop(prep$X %*% beta0) + prep$offset)
  M0 <- as.numeric(rowsum(mu0, prep$g, reorder = TRUE))
  p0 <- switch(family, NM2 = 2, NM1 = 1, NMP = 2)
  # Moment estimate from cluster totals: Var(Y_R) = M + alpha M^p
  alpha0 <- sum((prep$Y - M0)^2 - prep$Y) / sum(M0^p0)
  if (!is.finite(alpha0) || alpha0 <= 1e-3) alpha0 <- 0.1
  lnalpha0 <- log(alpha0)

  a_start <- rep(0, prep$n_alpha)
  int_a <- which(colnames(prep$Zc) == "(Intercept)")
  if (length(int_a) == 1L) a_start[int_a] <- lnalpha0

  p_start <- numeric(0)
  if (family == "NMP") {
    p_start <- rep(0, prep$n_p)
    int_p <- which(colnames(prep$Wc) == "(Intercept)")
    if (length(int_p) == 1L) p_start[int_p] <- p0
  }
  c(beta0, a_start, p_start)
}

.nm_fit_engine <- function(prep, start, family, method, max.iters,
                           print.level = 0L) {
  ll0 <- .nm_loglik(start, prep, family)
  if (anyNA(ll0)) {
    stop("The log-likelihood is not finite at the starting values. ",
         "Try supplying `start.vals`.")
  }
  maxLik::maxLik(
    logLik = function(par) .nm_loglik(par, prep, family),
    grad = function(par) .nm_grad(par, prep, family),
    start = start,
    method = method,
    control = maxLik::maxControl(iterlim = max.iters, printLevel = print.level)
  )
}

.nm_subset_prep <- function(prep, clusters) {
  rows <- unlist(lapply(clusters, function(cl) which(prep$g == cl)),
                 use.names = FALSE)
  new_g <- rep(seq_along(clusters), times = prep$n_g[clusters])
  out <- prep
  out$y <- prep$y[rows]
  out$X <- prep$X[rows, , drop = FALSE]
  out$offset <- prep$offset[rows]
  out$g <- new_g
  out$G <- length(clusters)
  out$n_g <- prep$n_g[clusters]
  out$Zc <- prep$Zc[clusters, , drop = FALSE]
  if (!is.null(prep$Wc)) out$Wc <- prep$Wc[clusters, , drop = FALSE]
  out$Y <- prep$Y[clusters]
  out$lgy <- prep$lgy[clusters]
  out
}

.nm_bootstrap_vcov <- function(prep, start, family, method, max.iters, B) {
  k <- length(start)
  est <- matrix(NA_real_, nrow = B, ncol = k)
  for (b in seq_len(B)) {
    cl <- sample.int(prep$G, prep$G, replace = TRUE)
    bp <- .nm_subset_prep(prep, cl)
    bf <- tryCatch(
      .nm_fit_engine(bp, start = start, family = family, method = method,
                     max.iters = max.iters),
      error = function(e) NULL)
    if (!is.null(bf) && all(is.finite(bf$estimate))) est[b, ] <- bf$estimate
  }
  est <- est[stats::complete.cases(est), , drop = FALSE]
  if (nrow(est) < 2L) {
    warning("Fewer than 2 bootstrap replications succeeded; ",
            "standard errors are NA.")
    return(matrix(NA_real_, k, k,
                  dimnames = list(prep$par_names, prep$par_names)))
  }
  if (nrow(est) < B) {
    warning(B - nrow(est), " of ", B, " bootstrap replications failed.")
  }
  vc <- stats::cov(est)
  dimnames(vc) <- list(prep$par_names, prep$par_names)
  attr(vc, "successful_bootstraps") <- nrow(est)
  vc
}

# Marginal means or EB expected frequencies for (new) data, used by
# predict.flexCountReg().
.nm_predict <- function(model, data, method = "Standard") {
  method <- toupper(method)
  coefs <- as.numeric(unlist(model$estimate))
  form_fixed <- stats::delete.response(stats::terms(model$formula))
  mf <- stats::model.frame(form_fixed, data, na.action = stats::na.pass)
  X <- stats::model.matrix(form_fixed, mf)
  off <- stats::model.offset(mf)
  if (is.null(off)) off <- rep(0, nrow(data))
  if (!is.null(model$offset) && all(model$offset %in% names(data))) {
    off <- off + rowSums(as.matrix(data[, model$offset, drop = FALSE]))
  }
  mu <- unname(exp(drop(X %*% coefs[seq_len(ncol(X))]) + off))

  if (method %in% c("STANDARD", "MARGINAL")) return(mu)
  if (method != "EB") {
    stop("For negative multinomial models, `method` must be \"Standard\" ",
         "or \"EB\".")
  }

  y_name <- all.vars(model$formula)[1L]
  need <- c(y_name, model$group_var)
  if (!all(need %in% names(data))) {
    stop("EB predictions require the outcome and grouping variable(s) (",
         paste(need, collapse = ", "), ") in `newdata`.")
  }
  prep <- .nm_prepare(formula = model$formula, group_var = model$group_var,
                      data = data, offset = model$offset,
                      dis_param_formula_1 = model$dis_param_formula_1,
                      dis_param_formula_2 = model$dis_param_formula_2,
                      family = model$family)
  cp <- .nm_components(coefs, prep, model$family)
  out <- rep(NA_real_, nrow(data))
  out[prep$rows_used] <- cp$mu * ((cp$r + cp$Y) / (cp$r + cp$M))[prep$g]
  out
}
