# Methods for countreg_hier() fits.

#' @importFrom nlme fixef
#' @export
nlme::fixef

#' @importFrom nlme ranef
#' @export
nlme::ranef

#' @importFrom nlme VarCorr
#' @export
nlme::VarCorr

# Which random-effect terms does `re.form` keep? NULL: all; NA or ~0: none;
# a formula: the terms whose label (e.g. "1 | site") it contains.
.fcr_hier_reform <- function(re.form, terms) {
  if (is.null(re.form)) return(rep(TRUE, length(terms)))
  if (length(re.form) == 1L && !inherits(re.form, "formula") &&
      is.na(re.form)) {
    return(rep(FALSE, length(terms)))
  }
  if (!inherits(re.form, "formula")) {
    stop("`re.form` must be NULL (all random effects), NA (none) or a ",
         "formula of random-effect terms, such as ~ (1 | site).",
         call. = FALSE)
  }
  bars <- reformulas::findbars(re.form)
  if (!length(bars)) return(rep(FALSE, length(terms)))
  labs <- vapply(bars, .fcr_hier_label, "")
  have <- vapply(terms, `[[`, "", "label")
  bad <- setdiff(labs, have)
  if (length(bad)) {
    stop("`re.form` term(s) ", toString(paste0("(", bad, ")")), " are not ",
         "in the model. Its terms: ", toString(paste0("(", have, ")")), ".",
         call. = FALSE)
  }
  have %in% labs
}

# Fixed part, distribution parameters and random-effect contributions for
# the rows of `d`.
.fcr_hier_parts <- function(object, d) {
  m <- object$model
  fd <- forge_new_data(m$blueprint, d)
  p <- m$par_internal
  n <- nrow(d)
  eta0 <- as.vector(fd$X %*% p[seq_len(m$k_beta)]) + fd$offset
  pos <- m$k_beta
  alpha <- sigma <- NULL
  if (m$n_alpha) {
    a <- p[pos + seq_len(m$n_alpha)]
    alpha <- if (is.null(fd$aux$dis1)) rep(exp(a), n) else
      exp(as.vector(fd$aux$dis1 %*% a))
    pos <- pos + m$n_alpha
  }
  if (m$n_sigma) {
    s <- p[pos + seq_len(m$n_sigma)]
    sigma <- if (is.null(fd$aux$dis2)) rep(exp(s), n) else
      exp(as.vector(fd$aux$dis2 %*% s))
  }
  env <- environment(object$formula) %||% globalenv()
  Sig <- .fcr_hier_sigmas(m$re_terms, p[m$k_fix + seq_len(m$n_theta)])
  re <- lapply(seq_along(m$re_terms), function(i) {
    t <- m$re_terms[[i]]
    Xl <- .fcr_hier_term_matrix(t, d)
    g <- .fcr_hier_group(t$gexpr, d, env)
    j <- match(as.character(g), t$levels)
    B <- matrix(m$b[t$uoff + seq_len(t$J * t$q)], t$J, t$q, byrow = TRUE)
    contrib <- rowSums(Xl * B[ifelse(is.na(j), 1L, j), , drop = FALSE])
    contrib[is.na(j)] <- 0
    contrib[is.na(g)] <- NA_real_
    list(X = Xl, j = j, missing_group = is.na(g), new = is.na(j) & !is.na(g),
         contrib = contrib, var = rowSums((Xl %*% Sig[[i]]) * Xl),
         Sigma = Sig[[i]])
  })
  list(eta0 = eta0, alpha = alpha, sigma = sigma, re = re)
}

#' Predictions from a hierarchical count model
#'
#' @param object a \code{\link{countreg_hier}()} fit.
#' @param newdata optional data frame (default: the estimation rows).
#' @param type \code{"response"} (expected count) or \code{"link"} (log of
#'   the mean parameter).
#' @param re.form which random effects to condition on: \code{NULL} (all, the
#'   default), \code{NA} (none, population level) or a formula of terms,
#'   such as \code{~ (1 | county)}.
#' @param allow.new.levels allow groups not seen in estimation; their
#'   effects are set to zero (or integrated out with \code{marginal}).
#' @param marginal integrate the random effects that are not conditioned on
#'   (those left out by \code{re.form}, and new groups) out of the expected
#'   count, instead of setting them to zero. For the log link the marginal
#'   mean is \eqn{\exp(\eta + z'\Sigma z / 2)}, larger than
#'   \eqn{\exp(\eta)}.
#' @param ... unused.
#'
#' @details Conditional predictions use the conditional modes of the random
#'   effects (as \code{lme4} does). With \code{marginal = FALSE} (the
#'   default, as in \code{lme4}) the effects that are not conditioned on are
#'   zero, so a population-level prediction is the count of a median group,
#'   not the mean over groups.
#' @returns A numeric vector with one value per row of \code{newdata}.
#' @export
predict.countreg_hier <- function(object, newdata = NULL,
                                  type = c("response", "link"),
                                  re.form = NULL, allow.new.levels = TRUE,
                                  marginal = FALSE, ...) {
  type <- match.arg(type)
  m <- object$model
  d <- if (is.null(newdata)) .fcr_training_rows(object) else
    as.data.frame(newdata)
  parts <- .fcr_hier_parts(object, d)
  incl <- .fcr_hier_reform(re.form, m$re_terms)
  eta <- parts$eta0
  extra <- numeric(length(eta))
  for (i in seq_along(parts$re)) {
    r <- parts$re[[i]]
    if (incl[i]) {
      if (any(r$new) && !allow.new.levels) {
        stop("`newdata` has groups of (", m$re_terms[[i]]$label, ") not ",
             "seen in estimation; set `allow.new.levels = TRUE`.",
             call. = FALSE)
      }
      eta <- eta + r$contrib
      extra[r$new] <- extra[r$new] + r$var[r$new]
    } else {
      extra <- extra + r$var
    }
  }
  if (isTRUE(marginal)) eta_mean <- eta + extra / 2 else eta_mean <- eta
  if (type == "link") return(unname(eta_mean))
  unname(.fcr_family_mean(m$family, exp(eta_mean), parts$alpha, parts$sigma))
}

#' @export
fitted.countreg_hier <- function(object, ...) {
  unname(object$model$fitted.values)
}

#' @export
residuals.countreg_hier <- function(object, type = c("response", "pearson"),
                                    ...) {
  type <- match.arg(type)
  m <- object$model
  r <- m$y - m$fitted.values
  if (type == "pearson") {
    v <- .fcr_family_variance(m$family, exp(m$eta), m$alpha, m$sigma)
    if (is.null(v)) {
      stop("Pearson residuals need the family's variance, which is not ",
           "available for family \"", m$family, "\".", call. = FALSE)
    }
    r <- r / sqrt(v)
  }
  unname(r)
}

#' Simulate responses from a hierarchical count model
#'
#' @param object a \code{\link{countreg_hier}()} fit.
#' @param nsim number of simulated response vectors.
#' @param seed optional random seed.
#' @param newdata optional data frame (default: the estimation rows).
#' @param use_re \code{"conditional"} keeps the estimated effects of the
#'   groups (conditional modes); \code{"new"} draws new effects for every
#'   group from their estimated distribution.
#' @param ... unused.
#' @returns A data frame with \code{nsim} columns.
#' @export
simulate.countreg_hier <- function(object, nsim = 1, seed = NULL,
                                   newdata = NULL,
                                   use_re = c("conditional", "new"), ...) {
  use_re <- match.arg(use_re)
  if (!is.null(seed)) set.seed(seed)
  m <- object$model
  d <- if (is.null(newdata)) .fcr_training_rows(object) else
    as.data.frame(newdata)
  parts <- .fcr_hier_parts(object, d)
  n <- nrow(d)
  out <- vapply(seq_len(nsim), function(s) {
    eta <- parts$eta0
    for (i in seq_along(parts$re)) {
      r <- parts$re[[i]]
      if (use_re == "conditional") {
        eta <- eta + r$contrib
      } else {
        t <- m$re_terms[[i]]
        g <- .fcr_hier_group(t$gexpr, d,
                             environment(object$formula) %||% globalenv())
        gi <- as.integer(factor(g))
        ng <- max(gi, na.rm = TRUE)
        L <- tryCatch(t(chol(r$Sigma)), error = function(e) {
          e <- eigen(r$Sigma, symmetric = TRUE)
          e$vectors %*% diag(sqrt(pmax(e$values, 0)), t$q)
        })
        B <- matrix(stats::rnorm(ng * t$q), ng, t$q) %*% t(L)
        eta <- eta + rowSums(r$X * B[gi, , drop = FALSE])
      }
    }
    .fcr_family_r(m$family, n, exp(eta), parts$alpha, parts$sigma)
  }, numeric(n))
  out <- as.data.frame(matrix(out, n, nsim))
  names(out) <- paste0("sim_", seq_len(nsim))
  out
}

#' Fixed effects, random effects and variance components
#'
#' \code{fixef()} returns the fixed coefficients of the mean,
#' \code{ranef()} the conditional modes of the random effects (one data
#' frame per grouping factor, with their posterior covariances in the
#' attribute \code{"postVar"}), and \code{VarCorr()} the estimated
#' covariance matrices of the random effects.
#'
#' @param object a \code{\link{countreg_hier}()} fit.
#' @param x a \code{VarCorr()} result.
#' @param condVar include the posterior covariances.
#' @param digits number of significant digits to print.
#' @param ... unused.
#' @name countreg_hier-effects
#' @returns \code{fixef()}: a named vector. \code{ranef()}: a named list of
#'   data frames. \code{VarCorr()}: a named list of covariance matrices with
#'   attributes \code{"stddev"} and \code{"correlation"}.
#' @examples
#' \donttest{
#' data("washington_roads")
#' h <- countreg_hier(Total_crashes ~ lnaadt + offset(lnlength) + (1 | ID),
#'                    data = washington_roads, family = "POISSON")
#' fixef(h)
#' head(ranef(h)$ID)
#' VarCorr(h)
#' }
NULL

#' @rdname countreg_hier-effects
#' @export
fixef.countreg_hier <- function(object, ...) {
  object$model$estimate[seq_len(object$model$k_beta)]
}

#' @rdname countreg_hier-effects
#' @export
ranef.countreg_hier <- function(object, condVar = TRUE, ...) {
  m <- object$model
  out <- list()
  for (i in seq_along(m$re_terms)) {
    t <- m$re_terms[[i]]
    B <- matrix(m$b[t$uoff + seq_len(t$J * t$q)], t$J, t$q, byrow = TRUE)
    df <- as.data.frame(B, row.names = t$levels)
    names(df) <- t$cnms
    if (condVar && !is.null(m$condvar)) attr(df, "postVar") <- m$condvar[[i]]
    nm <- t$group
    if (!is.null(out[[nm]]) && identical(rownames(out[[nm]]), t$levels)) {
      # Uncorrelated terms of the same factor, as from (x || site)
      pv <- attr(out[[nm]], "postVar")
      out[[nm]] <- cbind(out[[nm]], df)
      attr(out[[nm]], "postVar") <- NULL
    } else {
      out[[nm]] <- df
    }
  }
  out
}

#' @rdname countreg_hier-effects
#' @export
VarCorr.countreg_hier <- function(x, ...) {
  m <- x$model
  S <- .fcr_hier_sigmas(m$re_terms, m$par_internal[m$k_fix +
                                                      seq_len(m$n_theta)])
  out <- lapply(seq_along(S), function(i) {
    s <- sqrt(pmax(diag(S[[i]]), 0))
    R <- S[[i]] / outer(s, s)
    R[!is.finite(R)] <- NA_real_
    diag(R) <- 1
    structure(S[[i]], stddev = s, correlation = R)
  })
  names(out) <- vapply(m$re_terms, `[[`, "", "group")
  structure(out, class = "VarCorr.countreg_hier")
}

#' @rdname countreg_hier-effects
#' @export
print.VarCorr.countreg_hier <- function(x, digits = 4, ...) {
  rows <- list()
  for (i in seq_along(x)) {
    s <- attr(x[[i]], "stddev")
    R <- attr(x[[i]], "correlation")
    nm <- names(s)
    for (k in seq_along(s)) {
      cor <- if (k > 1L) paste(formatC(R[k, seq_len(k - 1L)], digits = 2,
                                       format = "f"), collapse = " ") else ""
      rows[[length(rows) + 1L]] <- data.frame(
        Groups = if (k == 1L) names(x)[i] else "", Name = nm[k],
        Variance = signif(s[k]^2, digits), Std.Dev. = signif(s[k], digits),
        Corr = cor, stringsAsFactors = FALSE)
    }
  }
  tab <- do.call(rbind, rows)
  if (all(tab$Corr == "")) tab$Corr <- NULL
  print(tab, row.names = FALSE, right = FALSE)
  invisible(x)
}

#' @export
print.countreg_hier <- function(x, ...) {
  m <- x$model
  cat("Hierarchical count model (", .fcr_family_label(m$family), ", ",
      .fcr_hier_engine_label(m), ")\n", sep = "")
  cat("Formula: ", paste(deparse(x$formula, width.cutoff = 500L),
                         collapse = " "), "\n", sep = "")
  cat("Log-likelihood: ", format(m$maximum, nsmall = 2), " on ",
      length(m$estimate), " parameters; ", m$nobs, " observations\n\n",
      sep = "")
  cat("Fixed coefficients and distribution parameters:\n")
  print(m$estimate[seq_len(m$k_fix)], digits = 4)
  cat("\nRandom effects:\n")
  print(VarCorr(x))
  if (isTRUE(m$singular)) cat("\nSingular fit: a variance is estimated at or",
                              "near zero.\n")
  invisible(x)
}

.fcr_hier_engine_label <- function(m) {
  if (m$engine == "aghq") {
    if (m$nAGQ == 1L) "Laplace approximation" else
      paste0("adaptive Gauss-Hermite quadrature, ", m$nAGQ, " points")
  } else {
    "Laplace approximation"
  }
}

#' @export
summary.countreg_hier <- function(object, level = 0.95, ...) {
  m <- object$model
  est <- m$estimate
  se <- .fcr_se(m$vcov, length(est))
  z <- est / se
  tab <- data.frame(Estimate = est, `Std. Error` = se, `z value` = z,
                    `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)),
                    check.names = FALSE)
  rownames(tab) <- names(est)
  fix <- seq_len(m$k_fix)
  vc <- tab[-fix, c("Estimate", "Std. Error"), drop = FALSE]
  structure(list(call = object$call, family = m$family,
                 engine = .fcr_hier_engine_label(m),
                 coefficients = tab[fix, , drop = FALSE], varcomp = vc,
                 VarCorr = VarCorr(object), logLik = m$maximum,
                 AIC = stats::AIC(object), BIC = stats::BIC(object),
                 nobs = m$nobs, n_groups = m$n_groups,
                 singular = isTRUE(m$singular),
                 converged = isTRUE(m$converged)),
            class = "summary.countreg_hier")
}

#' @export
print.summary.countreg_hier <- function(x, digits = 4, ...) {
  cat("Hierarchical count model: ", .fcr_family_label(x$family), "\n",
      sep = "")
  cat("Integration: ", x$engine, "\n", sep = "")
  cat("Call: ", paste(deparse(x$call, width.cutoff = 500L), collapse = " "),
      "\n\n", sep = "")
  cat("Random effects:\n")
  print(x$VarCorr, digits = digits)
  cat("Number of obs: ", x$nobs, "; groups: ",
      paste0(names(x$n_groups), ", ", x$n_groups, collapse = "; "), "\n\n",
      sep = "")
  cat("Fixed coefficients and distribution parameters:\n")
  stats::printCoefmat(x$coefficients, digits = digits, signif.legend = FALSE)
  cat("\nRandom-effect standard deviations and correlations (Wald SEs):\n")
  print(signif(as.matrix(x$varcomp), digits))
  cat("\nLog-likelihood: ", format(x$logLik, nsmall = 2), "   AIC: ",
      format(x$AIC, nsmall = 2), "   BIC: ", format(x$BIC, nsmall = 2), "\n",
      sep = "")
  if (x$singular) {
    cat("Singular fit: a random-effect variance is estimated at or near",
        "zero; consider removing that effect.\n")
  }
  if (!x$converged) cat("Warning: the optimizer may not have converged.\n")
  invisible(x)
}

#' Likelihood-ratio test of a random-effect term
#'
#' Refits a \code{\link{countreg_hier}()} model without one random-effect
#' term and compares the log-likelihoods. The null value of a variance is
#' on the boundary of its parameter space, so for a term with one effect
#' the p-value is that of an equal mixture of \eqn{\chi^2_0} and
#' \eqn{\chi^2_1} (Self and Liang 1987): half the naive p-value. For a term
#' with \eqn{q > 1} effects the mixture of \eqn{\chi^2_{k-q}} to
#' \eqn{\chi^2_k} (\eqn{k = q(q + 1)/2}) is approximated by an equal mixture
#' of \eqn{\chi^2_{k-1}} and \eqn{\chi^2_k}.
#'
#' @param object a \code{countreg_hier()} fit.
#' @param term the term to test, as its label (e.g. \code{"1 | site"});
#'   by default the last one.
#' @returns A data frame with the log-likelihoods, the statistic, its
#'   degrees of freedom and the p-value.
#' @references Self, S.G. and Liang, K.-Y. (1987). Asymptotic properties of
#'   maximum likelihood estimators and likelihood ratio tests under
#'   nonstandard conditions. \emph{Journal of the American Statistical
#'   Association}, 82(398), 605-610.
#' @export
re_lrtest <- function(object, term = NULL) {
  if (!inherits(object, "countreg_hier")) {
    stop("`object` must be a countreg_hier() fit.", call. = FALSE)
  }
  m <- object$model
  labs <- vapply(m$re_terms, `[[`, "", "label")
  if (is.null(term)) term <- labs[length(labs)]
  term <- gsub("^\\(|\\)$", "", trimws(term))
  i <- match(term, labs)
  if (is.na(i)) {
    stop("No random-effect term (", term, "). Terms: ",
         toString(paste0("(", labs, ")")), ".", call. = FALSE)
  }
  if (is.null(object$data)) {
    stop("The model was fit with keep_data = FALSE; refit it with the ",
         "data to test a term.", call. = FALSE)
  }
  t <- m$re_terms[[i]]
  rest <- labs[-i]
  fixed <- m$fixed_formula
  args <- list(data = object$data, family = m$family_spec %||% m$family,
               dis_param_formula_1 = m$dis_param_formula_1,
               dis_param_formula_2 = m$dis_param_formula_2,
               weights = if (any(m$weights != 1)) m$weights,
               ndraws = m$ndraws, keep_data = FALSE)
  start <- m$estimate[seq_len(m$k_fix)]
  if (length(rest)) {
    rhs <- paste(c(deparse(fixed[[3]], width.cutoff = 500L),
                   paste0("(", rest, ")")), collapse = " + ")
    f0 <- stats::as.formula(paste(deparse(fixed[[2]]), "~", rhs),
                            env = environment(object$formula))
    null <- suppressMessages(do.call(countreg_hier, c(list(
      formula = f0, engine = if (m$engine == "aghq" &&
                                   length(rest) == 1L) "aghq" else "laplace",
      start.vals = unname(start)), args)))
  } else {
    null <- suppressMessages(do.call(countreg, c(list(
      formula = fixed, method = "BFGS", start.vals = unname(start)),
      args[setdiff(names(args), "keep_data")], list(keep_data = FALSE))))
    if (m$engine == "laplace" || m$nAGQ == 1L) {
      message("Note: the Laplace approximation of the full model is ",
              "compared with the exact likelihood of the model without ",
              "random effects.")
    }
  }
  ll1 <- m$maximum
  ll0 <- as.numeric(null$model$maximum)
  stat <- max(0, 2 * (ll1 - ll0))
  k <- t$ntheta
  p <- if (t$q == 1L) {
    0.5 * stats::pchisq(stat, 1, lower.tail = FALSE)
  } else {
    0.5 * stats::pchisq(stat, k - 1, lower.tail = FALSE) +
      0.5 * stats::pchisq(stat, k, lower.tail = FALSE)
  }
  data.frame(term = paste0("(", term, ")"), logLik_without = ll0,
             logLik_with = ll1, statistic = stat, df = k, p.value = p,
             row.names = NULL)
}
