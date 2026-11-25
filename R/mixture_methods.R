# Methods for finite-mixture and latent-class models (countreg_mix objects).
# Statistical methods are implemented once for "countreg_mix"; countreg_lc
# objects inherit them.

# Settings needed to evaluate a fitted mixture on new rows.
.fcr_fm_pspec <- function(object) {
  m <- object$model
  haltons <- randtoolbox::halton(m$ndraws %||% 1500)
  list(blueprint = m$blueprint, map = m$map, K = m$K,
       mix_mode = m$mix_mode, fixed_mixing = m$fixed_mixing,
       haltons = haltons, normed_haltons = stats::qnorm(haltons),
       engine = m$engine %||% "quadrature", families = m$family,
       labels = m$labels,
       rp = if (!is.null(m$rpar)) utils::modifyList(
         m$rpar, list(n_units = m$rpar$n_units %||% m$n_units)))
}

# Component values and log prior probabilities for the rows of `data`.
# With random parameters, `train_unit` (the training unit of each row)
# gives the rows the estimation draws of their units instead of the shared
# prediction draws.
.fcr_fm_rows <- function(object, data, theta = object$model$estimate,
                         train_unit = NULL) {
  ps <- .fcr_fm_pspec(object)
  des <- .fcr_fm_design(ps, data)
  pars <- .fcr_fm_unpack(theta, ps$map)
  draws <- if (!is.null(ps$rp) && !is.null(train_unit)) {
    .fcr_fm_rp_train_draws(ps$rp, train_unit, ps$K)
  }
  cvs <- lapply(seq_len(ps$K), function(k) {
    cv <- .fcr_fm_comp_values(des, pars, k)
    # Random parameters: the component's expected counts at R draws, and
    # its exact moments
    if (!is.null(ps$rp)) {
      md <- .fcr_fm_rp_mu_draws(ps$rp, des, pars, k, draws)
      cv$mu_draws <- md$mu_draws
      cv$mu_fixed <- md$mu_fixed
      cv$draw_w <- md$draw_w
      cv$rp_moments <- .fcr_fm_rp_moments(ps$rp, des, pars, k)
    }
    cv
  })
  logpi <- .fcr_fm_logprior(ps, pars, des$Z, des$n)
  list(ps = ps, des = des, pars = pars, cvs = cvs, logpi = logpi,
       n = des$n)
}

# Row log-probabilities of observed counts y under each component
# (n x K; NA where y or the parameters are missing, and where a component
# is undefined at the row's parameters).
.fcr_fm_row_logpmf <- function(r, y) {
  LP <- matrix(NA_real_, r$n, r$ps$K)
  for (k in seq_len(r$ps$K)) {
    cv <- r$cvs[[k]]
    ok <- !is.na(y) & !is.na(cv$mu)
    if (!any(ok)) next
    if (!is.null(cv$mu_draws)) {
      LP[ok, k] <- .fcr_lmean_rows(.fcr_fm_logpmf_draws(
        r$ps, cv, r$ps$families[k], y, which(ok)), cv$draw_w)
      next
    }
    sv <- list(mu = cv$mu[ok], alpha = if (!is.null(cv$alpha)) cv$alpha[ok],
               sigma = if (!is.null(cv$sigma)) cv$sigma[ok])
    LP[ok, k] <- .fcr_fm_logpmf(r$ps, sv, r$ps$families[k], y[ok])
  }
  LP
}

# Rows that carry no outcome information: a missing outcome, or a missing
# covariate value of some component. These rows are left out of the class
# updates. They are told apart from rows whose component is undefined at
# valid inputs (NA in .fcr_fm_row_logpmf() too), which make their panel's
# probabilities NA instead of being dropped quietly.
.fcr_fm_row_missing <- function(r, y) {
  miss <- is.na(y)
  for (cv in r$cvs) {
    miss <- miss | is.na(cv$mu)
    if (!is.null(cv$alpha)) miss <- miss | is.na(cv$alpha)
    if (!is.null(cv$sigma)) miss <- miss | is.na(cv$sigma)
  }
  miss
}

# Panel identifiers of the rows of `data` (NULL for cross-sectional fits).
.fcr_fm_group_ids <- function(object, data) {
  gv <- object$model$group_var
  if (is.null(gv)) return(NULL)
  miss <- setdiff(gv, names(data))
  if (length(miss)) {
    stop("The data need the panel identifier column(s) ", toString(miss),
         " for this target.", call. = FALSE)
  }
  as.character(interaction(data[, gv, drop = FALSE], drop = TRUE,
                           sep = ":"))
}

# Log class probabilities (n x K) of the rows of `data` under a target:
#   prior     - membership probabilities from the covariates only;
#   posterior - updated with the outcomes of `data` itself (retrospective);
#   history   - updated with the outcomes in `history` (earlier rows of the
#               same panels), never with the rows being predicted.
# For panel models the membership probability belongs to the panel, so it
# is computed once per panel from the validated panel membership design
# (.fcr_fm_panel_design()), never from whichever row happens to come first:
# the result does not depend on the row order of `data` or `history`. The
# prior target of a panel model also works without the panel identifiers
# (each row is then treated as its own unit); with them, the membership
# variables are checked as in fitting.
#
# With random parameters, the posterior target evaluates the panels of
# training units (matched by the panel identifiers) with their estimation
# draws, so the training data give the stored posterior probabilities;
# other panels use the shared prediction draws. Without panels the rows
# cannot be matched to training units unless `train = TRUE` says that
# `data` is the stored training data (`object$data`).
.fcr_fm_log_class <- function(object, r, data, conditioning,
                              history = NULL, train = FALSE) {
  K <- r$ps$K
  if (K == 1L) return(r$logpi)
  gv <- object$model$group_var
  if (conditioning == "history" && is.null(gv)) {
    stop("`conditioning = \"history\"` needs a panel model (`group_var`): ",
         "without panels, earlier outcomes carry no information about the ",
         "class of a new row.", call. = FALSE)
  }
  if (conditioning == "prior" &&
      (is.null(gv) || !all(gv %in% names(data)))) {
    return(r$logpi)
  }
  update <- function(logpi_u, A_u) {
    lc <- logpi_u + A_u
    lc - .fcr_lse_rows(lc)
  }
  y <- NULL
  if (conditioning == "posterior") {
    y <- .fcr_response(object, data)
    if (is.null(y)) {
      stop("`conditioning = \"posterior\"` needs the outcome in the data.",
           call. = FALSE)
    }
    if (is.null(gv)) {
      rr <- if (train && !is.null(r$ps$rp)) {
        used <- object$model$rows_used %||% seq_len(r$n)
        .fcr_fm_rows(object, data, train_unit = match(seq_len(r$n), used))
      } else r
      return(update(r$logpi, .fcr_fm_row_logpmf(rr, y)))
    }
    miss_post <- .fcr_fm_row_missing(r, y)
  }
  g <- .fcr_fm_group_ids(object, data)
  if (conditioning != "history") {
    ug <- unique(g)
    gi <- match(g, ug)
    G <- length(ug)
    Zu <- .fcr_fm_panel_design(r$des$Z, gi, G, " (in the new data)")
    logpi_u <- .fcr_fm_logprior(r$ps, r$pars, Zu, G)
    if (conditioning == "prior") return(logpi_u[gi, , drop = FALSE])
    A <- .fcr_fm_post_A(object, r, data, y, g, gi, !miss_post)
    return(update(logpi_u, A)[gi, , drop = FALSE])
  }
  # history
  if (is.null(history)) {
    stop("`conditioning = \"history\"` needs `history`: the earlier rows ",
         "(with outcomes) of the panels being predicted.", call. = FALSE)
  }
  history <- as.data.frame(history)
  rh <- .fcr_fm_rows(object, history)
  yh <- .fcr_response(object, history)
  if (is.null(yh)) stop("`history` must contain the outcome.", call. = FALSE)
  gh <- .fcr_fm_group_ids(object, history)
  ug <- unique(c(g, gh))
  gi <- match(g, ug)
  hi <- match(gh, ug)
  G <- length(ug)
  # One membership design per panel: consistent within the target rows,
  # within the history rows, and between the two
  .fcr_fm_panel_design(r$des$Z, gi, G, " (in the new data)")
  .fcr_fm_panel_design(rh$des$Z, hi, G, " (in `history`)")
  Zu <- .fcr_fm_panel_design(rbind(r$des$Z, rh$des$Z), c(gi, hi), G,
                             " (between the new data and `history`)")
  logpi_u <- .fcr_fm_logprior(r$ps, r$pars, Zu, G)
  # Missing history: rows with a missing outcome or covariate are dropped
  # (with a warning); a panel with no usable history row keeps its prior,
  # like a panel absent from `history`.
  miss <- .fcr_fm_row_missing(rh, yh)
  if (any(miss)) {
    warning(sum(miss), " history row(s) with missing values are not ",
            "used.", call. = FALSE)
  }
  A <- matrix(0, G, K)
  if (any(!miss)) {
    As <- .fcr_fm_unit_A(rh, yh, hi, !miss)
    A[as.integer(rownames(As)), ] <- As
  }
  with_hist <- seq_len(G) %in% hi[!miss]
  post <- logpi_u
  post[with_hist, ] <- update(logpi_u[with_hist, , drop = FALSE],
                              A[with_hist, , drop = FALSE])
  out <- post[gi, , drop = FALSE]
  attr(out, "n_history_panels") <- length(unique(gi[with_hist[gi]]))
  out
}

# Posterior log-likelihoods of the panels (G x K) for the posterior target
# of a panel model. With random parameters, rows of training panels (by
# their identifiers) use the estimation draws of their units, the others
# the shared prediction draws.
.fcr_fm_post_A <- function(object, r, data, y, g, gi, use) {
  if (is.null(r$ps$rp)) return(.fcr_fm_unit_A(r, y, gi, use))
  tu <- match(g, object$model$unit_ids)
  tr <- !is.na(tu)
  A <- .fcr_fm_unit_A(r, y, gi, use & !tr)
  if (any(tr)) {
    rows <- which(tr)
    rt <- .fcr_fm_rows(object, data[rows, , drop = FALSE],
                       train_unit = tu[rows])
    At <- .fcr_fm_unit_A(rt, y[rows], gi[rows], use[rows])
    A[match(rownames(At), rownames(A)), ] <- At
  }
  A
}

# Component means (n x K), using each family's actual expectation. With
# random parameters, the mean over them: exact (from the moment generating
# functions) when the family's mean is proportional to mu, otherwise
# averaged over the draws; Inf, with a warning, where it does not exist.
.fcr_fm_comp_means <- function(r) {
  n_inf <- 0L
  M <- vapply(seq_len(r$ps$K), function(k) {
    cv <- r$cvs[[k]]
    if (!is.null(cv$mu_draws)) {
      fam <- r$ps$families[k]
      mom <- cv$rp_moments
      m <- if (.fcr_mean_is_linear(fam, cv$alpha)) {
        suppressWarnings(.fcr_family_mean(fam, cv$mu_fixed * exp(mom$logf),
                                          cv$alpha, cv$sigma))
      } else {
        .fcr_wmean_rows(.fcr_fm_draw_moment(r, k, "mean"), cv$draw_w)
      }
      m <- rep_len(m, r$n)
      no <- !is.na(mom$mean_ok) & !mom$mean_ok
      m[no] <- Inf
      n_inf <<- n_inf + sum(no)
      return(m)
    }
    suppressWarnings(.fcr_family_mean(r$ps$families[k], cv$mu, cv$alpha,
                                      cv$sigma))
  }, numeric(r$n))
  if (n_inf > 0L) .fcr_fm_rp_inf_warning(n_inf, "expected count")
  matrix(M, nrow = r$n, dimnames = list(NULL, r$ps$labels))
}

# Warning for component moments that do not exist because of the random
# parameters (as in countreg_rp() predictions).
.fcr_fm_rp_inf_warning <- function(n, what) {
  warning(n, " component ", what, "(s) are infinite: the ", what, " does ",
          "not exist for a lognormal random parameter with a positive x, ",
          "or a gamma random parameter with x * sd^2 / mean >= 1",
          if (what != "expected count") " (for the variance: at 2 * x)",
          ". See ?countreg_rp.", call. = FALSE)
}

.fcr_fm_comp_vars <- function(r) {
  n_inf <- 0L
  V <- vapply(seq_len(r$ps$K), function(k) {
    cv <- r$cvs[[k]]
    if (!is.null(cv$mu_draws)) {
      # Within-draw variance plus the variance of the means over draws;
      # infinite where the second moment over the random parameters does
      # not exist
      Md <- .fcr_fm_draw_moment(r, k, "mean")
      Vd <- .fcr_fm_draw_moment(r, k, "variance")
      m <- .fcr_wmean_rows(Md, cv$draw_w)
      v <- .fcr_wmean_rows(Vd, cv$draw_w) +
        .fcr_wmean_rows(Md^2, cv$draw_w) - m^2
      no <- !is.na(cv$rp_moments$var_ok) & !cv$rp_moments$var_ok
      v[no] <- Inf
      n_inf <<- n_inf + sum(no)
      return(v)
    }
    v <- suppressWarnings(.fcr_family_variance(r$ps$families[k], cv$mu,
                                               cv$alpha, cv$sigma))
    if (is.null(v)) {
      stop("A variance function is not available for family \"",
           r$ps$families[k], "\".", call. = FALSE)
    }
    rep_len(v, r$n)
  }, numeric(r$n))
  if (n_inf > 0L) .fcr_fm_rp_inf_warning(n_inf, "variance")
  matrix(V, nrow = r$n, dimnames = list(NULL, r$ps$labels))
}

# The family mean or variance of component k at each draw of the random
# parameters (n x R).
.fcr_fm_draw_moment <- function(r, k, what) {
  cv <- r$cvs[[k]]
  R <- ncol(cv$mu_draws)
  fun <- if (what == "mean") .fcr_family_mean else .fcr_family_variance
  v <- suppressWarnings(fun(r$ps$families[k], as.vector(cv$mu_draws),
                            if (!is.null(cv$alpha)) rep(cv$alpha, R),
                            if (!is.null(cv$sigma)) rep(cv$sigma, R)))
  if (is.null(v)) {
    stop("A variance function is not available for family \"",
         r$ps$families[k], "\".", call. = FALSE)
  }
  matrix(rep_len(v, r$n * R), r$n, R)
}

# Component k's log tail probability at q, averaged over the draws of the
# random parameters when there are any.
.fcr_fm_comp_logcdf <- function(r, k, q, lower.tail) {
  cv <- r$cvs[[k]]
  if (is.null(cv$mu_draws)) {
    return(.fcr_fm_logcdf(r$ps, cv, r$ps$families[k], q, lower.tail))
  }
  R <- ncol(cv$mu_draws)
  sv <- list(mu = as.vector(cv$mu_draws),
             alpha = if (!is.null(cv$alpha)) rep(cv$alpha, R),
             sigma = if (!is.null(cv$sigma)) rep(cv$sigma, R))
  L <- .fcr_fm_logcdf(r$ps, sv, r$ps$families[k], rep(q, R), lower.tail)
  .fcr_lmean_rows(matrix(L, r$n, R), cv$draw_w)
}

# Sum over components of P * M, with 0 * Inf = 0 for zero-probability
# components.
.fcr_fm_wsum <- function(P, M) {
  S <- ifelse(P > 0, P * M, 0)
  S[is.na(P) | is.na(M)] <- NA
  as.vector(rowSums(S))
}

#' Predictions for finite-mixture and latent-class models
#'
#' @param object a model from [countreg_mix()] or [countreg_lc()].
#' @param newdata,new_data optional data frame (give one or the other). When
#'   omitted, the data used to fit the model (rows removed for missing
#'   values give `NA`).
#' @param type the prediction:
#'   \describe{
#'   \item{`"response"`}{the expected count of the mixture,
#'     \eqn{\sum_k \pi_k m_k}, where \eqn{m_k} is the actual expected count
#'     of component \eqn{k} (for PLN and PLL it includes the lognormal
#'     factor).}
#'   \item{`"component"`}{a matrix of the component expected counts
#'     \eqn{m_k}, one column per component.}
#'   \item{`"variance"`}{the variance of the mixture,
#'     \eqn{\sum_k \pi_k (V_k + m_k^2) - m^2}, including the variation
#'     between components (`Inf` when a component with positive probability
#'     has infinite variance).}
#'   \item{`"prob"`}{a matrix of count probabilities
#'     \eqn{P(Y = y) = \sum_k \pi_k f_k(y)} for \eqn{y = 0, \dots,}
#'     `max_count` (passed in `...`; default the largest training count);
#'     `overflow = TRUE` adds a column with \eqn{P(Y >} `max_count`\eqn{)}.}
#'   \item{`"link"`}{a matrix of the component linear predictors
#'     \eqn{x'\beta_k} plus offset. There is no single linear predictor of
#'     the mixture. With random parameters (`rpar_formula`), it is the
#'     linear predictor at the class's coefficients of the random
#'     variables: the means of normal, uniform, triangular and gamma random
#'     parameters, but for lognormal ones the location parameter
#'     \eqn{\mu} (`Mean:`), whose coefficient is \eqn{e^{\mu + \sigma Z}}:
#'     the column then is not on the scale of the coefficient at all. It
#'     is not the log of the component's expected count either.}
#'   }
#' @param conditioning which class probabilities \eqn{\pi_k} are used:
#'   \describe{
#'   \item{`"prior"` (default)}{membership probabilities from the
#'     covariates only. The outcomes in `newdata` are never used, so these
#'     are predictions for new units and the basis of fitted values,
#'     residuals and diagnostics.}
#'   \item{`"history"`}{for panel models: the prior updated once with the
#'     likelihood of each panel's earlier outcomes in `history`,
#'     \eqn{\pi^H_{gk} \propto \pi_{gk} \prod_{i \in H_g} f_k(y_i)}. Rows of
#'     panels without history use the prior. The outcomes of `newdata` are
#'     not used.}
#'   \item{`"posterior"`}{probabilities conditional on the outcomes in
#'     `newdata` itself (whole panels for panel models). This is a
#'     retrospective summary, not an out-of-sample prediction: the outcomes
#'     being "predicted" were used.}
#'   }
#'   With random parameters within the classes, `"history"` and
#'   `"posterior"` update only the class probabilities: within each class
#'   the random parameters keep their population distribution (they are
#'   not updated to the unit's own coefficients), so the component means
#'   \eqn{m_k} are the same for every conditioning. The class likelihoods
#'   of the update integrate over the random parameters by simulation.
#'   For panels of the training data (identified by `group_var`) the
#'   `"posterior"` update uses each panel's estimation draws, so the
#'   training data reproduce the fitted posterior probabilities
#'   ([component_probs()]); other panels, `history` rows, and the rows of a
#'   cross-sectional model passed as `newdata` use one common set of
#'   `max(rpar_draws, 500)` Halton draws, so their probabilities differ
#'   from the fitted ones by simulation error. Without `newdata`, the
#'   training rows use their estimation draws.
#' @param history for `conditioning = "history"`: a data frame of earlier
#'   rows of the panels in `newdata`, with outcomes and the panel
#'   identifier. History rows with a missing outcome or covariate are
#'   dropped with a warning; a panel without usable history rows keeps its
#'   prior. For panel models the membership variables must be constant
#'   within each panel (in `newdata`, in `history` and between the two), as
#'   in fitting, so the class probabilities never depend on the row order.
#' @param interval only `"none"`: confidence, prediction and tolerance
#'   intervals are not yet available for mixture models.
#' @param level unused (for compatibility).
#' @param ... `max_count` and `overflow` for `type = "prob"`.
#'
#' @returns A numeric vector (`"response"`, `"variance"`) with one value per
#'   row of the data, or a matrix (`"component"`, `"prob"`, `"link"`).
#' @seealso [component_probs()] for the class probabilities themselves.
#' @examples
#' set.seed(1)
#' d <- data.frame(y = c(rpois(300, 1), rpois(200, 6)))
#' m <- countreg_mix(y ~ 1, data = d, components = 2, n_starts = 3,
#'                   seed = 1)
#' predict(m, newdata = d[1:3, , drop = FALSE])
#' predict(m, newdata = d[1:3, , drop = FALSE], type = "component")
#' predict(m, newdata = d[1:3, , drop = FALSE], type = "prob",
#'         max_count = 5, overflow = TRUE)
#' @export
predict.countreg_mix <- function(object, newdata = NULL, ...,
                                 type = c("response", "component",
                                          "variance", "prob", "link"),
                                 conditioning = c("prior", "history",
                                                  "posterior"),
                                 history = NULL, new_data = NULL,
                                 interval = c("none", "confidence",
                                              "prediction"),
                                 level = 0.95) {
  type <- match.arg(type)
  conditioning <- match.arg(conditioning)
  interval <- match.arg(interval)
  if (interval != "none") {
    stop("Confidence and prediction intervals are not yet available for ",
         "finite-mixture models: the interval methods of countreg() assume ",
         "a one-component model, and mixture intervals need their own ",
         "validation.", call. = FALSE)
  }
  dots <- list(...)
  if (!is.null(dots$method)) {
    stop("`method` does not apply to mixture models; use `conditioning`.",
         call. = FALSE)
  }
  if (is.null(newdata)) newdata <- new_data
  train <- is.null(newdata)
  data <- if (!train) as.data.frame(newdata) else {
    if (is.null(object$data)) {
      stop("The model was fit with keep_data = FALSE; pass `newdata`.",
           call. = FALSE)
    }
    as.data.frame(object$data)
  }
  r <- .fcr_fm_rows(object, data)
  if (type == "link") {
    L <- vapply(seq_len(r$ps$K), function(k) log(r$cvs[[k]]$mu),
                numeric(r$n))
    return(matrix(L, nrow = r$n, dimnames = list(NULL, r$ps$labels)))
  }
  M <- .fcr_fm_comp_means(r)
  if (type == "component") return(M)
  LC <- .fcr_fm_log_class(object, r, data, conditioning, history,
                          train = train)
  P <- exp(LC)
  if (type == "response") return(.fcr_fm_wsum(P, M))
  if (type == "variance") {
    V <- .fcr_fm_comp_vars(r)
    m <- .fcr_fm_wsum(P, M)
    out <- .fcr_fm_wsum(P, V + M^2) - m^2
    out[is.infinite(.fcr_fm_wsum(P, V))] <- Inf
    return(out)
  }
  # type == "prob"
  max_count <- dots$max_count
  if (is.null(max_count)) max_count <- max(object$model$y %||% 20)
  if (!is.numeric(max_count) || length(max_count) != 1L ||
      !is.finite(max_count) || max_count < 0 ||
      max_count != round(max_count)) {
    stop("`max_count` must be a single non-negative whole number.",
         call. = FALSE)
  }
  Pr <- .fcr_fm_prob_rows(r, P, max_count)
  if (isTRUE(dots$overflow)) {
    # P(Y > max_count) from the components' own upper tails, not as
    # 1 - sum(P): the subtraction loses all relative accuracy once the tail
    # is below the rounding error of 1 (Poisson(1): P(Y > 20) = 7.5e-21
    # would become 0).
    Pr <- cbind(Pr, exp(.fcr_fm_log_tail(r, LC, rep(max_count, r$n),
                                         lower.tail = FALSE)))
    colnames(Pr)[ncol(Pr)] <- paste0(">", max_count)
  }
  Pr
}

# log P(Y <= q) (lower.tail = TRUE) or log P(Y > q) of the mixture with log
# class probabilities LC (n x K): the log-sum-exp of log pi_k plus each
# component's direct log tail probability. Components with probability 0
# do not count, even where their tail is undefined.
.fcr_fm_log_tail <- function(r, LC, q, lower.tail) {
  L <- vapply(seq_len(r$ps$K), function(k)
    .fcr_fm_comp_logcdf(r, k, q, lower.tail), numeric(r$n))
  M <- LC + matrix(L, nrow = r$n)
  M[!is.na(LC) & LC == -Inf] <- -Inf
  .fcr_lse_rows(M)
}

# Count probabilities of the mixture for k = 0..max_count.
.fcr_fm_prob_rows <- function(r, P, max_count) {
  out <- matrix(0, r$n, max_count + 1,
                dimnames = list(NULL, as.character(0:max_count)))
  for (k in seq_len(r$ps$K)) {
    pk <- P[, k]
    for (c in 0:max_count) {
      lp <- .fcr_fm_row_logpmf_k(r, k, rep(c, r$n))
      out[, c + 1] <- out[, c + 1] + ifelse(pk > 0, pk * exp(lp), 0)
    }
  }
  out[is.na(P[, 1]), ] <- NA
  out
}

.fcr_fm_row_logpmf_k <- function(r, k, y) {
  cv <- r$cvs[[k]]
  out <- rep(NA_real_, r$n)
  ok <- !is.na(cv$mu) & !is.na(y)
  if (any(ok) && !is.null(cv$mu_draws)) {
    out[ok] <- .fcr_lmean_rows(.fcr_fm_logpmf_draws(
      r$ps, cv, r$ps$families[k], y, which(ok)), cv$draw_w)
  } else if (any(ok)) {
    sv <- list(mu = cv$mu[ok], alpha = if (!is.null(cv$alpha)) cv$alpha[ok],
               sigma = if (!is.null(cv$sigma)) cv$sigma[ok])
    out[ok] <- .fcr_fm_logpmf(r$ps, sv, r$ps$families[k], y[ok])
  }
  out
}

# Probability matrix used by rootogram()
.fcr_fm_prob_matrix <- function(object, data, max_count) {
  r <- .fcr_fm_rows(object, data)
  .fcr_fm_prob_rows(r, exp(r$logpi), max_count)
}

# Variance of the mixture at the prior probabilities (Pearson residuals)
.fcr_fm_variance <- function(object, data) {
  stats::predict(object, newdata = data, type = "variance")
}

#' Component (class) membership probabilities
#'
#' Prior or posterior probabilities of the mixture components (latent
#' classes). `class_probs()` is the same function under the latent-class
#' name.
#'
#' @param object a model from [countreg_mix()] or [countreg_lc()].
#' @param newdata optional data frame. Without it, the probabilities of the
#'   independent units of the training data are returned (one row per row,
#'   or per panel for panel models).
#' @param conditioning `"prior"` (default): membership probabilities from
#'   the covariates; `"posterior"`: conditional on the outcomes (each unit's
#'   own outcomes; a retrospective classification, not a prediction);
#'   `"history"`: for panel models, conditional on earlier outcomes in
#'   `history` (see [predict.countreg_mix()]).
#' @param history earlier rows of the panels, for `conditioning =
#'   "history"`.
#' @param ... passed from `class_probs()`.
#'
#' @returns A matrix with one column per component (named by the component
#'   labels) and one row per unit (training data) or per row of `newdata`.
#'   These are probabilities of the components, not of counts (for count
#'   probabilities use `predict(type = "prob")`).
#' @details The posterior probabilities give a descriptive classification
#'   (e.g. the most probable class of each unit), which is uncertain: units
#'   with probabilities near \eqn{1/K} are not well classified. The entropy
#'   in `summary()` summarizes that uncertainty.
#'
#'   With random parameters within the classes, the posterior of a unit
#'   integrates over its random parameters by simulation. Panels of the
#'   training data passed in `newdata` (matched by `group_var`) use their
#'   estimation draws and so reproduce the stored probabilities; new panels,
#'   and the rows of a cross-sectional model passed as `newdata`, use a
#'   common set of Halton draws (see [predict.countreg_mix()]).
#' @examples
#' set.seed(1)
#' d <- data.frame(y = c(rpois(300, 1), rpois(200, 6)))
#' m <- countreg_mix(y ~ 1, data = d, components = 2, n_starts = 3,
#'                   seed = 1)
#' head(component_probs(m, conditioning = "posterior"))
#' @export
component_probs <- function(object, newdata = NULL,
                            conditioning = c("prior", "posterior",
                                             "history"),
                            history = NULL) {
  if (!inherits(object, "countreg_mix")) {
    stop("component_probs() needs a countreg_mix() or countreg_lc() model.",
         call. = FALSE)
  }
  conditioning <- match.arg(conditioning)
  m <- object$model
  if (is.null(newdata) && conditioning != "history") {
    out <- if (conditioning == "prior") m$prior else m$posterior
    rownames(out) <- m$unit_ids
    attr(out, "conditioning") <- conditioning
    return(out)
  }
  if (is.null(newdata)) newdata <- .fcr_training_rows(object)
  data <- as.data.frame(newdata)
  r <- .fcr_fm_rows(object, data)
  out <- exp(.fcr_fm_log_class(object, r, data, conditioning, history))
  colnames(out) <- m$labels
  attr(out, "conditioning") <- conditioning
  out
}

#' @rdname component_probs
#' @export
class_probs <- function(object, ...) component_probs(object, ...)

#' @export
fitted.countreg_mix <- function(object, ...) {
  stats::predict(object, newdata = .fcr_training_rows(object),
                 type = "response")
}

#' @export
residuals.countreg_mix <- function(object,
                                   type = c("response", "pearson",
                                            "quantile"), ...) {
  type <- match.arg(type)
  d <- .fcr_training_rows(object)
  y <- as.numeric(object$model$y)
  if (type == "quantile") return(.fcr_fm_quantile_residuals(object, d, y))
  mu <- stats::predict(object, newdata = d, type = "response")
  r <- y - mu
  if (type == "pearson") {
    v <- .fcr_fm_variance(object, d)
    if (any(!is.finite(v) | v <= 0)) {
      stop("Pearson residuals are not available: the mixture variance is ",
           "infinite, undefined or zero for ", sum(!is.finite(v) | v <= 0),
           " row(s) (e.g. a Generalized Waring component with rho <= 2, or ",
           "a random parameter whose second moment does not exist). ",
           "Use quantile residuals.", call. = FALSE)
    }
    r <- r / sqrt(v)
  }
  r
}

# Randomized quantile residuals from the mixture distribution function
# F(y) = sum_k pi_k F_k(y), computed entirely on the log scale. With
# lp = log P(Y = y), U = F(y - 1) + u P(Y = y) and 1 - U = S(y) +
# (1 - u) P(Y = y), where S(y) = P(Y > y) comes from the components' own
# upper tails. Both are kept as logarithms, and the normal quantile is
# taken from the smaller of the two tails: an outcome far in the upper tail
# then gives a large finite residual instead of qnorm(1) = Inf.
.fcr_fm_quantile_residuals <- function(object, d, y, u = NULL) {
  r <- .fcr_fm_rows(object, d)
  n <- length(y)
  if (is.null(u)) u <- stats::runif(n)
  LC <- r$logpi
  log_lo <- .fcr_fm_log_tail(r, LC, y - 1, lower.tail = TRUE)
  log_s <- .fcr_fm_log_tail(r, LC, y, lower.tail = FALSE)
  M <- LC + .fcr_fm_row_logpmf(r, y)
  M[!is.na(LC) & LC == -Inf] <- -Inf
  lp <- .fcr_lse_rows(M)
  log_u <- pmin(.fcr_logadd(log_lo, log(u) + lp), 0)
  log_1mu <- pmin(.fcr_logadd(log_s, log1p(-u) + lp), 0)
  # F(y - 1) and P(Y > y) come from the same distribution function, so
  # F(y - 1) + P(Y > y) <= 1 (the old lower <= upper check)
  bad <- is.na(log_u) | is.na(log_1mu) |
    exp(log_lo) + exp(log_s) > 1 + 1e-12
  if (any(bad)) {
    warning(sum(bad), " row(s) have invalid fitted probabilities; their ",
            "quantile residuals are NA.", call. = FALSE)
  }
  low <- !bad & log_u <= log(0.5)
  out <- rep(NA_real_, n)
  out[low] <- stats::qnorm(log_u[low], log.p = TRUE)
  up <- !bad & !low
  out[up] <- stats::qnorm(log_1mu[up], lower.tail = FALSE, log.p = TRUE)
  out
}

#' @export
simulate.countreg_mix <- function(object, nsim = 1, seed = NULL, ...) {
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv()))
      get(".Random.seed", envir = globalenv())
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
      assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
  }
  rng_state <- if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv())
  m <- object$model
  d <- .fcr_training_rows(object)
  r <- .fcr_fm_rows(object, d)
  n <- r$n
  K <- r$ps$K
  unit <- m$unit %||% seq_len(n)
  prior_u <- m$prior
  n_u <- nrow(prior_u)
  classes <- matrix(NA_integer_, n_u, nsim)
  rp <- r$ps$rp
  sims <- lapply(seq_len(nsim), function(s) {
    cu <- if (K == 1L) rep(1L, n_u) else
      apply(prior_u, 1, function(p) sample.int(K, 1L, prob = p))
    classes[, s] <<- cu
    cr <- cu[unit]
    y <- numeric(n)
    # Random parameters: one new draw of each unit's coefficients
    hd <- if (!is.null(rp)) {
      Uu <- matrix(stats::runif(n_u * length(rp$names)), n_u)
      structure(lapply(seq_len(ncol(Uu)), function(j)
        matrix(Uu[unit, j], n, 1)), class = "fcr_unit_draws")
    }
    for (k in seq_len(K)) {
      idx <- which(cr == k)
      if (!length(idx)) next
      cv <- r$cvs[[k]]
      mu <- cv$mu
      if (!is.null(rp)) {
        cols <- match(rp$names, colnames(r$des$X))
        xb <- generate_random_draws(
          .fcr_draw_rows(hd, idx), r$pars$beta[[k]][cols], r$pars$rv[[k]],
          rp$dists, rp$names, r$des$X[idx, cols, drop = FALSE],
          correlated = rp$correlated)$xb_rand_mat[, 1]
        mu <- numeric(n)
        mu[idx] <- cv$mu_fixed[idx] * exp(xb)
      }
      y[idx] <- .fcr_family_r(r$ps$families[k], length(idx), mu[idx],
                              if (!is.null(cv$alpha)) cv$alpha[idx],
                              if (!is.null(cv$sigma)) cv$sigma[idx])
    }
    y
  })
  out <- as.data.frame(sims)
  names(out) <- paste0("sim_", seq_len(nsim))
  attr(out, "seed") <- rng_state
  attr(out, "classes") <- classes
  out
}

#' @export
print.countreg_mix <- function(x, digits = max(3L, getOption("digits") - 3L),
                               ...) {
  m <- x$model
  what <- if (inherits(x, "countreg_lc")) "Latent-class" else
    "Finite-mixture"
  cat(what, "count model:", m$K,
      if (inherits(x, "countreg_lc")) "class(es)" else "component(s)", "\n")
  cat("Families:", paste0(m$labels, " = ", m$family, collapse = ", "), "\n")
  if (!is.null(x$call)) {
    cat("Call: ")
    print(x$call)
  }
  cat("\nCoefficients:\n")
  print.default(format(stats::coef(x), digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\nLog-likelihood:", format(m$maximum, digits = digits + 3),
      "| N =", m$nobs)
  if (!is.null(m$n_groups)) cat(" | panels =", m$n_groups)
  cat("\nMixing proportions (prior means):",
      paste0(m$labels, " ", format(round(m$classes$prior_share, 3)),
             collapse = ", "), "\n")
  if (isFALSE(m$converged)) cat("Warning: the optimizer did not converge.\n")
  invisible(x)
}

#' Summaries of finite-mixture and latent-class models
#'
#' @param object a model from [countreg_mix()] or [countreg_lc()].
#' @param x a `summary.countreg_mix` object.
#' @param level confidence level of the coefficient intervals.
#' @param digits number of digits to print.
#' @param ... further arguments for [summary.flexCountReg()].
#' @details Besides the coefficient table (grouped by component, shared
#'   parameters and mixing model), the summary reports the prior mixing
#'   proportions, the posterior (effective) class sizes, the entropy of the
#'   posterior classification, the starts that reached the selected
#'   solution, and numerical checks of the observed information. No
#'   likelihood-ratio test of the number of components is reported: its
#'   usual chi-square reference distribution does not apply.
#'
#'   Entropy is \eqn{H = -\sum_u w_u \sum_k \tau_{uk}\log\tau_{uk}} over the
#'   independent units with posterior probabilities \eqn{\tau}; the relative
#'   entropy \eqn{1 - H / (N \log K)} is 1 for a perfect classification and
#'   0 when every unit has probabilities \eqn{1/K}. ICL is
#'   \eqn{BIC + 2H} (smaller is better), a classification-oriented
#'   criterion.
#' @returns An object of class `c("summary.countreg_mix",
#'   "summary.flexCountReg")`.
#' @export
summary.countreg_mix <- function(object, level = 0.95, ...) {
  s <- NextMethod()
  m <- object$model
  s$K <- m$K
  s$labels <- m$labels
  s$families <- m$family
  s$family_status <- m$family_status
  s$classes <- m$classes
  s$entropy <- m$entropy
  s$relative_entropy <- m$relative_entropy
  s$ICL <- m$ICL
  s$starts <- m$starts
  s$numerical <- m$numerical
  s$mix_mode <- m$mix_mode
  s$fixed_mixing <- m$fixed_mixing
  s$info <- m$map$info
  # Smooth terms are summarized by term (as in summary.flexCountReg())
  if (length(m$smooth_coefs)) s$info <- s$info[-m$smooth_coefs, , drop = FALSE]
  s$n_units <- m$n_units
  s$latent_class <- inherits(object, "countreg_lc")
  s$n_boot_ambiguous <- m$n_boot_ambiguous
  s$refined_from <- m$refined_from
  s$identical_pairs <- m$identical_pairs
  class(s) <- c("summary.countreg_mix", class(s))
  s
}

#' @rdname summary.countreg_mix
#' @export
print.summary.countreg_mix <- function(x, digits = x$digits, ...) {
  unit_word <- if (x$latent_class) "Class" else "Component"
  cat("Call:\n")
  if (!is.null(x$call)) print(x$call) else print(x$formula)
  word <- if (x$latent_class) {
    if (x$K > 1) "classes" else "class"
  } else {
    if (x$K > 1) "components" else "component"
  }
  cat("\n", if (x$latent_class) "Latent-class" else "Finite-mixture",
      " count model with ", x$K, " ", word, "\n", sep = "")
  fam <- data.frame(label = x$labels, family = x$families,
                    status = x$family_status)
  for (i in seq_len(nrow(fam))) {
    lab <- if (fam$label[i] == paste0(unit_word, i)) "" else
      paste0(" (", fam$label[i], ")")
    cat("  ", unit_word, " ", i, lab, ": ", fam$family[i],
        if (fam$status[i] != "validated") " [experimental]", "\n", sep = "")
  }
  mix <- switch(x$mix_mode,
                none = "none (one component)",
                constant = "constant proportions",
                covariate = "multinomial logit on covariates",
                fixed = paste0("fixed at ",
                               toString(format(x$fixed_mixing))))
  cat("Mixing:", mix)
  if (x$mix_mode %in% c("constant", "covariate")) {
    cat(" (reference: ", unit_word, " ", x$K, ")", sep = "")
  }
  cat("\nOptimizer:", x$method)
  if (!is.null(x$refined_from)) cat(" + Newton-Raphson refinement")
  cat(" |", trimws(x$message %||% ""), "\n")
  if (length(x$identical_pairs)) {
    cat("Warning: nearly identical components (",
        paste(x$identical_pairs, collapse = ", "), "): their shares and ",
        "means are arbitrary.\n", sep = "")
  }
  if (isFALSE(x$converged)) {
    cat("Warning: the optimizer did not report convergence.\n")
  }

  se_label <- switch(x$vcov_type %||% "hessian",
                     hessian = "observed information (Hessian)",
                     bayesian = "Bayesian (penalized Hessian)",
                     robust = "robust (sandwich, unit-level scores)",
                     bootstrap = paste0("bootstrapped (", x$n_boot_ok,
                                        " of ", x$bootstraps,
                                        " replicates)"),
                     x$vcov_type)
  cat("\nCoefficients (standard errors: ", se_label, "):\n", sep = "")
  tab <- x$coefficients[, 1:4, drop = FALSE]
  info <- x$info
  groups <- ifelse(info$block == "mixing",
                   "Mixing (multinomial logit)",
                   ifelse(is.na(info$component), "Shared parameters",
                          paste0(unit_word, " ", info$component, " (",
                                 ifelse(x$labels[pmax(1, info$component)] ==
                                          paste0(unit_word, info$component),
                                        "", paste0(x$labels[pmax(
                                          1, info$component)], ", ")),
                                 info$family, ")")))
  for (g in unique(groups)) {
    cat("\n", g, ":\n", sep = "")
    stats::printCoefmat(tab[groups == g, , drop = FALSE], digits = digits,
                        P.values = TRUE, has.Pvalue = TRUE, na.print = "NA")
  }
  if (length(x$offset)) {
    cat("Offset:", paste(x$offset, collapse = " + "),
        "(coefficient fixed at 1, all components)\n")
  }
  .fcr_print_smooth_table(x)

  cl <- x$classes
  cat("\n", unit_word, " sizes:\n", sep = "")
  show <- data.frame(
    label = cl$label, family = cl$family,
    prior = round(cl$prior_share, 4), posterior = round(cl$posterior_share, 4),
    effective_n = round(cl$effective_size, 1),
    modal_n = round(cl$modal_count, 1),
    mean_count = signif(cl$mean_count, 4))
  if (any(cl$small)) show$note <- ifelse(cl$small, "small", "")
  print(show, row.names = FALSE)
  cat("(prior: mean membership probability; posterior: mean posterior ",
      "probability;\n effective_n: sum of posterior probabilities; ",
      "modal_n: units whose most\n probable ", tolower(unit_word),
      " it is; mean_count: average expected count)\n", sep = "")

  cat("\nLog-likelihood:", format(round(x$logLik, digits)),
      "| AIC:", format(round(x$AIC, digits)),
      "| BIC:", format(round(x$BIC, digits)),
      "| ICL:", format(round(x$ICL, digits)), "\n")
  cat("N:", x$nobs)
  if (!is.null(x$n_groups)) cat(" | panels (independent units):", x$n_groups)
  if (!is.null(x$n_dropped) && x$n_dropped > 0) {
    cat(" | rows removed (missing values):", x$n_dropped)
  }
  cat("\n")
  if (x$K > 1) {
    cat("Entropy:", format(round(x$entropy, digits)),
        "| relative entropy:", format(round(x$relative_entropy, 3)), "\n")
  }
  st <- x$starts
  if (nrow(st) > 1) {
    cat("Starts:", nrow(st), "|", sum(st$reached_best),
        "reached the selected log-likelihood\n")
    if (sum(st$reached_best) < 2) {
      cat("Note: only one start reached this solution; try more starts.\n")
    }
  }
  num <- x$numerical
  cat("Observed information:", num$information)
  if (is.finite(num$condition)) {
    cat(" (condition number", format(signif(num$condition, 3)), ")")
  }
  cat(" | max |gradient|:", format(signif(num$max_abs_gradient, 3)), "\n")
  if (!is.null(x$n_boot_ambiguous) && x$n_boot_ambiguous > 0) {
    cat("Bootstrap refits with ambiguous label matching:",
        x$n_boot_ambiguous, "\n")
  }
  if (x$K > 1) {
    cat("No likelihood-ratio test of the number of ",
        tolower(unit_word), if (x$latent_class) "es" else "s",
        " is reported (non-standard distribution).\n", sep = "")
  }
  invisible(x)
}

#' @export
tidy.countreg_mix <- function(x, conf.int = FALSE, conf.level = 0.95,
                              exponentiate = FALSE, ...) {
  cf <- stats::coef(x)
  se <- .fcr_se(stats::vcov(x), length(cf))
  info <- x$model$map$info
  labels <- x$model$labels
  comp <- ifelse(info$block == "mean", "mean",
                 ifelse(info$block == "mixing", "mixing",
                        ifelse(info$block == "rpar", "random",
                               "distribution")))
  cls <- ifelse(is.na(info$component), "Shared",
                labels[pmax(1L, info$component)])
  out <- tibble::tibble(term = names(cf), component = comp, class = cls,
                        family = info$family, estimate = unname(cf),
                        std.error = unname(se))
  out$statistic <- out$estimate / out$std.error
  out$p.value <- 2 * stats::pnorm(-abs(out$statistic))
  if (conf.int) {
    z <- stats::qnorm((1 + conf.level) / 2)
    out$conf.low <- out$estimate - z * out$std.error
    out$conf.high <- out$estimate + z * out$std.error
  }
  if (exponentiate) {
    # Rate ratios for the mean coefficients only (mixing coefficients are
    # log odds ratios against the reference component)
    i <- out$component == "mean"
    for (col in intersect(c("estimate", "conf.low", "conf.high"),
                          names(out))) {
      out[[col]][i] <- exp(out[[col]][i])
    }
    out$std.error[i] <- out$estimate[i] * out$std.error[i]
  }
  out
}

#' @export
glance.countreg_mix <- function(x, ...) {
  m <- x$model
  tibble::tibble(
    model = class(x)[1],
    family = paste(m$family, collapse = "+"),
    components = m$K,
    logLik = as.numeric(stats::logLik(x)),
    AIC = stats::AIC(x),
    BIC = stats::BIC(x),
    ICL = m$ICL,
    df = length(stats::coef(x)),
    nobs = m$nobs,
    n_units = m$n_units,
    entropy = m$entropy,
    relative_entropy = m$relative_entropy,
    converged = isTRUE(m$converged),
    n_starts = m$n_starts,
    starts_reached_best = sum(m$starts$reached_best),
    information = m$numerical$information,
    vcov_type = m$vcov_type)
}
