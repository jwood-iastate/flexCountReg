#' Resampling count models, including models with several outcomes
#'
#' Fits a workflow or model specification to the analysis set of each
#' resample and evaluates it on the assessment set, for every outcome.
#' \pkg{tune}'s \code{fit_resamples()} computes metrics for one outcome
#' only, so models with several (\code{\link{countreg_mv}},
#' \code{\link{countreg_alloc}}) need this function. It also computes the
#' out-of-sample log score of the count distribution, which compares
#' distributions better than errors of the mean, for any engine.
#'
#' @param object a \pkg{workflows} workflow or a \code{\link{count_reg}}
#'   specification.
#' @param resamples an \pkg{rsample} resampling set, such as
#'   \code{rsample::group_vfold_cv()}.
#' @param formula the model formula, for a model specification.
#' @param metrics a \pkg{yardstick} metric set of numeric metrics (default
#'   RMSE and MAE).
#' @param log_score logical; also compute the mean log probability of the
#'   held-out counts (\code{"log_score"}, larger is better). For models with
#'   several outcomes it is the score of each outcome's marginal
#'   distribution.
#' @param summarize logical; average over the resamples (\code{TRUE}) or
#'   return the values of each?
#' @returns A tibble with columns \code{outcome}, \code{.metric},
#'   \code{mean}, \code{n} (resamples with a value) and \code{std_err}, or
#'   with \code{summarize = FALSE} \code{id}, \code{outcome}, \code{.metric}
#'   and \code{.estimate}.
#' @details Rows without a prediction (such as new groups of a fixed-effects
#'   model) are left out of the metrics. The options \code{marginal},
#'   \code{re.form}, \code{allow.new.levels} and \code{newgroup} of the
#'   engine argument \code{predict_args} also apply to the log score.
#' @examples
#' \donttest{
#' if (requireNamespace("parsnip", quietly = TRUE) &&
#'     requireNamespace("rsample", quietly = TRUE) &&
#'     requireNamespace("yardstick", quietly = TRUE)) {
#'   spec <- count_reg(family = "NB2") |>
#'     parsnip::set_engine("flexCountReg", method = "BFGS")
#'   set.seed(1)
#'   folds <- rsample::group_vfold_cv(washington_roads, group = ID, v = 3)
#'   fit_resamples_counts(spec, folds, Total_crashes ~ lnaadt + lnlength,
#'                        log_score = TRUE)
#' }
#' }
#' @export
fit_resamples_counts <- function(object, resamples, formula = NULL,
                                 metrics = NULL, log_score = FALSE,
                                 summarize = TRUE) {
  rlang::check_installed(c("parsnip", "rsample", "yardstick", "hardhat"),
                         reason = "to resample models.")
  if (!inherits(resamples, "rset")) {
    stop("`resamples` must be an rsample resampling set.", call. = FALSE)
  }
  is_spec <- inherits(object, "model_spec")
  if (!is_spec && !inherits(object, "workflow")) {
    stop("`object` must be a workflow or a model specification.",
         call. = FALSE)
  }
  if (is_spec && is.null(formula)) {
    stop("`formula` is needed with a model specification.", call. = FALSE)
  }
  if (!is_spec && !is.null(formula)) {
    stop("`formula` is for model specifications; a workflow has its own.",
         call. = FALSE)
  }
  if (is.null(metrics)) {
    metrics <- yardstick::metric_set(yardstick::rmse, yardstick::mae)
  }
  ids <- resamples[grep("^id", names(resamples))]
  id <- do.call(paste, c(as.list(ids), sep = " "))
  one <- function(split) {
    train <- rsample::analysis(split)
    test <- rsample::assessment(split)
    fit <- if (is_spec) parsnip::fit(object, formula, data = train) else
      parsnip::fit(object, data = train)
    p <- stats::predict(fit, new_data = test)
    engine <- hardhat::extract_fit_engine(fit)
    outs <- sub("^\\.pred_?", "", names(p))
    if (identical(outs, "")) outs <- .fcr_resample_outcome(fit, engine)
    miss <- setdiff(outs, names(test))
    if (length(miss)) {
      stop("The outcome(s) ", toString(miss), " are not columns of the ",
           "data.", call. = FALSE)
    }
    ls <- if (log_score) .fcr_resample_score(fit, engine, test, outs)
    do.call(rbind, lapply(seq_along(outs), function(j) {
      d <- data.frame(truth = as.numeric(test[[outs[j]]]),
                      estimate = as.numeric(p[[j]]))
      ok <- stats::complete.cases(d)
      r <- metrics(d[ok, , drop = FALSE], truth = truth, estimate = estimate)
      out <- tibble::tibble(outcome = outs[j], .metric = r$.metric,
                            .estimate = r$.estimate)
      if (log_score) {
        out <- rbind(out, tibble::tibble(
          outcome = outs[j], .metric = "log_score", .estimate = ls[[j]]))
      }
      out
    }))
  }
  res <- do.call(rbind, Map(function(s, i) {
    r <- one(s)
    tibble::add_column(r, id = i, .before = 1L)
  }, resamples$splits, id))
  if (!summarize) return(res)
  key <- unique(res[c("outcome", ".metric")])
  stat <- lapply(seq_len(nrow(key)), function(i) {
    v <- res$.estimate[res$outcome == key$outcome[i] &
                         res$.metric == key$.metric[i]]
    v <- v[!is.na(v)]
    c(mean(v), length(v), stats::sd(v) / sqrt(length(v)))
  })
  key$mean <- vapply(stat, `[`, 0, 1L)
  key$n <- as.integer(vapply(stat, `[`, 0, 2L))
  key$std_err <- vapply(stat, `[`, 0, 3L)
  key
}

# Name of the single outcome column of a fit.
.fcr_resample_outcome <- function(fit, engine) {
  if (inherits(fit, "workflow")) {
    o <- names(hardhat::extract_mold(fit)$outcomes)
    if (length(o) == 1L) return(o)
  }
  all.vars(engine$formula[[2L]])[1L]
}

# Mean log probability of the held-out counts of each outcome.
.fcr_resample_score <- function(fit, engine, test, outs) {
  extra <- engine$parsnip$predict_args
  extra <- extra[intersect(names(extra), c("marginal", "re.form", "newgroup",
                                           "allow.new.levels"))]
  Y <- lapply(outs, function(o) as.numeric(test[[o]]))
  cs <- sort(unique(unlist(Y)))
  type <- if (inherits(engine, "countreg_alloc")) "count_prob" else "prob"
  P <- stats::predict(fit, new_data = test, type = "raw",
                      opts = c(list(type = type, counts = cs), extra))
  if (!is.list(P)) P <- list(P)
  Map(function(p, y) mean(log(p[cbind(seq_along(y), match(y, cs))]),
                          na.rm = TRUE), P, Y)
}
