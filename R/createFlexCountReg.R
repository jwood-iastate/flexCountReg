# Internal constructor for flexCountReg objects.
#
# Every estimator returns an S3 object of class c(<subclass>, "flexCountReg"),
# a list with `call`, `model` (the fitted maxLik object with the package's
# fields added), `data` (the training data, unless keep_data = FALSE) and
# `formula`. The subclass is derived from model$modelType:
#   countreg     -> "countreg"       countreg_rp -> "countreg_rp"
#   RENB         -> "renb"           poisLindRE  -> "poislind_re"
#   REPoisson    -> "repois"         NegMultinom -> "negmultinom"
#   countreg_mix -> "countreg_mix"   countreg_lc -> c("countreg_lc",
#                                                     "countreg_mix")
# The constructor also fills in the fields that the S3 methods use, so the
# methods work the same way for every estimator: a named `estimate`, a
# covariance matrix `vcov` with `vcov_type`, `nobs`, `y`, and `converged`.
.createFlexCountReg <- function(model = NULL, data = NULL, call = NULL, ...,
                                keep_data = TRUE) {
  formula <- formula(stats::as.formula(model$formula))
  model <- .fcr_normalize_model(model)
  model$data <- NULL   # the training data is stored once, in obj$data
  # maxLik keeps the log-likelihood (objectiveFn) and any gradient or Hessian
  # functions. Their closures hold the training data, design matrices and
  # draws, so with them a keep_data = FALSE fit still carried (and shared)
  # the data, and every fit was several times larger than needed. No method
  # uses them after estimation.
  model[vapply(model, is.function, logical(1))] <- NULL
  # Without the data, the LR sample evidence keeps only fingerprints of the
  # model's variables, not their values
  if (!isTRUE(keep_data)) model$lr_sample <- .fcr_lr_sample_strip(model$lr_sample)
  # Package version that created the object (for migrating saved fits)
  model$fcr_version <- tryCatch(
    as.character(utils::packageVersion("flexCountReg")),
    error = function(e) NA_character_)
  # Schema of the fitted object: the fields that prediction, comparison and
  # refitting rely on. See .fcr_schema_version.
  model$fcr_schema <- .fcr_schema_version

  subclass <- switch(model$modelType %||% "countreg",
                     countreg = "countreg",
                     countreg.rp = "countreg_rp",
                     countreg_fe = "countreg_fe",
                     countreg_hier = "countreg_hier",
                     RENB = "renb",
                     poisLindRE = "poislind_re",
                     REPoisson = "repois",
                     NegMultinom = "negmultinom",
                     countreg_mix = "countreg_mix",
                     countreg_lc = c("countreg_lc", "countreg_mix"),
                     "countreg")

  structure(list(
    call = call,
    model = model,
    data = if (isTRUE(keep_data)) data else NULL,
    formula = formula
  ), class = c(subclass, "flexCountReg"))
}

# Fill in the standard fields of a fitted model.
.fcr_normalize_model <- function(model) {
  est <- model[["estimate"]]
  if (is.list(est)) est <- unlist(est)
  if (is.null(names(est)) && !is.null(model[["x_names"]])) names(est) <- model[["x_names"]]
  model$estimate <- est
  model$coefficients <- est
  nm <- names(est)

  if (is.null(model[["vcov"]]) || !is.matrix(model[["vcov"]])) {
    V <- tryCatch(suppressWarnings(.fcr_hessian_vcov(model)),
                  error = function(e) NULL)
    if (is.null(V)) V <- matrix(NA_real_, length(est), length(est))
    model$vcov <- V
    if (is.null(model[["vcov_type"]])) model$vcov_type <- "hessian"
  }
  if (!is.null(nm) && all(dim(model[["vcov"]]) == length(nm))) {
    dimnames(model[["vcov"]]) <- list(nm, nm)
  }
  if (is.null(model[["vcov_type"]])) {
    model$vcov_type <- switch(model[["stderr"]] %||% "normal",
                              robust = "robust", boot = "bootstrap",
                              "hessian")
  }

  if (is.null(model[["y"]])) model$y <- model[["observed"]]
  if (is.null(model[["nobs"]]) && !is.null(model[["y"]])) model$nobs <- length(model[["y"]])
  if (is.null(model[["logLik"]]) || !is.numeric(model[["logLik"]])) {
    model$logLik <- model[["maximum"]] %||% model[["logLik_value"]]
  }
  if (is.null(model[["maximum"]])) model$maximum <- model[["logLik"]]
  if (is.null(model[["converged"]])) {
    ml <- model[["maxLik_object"]]
    if (!is.null(model[["code"]])) {
      model$converged <- .fcr_converged(model, model$method %||% "")
    } else if (!is.null(ml) && !is.null(ml$code)) {
      model$converged <- .fcr_converged(ml, toupper(model$method %||% ""))
    }
  }
  model
}

# Fitted-object schema.
#   1 (no `fcr_schema` field): flexCountReg 0.1.1 and earlier, and
#     development builds of 0.1.2 before the follow-up review.
#   2: adds `lr_sample` (sample evidence for likelihood-ratio tests),
#      `panel_index` and row-aligned `weights` (bootstrap refits), and
#      `engine` "quadrature" for the PLN, PGE, PW and PLL families.
# Older objects keep working where that is reliable: predictions and
# probabilities use the engine they were fit with ("halton" when none is
# stored). LR tests between fitted models and marginal-effect bootstraps of
# panel models need schema 2 and ask for a refit otherwise.
.fcr_schema_version <- 2L

# Schema of a fitted object (1 when it predates the field).
.fcr_schema <- function(object) object$model$fcr_schema %||% 1L
