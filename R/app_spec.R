# The versioned analysis specification. Every control in the app edits a
# draft specification; a fit freezes a copy of it, and the code that runs in
# the worker and the code that is exported are both compiled from it.

.spec_schema <- 1L

# Create an analysis specification
#
app_spec <- function(data = list(), recipe = list(), model = list()) {
  d <- list(source = "bundled", name = "abalone_rings")
  r <- list(response = NULL, predictors = character(), factors = character(),
            exposure = NULL, log_offset = NULL, weights = NULL,
            group = NULL)
  m <- list(estimator = "countreg", family = "NB2", components = 2L,
            dis1 = character(), mixing = character(),
            shared_mean = character(), shared_dis1 = FALSE,
            n_starts = 10L, seed = 1L, method = "BFGS", max_iters = 1000L,
            stderr = "normal")
  merge <- function(base, new) {
    unknown <- setdiff(names(new), names(base))
    if (length(unknown)) {
      stop("Unknown specification field(s): ", toString(unknown),
           call. = FALSE)
    }
    base[names(new)] <- new
    base
  }
  structure(list(schema = .spec_schema, data = merge(d, data),
                 recipe = merge(r, recipe), model = merge(m, model)),
            class = "fcr_app_spec")
}

# Problems that prevent fitting, as readable sentences (empty when the
# specification can be estimated). `data` is the source data frame.
validate_spec <- function(spec, data) {
  p <- character(0)
  add <- function(...) p <<- c(p, paste0(...))
  r <- spec$recipe
  m <- spec$model
  cols <- names(data)
  if (is.null(r$response) || !nzchar(r$response)) {
    add("Choose a response (count) column.")
  } else if (!r$response %in% cols) {
    add("The response column \"", r$response, "\" is not in the data.")
  } else {
    y <- data[[r$response]]
    yy <- y[!is.na(y)]
    if (!is.numeric(y)) {
      add("The response must be numeric counts.")
    } else if (any(yy < 0) || any(abs(yy - round(yy)) > 1e-8)) {
      add("The response must contain nonnegative whole numbers.")
    }
  }
  roles <- c(r$predictors, r$factors, r$exposure, r$log_offset, r$weights,
             r$group, m$dis1, m$mixing)
  miss <- setdiff(roles, cols)
  if (length(miss)) add("Column(s) not in the data: ", toString(miss), ".")
  if (!is.null(r$response) && r$response %in%
      c(r$predictors, m$dis1, m$mixing)) {
    add("The response cannot also be a predictor.")
  }
  if (!is.null(r$exposure) && !is.null(r$log_offset)) {
    add("Use either an exposure (logged by the model) or an offset that is ",
        "already on the log scale, not both (that would log it twice).")
  }
  if (!is.null(r$exposure) && r$exposure %in% cols) {
    e <- data[[r$exposure]]
    if (!is.numeric(e) || any(e[!is.na(e)] <= 0)) {
      add("The exposure must be positive numbers (it enters as ",
          "log(exposure)); zero exposure has no defined offset.")
    }
  }
  if (!is.null(r$weights) && r$weights %in% cols) {
    w <- data[[r$weights]]
    if (!is.numeric(w) || any(w[!is.na(w)] < 0)) {
      add("Frequency weights must be nonnegative numbers.")
    }
  }
  if (!is.null(r$group) && r$group %in% c(r$predictors, m$dis1, m$mixing)) {
    add("The panel identifier should not also be a predictor.")
  }
  if (!m$estimator %in% .app_estimators) {
    add("Unknown estimator \"", m$estimator, "\".")
  }
  ft <- fcr_family_table()
  if (!length(m$family) || any(!m$family %in% ft$id)) {
    add("Unknown family: ", toString(setdiff(m$family, ft$id)), ".")
  }
  mixture <- m$estimator %in% c("countreg_mix", "countreg_lc")
  if (mixture) {
    if (!is.numeric(m$components) || m$components < 1 ||
        m$components != round(m$components)) {
      add("The number of classes must be a whole number of at least 1.")
    }
    if (length(m$family) > 1 && length(m$family) != m$components) {
      add("Give one family, or one per class.")
    }
    if (!is.numeric(m$n_starts) || m$n_starts < 1) {
      add("The number of starts must be at least 1.")
    }
  } else {
    if (length(m$family) != 1) add("Count regression uses one family.")
    if (!is.null(r$group)) {
      add("A panel identifier is used by the mixture estimators (one class ",
          "per panel); for count regression leave it empty.")
    }
    if (length(m$mixing) || length(m$shared_mean)) {
      add("Membership predictors and shared terms apply to mixture models.")
    }
  }
  if (length(m$dis1) &&
      all(is.na(ft$param1[ft$id %in% m$family]))) {
    add("The ", toString(unique(m$family)), " family has no distribution ",
        "parameter to model; clear the distribution-parameter predictors.")
  }
  bad_shared <- setdiff(m$shared_mean, c("(Intercept)", r$predictors))
  if (length(bad_shared)) {
    add("Shared terms must be model terms: ", toString(bad_shared), ".")
  }
  if (!is.numeric(m$max_iters) || m$max_iters < 1) {
    add("The iteration limit must be at least 1.")
  }
  if (!is.numeric(m$seed) || length(m$seed) != 1) add("The seed must be a number.")
  p
}

# ---- Settings files -------------------------------------------------------------

# Write and read app settings
#
# Settings are stored as versioned JSON. Reading validates the schema and
# rejects unknown fields instead of ignoring them.
spec_to_json <- function(spec, path = NULL) {
  x <- unclass(spec)
  txt <- jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", pretty = TRUE,
                          digits = NA)
  if (!is.null(path)) writeLines(txt, path, useBytes = TRUE)
  txt
}

spec_from_json <- function(txt = NULL, path = NULL) {
  if (!is.null(path)) txt <- paste(readLines(path, warn = FALSE),
                                   collapse = "\n")
  x <- jsonlite::fromJSON(txt, simplifyVector = TRUE)
  if (!identical(as.integer(x$schema), .spec_schema)) {
    stop("Unsupported settings schema ", format(x$schema), " (this app ",
         "reads schema ", .spec_schema, ").", call. = FALSE)
  }
  extra <- setdiff(names(x), c("schema", "data", "recipe", "model"))
  if (length(extra)) {
    stop("Unknown top-level settings field(s): ", toString(extra), ".",
         call. = FALSE)
  }
  vec <- function(v) if (is.null(v) || (is.list(v) && !length(v)))
    character() else unlist(v)
  rec <- x$recipe
  for (nm in c("predictors", "factors")) rec[[nm]] <- vec(rec[[nm]])
  mod <- x$model
  for (nm in c("dis1", "mixing", "shared_mean")) mod[[nm]] <- vec(mod[[nm]])
  mod$family <- vec(mod$family)
  chk <- function(v, nm) {
    if (!is.null(v) && !(is.character(v) && length(v) == 1)) {
      stop("Settings field ", nm, " must be one name.", call. = FALSE)
    }
  }
  for (nm in c("response", "exposure", "log_offset", "weights", "group")) {
    # An empty JSON object or array means "not set"
    if (is.list(rec[[nm]]) && !length(rec[[nm]])) rec[nm] <- list(NULL)
    chk(rec[[nm]], nm)
  }
  app_spec(data = x$data, recipe = rec, model = mod)
}
