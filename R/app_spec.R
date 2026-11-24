# The versioned model specification of the app. Every Build control edits a
# draft specification; a fit freezes a copy, and the code that runs and the
# code that is exported are both compiled from it.

# Schema 3 added recipe$interactions; schema 4 added smooth terms, the
# hierarchical and fixed-effects models, correlated random effects,
# cluster-robust standard errors and sampling weights. Older settings are
# migrated on reading (see .spec_migrate()).
.spec_schema <- 4L

# Create a model specification.
#   data   - source ("bundled" or "upload"), name (dataset or file name) and
#            object (the name of the data object in the session script,
#            e.g. "data_1").
#   recipe - response, predictors, factors (columns converted with
#            factor()), exposure (entered as offset(log(exposure))),
#            log_offset (entered as offset(column)), weights (frequency
#            weights), group (panel identifier), interactions (a list
#            of pairs of predictor names; each pair enters the mean formula
#            as a:b next to the main effects), smooth (predictors entered
#            as smooth terms s(x)), group2 (a second grouping column:
#            crossed or nested random effects, or a second fixed-effect
#            dimension) and cluster (the cluster column of cluster-robust
#            standard errors).
#   model  - estimator (see .app_estimators), family (one code, or one per
#            mixture component), components, dis1 (predictors of the first
#            distribution parameter), mixing (membership predictors),
#            shared_mean, shared_dis1, rpar (random parameters),
#            rpar_intercept (random intercept), rpar_dist, correlated,
#            ndraws, nm_family (negative multinomial), n_starts, seed,
#            method, max_iters, stderr ("normal", "robust" or "cluster"),
#            weight_type ("frequency" or "sampling"), cre (predictors
#            whose group means are added), re_slopes (random slopes of a
#            hierarchical model), re_correlated, group2_type ("crossed"
#            or "nested"), re_dist (distribution of hierarchical random
#            effects) and fe_jackknife.
app_spec <- function(data = list(), recipe = list(), model = list()) {
  d <- list(source = "bundled", name = "abalone_rings", object = "data_1")
  r <- list(response = NULL, predictors = character(), factors = character(),
            exposure = NULL, log_offset = NULL, weights = NULL,
            group = NULL, interactions = list(), smooth = character(),
            group2 = NULL, cluster = NULL)
  m <- list(estimator = "countreg", family = "NB2", components = 2L,
            dis1 = character(), mixing = character(),
            shared_mean = character(), shared_dis1 = FALSE,
            rpar = character(), rpar_intercept = FALSE, rpar_dist = "n",
            correlated = FALSE, ndraws = 500L, nm_family = "NM2",
            n_starts = 10L, seed = 1L, method = "BFGS", max_iters = 1000L,
            stderr = "normal", weight_type = "frequency", cre = character(),
            re_slopes = character(), re_correlated = TRUE,
            group2_type = "crossed", re_dist = "normal",
            fe_jackknife = FALSE)
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
  est <- .app_estimator(m$estimator)
  if (nrow(est) != 1 || is.na(est$id)) {
    add("Unknown model type \"", m$estimator, "\".")
    return(p)
  }
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
             r$group, r$group2, r$cluster, m$dis1, m$mixing, m$rpar)
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
  # Panels
  if (est$group == "required" && is.null(r$group)) {
    add("This model needs a panel identifier (the column that groups the ",
        "repeated observations).")
  }
  if (est$group == "none" && !is.null(r$group)) {
    add("Count regression does not use a panel identifier; choose a ",
        "random-parameter, random-effects or mixture model for panels.")
  }
  if (!is.null(r$group) && r$group %in% c(r$predictors, m$dis1, m$mixing,
                                          m$rpar, r$group2)) {
    add("The panel identifier should not also be a predictor.")
  }
  # Families
  ft <- fcr_family_table()
  if (est$uses_family) {
    if (!length(m$family) || any(!m$family %in% ft$id)) {
      add("Unknown family: ", toString(setdiff(m$family, ft$id)), ".")
    }
  }
  if (m$estimator == "negmultinom" &&
      !m$nm_family %in% c("NM1", "NM2", "NMP")) {
    add("The negative multinomial family must be NM1, NM2 or NMP.")
  }
  mixture <- .app_is_mixture(m$estimator)
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
    if (est$uses_family && length(m$family) != 1) {
      add("This model uses one family.")
    }
    if (length(m$mixing) || length(m$shared_mean)) {
      add("Membership predictors and shared terms apply to mixture models.")
    }
  }
  # Random parameters
  if (m$estimator == "countreg_rp") {
    if (!length(m$rpar) && !isTRUE(m$rpar_intercept)) {
      add("Choose at least one random parameter (or a random intercept).")
    }
    if (!m$rpar_dist %in% c("n", "ln", "t", "u", "g")) {
      add("Unknown random-parameter distribution \"", m$rpar_dist, "\".")
    }
    if (isTRUE(m$correlated) && m$rpar_dist != "n") {
      add("Correlated random parameters must be normal.")
    }
    if (!is.numeric(m$ndraws) || m$ndraws < 10) {
      add("Use at least 10 Halton draws.")
    }
  } else if (length(m$rpar) || isTRUE(m$rpar_intercept)) {
    add("Random parameters apply to the random-parameters model.")
  }
  if (length(m$dis1)) {
    if (!est$dis1) {
      add("Distribution-parameter predictors are not available for this ",
          "model type.")
    } else if (all(is.na(ft$param1[ft$id %in% m$family]))) {
      add("The ", toString(unique(m$family)), " family has no distribution ",
          "parameter to model; clear the distribution-parameter predictors.")
    }
  }
  # Smooth terms: numeric predictors entered as s(x)
  if (length(r$smooth)) {
    if (!isTRUE(est$smooth)) {
      add("Smooth effects are not available for this model type.")
    }
    bad <- setdiff(r$smooth, r$predictors)
    if (length(bad)) {
      add("Smooth effects must be predictors: ", toString(bad), ".")
    }
    sm <- intersect(r$smooth, cols)
    nonnum <- sm[!vapply(data[sm], is.numeric, logical(1)) |
                   sm %in% r$factors]
    if (length(nonnum)) {
      add("Smooth effects need numeric predictors (not categorical): ",
          toString(nonnum), ".")
    }
    if (length(intersect(r$smooth, m$rpar))) {
      add("A predictor cannot be both a smooth effect and a random ",
          "parameter: ", toString(intersect(r$smooth, m$rpar)), ".")
    }
    few <- sm[vapply(data[setdiff(sm, nonnum)], function(v)
      length(unique(stats::na.omit(v))) < 10, logical(1))]
    if (length(few)) {
      add("A smooth effect needs a predictor with at least 10 distinct ",
          "values: ", toString(few), ".")
    }
  }
  # Hierarchical and fixed-effects models
  hier <- m$estimator == "countreg_hier"
  fe <- m$estimator == "countreg_fe"
  if (!is.null(r$group2)) {
    if (!hier && !fe) {
      add("A second grouping column applies to hierarchical and ",
          "fixed-effects models.")
    } else if (identical(r$group2, r$group)) {
      add("The second grouping column must differ from the first.")
    }
  }
  if (length(m$re_slopes)) {
    if (!hier) {
      add("Random slopes apply to the hierarchical model.")
    } else if (length(setdiff(m$re_slopes, r$predictors))) {
      add("Random slopes must be predictors: ",
          toString(setdiff(m$re_slopes, r$predictors)), ".")
    }
  }
  if (!m$group2_type %in% c("crossed", "nested")) {
    add("The second grouping must be crossed or nested.")
  }
  if (!m$re_dist %in% .fcr_re_dists) {
    add("Unknown random-effect distribution \"", m$re_dist, "\".")
  } else if (hier && m$re_dist != "normal") {
    if (!is.null(r$group2)) {
      add("Non-normal random effects need a single grouping column.")
    }
    if (m$re_dist == "gamma" && length(m$re_slopes)) {
      add("Gamma random effects are available for a random intercept ",
          "only (no random slopes).")
    }
  }
  # Correlated random effects
  if (length(m$cre)) {
    if (!isTRUE(est$cre)) {
      add("Correlated random effects (group means) are not available for ",
          "this model type.")
    } else if (is.null(r$group)) {
      add("Correlated random effects need the panel identifier.")
    }
    bad <- setdiff(m$cre, r$predictors)
    if (length(bad)) {
      add("Group-mean (correlated random effects) variables must be ",
          "predictors: ", toString(bad), ".")
    }
    cr <- intersect(m$cre, cols)
    if (any(!vapply(data[cr], is.numeric, logical(1)) | cr %in% r$factors)) {
      add("Group means need numeric predictors (not categorical).")
    }
  }
  # Standard errors and weights
  if (!m$stderr %in% c("normal", "robust", "cluster")) {
    add("Unknown standard-error type \"", m$stderr, "\".")
  } else if (m$stderr == "robust" && !isTRUE(est$robust)) {
    add("Robust (sandwich) standard errors are not available for this ",
        "model type.")
  } else if (m$stderr == "cluster") {
    if (!isTRUE(est$cluster)) {
      add("Cluster-robust standard errors are not available for this ",
          "model type.")
    } else if (is.null(r$cluster) && !fe) {
      add("Choose the cluster column for cluster-robust standard errors.")
    }
  }
  if (!is.null(r$cluster) && m$stderr != "cluster") {
    add("A cluster column is used only with cluster-robust standard ",
        "errors.")
  }
  if (!m$weight_type %in% c("frequency", "sampling")) {
    add("Weights must be frequency or sampling weights.")
  } else if (m$weight_type == "sampling") {
    if (!isTRUE(est$sampling)) {
      add("Sampling weights are available for count regression and ",
          "random-parameter models.")
    } else if (is.null(r$weights)) {
      add("Sampling weights need a weights column.")
    }
  }
  # Interactions: pairs of two different chosen predictors (random
  # parameters of countreg_rp are among the predictors too)
  for (pr in r$interactions) {
    pr <- unlist(pr)
    if (!is.character(pr) || length(pr) != 2 || anyNA(pr)) {
      add("Each interaction must name two predictors.")
    } else if (pr[1] == pr[2]) {
      add("An interaction needs two different predictors (\"", pr[1],
          "\" is given twice).")
    } else if (!all(pr %in% r$predictors)) {
      add("The interaction ", pr[1], " \u00d7 ", pr[2], " uses ",
          toString(paste0("\"", setdiff(pr, r$predictors), "\"")),
          ", which is not among the predictors; choose it as a predictor ",
          "first.")
    }
  }
  keys <- vapply(r$interactions, function(pr)
    paste(sort(as.character(unlist(pr))), collapse = "\r"), "")
  if (anyDuplicated(keys)) add("An interaction is given more than once.")
  bad_shared <- setdiff(m$shared_mean, c("(Intercept)", r$predictors))
  if (length(bad_shared)) {
    add("Shared terms must be model terms: ", toString(bad_shared), ".")
  }
  if (!is.numeric(m$max_iters) || is.na(m$max_iters) || m$max_iters < 1) {
    add("The iteration limit must be at least 1.")
  }
  if (!is.numeric(m$seed) || length(m$seed) != 1 || is.na(m$seed)) {
    add("The seed must be a number.")
  }
  p
}

# ---- Settings files -------------------------------------------------------------
#
# Specifications are stored as versioned JSON. Reading validates the schema
# and rejects unknown fields instead of ignoring them.

spec_to_json <- function(spec, path = NULL) {
  txt <- jsonlite::toJSON(unclass(spec), auto_unbox = TRUE, null = "null",
                          pretty = TRUE, digits = NA)
  if (!is.null(path)) writeLines(txt, path, useBytes = TRUE)
  txt
}

spec_from_json <- function(txt = NULL, path = NULL) {
  if (!is.null(path)) txt <- paste(readLines(path, warn = FALSE),
                                   collapse = "\n")
  x <- jsonlite::fromJSON(txt, simplifyVector = TRUE)
  .spec_from_list(x)
}

# Bring settings of an older schema up to the current one. Schema 2 had no
# interactions.
.spec_migrate <- function(x) {
  sc <- suppressWarnings(as.integer(x$schema))
  if (length(sc) != 1 || is.na(sc) || !sc %in% 2:.spec_schema) {
    stop("Unsupported settings schema ", format(x$schema), " (this app ",
         "reads schemas 2 to ", .spec_schema, ").", call. = FALSE)
  }
  if (sc == 2L) x$recipe["interactions"] <- list(list())
  # Fields added by later schemas take their defaults (app_spec())
  x$schema <- .spec_schema
  x
}

# Interactions read from JSON: jsonlite simplifies a list of pairs to a
# two-column character matrix (or keeps a list); either becomes a list of
# character pairs
.spec_pairs <- function(v) {
  if (is.null(v) || !length(v)) return(list())
  if (is.matrix(v)) {
    return(lapply(seq_len(nrow(v)), function(i) as.character(v[i, ])))
  }
  if (is.list(v)) return(lapply(v, function(p) as.character(unlist(p))))
  stop("Settings field interactions must be a list of pairs of names.",
       call. = FALSE)
}

.spec_from_list <- function(x) {
  x <- .spec_migrate(x)
  extra <- setdiff(names(x), c("schema", "data", "recipe", "model"))
  if (length(extra)) {
    stop("Unknown top-level settings field(s): ", toString(extra), ".",
         call. = FALSE)
  }
  vec <- function(v) if (is.null(v) || (is.list(v) && !length(v)))
    character() else unlist(v)
  rec <- x$recipe
  for (nm in c("predictors", "factors", "smooth")) {
    if (!is.null(rec[[nm]])) rec[[nm]] <- vec(rec[[nm]])
  }
  rec["interactions"] <- list(.spec_pairs(rec$interactions))
  mod <- x$model
  for (nm in c("dis1", "mixing", "shared_mean", "rpar", "family", "cre",
               "re_slopes")) {
    if (!is.null(mod[[nm]]) || nm %in% c("dis1", "mixing", "shared_mean",
                                         "rpar", "family")) {
      mod[[nm]] <- vec(mod[[nm]])
    }
  }
  chk <- function(v, nm) {
    if (!is.null(v) && !(is.character(v) && length(v) == 1)) {
      stop("Settings field ", nm, " must be one name.", call. = FALSE)
    }
  }
  for (nm in c("response", "exposure", "log_offset", "weights", "group",
              "group2", "cluster")) {
    # An empty JSON object or array means "not set"
    if (is.list(rec[[nm]]) && !length(rec[[nm]])) rec[nm] <- list(NULL)
    chk(rec[[nm]], nm)
  }
  app_spec(data = x$data, recipe = rec, model = mod)
}
