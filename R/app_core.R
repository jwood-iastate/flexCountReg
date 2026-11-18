# The app's view of the package: families, estimators and which operations
# each fitted model supports. Everything the interface needs to know about
# the statistical functions comes from here (and from the family registry
# in families.R), so a change there is adapted in one file and caught by
# the app tests.

# Display names of the count families (codes as in ?count_families)
.fcr_family_names <- c(
  POISSON = "Poisson", GP1 = "Generalized Poisson 1",
  GP2 = "Generalized Poisson 2", NB1 = "Negative binomial 1",
  NB2 = "Negative binomial 2", NBP = "Negative binomial P",
  PLN = "Poisson-lognormal", PGE = "Poisson-generalized exponential",
  PIG1 = "Poisson-inverse Gaussian 1", PIG2 = "Poisson-inverse Gaussian 2",
  PIG = "Poisson-inverse gamma", PL = "Poisson-Lindley",
  PLG = "Poisson-Lindley-gamma", PLL = "Poisson-Lindley-lognormal",
  PW = "Poisson-Weibull", SI = "Sichel", GW = "Generalized Waring",
  COM = "Conway-Maxwell-Poisson", DW = "Discrete Weibull",
  DPO = "Double Poisson", DELAPORTE = "Delaporte",
  HERMITE = "Hermite (order 2)", HERMITE3 = "Hermite (order 3)",
  HERMITE4 = "Hermite (order 4)", HYPERPOIS = "Hyper-Poisson",
  GAMMACOUNT = "Gamma count", PA = "Polya-Aeppli (Poisson-Tweedie 1.5)",
  NBGE = "Negative binomial-generalized exponential")

# Count families: code, display name, distribution-parameter labels and the
# status of the family as a mixture component.
fcr_family_table <- function() {
  ids <- .fcr_families
  par <- lapply(ids, get_params)
  lab <- function(i) vapply(par, function(p) p[[i]] %||% NA_character_,
                            character(1))
  data.frame(id = ids,
             name = unname(ifelse(ids %in% names(.fcr_family_names),
                                  .fcr_family_names[ids], ids)),
             param1 = lab(1), param2 = lab(2),
             mixture_status = ifelse(ids %in% .fcr_fm_validated, "validated",
                                     "experimental"),
             stringsAsFactors = FALSE)
}

# Estimators offered by the app: the function that fits them, whether they
# use a count family, whether a panel identifier is required, allowed or
# not used, and a short description.
.app_estimators <- data.frame(
  id = c("countreg", "countreg_rp", "renb", "poisLind.re", "repois",
         "negmultinom", "countreg_mix", "countreg_lc"),
  label = c("Count regression", "Random parameters",
            "Random effects: negative binomial (panel)",
            "Random effects: Poisson-Lindley (panel)",
            "Random effects: Poisson (panel)",
            "Negative multinomial (panel)", "Finite mixture",
            "Latent classes"),
  uses_family = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
  group = c("none", "optional", "required", "required", "required",
            "required", "optional", "optional"),
  dis1 = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
  stringsAsFactors = FALSE)

.app_estimator_choices <- function() {
  stats::setNames(.app_estimators$id, .app_estimators$label)
}

.app_estimator <- function(id) {
  .app_estimators[match(id, .app_estimators$id), , drop = FALSE]
}

.app_is_mixture <- function(estimator) {
  estimator %in% c("countreg_mix", "countreg_lc")
}

# Capabilities of a family/estimator combination while building a model
fcr_capabilities <- function(estimator, family) {
  est <- .app_estimator(estimator)
  ft <- fcr_family_table()
  fam <- ft[ft$id %in% family, , drop = FALSE]
  uses_family <- isTRUE(est$uses_family)
  mixture <- .app_is_mixture(estimator)
  list(ok = nrow(est) == 1,
       uses_family = uses_family,
       dis1 = uses_family && isTRUE(est$dis1) && any(!is.na(fam$param1)),
       dis1_label = paste(unique(stats::na.omit(fam$param1)), collapse = ", "),
       mixture = mixture,
       group = est$group,
       experimental = mixture && any(fam$mixture_status == "experimental"))
}

# What a fitted model supports in the app, with the reason when it does not.
# Returns list(ok, reason).
app_supports <- function(fit, operation) {
  cls <- class(fit)
  yes <- list(ok = TRUE, reason = NULL)
  no <- function(r) list(ok = FALSE, reason = r)
  countreg <- inherits(fit, "countreg")
  mixture <- inherits(fit, "countreg_mix")
  switch(operation,
    rootogram = , quantile = if (countreg || mixture) yes else
      no(paste("Rootograms and quantile residuals need the fitted count",
               "distribution of each row; they are available for count",
               "regression and mixture models, not for", cls[1], "models.")),
    cure = yes,
    interval = if (countreg) yes else
      no(paste("Confidence and prediction intervals are available for",
               "count regression (countreg) models; flexCountReg does not",
               "yet compute them for", cls[1], "models.")),
    tolerance = if (countreg) yes else
      no(paste("Tolerance intervals are available for count regression",
               "(countreg) models; flexCountReg does not yet compute them",
               "for", cls[1], "models.")),
    no("Unknown operation."))
}

# Prediction methods of a fitted model (predict(..., method =)), if any
app_predict_methods <- function(fit) {
  if (inherits(fit, "countreg_rp")) {
    return(c("Exact (moment-generating functions)" = "Exact",
             "Simulated (Halton draws)" = "Simulated"))
  }
  if (inherits(fit, "negmultinom")) return(c("Standard" = "Standard"))
  NULL
}

# Example datasets bundled with flexCountReg
fcr_bundled_datasets <- function() {
  data.frame(
    name = c("abalone_rings", "stream_salamanders", "washington_roads"),
    title = c("Abalone shell-ring counts (UCI; CC BY 4.0)",
              "Stream salamander counts by site and species (CC0)",
              "Washington State road-segment crash counts"),
    response = c("rings", "count", "Total_crashes"),
    group = c(NA, "site", "ID"),
    stringsAsFactors = FALSE)
}

# Load a bundled dataset (documented data access, into a fresh environment)
load_bundled_data <- function(name) {
  if (!name %in% fcr_bundled_datasets()$name) {
    stop("Unknown bundled dataset: ", name, call. = FALSE)
  }
  env <- new.env(parent = emptyenv())
  utils::data(list = name, package = "flexCountReg", envir = env)
  as.data.frame(env[[name]])
}

# Versions recorded with every fit and export
fcr_versions <- function() {
  pv <- function(p) tryCatch(as.character(utils::packageVersion(p)),
                             error = function(e) NA_character_)
  c(R = paste(R.version$major, R.version$minor, sep = "."),
    flexCountReg = pv("flexCountReg"), shiny = pv("shiny"))
}
