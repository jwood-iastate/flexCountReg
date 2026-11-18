# The app's view of the package: families, estimators and capabilities.
# Everything the interface needs to know about families and estimators
# comes from here (and from the registry in families.R), so a change in the
# statistical code is adapted in one file and caught by the app tests.

# Display names of the count families (the codes are those of
# ?flexCountReg::count_families).
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

# Count families available in the app
#
# One row per family of the installed 'flexCountReg': code, display name,
# the labels of its distribution parameters (`NA` when the family has no
# such parameter) and its status as a mixture component.
#
fcr_family_table <- function() {
  ids <- .fcr_families
  validated <- .fcr_fm_validated
  par <- lapply(ids, get_params)
  lab <- function(i) vapply(par, function(p) p[[i]] %||% NA_character_,
                            character(1))
  data.frame(id = ids,
             name = unname(ifelse(ids %in% names(.fcr_family_names),
                                  .fcr_family_names[ids], ids)),
             param1 = lab(1), param2 = lab(2),
             mixture_status = ifelse(ids %in% validated, "validated",
                                     "experimental"),
             stringsAsFactors = FALSE)
}

# Estimators offered by the app, with the core function that fits them.
.app_estimators <- c(
  "Count regression" = "countreg",
  "Finite mixture" = "countreg_mix",
  "Latent classes" = "countreg_lc")

# Capabilities of a family/estimator combination: which distribution-
# parameter formulas apply, and why an option is unavailable.
fcr_capabilities <- function(estimator, family) {
  ft <- fcr_family_table()
  fam <- ft[ft$id %in% family, , drop = FALSE]
  if (!nrow(fam)) {
    return(list(ok = FALSE, reason = paste("Unknown family", family)))
  }
  mixture <- estimator %in% c("countreg_mix", "countreg_lc")
  list(ok = TRUE,
       dis1 = any(!is.na(fam$param1)),
       dis1_label = paste(unique(stats::na.omit(fam$param1)), collapse = ", "),
       mixture = mixture,
       experimental = mixture && any(fam$mixture_status == "experimental"),
       intervals = !mixture,
       intervals_reason = if (mixture)
         "Intervals for mixture predictions are not yet available in flexCountReg.")
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

# Load a bundled dataset into a fresh environment (documented data access,
# never the global environment).
load_bundled_data <- function(name) {
  if (!name %in% fcr_bundled_datasets()$name) {
    stop("Unknown bundled dataset: ", name, call. = FALSE)
  }
  env <- new.env(parent = emptyenv())
  utils::data(list = name, package = "flexCountReg", envir = env)
  as.data.frame(env[[name]])
}

# Versions recorded with every fit and export.
fcr_versions <- function() {
  pv <- function(p) tryCatch(as.character(utils::packageVersion(p)),
                             error = function(e) NA_character_)
  c(R = paste(R.version$major, R.version$minor, sep = "."),
    flexCountReg = pv("flexCountReg"),
    shiny = pv("shiny"))
}
