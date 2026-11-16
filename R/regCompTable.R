#' Create a Table Comparing Regression Models with AIC, BIC, and McFadden's
#' Pseudo-R-Squared
#'
#' This function creates tables comparing the flexCountReg package models
#' supplied to the function.
#'
#' @name regCompTable
#' @param models A named list of fitted flexCountReg model objects. This must
#'   include 2 or more models.
#' @param coefs A logical. The default value `TRUE` indicates that the
#'   coefficients from the models should be included in the table of
#'   comparisons.
#' @param AIC A logical. The default value `TRUE` indicates that AIC values for
#'   the models should be included.
#' @param BIC A logical. The default value `TRUE` indicates that BIC values for
#'   the models should be included.
#' @param RSquare A logical. The default value `TRUE` indicates that the
#'   McFadden's Pseudo-R-Squared statistic should be included. By default it
#'   compares every model against the same intercept-only Poisson model (the
#'   base count model), so the values are comparable across models; see
#'   \code{basemodel} and \code{variables}.
#' @param basemodel the family of the baseline model for McFadden's
#'   Pseudo-R-Squared: "Poisson" (default) or a family of
#'   \code{\link{countreg}}.
#' @param variables logical; if \code{TRUE}, the baseline model includes the
#'   mean-model covariates of each model; if \code{FALSE} (default), it is an
#'   intercept-only model (with the model's offset).
#' @param tableType The type of table format to return. Options include "tibble"
#'   for returning the table as a tibble, "gt" for a \link[gt]{gt} table object,
#'   or "latex" for a latex table. The default is "tibble".
#' @param digits An integer value indicating the number of decimals to round the
#'   table values to.
#'
#' @include regCompTest.R
#' @importFrom knitr kable
#' @importFrom tibble tibble
#' @importFrom stats as.formula model.frame model.response dpois glm
#' @importFrom gt gt tab_header fmt_number
#' 
#' @returns A table comparing the models supplied to the function in the format 
#'  specified The table includes coefficients, standard errors, statistical 
#'  significance, AIC, BIC, and McFadden's Pseudo-R-Squared. Table formats 
#'  include: "tibble", "gt", and "latex".
#' 
#' @examples
#' 
#' # Comparing the NBP model with the NB2 and NB1 models
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
#'
#' nb.1 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                     ShouldWidth04 + AADTover10k,
#'                     data=washington_roads, family = 'NB1', method = 'NM')
#'                     
#' nb.2 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                     ShouldWidth04 + AADTover10k,
#'                     data=washington_roads, family = 'NB2', method = 'NM')
#'                     
#'                     
#' nb.p <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                     ShouldWidth04 + AADTover10k,
#'                     data=washington_roads, family = 'NBP', method = 'NM')
#'                     
#'                     
#' comptable <- 
#'  regCompTable(list("NB-1"=nb.1, "NB-2"=nb.2, "NB-P"=nb.p), tableType="latex")
#' print(comptable)
#' @export
regCompTable <- function(models, coefs=TRUE, AIC=TRUE, BIC=TRUE, RSquare=TRUE, 
                         tableType="tibble", digits=3,
                         basemodel = "Poisson", variables = FALSE){
  
  Nmodels <- length(models)
  if(Nmodels < 2) {
    msg <- paste0("Provide a list of 2 or more models for the `models` input ", 
                  "to call the `regCompTable` function")
    warning(msg)
  }
  
  # Extract the names of all coefficients in any of the models supplied and
  # create an initial tibble with a column of the coefficient names and other
  # stats
  vars <- c() # names of estimated coefficients
  modNames <- names(models)
  for (i in seq_along(models)){
    coef_names_i <- names(models[[i]]$model$estimate) # Coefficient names
    for (item in coef_names_i) { # Add any new, unique, variables/coefficient names
      if (!(item %in% vars)) {
        vars <- append(vars, item)
      }
    }
  }
  vars <- append(vars, "N Obs.") # Add the number of observed datapoints
  vars <- append(vars, "LL") # Add the Log Likelihood
  
  if(AIC){
    vars <- append(vars, "AIC")
  }
  if(BIC){
    vars <- append(vars, "BIC")
  }
  if(RSquare){
    vars <- append(vars, "Pseudo-R-Sq.")
  }
  
  summaryTibble <- tibble(`Parameter` = vars) # Initial tibble
  
  for (i in seq_along(modNames)){ 
    m_i <- models[[i]]
    # Log-likelihood, AIC and BIC from the model's methods; McFadden's R^2
    # against the same null model as regCompTest()
    LL <- as.numeric(stats::logLik(m_i))
    LLbase <- .fit_null_ll(m_i, .fcr_fit_pieces(m_i), basemodel = basemodel,
                           variables = variables)
    aic <- round(stats::AIC(m_i), digits)
    bic <- round(stats::BIC(m_i), digits)
    rsq <- round(.mcfadden_r2(LL, LLbase), digits)
    Nobs <- stats::nobs(m_i)
    
    coefs <- stats::coef(m_i)
    coefNames <- names(coefs)
    se <- .fcr_se(stats::vcov(m_i))
    t <- abs(coefs / se)
    coefsRounded <- round(coefs, digits)
    seRounded <- round(se, digits)
    tablecoefs <- c()
    
    for (j in seq_along(vars)){
      if (vars[j] %in% coefNames){
        for (k in seq_along(coefNames)){
          if (coefNames[k] == vars[j]){
            # Standard errors are NA when the estimated variance is
            # negative (weakly identified parameters): no stars then
            coefval <- paste0(coefsRounded[k], " (", seRounded[k], ")",
                              .fcr_stars(t[k]))
            tablecoefs <- append(tablecoefs, coefval)
          }
        }
      } else if (vars[j] == "N Obs.") {
        tablecoefs <- append(tablecoefs, Nobs)
      } else if (vars[j]  == "LL") {
        tablecoefs <- append(tablecoefs, round(as.numeric(LL), digits))
      } else if (vars[j]  == "AIC") {
        tablecoefs <- append(tablecoefs, aic)
      } else if (vars[j]  == "BIC") {
        tablecoefs <- append(tablecoefs, bic)
      } else if (vars[j]  == "Pseudo-R-Sq.") {
        tablecoefs <- append(tablecoefs, rsq)
      } else{
        tablecoefs <- append(tablecoefs, "---")
      }
    }
    summaryTibble[[modNames[i]]] <- tablecoefs
  }
  
  # Adding a row with the notation
  # notation_row <- tibble(
  #   Parameter = "Note",
  #   !!!setNames(
  #     rep("p-value codes: * (p<=0.05), ** (p<=0.01), *** (p<=0.001)",
  #         length(modNames)), modNames)
  # )
  
  
  if(tableType == "tibble"){
    return(summaryTibble)
  } 
  else if (tableType == "gt"){
    gtTable <- gt(summaryTibble) |>
      tab_header(
        title = "Model Comparisons",
        subtitle = "Mean (Standard Error)"
      ) |>
      gt::tab_source_note(
        source_note = "p-value codes: *=(p<=0.05), **=(p<=0.01), ***=(p<=0.001)"
      )
    return(gtTable)
  }
  else if (tableType == "latex"){
    latexTable <- knitr::kable(summaryTibble, 
                               format = "latex", 
                               booktabs = TRUE, 
                               caption = "Model Comparison Statistics")
    return(latexTable)
  }
}

# Significance stars for an absolute z statistic ("" when NA).
.fcr_stars <- function(t) {
  if (is.na(t)) return("")
  if (t >= 3.29) "***" else if (t >= 2.5758) "**" else if (t >= 1.96) "*" else ""
}
