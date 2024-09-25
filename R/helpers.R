# Helper functions for the flexCountReg package

# create function to clean data and run maxLik for bootstrapping
mod.boot <- function(data){
  # Prepare model matrices
  mod1_frame <- stats::model.frame(formula, data)
  X <- stats::model.matrix(formula, data)
  
  # If an offset is specified, create a vector for the offset
  if (!is.null(offset)){
    X_offset <- data[, offset]
  }
  
  if (!is.null(dis_param_formula_1)) {
    mod_alpha_frame <- as.matrix(modelr::model_matrix(data, dis_param_formula_1))
  } 
  
  if (!is.null(dis_param_formula_2)) {
    mod_sigma_frame <- as.matrix(modelr::model_matrix(data, dis_param_formula_2))
    
  } 
  
  if (!is.null(underreport_formula)){
    X_underreport <- as.matrix(modelr::model_matrix(data, underreport_formula))
  }
  y <- stats::model.response(mod1_frame)
  
  int_res <- maxLik::maxLik(logLikFunc, 
                            start = fit$estimate,
                            method = method, 
                            control = list(iterlim = max.iters, 
                                           printLevel = print.level))
  return(int_res)
}


# Determine model to estimate, probability distribution to use, and parameters
get_params <- function(family) {
  switch(
    family,
    "Poisson" = list(NULL, NULL),
    "NB1" = list("ln(alpha)", NULL),
    "NB2" = list("ln(alpha)", NULL),
    "NBP" = list("ln(alpha)", "ln(p)"),
    "PLN" = list("ln(sigma)", NULL),
    "PGE" = list("ln(shape)", "ln(scale)"),
    "PIG1" = list("ln(eta)", NULL),
    "PIG2" = list("ln(eta)", NULL),
    "PL" = list("ln(theta)", NULL),
    "PLG" = list("ln(theta)", "ln(alpha)"),
    "PLL" = list("ln(theta)", "ln(sigma)"),
    "PW" = list("ln(alpha)", "ln(sigma)"),
    "SI" = list("gamma", "ln(sigma)"),
    "GW" = list("ln(k)", "ln(rho)"),
    "COM" = list("ln(nu)", NULL),
  )
}

get_probFunc <- function(family){
  switch( # get the probability function for the specified distribution
  family,
  "Poisson" = function(y, predicted, alpha, sigma, haltons, normed_haltons) {
    return(stats::dpois(y, predicted))
  }
  "NB1" = function(y, predicted, alpha, sigma, haltons, normed_haltons) {
    mu <- predicted
    return(stats::dnbinom(y, size = mu/alpha, mu = mu))
  },
  "NB2" = function(y, predicted, alpha, sigma, haltons, normed_haltons) {
    mu <- predicted
    return(stats::dnbinom(y, size = alpha, mu = mu))
  },
  "NBP" = function(y, predicted, alpha, sigma, haltons, normed_haltons){
    mu <- predicted
    return(stats::dnbinom(y, size = (mu^(2-sigma))/alpha, mu = mu))
  },
  "PLN" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dpLnorm_cpp(x=y, mean=predicted, sigma=alpha, h=normed_haltons),
  "PGE" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dpge(y, mean=predicted, shape=alpha, scale=sigma, haltons=haltons),
  "PIG1" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dpinvgaus(y, mu=predicted, eta=alpha),
  "PIG2" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dpinvgaus(y, mu=predicted, eta=alpha, form="Type 2"),
  "PL" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dplind(y, mean=predicted, theta=alpha),
  "PLG" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dplindGamma(x=y, mean=predicted, theta=alpha, alpha=sigma, hdraws=haltons),
  "PLL" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dplindLnorm(x=y, mean=predicted, theta=alpha, sigma=sigma, hdraws=normed_haltons),
  "PW" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dpWeib_cpp(y, mean=predicted, alpha=alpha, sigma=sigma, h=haltons),
  "SI" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dsichel(x=y, mu= predicted, sigma=sigma, gamma=log(alpha)),
  "GW" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dgwar(y, mu= predicted, k=alpha, rho=sigma),
  "COM" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dcom(x=y, mu=predicted, nu=alpha)
  )
}