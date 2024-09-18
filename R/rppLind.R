#' Function for estimating a random parameter Poisson-Lindley model
#'
#' @name rppLind
#' @param formula an R formula.. This formula should specify the outcome and the independent variables that have fixed parameters.
#' @param rpar_formula a symbolic description of the model related specifically to the random parameters. This should not include an outcome variable. If the intercept is random, include it in this formula. If the intercept is fixed, include it in \code{formula} but not in \code{rpar_formula}. To remove the intercept, use \code{0 + vars} or \code{-1 + vars},
#' @param data a dataframe that has all of the variables in the \code{formula} and \code{rpar_formula},
#' @param rpardists an optional named vector whose names are the random parameters and values the distribution. The distribution options include normal ("n"), lognormal ("ln"), triangular ("t"), uniform ("u"), and gamma ("g"). If this is not provided, normal distributions are used for all random coefficients,
#' @param ndraws the number of Halton draws to use for estimating the random parameters,
#' @param scrambled if the Halton draws should be scrambled or not. \code{scrambled = FALSE} results in standard Halton draws while \code{scrambled = TRUE} results in scrambled Halton draws,
#' @param correlated if the random parameters should be correlated (\code{correlated = FALSE} results in uncorrelated random coefficients, \code{correlated = TRUE} results in correlated random coefficients). If the random parameters are correlated, only the normal distribution is used for the random coefficients,
#' @param method a method to use for optimization in the maximum likelihood estimation. For options, see \code{\link[maxLik]{maxLik}},
#' @param max.iters the maximum number of iterations to allow the optimization method to perform,
#' @param start.vals an optional vector of starting values for the regression coefficients
#' @param print.level Integer specifying the verbosity of output during optimization.
#' @param bootstraps Optional integer specifying the number of bootstrap samples to be used
#'        for estimating standard errors. If not specified, no bootstrapping is performed.
#' @param weights the name of the weighting variable. This is an optional 
#'   parameter for weighted regression.
#'    
#' @import randtoolbox maxLik stats modelr
#' @importFrom MASS glm.nb
#' @importFrom utils head  tail
#' @importFrom purrr map map_df
#' @importFrom broom tidy
#' @importFrom dplyr group_by %>% summarise
#' @importFrom tibble deframe
#' @include tri.R plind.R
#'
#' @examples
#' \donttest{
#' ## Random Parameters Poisson-Lindley
#' data("washington_roads")
#' plind.rp <- rppLind(Animal ~ lnlength + lnaadt,
#'                              rpar_formula = ~ -1 + speed50,
#'                              data = washington_roads,
#'                              ndraws = 10,
#'                              correlated = FALSE,
#'                              rpardists = c(speed50="n"),
#'                              method = "nm",
#'                              print.level = 2)
#'
#' summary(plind.rp)}
#' @export
rppLind <- function(formula, rpar_formula, data, 
                 rpardists = NULL,
                 ndraws = 1500, scrambled = FALSE,
                 correlated = FALSE, method = 'BHHH', max.iters = 1000,
                 start.vals = NULL, print.level = 0,
                 weights=NULL, bootstraps=NULL) {
  # start.vals can be a vector or a named vector with the starting values for the parameters
  # print.level is used to determine the level of details for the optimization to print (for the maxLik function call)

  sd.start <- 0.1 # starting value for each of the standard deviations of random parameters

  # Check and correct the rpardists if the random parameters are correlated
  if(correlated && !is.null(rpardists) && sum(ifelse(rpardists=="n",0,1))>0){
    print("When the random parameters are correlated, only the normal distribution is used.")
    rpardists <- NULL
  }

  # Function to generate Halton draws
  halton_draws <- function(ndraws, rpar, scrambled) {
    num_params <- length(rpar)
    halton_seq <- randtoolbox::halton(ndraws, num_params, mixed = scrambled)
    return(halton_seq)
  }

  # Generating model matrices
  mod1_frame <- stats::model.frame(formula, data)
  X_Fixed <- modelr::model_matrix(data, formula)
  X_rand <- modelr::model_matrix(data, rpar_formula)
  #X_rand <-stats::model.matrix(rpar_formula, data)
  y_name <- all.vars(formula)[1]
  y <- stats::model.response(mod1_frame)
  
  
  # Preparing the weights
  if (!is.null(weights)){
    if (is.character(weights) && weights %in% colnames(data)) {
      # Use weights variable from data
      weights <- as.numeric(data[[weights]])
    } else {
      stop("Weights should be the name of a variable in the data.")
    }
  } else {
    weights <- rep(1, length(y)) # Default weights of 1
  }

  # Create named vectors
  fixed_terms <- names(X_Fixed)
  rand_terms <- names(X_rand)
  predictor_terms <- c(fixed_terms, rand_terms)
  nb_vars <- predictor_terms[!grepl("ntercept", predictor_terms)] # remove the intercept

  X_Fixed <- as.matrix(X_Fixed)
  X_rand <- as.matrix(X_rand)

  # check if the intercept is included in both the fixed and random parameters
  if("\\(Intercept\\)" %in% colnames(X_Fixed) && "\\(Intercept\\)" %in% colnames(X_rand)){
    print("Do not include the intercept in both the fixed parameters (in `formula`) and random parameters (in `rpar_formula`). Use `- 1` in the formula you want to remove the intercept from.")
    stop()
  }


  if(!is.null(rpardists)){ # check the random parameter distributions for the intercept (and correct, if needed)
    names(rpardists) <- gsub("intercept", "\\(Intercept\\)", names(rpardists)) # correct intercept name
    names(rpardists) <- gsub("constant", "\\(Intercept\\)", names(rpardists)) # correct intercept name
  }

  # Check if the specified distributions include all specified random parameters and no extras
  if (!correlated && !is.null(rpardists)){
    prob = FALSE
    for (i in rand_terms){
      if(i %in% names(rpardists)){
        next
      }
      else{
        print(paste("Variable ", i, " is specified as random but is not included in the `rpardists`. Please include it or do not specify the distributions."))
        prob = TRUE
      }
    }

    for (i in names(rpardists)){
      if (i %in% rand_terms){
        next
      }
      else{
        print(paste("Variable ", i, " is included in `rpardists` but is not in the specified random parameters. Please include it as a random parameter or remove it from `rpardists`."))
        prob = TRUE
      }
    }
    if (prob) stop() # if a problem was identified, the instructions were printed and the function will now stop executing
  }

  nb_formula <- reformulate(nb_vars, response = y_name)

  x_fixed_names <- colnames(X_Fixed)
  rpar <- colnames(X_rand)

  hdraws <- halton_draws(ndraws, rpar, scrambled)

  if(!is.null(start.vals)){
    params <- unname(start.vals)
    Lparams <- length(params)
    Lrpar = length(rpar)
    start <- params
    x_names <- names(start.vals)
  }
  else{
    nb_model <- glm.nb(nb_formula, data)
    params <- coef(nb_model)
    Lparams <- length(params)
    Lrpar = length(rpar)

    if (length(rpar)<2){
      correlated = FALSE
    }

    start <- params
    # Make sure the start parameters are ordered correctly, including the location of the intercept
    varnameorder <- c(x_fixed_names, rpar)
    match_indices <- match(varnameorder, names(start))
    start <- start[match_indices] # correctly ordered

    x_rand_names_mean <- paste0(rpar, ':Mean')
    x_names = c(x_fixed_names, x_rand_names_mean)

    if (correlated){
      randparam_means = tail(start, Lrpar)
      rparam_var <- abs(randparam_means)/2
      rparam_var <- diag(rparam_var)
      Chl <- chol(rparam_var)

      for (i in 1:length(rpar)){
        for (j in 1:length(rpar)){
          if (i >= j){
            start <- append(start, Chl[j,i])
            x_names <- append(x_names, paste('Cholesky Value for' ,paste(rpar[j], rpar[i], sep=":")))
          }
        }
      }
    }
    else{
      start <- append(start, rep(sd.start, length(rpar)))
      x_rand_names_sd <- paste0(rpar, ':St. Dev.')
      x_names <- c(x_names, x_rand_names_sd)
    }
    
    start <- append(start, log(0.1)) # initial log of theta parameter
    x_names <- append(x_names, 'ln(alpha)')
  
    t <- 116761/exp(-9.04*nb_model$theta)
    theta <- ifelse(t<100,t,1) # Approximate Initial Theta
      
    start <- append(start, log(theta))
    x_names <- append(x_names, 'ln(theta)')
    names(start) <- x_names
  }
  
  lind_dens<-function(mu, theta){
    return(dplind(x=y, mean=mu, theta=theta))
  }
    
  # main function for estimating log-likelihoods
  p_plind_rp <- function(p, y, X_Fixed, X_rand, ndraws, rpar, correlated){
    est_method <- method
    if (!correlated) exact.gradient=TRUE else exact.gradient=FALSE # use numerical gradient if using correlated random parameters
    N_fixed = length(x_fixed_names)
    N_rand = length(rpar)
    coefs <- as.array(p)
    fixed_coefs <- head(coefs,N_fixed)
    h <- head(coefs, (N_fixed + N_rand))
    random_coefs_means <- tail(h, N_rand)
    
    dist_params <- tail(coefs,1)
    log_theta <- dist_params[1]


    theta <- exp(log_theta)
    mu_fixed <- exp(X_Fixed %*% fixed_coefs)

    if(length(rpar)==1){
      rand_sdevs <- coefs[length(coefs)-1]
    }else{
      t <- head(t,-1)
    }

    # generate and scale random draws
    if (length(rpar)>1){
      if (correlated){ # Generate correlated random draws
        chol_vals <- t
        Ch <- matrix(0, N_rand, N_rand)
        counter = 1
        for (i in 1:N_rand){
          for (j in 1:N_rand){
            if (j<=i){
              Ch[j,i] <- chol_vals[counter]
              counter <- counter + 1
            }
          }
        }
        scaled_draws <- qnorm(hdraws) %*% Ch
        draws <- apply(scaled_draws, 1, function(x) x + random_coefs_means)
      }
      else{
        rand_sdevs <- t
        draws <- hdraws #initialize the matrix

        if (is.null(rpardists)){
          draws <- apply(hdraws, 1, function(x) stats::qnorm(x, random_coefs_means, rand_sdevs))
        }
        else{
          for (i in 1:length(rpar)){
            if (rpardists[i]=="ln"){
              draws[,i] <- stats::qlnorm(hdraws[,i], random_coefs_means[i], abs(rand_sdevs[i]))
            }
            else if (rpardists[i]=="t"){
              draws[,i] <- qtri(hdraws[,i], mode=random_coefs_means[i], sigma=abs(rand_sdevs[i]))
            }
            else if (rpardists[i]=="u"){
              draws[,i] <- random_coefs_means[i] + (hdraws[,i] - 0.5)*abs(rand_sdevs[i])
            }
            else if (rpardists[i]=="g"){
              draws[,i] <- stats::qgamma(hdraws[,i], shape=random_coefs_means[i]^2/(rand_sdevs[i]^2), 
                                               rate=random_coefs_means[i]/(rand_sdevs[i]^2))
            }
            else{
              draws[,i] <- stats::qnorm(hdraws[,i], random_coefs_means[i], abs(rand_sdevs[i]))
            }
          }
        }
        draws <- t(draws)
      }
      xb_rand_mat <- crossprod(t(X_rand), draws)
    }
    else{
      if (is.null(rpardists)){
        scaled_draws <- hdraws * rand_sdevs
        draws <- scaled_draws + random_coefs_means[1]
      }
      else{
        if (rpardists[1]=="ln"){
          draws <- stats::qlnorm(hdraws, random_coefs_means, abs(rand_sdevs))
        }
        else if (rpardists[1]=="t"){
          draws <- qtri(hdraws, random_coefs_means, abs(rand_sdevs))
        }
        else if (rpardists[1]=="u"){
          draws <- random_coefs_means + (hdraws-0.5)*abs(rand_sdevs)
        }
        else if (rpardists[1]=="g"){
          draws <- stats::qgamma(hdraws, shape=random_coefs_means^2/(rand_sdevs^2), 
                                 rate=random_coefs_means/(rand_sdevs^2))
        }
        else{
          draws <- stats::qnorm(hdraws, random_coefs_means, abs(rand_sdevs))
        }

      }
      xb_rand_mat <- sapply(draws, function(x) X_rand * x)
    }

    rpar_mat <- exp(xb_rand_mat)
    pred_mat <- apply(rpar_mat, 2, function(x) x * mu_fixed)
    
    if (sum(is.na(rpar_mat))>0) print('NA in rpar_mat')
    if (sum(is.na(pred_mat))>0) print('NA in pred_mat')
    prob_mat <- apply(pred_mat, 2, lind_dens, theta = theta) # Pitr
    probs <- rowMeans(prob_mat) # Pi
    
    ll <- sum(log(probs))
    
    if (est_method == 'bhhh' || method == 'BHHH'){
      return(log(probs))
    } else{return(ll)}
  }

  
  fit <- maxLik::maxLik(p_plind_rp,
                start = start,
                y = y,
                X_Fixed = X_Fixed,
                X_rand = X_rand,
                ndraws = ndraws,
                rpar = rpar,
                correlated = correlated,
                method = method,
                control = list(iterlim = max.iters, printLevel = print.level))
  
  # Optionally, compute bootstrapped standard errors
  # create function to clean data and run maxLik
  plind.boot <- function(data){
    mod1_frame <- stats::model.frame(formula, data)
    X_Fixed <- stats::model.matrix(formula, data)
    y <- stats::model.response(mod1_frame)
    
    int_res <- maxLik::maxLik(p_plind_rp,
                              start = fit$estimate,
                              y = y,
                              X_Fixed = X_Fixed,
                              X_rand = X_rand,
                              ndraws = ndraws,
                              rpar = rpar,
                              correlated = correlated,
                              method = method,
                              control = list(iterlim = max.iters, printLevel = print.level))
    return(int_res)
  }
  
  
  if (!is.null(bootstraps) & is.numeric(bootstraps)) {
    bs.data <- modelr::bootstrap(data, n = bootstraps)
    
    
    mod1_frame <- stats::model.frame(formula, data)
    X_Fixed <- stats::model.matrix(formula, data)
    y <- stats::model.response(mod1_frame)
    
    models <- map(bs.data$strap, ~ plind.boot(data = .))
    tidied <- map_df(models, broom::tidy, .id = "id")
    
    SE <- tidied %>%
      group_by(term) %>%
      summarise(sd = sd(estimate)) %>% deframe()
    
    fit$bootstrapped_se <- SE
  }
  
  fit$bootstraps = if (!is.null(bootstraps)) bootstraps else NULL

  N_fixed = length(x_fixed_names)
  N_rand = length(rpar)
  
  # names(fit$estimate) <- x_names

  param.splits <- as.factor(ifelse((grepl("St. Dev", x_names) + grepl("Cholesky", x_names)==1), "rpr", "coef"))

  coefs <- as.array(fit$estimate)

  # split.coefs <- split(coefs, param.splits)

  betas <- coefs[1:(N_fixed+N_rand)] # includes means of random parameters
  t <- coefs[(N_fixed+N_rand+1):length(coefs)]

  if (correlated){
    chol_vals <-coefs[(N_fixed+N_rand+1):(length(coefs)-1)]
    Cholesky <- matrix(0, N_rand, N_rand)
    counter = 1
    for (i in 1:N_rand){
      for (j in 1:N_rand){
        if (j<=i){
          Cholesky[j,i] <- chol_vals[counter]
          counter <- counter + 1
        }
      }
    }

    Covariance <- t(Cholesky) %*% Cholesky
    Correlation <- cov2cor(Covariance)
  }
  else{
    sd <- abs(coefs[(N_fixed+N_rand+1):(length(coefs)-1)])
    fit$estimate[(N_fixed+N_rand+1):(length(coefs)-1)] <- sd
  }
  
  theta <- exp(tail(fit$estimate,1))
  
  if(correlated){
    fit$Cholesky <- Cholesky
    fit$Covariance <- Covariance
    fit$Correlation <- Correlation
    fit$sd <- as.vector(sqrt(diag(Covariance)))
  }
  else{
    fit$sd <- abs(sd)
  }
  fit$coefs <- betas
  fit$theta = theta
  names(fit$estimate) <- x_names # ensure it retains the names
  fit$formula <- formula
  fit$x_names <- x_names
  fit$rpar_formula <- rpar_formula
  fit$scrambled <- scrambled
  fit$numdraws <- ndraws
  fit$correlated <- correlated
  fit$bootstraps <- NULL
  if (!correlated){
    fit$rpardists = rpardists
    names(fit$sd) <- colnames(X_rand)
  }
  else{
    dst <- as.list(rep("n", length(colnames(X_rand))))
    names(dst) <- colnames(X_rand)
    fit$rpardists = dst
  }

  fit$modelType <- "rppLind"
  
  obj = .createFlexCountReg(model = fit, data = data, call = match.call(), formula = formula)
  return(obj)
}
