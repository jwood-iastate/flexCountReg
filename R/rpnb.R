#' Function for estimating a random parameter negative binomial with the ability to specify if the NB-1, NB-2, or NB-P should be used
#'
#' @name rpnb
#' @param formula an R formula.. This formula should specify the outcome and the independent variables that have fixed parameters.
#' @param rpar_formula a symbolic description of the model related specifically to the random parameters. This should not include an outcome variable. If the intercept is random, include it in this formula. If the intercept is fixed, include it in \code{formula} but not in \code{rpar_formula}. To remove the intercept, use \code{0 + vars} or \code{-1 + vars},
#' @param data a dataframe that has all of the variables in the \code{formula} and \code{rpar_formula},
#' @param form the version of the negative binomial to estimate (\code{"nb2"} estimates the NB-2, \code{"nb1"} estimates the NB-1, \code{"nbp"} estimates the NB-P)
#' @param rpardists an optional named vector whose names are the random parameters and values the distribution. The distribution options include normal ("n"), lognormal ("ln"), triangular ("t"), uniform ("u"), and gamma ("g"). If this is not provided, normal distributions are used for all random coefficients,
#' @param ndraws the number of Halton draws to use for estimating the random parameters,
#' @param scrambled if the Halton draws should be scrambled or not. \code{scrambled = FALSE} results in standard Halton draws while \code{scrambled = FALSE} results in scrambled Halton draws,
#' @param correlated if the random parameters should be correlated (\code{correlated = FALSE} results in uncorrelated random coefficients, \code{correlated = TRUE} results in correlated random coefficients). If the random parameters are correlated, only the normal distribution is used for the random coefficients,
#' @param method a method to use for optimization in the maximum likelihood estimation. For options, see \code{\link[maxLik]{maxLik}},
#' @param max.iters the maximum number of iterations to allow the optimization method to perform,
#' @param start.vals an optional vector of starting values for the regression coefficients
#' @param print.level determines the level of verbosity for printing details of the optimization as it is computed. A value of 0 does not print out any information, a value of 1 prints minimal information, and a value of 2 prints the most information.
#' @import randtoolbox maxLik stats modelr
#' @importFrom MASS glm.nb
#' @importFrom utils head  tail
#' @include tri.R
#'
#' @examples
#' \donttest{
#'
#' ## Random Parameters Negative Binomial model (NB-1)
#' data("washington_roads")
#' nb1.rp <- rpnb(Total_crashes ~ - 1 + lnlength + lnaadt,
#'                rpar_formula = ~ speed50,
#'                data = washington_roads,
#'                ndraws = 100,
#'                correlated = FALSE,
#'                rpardists = c(intercept="u", speed50="t"),
#'                form = 'nb1',
#'                method = "bfgs",
#'                print.level = 2)
#'
#' summary(nb1.rp)
#'
#' ## Random Parameters Negative Binomial model (NB-2)
#' nb2.rp <- rpnb(Total_crashes ~ - 1 + lnlength + lnaadt,
#'                rpar_formula = ~ speed50,
#'                data = washington_roads,
#'                ndraws = 100,
#'                correlated = TRUE,
#'                form = 'nb2',
#'                method = "bfgs",
#'                print.level = 1)
#'
#' summary(nb2.rp)
#'
#' ## Random Parameters Negative Binomial model (NB-P)
#' nbp.rp <- rpnb(Total_crashes ~ - 1 + lnlength + lnaadt,
#'                rpar_formula = ~ speed50,
#'                data = washington_roads,
#'                ndraws = 100,
#'                correlated = TRUE,
#'                form = 'nbp',
#'                method = "bfgs",
#'                print.level = 1)
#'
#' summary(nbp.rp)}
#' @export
rpnb <- function(formula, rpar_formula, data, form = 'nb2',
                 rpardists = NULL,
                 ndraws = 1500, scrambled = FALSE,
                 correlated = FALSE, method = 'BHHH', max.iters = 1000,
                 start.vals = NULL, print.level = 0) {
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

  # function to compute probabilities
  nb_prob <- function(y, mu, alpha, p) {
    if (form=='nb2'){
      return(stats::dnbinom(y, size = alpha, mu = mu))
    } else if (form=='nb1'){
      return(stats::dnbinom(y, size = mu/alpha, mu = mu))
    } else{
      return(stats::dnbinom(y, size = (mu^(2-p))/alpha, mu = mu))
    }
  }

  # Generating model matrices
  mod1_frame <- stats::model.frame(formula, data)
  X_Fixed <- modelr::model_matrix(data, formula)
  X_rand <- modelr::model_matrix(data, rpar_formula)
  #X_rand <-stats::model.matrix(rpar_formula, data)
  y_name <- all.vars(formula)[1]
  y <- stats::model.response(mod1_frame)
  
  X_expanded <- cbind(X_Fixed, X_rand, X_rand) # for use in the gradient

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
    start <- append(start, log(0.1)) # initial log of overdispersion parameter
    x_names <- append(x_names, 'ln(alpha)')

    if (form=='nbp'){
      start <- append(start, 1.5) # initial value for parameter P
      x_names <- append(x_names, 'P')
    }
    names(start) <- x_names
  }
  
  
  # main function for estimating log-likelihoods
  p_nb_rp <- function(p, y, X_Fixed, X_rand, ndraws, rpar, correlated){
    est_method <- method
    if (!correlated) exact.gradient=TRUE else exact.gradient=FALSE # use numerical gradient if using correlated random parameters
    N_fixed = length(x_fixed_names)
    N_rand = length(rpar)
    coefs <- as.array(p)
    fixed_coefs <- head(coefs,N_fixed)
    h <- head(coefs, (N_fixed + N_rand))
    random_coefs_means <- tail(h, N_rand)

    if (form=='nbp'){
      dist_params <- tail(coefs,2)
      log_alpha <- dist_params[1]
      p <- dist_params[2]
    } else{
      dist_params <- tail(coefs,1)
      log_alpha <- dist_params[1]
      p <- NULL
    }

    alpha <- exp(log_alpha)
    mu_fixed <- exp(X_Fixed %*% fixed_coefs)

    if(length(rpar)==1){
      rand_sdevs <- coefs[length(coefs)-1]
    }else{
      numtail <- length(coefs) - N_fixed - N_rand
      t <- tail(coefs, numtail)
      if (form=='nbp') t <- head(t,-2)
      else t <- head(t,-1)
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
              draws[,i] <- qtri(hdraws[,i], random_coefs_means[i], abs(rand_sdevs[i]))
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
    prob_mat <- apply(pred_mat, 2, nb_prob, y = y, alpha = alpha, p=p) # Pitr
    probs <- rowMeans(prob_mat) # Pi
    
    ll <- sum(log(probs))
    
    if (est_method == 'bhhh' || method == 'BHHH'){
      return(log(probs))
    } else{return(ll)}
  }
  
  gradFun <- function(p, y, X_Fixed, X_rand, ndraws, rpar, correlated) {
    params <- p
    N_fixed <- ncol(X_Fixed)
    N_rand <- ncol(X_rand)
    
    fixed_coefs <- params[1:N_fixed]
    
    if (N_rand > 1) {
      rand_means <- params[(N_fixed + 1):(N_fixed + N_rand)]
      rand_sdevs <- params[(N_fixed + N_rand + 1):(N_fixed + 2 * N_rand)]
    } else {
      rand_means <- params[(N_fixed + 1)]
      rand_sdevs <- params[(N_fixed + 2)]
    }
    
    log_alpha <- params[length(params) - ifelse(form == 'nbp', 2, 1)]
    alpha <- exp(log_alpha)
    
    mu_fixed <- exp(X_Fixed %*% fixed_coefs)
    
    if (correlated && N_rand > 1) {
      chol_vals <- p[(N_fixed + N_rand + 1):(length(p) - ifelse(form == 'nbp', 2, 1))]
      Ch <- matrix(0, N_rand, N_rand)
      Ch[lower.tri(Ch, diag = TRUE)] <- chol_vals
      scaled_draws <- qnorm(hdraws) %*% Ch
    } else if (N_rand > 1) {
      rpardists <- rep("n", length(rpar))
      names(rpardists) <- rpar
      scaled_draws = matrix(0, nrow = nrow(hdraws), ncol = ncol(hdraws))
      for (i in 1:N_rand) {
        dist_type <- rpardists[rpar[i]]
        scaled_draws[, i] <- switch(dist_type,
                                    ln = stats::qlnorm(hdraws[, i], rand_means[i], abs(rand_sdevs[i])),
                                    t = qtri(hdraws[, i], rand_means[i], abs(rand_sdevs[i])),
                                    u = rand_means[i] + (hdraws[, i] - 0.5) * abs(rand_sdevs[i]),
                                    g = stats::qgamma(hdraws[, i], shape = rand_means[i]^2 / (rand_sdevs[i]^2), rate = rand_means[i] / (rand_sdevs[i]^2)),
                                    n = stats::qnorm(hdraws[, i], rand_means[i], abs(rand_sdevs[i])))
      }
    } else {
      dist_type <- rpardists <- "n"
      names(rpardists) <- rpar
      
      scaled_draws <- switch(dist_type,
                             ln = stats::qlnorm(hdraws, rand_means, abs(rand_sdevs)),
                             t = qtri(hdraws, rand_means, abs(rand_sdevs)),
                             u = rand_means + (hdraws - 0.5) * abs(rand_sdevs),
                             g = stats::qgamma(hdraws, shape = rand_means^2 / (rand_sdevs^2), rate = rand_means / (rand_sdevs^2)),
                             n = stats::qnorm(hdraws, rand_means, abs(rand_sdevs)))
    }
    
    xb_rand_mat <- crossprod(t(X_rand), t(scaled_draws))
    
    rpar_mat <- exp(xb_rand_mat)
    pred_mat <- apply(rpar_mat, 2, function(x) x * mu_fixed)
    
    mu <- rowMeans(pred_mat)
    
    grad <- matrix(0, nrow = nrow(X_rand), ncol = length(params))
    
    if (form == 'nb1') {
      p_val = 1
    } else if (form == 'nb2') {
      p_val = 2
    } else if (form == 'nbp') {
      p_val <- params[length(params)]
    }
    
    R <- mu^(2 - p_val) / alpha
    S_sum <- y + R
    M_sum <- mu + R
    
    if (form != "nb2") {
      ln_alpha_grad_cons <- R * (alpha * y - alpha * mu + alpha * R * (log(M_sum) - log(R) - digamma(S_sum) + digamma(R))) / (alpha * M_sum)
      beta_grad_cons <- R * (alpha * y * (p_val - 1) - alpha * mu + (p_val - 2) * (alpha * M_sum * (log(M_sum) - log(R) - digamma(S_sum) + digamma(R)) + mu^(2 - p_val))) / (alpha * M_sum)
    } else {
      ln_alpha_grad_cons <- ((y - mu) + (alpha * mu + 1) * (log(alpha * mu + 1) - digamma(y + 1 / alpha) + digamma(1 / alpha))) / (mu + 1 / alpha)
      beta_grad_cons <- (y - mu) / (alpha * mu + 1)
    }
    
    # Gradient for fixed coefficients
    g_beta <- apply(X_Fixed, 2, function(x) x * beta_grad_cons)
    
    # Start of gradient for random coefficients
    g_beta_r <- apply(X_rand, 2, function(x) x * beta_grad_cons)
    
    ## Get individual-specific coefficients
    n_obs <- nrow(X_rand)
    num_vars_rand <- ncol(X_rand)
    ind_coefs <- matrix(0, nrow = n_obs, ncol = num_vars_rand)
    prob_mat <- apply(pred_mat, 2, function(mu_i) {
      if (form == 'nb2') {
        return(stats::dnbinom(y, size = 1 / alpha, mu = mu_i))
      } else {
        return(stats::dpois(y, lambda = mu_i))
      }
    })
    
    draws <- scaled_draws
    
    if (num_vars_rand > 1) {
      ind_coefs_pred <- matrix(0, nrow = n_obs, ncol = num_vars_rand)
      for (i in 1:num_vars_rand) {
        hals <- draws[, i]
        b_i <- as.vector(crossprod(t(prob_mat), hals))
        ind_coefs[, i] <- b_i / rowSums(prob_mat)
      }
    } else {
      hals <- draws
      b_i <- as.vector(crossprod(t(prob_mat), hals))
      ind_coefs <- b_i / rowSums(prob_mat)
    }
    
    d_dist <- function(dist, mu, sigma, b) {
      if (dist == "t") {
        if (b <= mu) {
          d_mu <- -1 / sigma^2      
          d_sigma <- 2 * (b - mu + sigma) / sigma^2  
        } else {
          d_mu <- 1 / sigma^2      
          d_sigma <- -2 * (mu + sigma - b) / sigma^2
        }
      } else if (dist == "u") {
        d_mu <- 0    
        d_sigma <- -1 / sigma^2
      } else if (dist == "g") {
        shape <- mu^2 / sigma^2
        rate <- mu / sigma^2
        d_mu <- (rate - log(b) / sigma^2) - 2 * (rate * (log(rate) - log(b) - digamma(shape)))    
        d_sigma <- shape / sigma * (log(b) - log(rate) + digamma(shape)) - 2 * (rate / sigma * (mu - b))
      } else if (dist == "ln") {
        dens <- (1 / (b * sigma * sqrt(2 * pi))) * exp(-(log(b) - mu)^2 / (2 * sigma^2))
        d_mu <- dens * (log(b) - mu) / sigma^2    
        d_sigma <- dens * (-1 / sigma + (log(b) - mu)^2 / sigma^3)
      } else if (dist == "n") {
        dens <- dnorm(b, abs(mu), abs(sigma))
        dens <- ifelse(is.na(dens), 1e-200,dens)
        d_mu <- dens * (b - mu) / sigma^2    
        d_sigma <- dens * (-1 / sigma + (b - mu)^2 / sigma^3)
      }
      return(cbind(d_mu, d_sigma))
    }
    
    # Compute gradients for the random parameters using the chain rule
    
    g_beta_r_mu <- matrix(0, ncol = num_vars_rand, nrow = n_obs)
    g_beta_r_sd <- matrix(0, ncol = num_vars_rand, nrow = n_obs)
    
    if (num_vars_rand>1){
      for (i in 1:num_vars_rand) {
        grads_rp <- d_dist(dist = rpardists[i], mu = rand_means[i], sigma = rand_sdevs[i], b = ind_coefs[, i])
        g_beta_r_mu[, i] <- grads_rp[, 1] * g_beta_r[, i]
        g_beta_r_sd[, i] <- grads_rp[, 2] * g_beta_r[, i]
      }
    }
    else{
      grads_rp <- d_dist(dist = rpardists, mu = rand_means, sigma = rand_sdevs, b = ind_coefs)
      g_beta_r_mu <- grads_rp[, 1] * g_beta_r
      g_beta_r_sd <- grads_rp[, 2] * g_beta_r
    }
    
    
    # Gradient for p if form is nbp
    if (form == 'nbp') {
      grad <- cbind(g_beta, g_beta_r_mu, g_beta_r_sd, ln_alpha_grad_cons, p_val * log(mu) * ln_alpha_grad_cons)
    } else {
      grad <- cbind(g_beta, g_beta_r_mu, g_beta_r_sd, ln_alpha_grad_cons)
    }
    return(grad)
  }
  
  
  fit <- maxLik::maxLik(p_nb_rp,
                start = start,
                y = y,
                X_Fixed = X_Fixed,
                X_rand = X_rand,
                ndraws = ndraws,
                rpar = rpar,
                correlated = correlated,
                # grad = if (correlated) NULL else if (method == 'BHHH') gradFun else function(...) colSums(gradFun(...)),
                method = method,
                control = list(iterlim = max.iters, printLevel = print.level))

  N_fixed = length(x_fixed_names)
  N_rand = length(rpar)
  
  # names(fit$estimate) <- x_names

  param.splits <- as.factor(ifelse((grepl("St. Dev", x_names) + grepl("Cholesky", x_names)==1), "rpr", "coef"))

  coefs <- as.array(fit$estimate)

  # split.coefs <- split(coefs, param.splits)

  betas <- coefs[1:(N_fixed+N_rand)] # includes means of random parameters
  t <- coefs[(N_fixed+N_rand+1):length(coefs)]

  if (correlated){
    if(form=='nbp'){
      chol_vals <- coefs[(N_fixed+N_rand+1):(length(coefs)-2)]
    }
    else{
      chol_vals <-coefs[(N_fixed+N_rand+1):(length(coefs)-1)]
    }
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
    if(form=='nbp'){
      sd <- abs(coefs[(N_fixed+N_rand+1):(length(coefs)-2)])
      fit$estimate[(N_fixed+N_rand+1):(length(coefs)-2)] <- sd
    }
    else{
      sd <- abs(coefs[(N_fixed+N_rand+1):(length(coefs)-1)])
      fit$estimate[(N_fixed+N_rand+1):(length(coefs)-1)] <- sd
    }
  }
  

  if(form =='nbp'){
    t <- tail(coefs,2)
    alpha <- exp(t[1])
    P <- t[2]
  }
  else{
    alpha <- exp(tail(fit$estimate,1))
    P <- NULL
  }
  #names(fit$estimate) <- x_names
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
  fit$alpha = alpha
  fit$P <- P
  #fit$estimate <- ifelse(grepl("St.", names(fit$estimate)), abs(fit$estimate ), fit$estimate)
  names(fit$estimate) <- x_names # ensure it retains the names
  fit$formula <- formula
  fit$x_names <- x_names
  fit$rpar_formula <- rpar_formula
  fit$scrambled <- scrambled
  fit$numdraws <- ndraws
  fit$correlated <- correlated
  fit$bootstraps <- NULL
  fit$form = form
  if (!correlated){
    fit$rpardists = rpardists
    names(fit$sd) <- colnames(X_rand)
  }
  else{
    dst <- as.list(rep("n", length(colnames(X_rand))))
    names(dst) <- colnames(X_rand)
    fit$rpardists = dst
  }

  fit$modelType <- "rpnb"
  
  obj = .createFlexCountReg(model = fit, data = data, call = match.call(), formula = formula)
  return(obj)
}
