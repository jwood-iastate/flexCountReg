#' Predictions for non-random parameters count models
#'
#' @name predict.flexCountReg
#' @param object a model object estimated using this R package.
#' @param ... optional arguments passed to the function. These include `data` 
#'        and `method`.
#' 
#' @note optional parameter `data`: a dataframe that has all of the variables in 
#'        the \code{formula} and \code{rpar_formula}. This can be the data used 
#'        for estimating the model or another dataframe.
#' @note optional parameter `method`: Only valid for random parameters models. 
#'        The method to be used in generating the predictions (options include 
#'        \code{Simulated}, \code{Exact}, or \code{Individual}).
#' @note the method option \code{Individual} requires that the outcome be 
#'        observed for all observations. This is due to the use of Bayesian 
#'        methods for computing the individual observation coefficients.
#'
#' @import randtoolbox stats lamW modelr rlang 
#' @importFrom utils head  tail
#' @include corr_haltons.R halton_dists.R get_chol.R
#'
#' @examples
#'
#' # Standard Negative Binomial
#' data("washington_roads")
#' washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>=10000,1,0)
#' nb2 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
#'                data = washington_roads, family = "NB2",
#'                dis_param_formula_1 = ~ speed50, method='BFGS')
#'
#' predict(nb2, data=washington_roads)
#' 
#' # Poisson-Lognormal
#' pln <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
#'               data = washington_roads, family = "PLN", ndraws=10)
#' predict(pln, data=washington_roads)
#' 
#' # Random Parameter NB2
#' nb2.rp <- rpnb(Total_crashes ~ - 1 + lnlength + lnaadt,
#'                               rpar_formula = ~ speed50,
#'                               data = washington_roads,
#'                               ndraws = 100,
#'                               correlated = FALSE,
#'                               rpardists = c(intercept="u", speed50="t"),
#'                               form = 'nb2',
#'                               method = "bfgs")
#'                
#' predict(nb2.rp, list(data=washington_roads, method="Simulated"))
#' 
#' # NB2 with underreporting (logit)
#' nb2_underreport <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
#'               data = washington_roads, family = "NB2",
#'               underreport_formula = ~ speed50 + AADT10kplus)
#' predict(nb2_underreport, data=washington_roads)
#' 
#' # Poisson-Lognormal with underreporting (probit)
#' plogn_underreport <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
#'               data = washington_roads, family = "NB2",
#'               underreport_formula = ~ speed50 + AADT10kplus, underreport_family = "probit")
#' predict(plogn_underreport, data=washington_roads)
#' 
#' @export
predict.flexCountReg <- function(object, ...){
  # Extract optional parameters from '...'
  additional_args  <- list(...)
  if (is.null(additional_args$data)) data <- object$data else data <- as.data.frame(additional_args$data)
  if (is.null(additional_args$method)) method <- "Exact" else data <- additional_args$method

  model <- object$model
  
  modtype <- model$modelType 
  if(modtype=="rpnb" || modtype=="rppois"){
    # The Individual method uses a Bayesian approach to get individual coefficients and variance of the coefficients
    
    form <- model$form
    
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
    
    rpar_formula <- model$rpar_formula
    rpardists <- model$rpardists
    formula <- delete.response(model$formula)
    data <- as.data.frame(data)
    
    print(formula)
    print(head(data))
    
    
    X_Fixed <- as.matrix(modelr::model_matrix(data, formula))
    X_rand <- as.matrix(modelr::model_matrix(data, model$rpar_formula))
    
    X <- cbind(X_Fixed, X_rand)
    
    x_fixed_names <- colnames(X_Fixed)
    rpar <- colnames(X_rand)
    
    N_fixed <- num_vars_fixed <- length(x_fixed_names)
    Nrand <- N_rand <- num_vars_rand <- length(unname(rpar))
    total_vars <- num_vars_fixed + num_vars_rand
    
    coefs <- unlist(model$estimate, recursive = TRUE, use.names = FALSE)
    fixed_coefs <- head(coefs,num_vars_fixed)
    h <- head(coefs, total_vars)
    random_coefs_means <- tail(h, num_vars_rand)
    
    
    if (N_rand > 1) {
      rand_sdevs <- abs(coefs[(total_vars+ 1):(total_vars+N_rand)])
    } 
    else {
      rand_sdevs <- abs(coefs[(N_fixed + 2)])
    }
    
    sd <- rpar_sd <- rand_sdevs

    if (modtype=="rpnb" ){
      alpha <- model$alpha
      p <- model$P
    }
    
    mu_fixed <- exp(X_Fixed %*% fixed_coefs)
    correlated <- model$correlated
    scrambled <- model$scrambled
    ndraws <- max(model$numdraws,2000)
    dists <- model$rpardists
    
    hdraws <- randtoolbox::halton(ndraws, num_vars_rand, mixed = scrambled)
    
    # function to adjust for random distributions
    rpar.adjust <- function(dist, mu, sigma, xrand){ 
      if(dist=="n"){
        adj <- exp(xrand*mu + xrand^2*sigma^2/2)
        return(adj)
      }
      else if (dist=="ln"){
        if (sigma<=sqrt(exp(-mu-1))){
          W <- lamW::lambertW0(-xrand*sigma^2*exp(mu))
          W <- ifelse(is.na(W), lamW::lambertWm1(-xrand*sigma^2*exp(mu)), W)
          adj <- ifelse(xrand*exp(mu-W)-1/(sigma^2)==0,1,(exp(xrand*exp(mu-W))-W^2/(2*sigma^2))/(sigma*sqrt(abs(xrand*exp(mu-W)-1/(sigma^2)))))
        }
        else{
          draws <- stats::qlnorm(randtoolbox::halton(ndraws, 1), random_coefs_means, rand_sdevs)
          xs <- exp(crossprod(xrand, draws))
          adj <- rowMeans(xs)
        }
        return(adj)
      }
      else if (dist=="t"){
        adj <- ifelse(xrand!=0, (exp(2*sigma*xrand)-2*exp(sigma*xrand)+1)*exp(xrand*(mu-sigma))/(sigma^2*xrand^2), 1)
        return(adj)
      }
      else if (dist=="u"){
        adj <- ifelse(xrand!=0,(-exp(xrand*(mu-sigma))+exp(xrand*(mu+sigma)))/(2*xrand*sigma),1)
        return(adj)
      }
      else if (dist=="g"){
        adj <- ifelse(xrand*sigma^2/(mu)==1,1,((mu-sigma^2*xrand)/mu)^(-1*(mu^2/(sigma^2))))
        return(adj)
      }
    }
    
    if (method == 'Exact'){
      rand_pred<- rep(1, length(mu_fixed))
      
      if (length(dists)>1){
        for (i in 1:length(dists)){
          
          #print(paste("N Rand Perd =", length(rand_pred), ", dists[i] =", dists[i], ", Coef Mean =", random_coefs_means[i], ", Coef SD =", sd[i], ", N X[,i] =", length(X_rand[,i])) )
          
          rand_pred <- rand_pred*rpar.adjust(dists[i], random_coefs_means[i], sd[i], X_rand[,i])
        }
      }
      else{
        print(paste(dists, random_coefs_means, sd, length(X_rand)))
        rand_pred <- rpar.adjust(dists, random_coefs_means, sd, X_rand)
      }
      
      fixed_predictions <- exp(X_Fixed %*% fixed_coefs)
      predictions <- fixed_predictions*rand_pred
      return(as.vector(predictions))
    }
    else if (method=='Simulated'){
      
      # generate and scale random draws
      if(length(rpar)==1){
        if(is.null(rpardists)){
          draws <- halton_dists(dist="n", hdraw=hdraws, mean=random_coefs_means, sdev=rand_sdevs)
        }
        else{
          draws <- halton_dists(dist=rpardists, hdraw=hdraws, mean=random_coefs_means, sdev=rand_sdevs)
        }
        xb_rand_mat <- sapply(draws, function(x) X_rand * x)
      }else{
        if (correlated){ # Generate correlated random draws
          Ch <- get_chol(t, Nrand)
          draws <- corr_haltons(random_coefs_means, cholesky = Ch, hdraws=hdraws)
        }else{
          draws <- hdraws # initialize the matrix
          for (i in 1:Nrand){
            draws[,i] <- halton_dists(dist=rpardists[i], hdraw=hdraws[,i], mean=random_coefs_means[i], sdev=rand_sdevs[i])
          }
        }
        draws <- t(draws)
        xb_rand_mat <- crossprod(t(X_rand), draws)
      }
 
      rpar_mat <- exp(xb_rand_mat)
      pred_mat <- apply(rpar_mat, 2, function(x) x * mu_fixed)
      mui <- rowMeans(pred_mat)
      
      return(mui)
    }
    else if (method=='Individual'){
      # generate and scale random draws
      if(length(rpar)==1){
        if(is.null(rpardists)){
          draws <- halton_dists(dist="n", hdraw=hdraws, mean=random_coefs_means, sdev=rand_sdevs)
        }
        else{
          draws <- halton_dists(dist=rpardists, hdraw=hdraws, mean=random_coefs_means, sdev=rand_sdevs)
        }
        xb_rand_mat <- sapply(draws, function(x) X_rand * x)
      }else{
        if (correlated){ # Generate correlated random draws
          Ch <- get_chol(t, Nrand)
          draws <- corr_haltons(random_coefs_means, cholesky = Ch, hdraws=hdraws)
        }else{
          draws <- hdraws # initialize the matrix
          for (i in 1:Nrand){
            draws[,i] <- halton_dists(dist=rpardists[i], hdraw=hdraws[,i], mean=random_coefs_means[i], sdev=rand_sdevs[i])
          }
        }
        draws <- t(draws)
        xb_rand_mat <- crossprod(t(X_rand), draws)
      }
      
      rpar_mat <- exp(xb_rand_mat)
      
      pred_mat <- apply(rpar_mat, 2, function(x) x * mu_fixed)
      
      variable_names <- all.vars(model$formula)
      y <- data[,variable_names[1]]
      n_obs <- length(mu_fixed)
      ind_coefs <- matrix(0, nrow=n_obs, ncol=num_vars_rand)
      ind_var <- matrix(0, nrow=n_obs, ncol=num_vars_rand)
      if (modtype=="rpnb"){
        prob_mat <- apply(pred_mat, 2, nb_prob, y=y, alpha=alpha, p=p)
      }
      else{
        prob_mat <- apply(pred_mat, 2, dpois, x=y)
      }
      
      
      pred_i <- rep(0,n_obs)
      
      if(num_vars_rand>1){
        ind_coefs_pred <- matrix(0, nrow=n_obs, ncol=num_vars_rand)
        for (i in 1:num_vars_rand){
          hals <- draws[i,]
          b_i <- t(apply(prob_mat, 1, function(x) hals * x))
          bb_i <- t(apply(prob_mat, 1, function(x) hals * hals * x))
          bi <- rowSums(b_i)/rowSums(prob_mat)
          bbi <- rowSums(bb_i)/rowSums(prob_mat)
          var <- bbi - bi^2
          ind_coefs[,i] <- bi
          ind_var[,i] <- var
          ind_coefs_pred[,i] <- bi + var/2
        }
      }
      
      xr <- exp(rowSums(ind_coefs_pred * X_rand))
      pred_i <- diag(outer(xr, as.vector(mu_fixed)))
      return(pred_i)
    }
    
    else{print('Please use one of the following methods: Approximate, Simulated, or Individual')}
  }
  else{
    #mod_df <- stats::model.frame(model$formula, data)
    
    X <- as.matrix(modelr::model_matrix(data, model$formula))
    #y <- as.numeric(stats::model.response(mod_df))
    coefs <- unlist(model$estimate, recursive = TRUE, use.names = FALSE)
    N_x <- ncol(X)
    N_pars <- length(coefs)
    
    
    
    beta_pred <- as.vector(coefs[1:ncol(X)])
    
    if (!is.null(model$dis_param_formula_1)){
      alpha_X <- as.matrix(modelr::model_matrix(data, model$dis_param_formula_1))
      N_alpha <- ncol(alpha_X)
      alpha_pars <- coefs[(N_x+1):(N_x+N_alpha)]
      alpha <- exp(alpha_X %*% alpha_pars)
    }else{
      N_alpha <- 1
      alpha <- exp(coefs[N_x+1])
    }
    
    if(!is.null(model$underreport_formula)){
      X_underreport <- as.matrix(modelr::model_matrix(data, model$underreport_formula))
      N_underreport <- ncol(X_underreport)
      underrep_pars <- coefs[(N_pars-N_underreport+1):N_pars]
      underrep_lin <- X_underreport%*%underrep_pars
      
      if(model$underreport_family=="logit"){
        mu_adj <- 1/(1+exp(-underrep_lin))
      }
      else{
        mu_adj <- pnorm(underrep_lin, lower.tail = FALSE)
        
      }
    }else{mu_adj=1} # If no underreporting model, set the probability to 1
    
    if (modtype=="countreg" & model$family=="PLN"){
      predictions <- exp(X %*% beta_pred + (alpha^2)/2)*mu_adj
    }
    else{
      predictions <- exp(X %*% beta_pred)*mu_adj
    }
    return(c(predictions))
  }
}
