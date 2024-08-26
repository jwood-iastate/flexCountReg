#' Predictions for non-random parameters count models
#'
#' @name predict.flexCountReg
#' @param object a model object estimated using this R package.
#' @param ... optional arguments passed to the function. These include `data` and `method`.
#' 
#' @note optional parameter `data`: a dataframe that has all of the variables in the \code{formula} and \code{rpar_formula}. This can be the data used for estimating the model or another dataframe.
#' @note optional parameter `method`: Only valid for random parameters models. The method to be used in generating the predictions (options include \code{Simulated}, \code{Exact}, or \code{Individual}).
#' @note the method option \code{Individual} requires that the outcome be observed for all observations. This is due to the use of Bayesian methods for computing the individual observation coefficients.
#'
#' @import randtoolbox stats lamW modelr rlang 
#' @importFrom utils head  tail
#' @include corr_haltons.R halton_dists.R get_chol.R
#'
#' @examples
#'
#' # Standard Negative Binomial
#' data("washington_roads")
#' nbg <- nbg(Total_crashes ~ lnaadt + lnlength + speed50,
#'            ln.alpha.formula = ~ lnaadt,
#'            data = washington_roads)
#' summary(nbg)
#'
#' predict(nbg, data=washington_roads)
#' @export
predict.flexCountReg <- function(object, ...){
  # Extract optional parameters from '...'
  args <- list(...)
  if (is.null(args$data)) data <- object$data else data <- args$data
  if (is.null(args$method)) method <- "Exact" else data <- args$method

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
    formula <- model$formula
    
    mod1_frame <- stats::model.frame(formula, data)
    X_Fixed <- as.matrix(modelr::model_matrix(data, formula))
    X_rand <- as.matrix(modelr::model_matrix(data, rpar_formula))
    
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
      rand_pred<- rep(1, length(data))
      
      if (length(dists)>1){
        for (i in 1:length(dists)){
          
          rand_pred <- rand_pred*rpar.adjust(dists[i], random_coefs_means[i], sd[i], X_rand[,i])
        }
      }
      else{
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
      
      
      y <- model.response(mod1_frame)
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
    print(dim(data))
    #mod_df <- stats::model.frame(model$formula, data)
    
    X <- as.matrix(modelr::model_matrix(data, model$formula))
    #y <- as.numeric(stats::model.response(mod_df))
    
    
    beta_pred <- model$beta_pred
    
    if (modtype=="posLogn"){
      if (!is.null(model$ln.sigma.formula)){
        sigma_X <- stats::model.matrix(model$ln.sigma.formula, data)
        sigma_pars <- model$sigma_coefs
        sigma <- exp(sigma_X %*% sigma_pars)
        predictions <- exp(X %*% beta_pred + (sigma^2)/2)
      }
      else{
        predictions <- exp(X %*% beta_pred + (model$sigma^2)/2)
      }
    }
    else{
      predictions <- exp(X %*% beta_pred)
    }
    return(c(predictions))
  }
}
