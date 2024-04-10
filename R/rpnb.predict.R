#' Function for generating predictions based on the random parameters negative binomial with multiple optional methods
#'
#' @name rpnb.predict
#' @param model a model object estimated using the \code{rpnb} function,
#' @param data a dataframe that has all of the variables in the \code{formula} and \code{rpar_formula}. This can be the data used for estimating the model or another dataframe,
#' @param method the method to be used in generating the predictions (options include \code{Simulated}, \code{Exact}, or \code{Individual}).
#' @note the method option \code{Individual} requires that the outcome be observed for all observations. This is due to the use of Bayesian methods for computing the individual observation coefficients.
#'
#' @import nlme randtoolbox stats lamW modelr
#' @importFrom utils head  tail
#' @include tri.R
#'
#' @examples
#' \donttest{
#'
#' ## Random Parameters Negative Binomial model (NB-2)
#'
#' data("washington_roads")
#'
#' nb2.rp <- rpnb(Total_crashes ~ - 1 + lnlength + lnaadt,
#'                rpar_formula = ~ speed50,
#'                data = washington_roads,
#'                ndraws = 100,
#'                correlated = TRUE,
#'                form = 'nb2',
#'                method = "bfgs",
#'                print.level = 1)
#'
#' ## Exact Prediction
#' hist(rpnb.predict(nb2.rp, washington_roads))
#'
#' ## Simulated Prediction
#' hist(rpnb.predict(nb2.rp, washington_roads, method="Simulated"))
#'
#' ## Exact Prediction
#' hist(rpnb.predict(nb2.rp, washington_roads, method="Individual"))
#' }
#' @export
rpnb.predict <- function(model, data, method='Exact') {
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

  num_vars_fixed <- length(x_fixed_names)
  N_rand <- num_vars_rand <- length(unname(rpar))
  total_vars <- num_vars_fixed + num_vars_rand

  coefs <- unlist(model$estimate, recursive = TRUE, use.names = FALSE)
  fixed_coefs <- head(coefs,num_vars_fixed)
  h <- head(coefs, total_vars)
  random_coefs_means <- tail(h, num_vars_rand)
  alpha <- model$alpha
  r <- 1/alpha # inverse of alpha
  mu_fixed <- exp(X_Fixed %*% fixed_coefs)
  correlated <- model$correlated
  scrambled <- model$scrambled
  ndraws <- max(model$numdraws,2000)
  p <- model$P
  dists <- model$rpardists

  hdraws <- randtoolbox::halton(ndraws, num_vars_rand, mixed = scrambled)

  # function to adjust for random distributions
  rpar.adjust <- function(dist, mu, sigma){ # adjusted coefficients for random parameters
    if(dist=="n"){
      adj <- mu + sigma^2/2
      return(adj)
    }
    else if (dist=="ln"){
      if (sigma<=sqrt(exp(-mu-1))){
        W <- lamW::lambertW0(-sigma^2*exp(mu))
        W <- ifelse(is.na(W), lamW::lambertWm1(-sigma^2*exp(mu)), W)
        adj <- exp(mu-W)-W^2/(2*sigma^2)-log(sigma)-0.5*log(abs(exp(mu-W)-1/(sigma^2)))
      }
      else{
        h <- randtoolbox::halton(1)
        adj <- log(mean(exp((stats::qlnorm(h, mu, abs(sigma))))))
      }
      
      return(adj)
    }
    else if (dist=="t"){
      adj <- mu-2*log(sigma)+log(2*cosh(sigma)-2)
      return(adj)
    }
    else if (dist=="u"){
      adj <- mu-log(sigma)+log(sinh(sigma))
      return(adj)
    }
  }

  if(length(rpar)==1){
    rpar_sd = rand_sdevs <- coefs[length(coefs)]
  }else{
    numtail <- length(coefs) - total_vars
    t <- tail(coefs, numtail)
  }

  if (method == 'Exact'){
    sd <- abs(model$sd)
    rand_coefs <- rep(0, length(random_coefs_means))

    for (i in 1:length(dists)){
      rand_coefs[i] <- rpar.adjust(dists[rpar[i]], random_coefs_means[i], sd[rpar[i]])
    }

    X <- cbind(X_Fixed, X_rand)

    betas <- c(fixed_coefs, rand_coefs)
    predictions <- exp(X %*% betas)
    return(as.vector(predictions))
  }
  else if (method=='Simulated'){
    # generate and scale random draws
    if (length(rpar)>1){
      if (correlated){ # Generate correlated random draws
        Ch <- model$Cholesky
        scaled_draws <- qnorm(hdraws) %*% Ch
        draws <- apply(scaled_draws, 1, function(x) x + random_coefs_means)
      }
      else{
        rand_sdevs <- t
        draws <- hdraws #initialize the matrix

        if (is.null(rpardists)){
          print(random_coefs_means)
          print(rand_sdevs)
          draws <- apply(hdraws, 1, function(x) stats::qnorm(x, random_coefs_means, rand_sdevs))
        }
        else{
          for (i in rpar){
            counter=1
            if (rpardists[i]=="ln"){
              draws[,counter] <- stats::qlnorm(hdraws[,counter], random_coefs_means[counter], abs(rand_sdevs[counter]))
              counter=counter+1
            }
            else if (rpardists[i]=="t"){
              draws[,counter] <- qtri(hdraws[,counter], random_coefs_means[counter], abs(rand_sdevs[counter]))
              counter=counter+1
            }
            else if (rpardists[i]=="u"){
              draws[,counter] <- random_coefs_means[counter] + (hdraws[,counter] - 0.5)*abs(rand_sdevs[counter])
              counter=counter+1
            }
            else{
              draws[,counter] <- stats::qnorm(hdraws[,counter], random_coefs_means[counter], abs(rand_sdevs[counter]))
              counter=counter+1
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
        else{
          draws <- stats::qnorm(hdraws, random_coefs_means, abs(rand_sdevs))
        }

      }
      xb_rand_mat <- sapply(draws, function(x) X_rand * x)
    }

    rpar_mat <- exp(xb_rand_mat)

    # Efficient computation of predicted matrix
    pred_mat <- apply(rpar_mat, 2, function(x) x * mu_fixed)

    mui = rowMeans(pred_mat)
    return(mui)
  }
  else if (method=='Individual'){
    # generate and scale random draws
    if (length(rpar)>1){
      if (correlated){ # Generate correlated random draws
        Ch <- model$Cholesky
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
          for (i in rpar){
            counter=1
            if (rpardists[i]=="ln"){
              draws[,counter] <- stats::qlnorm(hdraws[,counter], random_coefs_means[counter], abs(rand_sdevs[counter]))
              counter=counter+1
            }
            else if (rpardists[i]=="t"){
              draws[,counter] <- qtri(hdraws[,counter], random_coefs_means[counter], abs(rand_sdevs[counter]))
              counter=counter+1
            }
            else if (rpardists[i]=="u"){
              draws[,counter] <- random_coefs_means[counter] + (hdraws[,counter] - 0.5)*abs(rand_sdevs[counter])
              counter=counter+1
            }
            else{
              draws[,counter] <- stats::qnorm(hdraws[,counter], random_coefs_means[counter], abs(rand_sdevs[counter]))
              counter=counter+1
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
        else{
          draws <- stats::qnorm(hdraws, random_coefs_means, abs(rand_sdevs))
        }

      }
      xb_rand_mat <- sapply(draws, function(x) X_rand * x)
    }

    rpar_mat <- exp(xb_rand_mat)

    pred_mat <- apply(rpar_mat, 2, function(x) x * mu_fixed)


    y <- model.response(mod1_frame)
    n_obs <- length(mu_fixed)
    ind_coefs <- matrix(0, nrow=n_obs, ncol=num_vars_rand)
    ind_var <- matrix(0, nrow=n_obs, ncol=num_vars_rand)
    prob_mat <- apply(pred_mat, 2, nb_prob, y=y, alpha=alpha, p=p)
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
