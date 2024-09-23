get_halton_draws <-function(ndraws=500, dim=1, scrambled=FALSE){

  draws <- halton(ndraws, dims, mixed=scrambled)
  
  rpardraws <- draws[,1:dim]
  
  distdraws <- draws[,draws] # draws used for count distributions that require numerical integration
  
  return(c(rpardraws, distdraws))
}

data_prep <- function(formula,
                      data,
                      rpar_formula=NULL, 
                      panel_rpar_formula=NULL, 
                      underrport_formula=NULL, 
                      parmam1_formula=NULL, 
                      param2_formula=NULL, 
                      underreport_rpar_formula=NULL, 
                      het_means_formula=NULL, 
                      het_variance_formula=NULL, 
                      panel_ids=NULL){
  
  require(modelr)
  require(dplyr)
  require(stats)
  require(plm)
  
  y_name <- all.vars(formula)[1]
  y <- data[[y_name]]
  
  if (!is.null(panel_ids)){
    data <- plm::pdata.frame(data, panel_ids)
    x_matrix <- function(data, formula){ # model.matrix function for panel data
      plm::model.matrix(data, formula)
    }
    
  }else{
    x_matrix <- function(data, formula){ # model.matrix function for cross-sectional data
      as.matrix(modelr::model_matrix(data, formula))
    }
  }
  
  X_fixed <- x_matrix(data, formula)
  if (!is.null(rpar_formula)) X_rand_crosssectional <- x_matrix(data, rpar_formula) else X_rand_crosssectional <- NULL
  if (!is.null(panel_rpar_formula)) X_rand_panel <- x_matrix(data, panel_rpar_formula) else X_rand_panel <- NULL
  if (!is.null(underrport_formula)) X_underreport <- x_matrix(data, underrport_formula) else X_underreport <- NULL
  if (!is.null(parmam1_formula)) X_param1 <- x_matrix(data, parmam1_formula) else X_param1 <- NULL
  if (!is.null(param2_formula)) X_param2 <- x_matrix(data, param2_formula) else X_param2 <- NULL
  if (!is.null(underreport_rpar_formula)) X_underreport_rpar <- x_matrix(data, underreport_rpar_formula)
  if (!is.null(het_means_formula)) X_het_means <- x_matrix(data, het_means_formula) else X_het_means <- NULL
  if (!is.null(het_variance_formula)) X_het_variance <- x_matrix(data, het_variance_formula) else X_het_variance <- NULL
  
  return(list(y=y, 
              X_fixed=X_fixed, 
              X_rand_crosssectional=X_rand_crosssectional,
              X_rand_panel=X_rand_panel, 
              X_underreport=X_underreport,
              X_param1=X_param1, 
              X_param2=X_param2, 
              X_underreport_rpar=X_underreport_rpar, 
              X_het_means=X_het_means, 
              X_het_variance=X_het_variance))
  
}

# Get number of parameters required for each part of the model
get_Nparams <- function(prepped_data, family="NB2"){ # inputs are the output from `data_prep` and the family name
  
  N_fixed <- ncol(prepped_data$X_fixed)
  N_rand_crosssectional <- ifelse(is.null(prepped_data$X_rand_crosssectional), 0, ncol(prepped_data$X_rand_crosssectional))
  N_rand_panel <- ifelse(is.null(prepped_data$X_rand_panel), 0, ncol(prepped_data$X_rand_panel))
  N_underreport <- ifelse(is.null(prepped_data$X_underreport), 0, ncol(prepped_data$X_underreport))
  N_param1 <- ifelse(is.null(prepped_data$X_param1), 0, ncol(prepped_data$X_param1))
  N_param2 <- ifelse(is.null(prepped_data$X_param2), 0, ncol(prepped_data$X_param2))
  N_underreport_rpar <- ifelse(is.null(prepped_data$X_underreport_rpar), 0, ncol(prepped_data$X_underreport_rpar))
  N_het_means <- ifelse(is.null(prepped_data$X_het_means), 0, ncol(prepped_data$X_het_means))
  N_het_variance <- ifelse(is.null(prepped_data$X_het_variance), 0, ncol(prepped_data$X_het_variance))
  
  family_params <- length(Filter(Negate(is.null),  get_params(family))) # No. of params for the distribution
  
  return(list(N_fixed=N_fixed, 
              N_rand_crosssectional=N_rand_crosssectional, 
              N_rand_panel=N_rand_panel, 
              N_underreport=N_underreport, 
              N_param1=N_param1, 
              N_param2=N_param2, 
              N_underreport_rpar=N_underreport_rpar, 
              N_het_means=N_het_means, 
              N_het_variance=N_het_variance, 
              N_family_params=family_params,
              N_rand_total = N_rand_crosssectional + N_rand_panel + N_underreport_rpar))
  
}

# Create starting values for the optimization
get_start_vals <- function(formula, 
                           prepped_data, # use prepped_data from data_prep and the original data
                           data,
                           family="NB2",
                           rpar_formula=NULL, 
                           panel_rpar_formula=NULL,
                           correlated=FALSE){ # correlated is a boolean value
  
  X_matrix_counts <- get_Nparams(prepped_data, family)
  family_params <- get_params(family)
  
  X_fixed <- prepped_data$X_fixed
  X_rand_crosssectional <- prepped_data$X_rand_crosssectional
  X_rand_panel <- prepped_data$X_rand_panel 
  X_underreport <- prepped_data$X_underreport
  X_param1 <- prepped_data$X_param1 
  X_param2 <- prepped_data$X_param2 
  X_underreport_rpar <- prepped_data$X_underreport_rpar 
  X_het_means <- prepped_data$X_het_means 
  X_het_variance <- prepped_data$X_het_variance
  
  y_name <- all.vars(formula)[1]
  
  # Create named vectors
  modelterms <- colnames(X_fixed)
  if (!is.null(X_rand_crosssectional)) modelterms <- c(modelterms, colnames(X_rand_crosssectional))
  if (!is.null(X_rand_panel)) modelterms <- c(modelterms, colnames(X_rand_panel))
  
  nb_vars <- predictor_terms[!grepl("ntercept", modelterms)] # remove the intercept
  nb_formula <- reformulate(nb_vars, response = y_name, intercept = TRUE)
  
  nb_model <- glm.nb(nb_formula, data)
  params <- coef(nb_model)
  
  varnameorder <- colnames(X_fixed)
  
  if (!is.null(X_rand_crosssectional)) varnameorder <- c(varnameorder, colnames(X_rand_crosssectional))
  if (!is.null(X_rand_panel)) varnameorder <- c(varnameorder, colnames(X_rand_panel))
  
  varnameorder <- c(colnames(X_fixed), colnames(X_rand))
  
  match_indices <- match(varnameorder, names(params))
  start_init <- params[match_indices] # correctly ordered starting values
  
  start <- start_init[1:N_fixed]
  if (!is.null(X_rand_crosssectional)) start_rand_crosssectional <- start_init[(N_fixed+1):(N_fixed+N_rand_crosssectional)] else start_rand_crosssectional <- NULL
  if (!is.null(X_rand_panel)) start_rand_panel <- start_init[(N_fixed+N_rand_crosssectional+1):(N_fixed+N_rand_crosssectional+N_rand_panel)] else start_rand_panel <- NULL
  
  x_names <- colnames(X_fixed)
  
  if (X_matrix_counts$N_underreport > 0){
    x_names <- c(x_names, paste0("Underreporting:",colnames(X_underreport)))
    start <- c(start, rep(0.1, N_underreport))
  }
  
  if (correlated){
    randnames_chol <- c()
    randnames <- c()
    # Add mean values
    if (!is.null(X_rand_crosssectional)) {
      randnames <- c(randnames, paste0("Mean:", colnames(X_rand_crosssectional)))
      start <- c(start, start_rand_crosssectional)
    }
    if (!is.null(X_rand_panel)) {
      randnames <- c(randnames, paste0("Mean:", colnames(X_rand_panel)))
      start <- c(start, start_rand_panel)
    }
    if (!is.null(X_underreport_rpar)) {
      randnames <- c(randnames, paste0("Underreport:Mean:",colnames(X_underreport_rpar)))
      start <- c(start, rep(0, N_underreport_rpar))
    }

    if (!is.null(X_rand_crosssectional)) randnames_chol <- c(randnames_chol, colnames(X_rand_crosssectional))
    if (!is.null(X_rand_panel)) randnames_chol <- c(randnames_chol, colnames(X_rand_panel))
    if (!is.null(X_underreport_rpar)) randnames_chol <- c(randnames_chol, paste0("Underreport:",colnames(X_underreport_rpar)))
    
    N_Rand_Total_chol <- length(randnames_chol)
    
    rparam_var <- rep(0.1, N_Rand_Total_chol)
    rparam_var <- diag(rparam_var)
    Chl <- chol(rparam_var) # initial cholesky matrix
    
    for (i in 1:N_Rand_Total_chol){
      for (j in 1:N_Rand_Total_chol){
        if (i >= j){
          start <- append(start, Chl[j,i])
          x_names <- append(x_names, paste('Cholesky Value for' ,paste(randnames_chol[j], randnames_chol[i], sep=":")))
        }
      }
    }
    
  }else{
    if (X_matrix_counts$N_rand_crosssectional > 0){
      x_names <- c(x_names, paste0("Mean:",colnames(X_rand_crosssectional)))
      start <- c(start, start_rand_crosssectional)
    }
    if (X_matrix_counts$N_rand_panel > 0) {
      x_names <- c(x_names, paste0("Mean:",colnames(X_rand_panel)))
      start <- c(start, start_rand_panel)
    }
    if (X_matrix_counts$N_underreport_rpar > 0) {
      x_names <- c(x_names, paste0("Underreporting:Mean:",colnames(X_underreport_rpar))) 
      start <- c(start, rep(0, N_underreport_rpar))
    }
    
    if (X_matrix_counts$N_rand_crosssectional > 0){
      x_names <- c(x_names, paste0("Std:",colnames(X_rand_crosssectional)))
      start <- c(start, rep(0.1, N_rand_crosssectional))
    }
    if (X_matrix_counts$N_rand_panel > 0) {
      x_names <- c(x_names, paste0("Std:",colnames(X_rand_panel)))
      start <- c(start, rep(0.1, N_rand_panel))
    }
    if (X_matrix_counts$N_underreport_rpar > 0) {
      x_names <- c(x_names, paste0("Underreporting:Std:",colnames(X_underreport_rpar)))
      start <- c(start, rep(0.1, N_underreport_rpar))
    }
  }
  
  if (X_matrix_counts$N_het_means > 0) {
    x_names <- c(x_names, paste0("HetMeans:",colnames(X_het_means)))
    start <- c(start, rep(0, N_het_means))
  }
  if (X_matrix_counts$N_het_variance > 0) {
    x_names <- c(x_names, paste0("HetVariance:",colnames(X_het_variance)))
    start <- c(start, rep(0.1, N_het_variance))
  }
  
  if (X_matrix_counts$N_param1 > 0) {
    x_names <- c(x_names, paste0(family_params[1],colnames(X_param1)))
    start <- c(start, rep(0.1, N_param1))
  }else {
    x_names <- c(x_names, family_params[1])
    start <- append(start, 0.1)
  }
  if (X_matrix_counts$N_param2 > 0) {
    x_names <- c(x_names, paste0(family_params[2],colnames(X_param2)))
    start <- c(start, rep(0.1, N_param2))
  }
  else if (!is.null(family_params[2])) {
    x_names <- c(x_names, family_params[2])
    start <- append(start, 0.1)
  }
  names(start) <- x_names
  
  return(start) 
}

# Extract coefficients for different parts of the model
coef_vals <- function(coeffs, 
                      formula, 
                      prepped_data, # use prepped_data from data_prep and the original data
                      data,
                      family="NB2",
                      correlated=FALSE){
  
  X_matrix_counts <- get_Nparams(prepped_data, family)
  
  N_fixed <- X_matrix_counts$N_fixed 
  N_rand_crosssectional <- X_matrix_counts$N_rand_crosssectional 
  N_rand_panel <- X_matrix_counts$N_rand_panel 
  N_underreport <- X_matrix_counts$N_underreport 
  N_param1 <- X_matrix_counts$N_param1 
  N_param2 <- X_matrix_counts$N_param2 
  N_underreport_rpar <- X_matrix_counts$N_underreport_rpar 
  N_het_means <- X_matrix_counts$N_het_means 
  N_het_variance <- X_matrix_counts$N_het_variance
  N_rand_total <- X_matrix_counts$N_rand_total
  
  if (correlated) N_rpar_std <- (N_rand_total*(N_rand_total+1))/2  else  N_rpar_std <- N_rand_total
  
  fixed_indices <- 1:N_fixed
  underrep_indices <- (N_fixed+1):(N_fixed+N_underreport)
  rand_crssc_indices <- (N_fixed+N_underreport+1):(N_fixed+N_underreport+N_rand_crosssectional)
  rand_panel_indices <- (N_fixed+N_underreport+N_rand_crosssectional+1):(N_fixed+N_underreport+N_rand_crosssectional+N_rand_panel)
  underrep_rpar_indices <- (N_fixed+N_underreport+N_rand_crosssectional+N_rand_panel+1):(N_fixed+N_underreport+N_rand_crosssectional+N_rand_panel+N_underreport_rpar)
  rand_crssc_std_indices <- (N_fixed+N_underreport+N_rand_crosssectional+N_rand_panel+N_underreport_rpar+1):(N_fixed+N_underreport+N_rand_crosssectional+N_rand_panel+N_underreport_rpar+N_rand_crosssectional)
  rand_panel_std_indices <- (N_fixed+N_underreport+N_rand_crosssectional+N_rand_panel+N_underreport_rpar+N_rand_crosssectional+1):(N_fixed+N_underreport+N_rand_crosssectional+N_rand_panel+N_underreport_rpar+N_rand_crosssectional+N_rand_panel)
  underrep_rpar_std_indices <- (N_fixed+N_underreport+N_rand_crosssectional+N_rand_panel+N_underreport_rpar+N_rand_crosssectional+N_rand_panel+1):(N_fixed+N_underreport+N_rand_crosssectional+N_rand_panel+N_underreport_rpar+N_rand_crosssectional+N_rand_panel+N_underreport_rpar)
  het_means_indices <- (N_fixed+N_underreport+N_rand_total+N_rpar_std+1):(N_fixed+N_underreport+N_rand_total+N_rpar_std+N_het_means)
  het_variance_indices <- (N_fixed+N_underreport+N_rand_total+N_rpar_std+N_het_means+1):(N_fixed+N_underreport+N_rand_total+N_rpar_std+N_het_means+N_het_variance)
  family_param1_indices <- (N_fixed+N_underreport+N_rand_total+N_rpar_std+N_het_means+N_het_variance+1):(N_fixed+N_underreport+N_rand_total+N_rpar_std+N_het_means+N_het_variance+N_param1)
  family_param2_indices <- (N_fixed+N_underreport+N_rand_total+N_rpar_std+N_het_means+N_het_variance+N_param1+1):(N_fixed+N_underreport+N_rand_total+N_rpar_std+N_het_means+N_het_variance+N_param1+N_param2)
  if (correlated) chol_indices <- (N_fixed+N_underreport+N_rand_total+1):(N_fixed+N_underreport+N_rand_total+N_rpar_std)
  
  
  beta_fixed <- coeffs[fixed_indices]
  if (X_matrix_counts$N_underreport>0) beta_underreport <- coeffs[underrep_indices] else beta_underreport <- NULL
  if (!is.null(X_rand_crosssectional)) beta_crss_mean <- coeffs[rand_crssc_indices] else beta_rand_crosssectional <- NULL
  if (!is.null(X_rand_panel)) beta_panel_mean <- coeffs[rand_panel_indices] else beta_rand_panel <- NULL
  if (!is.null(X_underreport_rpar)) beta_underreport_rpar <- coeffs[underrep_rpar_indices] else beta_underreport_rpar <- NULL
  if (correlated) beta_chol <- coeffs[chol_indices] else beta_chol <- NULL
  if (!correlated) {
    if (!is.null(X_rand_crosssectional)) beta_crss_std <- coeffs[rand_crssc_std_indices] else beta_rand_crosssectional <- NULL
    if (!is.null(X_rand_panel)) beta_panel_std <- coeffs[rand_panel_std_indices] else beta_rand_panel <- NULL
    if (!is.null(X_underreport_rpar)) beta_underreport_rpar_std <- coeffs[underrep_rpar_std_indices] else beta_underreport_rpar_std <- NULL
  }
  if (N_het_means>0) beta_het_means <- coeffs[het_means_indices] else beta_het_means <- NULL
  if (N_het_variance>0) beta_het_variance <- coeffs[het_variance_indices] else beta_het_variance <- NULL
  beta_param1 <- coeffs[family_param1_indices]
  if (N_param2>0) beta_param2 <- coeffs[family_param2_indices] else beta_param2 <- NULL
  
  return(list(beta_fixed=as.vector(beta_fixed), 
              beta_underreport=as.vector(beta_underreport), 
              beta_crss_mean=as.vector(beta_crss_mean), 
              beta_panel_mean=as.vector(beta_panel_mean), 
              beta_underreport_rpar=as.vector(beta_underreport_rpar), 
              beta_crss_std=as.vector(beta_crss_std), 
              beta_panel_std=as.vector(beta_panel_std), 
              beta_underreport_rpar_std=as.vector(beta_underreport_rpar_std), 
              beta_het_means=as.vector(beta_het_means), 
              beta_het_variance=as.vector(beta_het_variance), 
              beta_param1=as.vector(beta_param1), 
              beta_param2=as.vector(beta_param2), 
              beta_chol=as.vector(beta_chol)))
          
}