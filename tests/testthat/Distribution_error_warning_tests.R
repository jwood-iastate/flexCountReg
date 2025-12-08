library(testthat)
library(dplyr)
library(flexCountReg)

# expect_warning

test_that("log dqwar and others", {
  expect_equal(dgwar(0, mu=1, k=2, rho=3, log=TRUE), -0.5108256)
  expect_equal(pgwar(1, mu=1, k=2, rho=3, log=TRUE, lower.tail=FALSE), -1.609437912)
  expect_equal(plindley(0.5, 1.5, log=TRUE, lower.tail=FALSE), -0.4876357)
  expect_equal(plindley(0.5, 1.5, log=TRUE), -0.9521162)
  expect_equal(qlindley(-1, theta = 1.5, log.p = TRUE), 0.4720311)
})

test_that("qgwar errors and similar", {
  expect_error(qgwar(1.8, mu=1, k=2, rho=3))
})

test_that("various distributions warnings", {
  expect_error(dlindley(0, -1.5))
  expect_error(plindley (0.1, -1))
  expect_warning(qlindley(c(0.1,-0.1,0.2), theta=1))
  expect_warning(qlindley(c(1.1,-0.1,0.2), theta=1))
  expect_warning(qlindley(c(0.1,0.3,0.2), theta=-1))
  expect_warning(rlindley(1.5,1))
})

test_that("cor2cov", {
  C <- matrix(c(1,-0.3,0.7,-0.3,1,-0.2,0.7,-0.2,1), 3, 3)
  C1 <- c(1,-0.3,0.7,-0.3,1,-0.2,0.7,-0.2,1)
  S <- c(0.5, 2, 1.25, 3)
  S1 <- matrix(c(0.5, 2, 1.25))
  expect_warning(cor2cov(C1,S))
  expect_warning(cor2cov(C,S1))
  expect_warning(cor2cov(C,S))
})

test_that("Correlated Halton Draws",{
  means <- c(3, 2, 0.9)
  sdevs <- c(0.25,1.5,0.8)
  CORR <- matrix(c(1, -0.3, 0.5, -0.3, 1, -0.2, 0.5, -0.2, 1), 3, 3)
  
  # Create the Cholesky decomposition matrix and set values for ndraws, etc.
  ndraws <- 5000
  scrambled <- TRUE
  dist <- "normal"
  dist2 <- "lognormal"
  
  expect_error(corr_haltons(means, stdev=sdevs, correlations=CORR,
                            ndraws=ndraws, scrambled=scrambled,
                            dist=dist2))
  
  expect_error(corr_haltons(means, stdev=sdevs,hdraws=matrix(c(0.1,0.3,0.9)), correlations=CORR,
                            dist=dist))
  
  expect_error(corr_haltons(means, stdev=sdevs,hdraws=matrix(c(0.1,0.3,0.9, 0.5,0.6,0.2), ncol=3), 
                            dist=dist))
  
  corr_haltons(means, stdev=sdevs,hdraws=matrix(c(0.1,0.3,0.9, 0.5,0.6,0.2), ncol=3), correlations=CORR,
               dist=dist)
})

test_that("Countreg Bootstrapping",{
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)
  # Estimate an NB2 model with a dispersion parameter as a function of the
  # variable `speed50` (i.e., generalized NB2), verbose output, and use the
  # BFGS optimization method
  mod <- countreg(Total_crashes ~ lnaadt + offset(lnlength) + speed50 + AADT10kplus,
                 data = washington_roads, family = "PIG", verbose = TRUE,
                 method='SN',
                 bootstraps = 3)
  
  expect_s3_class(mod, "flexCountReg")
  
  expect_error(countreg(Total_crashes ~ lnaadt + offset(lnlength) + speed50 + AADT10kplus,
                          data = washington_roads, family = "Lindley"))
  mod.rp <- countreg.rp(Total_crashes ~ lnaadt + lnlength,
                                    rpar_formula = ~ -1 + speed50,
                                    data = washington_roads,
                                    family = "NB2",
                                    rpardists = c(speed50 = "g"),
                                    ndraws = 100,
                                    method = "BHHH")
  
  predict(mod.rp, method == 'Exact')
  
  expect_s3_class(mod, "mod.rp")
  
})
