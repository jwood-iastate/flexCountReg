test_that("log dqwar and others", {
  expect_equal(dgwar(0, mu=1, k=2, rho=3, log=TRUE), -0.510825624)
  expect_equal(
    pgwar(1, mu=1, k=2, rho=3, log=TRUE, lower.tail=FALSE), -1.609437912)
  expect_equal(
    plindley(0.5, 1.5, log=TRUE, lower.tail=FALSE), -0.487635736)
  expect_equal(plindley(0.5, 1.5, log=TRUE), -0.952116164)
  skip_if_not_installed("lamW")
  expect_equal(qlindley(-1, theta = 1.5, log.p = TRUE), 0.4720311)
})

test_that("various distributions warnings", {
  expect_warning(rlindley(1.5,1))
  skip_if_not_installed("lamW")
  expect_warning(qlindley(c(0.1,-0.1,0.2), theta=1))
  expect_warning(qlindley(c(1.1,-0.1,0.2), theta=1))
  expect_warning(qlindley(c(0.1,0.3,0.2), theta=-1))
})



test_that("cor2cov", {
  C <- matrix(c(1,-0.3,0.7,-0.3,1,-0.2,0.7,-0.2,1), 3, 3)
  C1 <- c(1,-0.3,0.7,-0.3,1,-0.2,0.7,-0.2,1)
  S <- c(0.5, 2, 1.25, 3)
  S1 <- matrix(c(0.5, 2, 1.25))
  expect_error(cor2cov(C1,S))
  expect_error(cor2cov(C,S1))
  expect_error(cor2cov(C,S))
})

test_that("Correlated Halton Draws",{
  
  expect_error(corr_haltons(
    c(3, 2, 0.9), 
    c(0.25,1.5,0.8), 
    correlations = matrix(c( 1,    -0.3, 0.5, 
                            -0.3,   1,  -0.2, 
                             0.5,  -0.2, 1), 3, 3),
    ndraws=5000, scrambled=TRUE,
    dist="lognormal"), "Unsupported")
  
  expect_error(corr_haltons(
    c(3, 2, 0.9), 
    c(0.25,1.5,0.8),
    hdraws = matrix(c(0.1,0.3,0.9, 0.5,0.6,0.2), ncol = 3), 
    dist="normal"), "must be a 3 x 3 matrix")
  
})



test_that("Countreg Bootstrapping",{
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)

  mod <- countreg(
    Total_crashes ~ lnaadt + offset(lnlength) + speed50 + AADT10kplus,
    data = washington_roads, family = "PIG",
    method='BFGS',
    stderr = "boot",
    bootstraps = 3)
  
  expect_s3_class(mod, "flexCountReg")
  
  expect_error(countreg(
    Total_crashes ~ lnaadt + offset(lnlength) + speed50 + AADT10kplus,
    data = washington_roads, family = "Lindley"))
  # A short Nelder-Mead run stops at the iteration limit
  mod.rp <- suppressWarnings(countreg.rp(Total_crashes ~ lnlength + speed50,
                        rpar_formula = ~ -1 + lnaadt,
                        data = washington_roads,
                        family = "NB2",
                        rpardists = c( lnaadt = "g"),
                        ndraws = 10,
                        method = "NM",
                        max.iters = 100))
  
  pr <- predict(mod.rp, method = 'Exact')
  expect_length(pr, nrow(washington_roads))
  expect_true(all(is.finite(pr) & pr > 0))
  expect_s3_class(mod.rp, "flexCountReg")
  
})


