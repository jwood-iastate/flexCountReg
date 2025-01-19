test_that("NB-2 model runs and returns predictions", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>=10000,1,0)
  nb2 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                  data = washington_roads, family = "NB2",
                  dis_param_formula_1 = ~ speed50, method='BFGS')
  
  pred <- predict(nb2, data=washington_roads)
  
  expect_true(length(pred) > 0)  # Ensure predictions are returned
})

test_that("NB-2 model runs and returns predictions", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>=10000,1,0)
  nb2 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                  data = washington_roads, family = "NB2",
                  dis_param_formula_1 = ~ speed50, method='BFGS')
  
  pred <- predict(nb2, data=washington_roads)
  
  expect_true(length(pred) > 0)  # Ensure predictions are returned
})

test_that("Poisson-Lognormal model runs and returns predictions", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>=10000,1,0)
  pln <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                  data = washington_roads, family = "PLN", ndraws=10)
  
  pred <- predict(pln, data=washington_roads)
  
  expect_true(length(pred) > 0)  # Ensure predictions are returned
})

test_that("Random Parameter NB2 model runs and returns predictions", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>=10000,1,0)
  nb2.rp <- rpnb(Total_crashes ~ - 1 + lnlength + lnaadt,
                 rpar_formula = ~ speed50,
                 data = washington_roads,
                 ndraws = 100,
                 correlated = FALSE,
                 rpardists = c(intercept="u", speed50="t"),
                 form = 'nb2',
                 method = "bfgs")
  
  pred <- predict(nb2.rp, list(data=washington_roads, method="Simulated"))
  pred2 <- predict(nb2.rp, list(data=washington_roads, method="Exact"))
  pred3 <- predict(nb2.rp, list(data=washington_roads, method="Individual"))
  
  expect_true(length(pred) > 0)  # Ensure predictions are returned
  expect_true(length(pred2) > 0)  # Ensure predictions are returned
  expect_true(length(pred3) > 0)  # Ensure predictions are returned
})

# test_that("Random Parameter NB1 with correlation model runs and returns predictions", {
#   data("washington_roads")
#   washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>=10000,1,0)
#   nb1.rp <- rpnb(Total_crashes ~ - 1 + lnlength + lnaadt,
#                  rpar_formula = ~ speed50,
#                  data = washington_roads,
#                  ndraws = 10,
#                  correlated = TRUE,
#                  rpardists = c(intercept="n", speed50="n"),
#                  form = 'nb1',
#                  method = "bfgs")
#   
#   pred <- predict(nb1.rp, list(data=washington_roads, method="Simulated"))
#   pred2 <- predict(nb1.rp, list(data=washington_roads, method="Exact"))
#   pred3 <- predict(nb1.rp, list(data=washington_roads, method="Individual"))
#   
#   expect_true(length(pred) > 0)  # Ensure predictions are returned
#   expect_true(length(pred2) > 0)  # Ensure predictions are returned
#   expect_true(length(pred3) > 0)  # Ensure predictions are returned
# })

test_that("Random Parameter NBP model runs and returns predictions", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>=10000,1,0)
  nbp.rp <- rpnb(Total_crashes ~ lnlength + speed50,
                 rpar_formula = ~ - 1 + lnaadt,
                 data = washington_roads,
                 ndraws = 10,
                 correlated = FALSE,
                 rpardists = c(lnaadt="g"), # trying gamma distribution
                 form = 'nbp',
                 method = "nm")
  
  pred <- predict(nbp.rp, list(data=washington_roads, method="Simulated"))
  pred2 <- predict(nbp.rp, list(data=washington_roads, method="Exact"))
  pred3 <- predict(nbp.rp, list(data=washington_roads, method="Individual"))
  
  expect_true(length(pred) > 0)  # Ensure predictions are returned
  expect_true(length(pred2) > 0)  # Ensure predictions are returned
  expect_true(length(pred3) > 0)  # Ensure predictions are returned
})

test_that("Random Parameter NB2 model runs and returns predictions - 2", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>=10000,1,0)
  nb2.rp <- rpnb(Total_crashes ~ speed50 ,
                 rpar_formula = ~ -1 + lnlength+ lnaadt,
                 data = washington_roads,
                 ndraws = 10,
                 correlated = FALSE,
                 rpardists = c(lnlength="u", lnaadt="ln"), # trying lognormal distribution
                 form = 'nb2',
                 method = "NM")
  
  pred <- predict(nb2.rp, list(data=washington_roads, method="Simulated"))
  pred2 <- predict(nb2.rp, list(data=washington_roads, method="Exact"))
  pred3 <- predict(nb2.rp, list(data=washington_roads, method="Individual"))
  
  expect_true(length(pred) > 0)  # Ensure predictions are returned
  expect_true(length(pred2) > 0)  # Ensure predictions are returned
  expect_true(length(pred3) > 0)  # Ensure predictions are returned
})



test_that("NB2 with underreporting (logit) model runs and returns predictions", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>=10000,1,0)
  nb2_underreport <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                              data = washington_roads, family = "NB2",
                              underreport_formula = ~ speed50 + AADT10kplus)
  
  pred <- predict(nb2_underreport, data=washington_roads)
  
  expect_true(length(pred) > 0)  # Ensure predictions are returned
})

test_that("Poisson-Lognormal with underreporting (probit) model runs and returns predictions", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>=10000,1,0)
  plogn_underreport <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                                data = washington_roads, family = "NB2",
                                underreport_formula = ~ speed50 + AADT10kplus, underreport_family = "probit")
  
  pred <- predict(plogn_underreport, data=washington_roads)
  
  expect_true(length(pred) > 0)  # Ensure predictions are returned
})