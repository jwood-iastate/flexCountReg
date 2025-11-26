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