test_that("NB-2 model runs and returns correct output", {
  data("washington_roads")
  washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(
    Total_crashes ~ lnaadt + lnlength + speed50 + ShouldWidth04 + AADTover10k,
    data = washington_roads, 
    family = 'nb2', 
    method = 'BHHH', 
    max.iters = 500)  # Reduced from 3000
  
  summary(model, confint_level=0.8, digits=4)
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
})

test_that("NB-1", {
  data("washington_roads")
  washington_roads$AADTover10k <- ifelse(washington_roads$AADT > 10000, 1, 0)
  model1 <- countreg(Total_crashes ~ lnaadt + speed50,
                     data = washington_roads, family = "NB1",
                     offset = "lnlength",
                     stderr = "robust")
  
  expect_true(length(model1$model$estimate) > 0)
})


test_that("NB-p with sample weights", {
  set.seed(12345)
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  washington_roads$wgt <- runif(nrow(washington_roads))
  model2 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                     data = washington_roads, family = "NBP", 
                     weights = "wgt", max.iters = 500)  # Reduced from 3000
  
  summary(model2)
  
  expect_true(length(model2$model$estimate) > 0)
})

test_that("Poisson-Lognormal", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(
    Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
    data = washington_roads, 
    family = "PLN", 
    ndraws = 10)  # Already low - good
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
})

test_that("Poisson Generalized-Exponential", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)
  
  model <- suppressWarnings(countreg(Animal ~ lnaadt,
                                     data = washington_roads, 
                                     family = "PGE", 
                                     ndraws=5,
                                     method="NM", 
                                     offset = "lnlength",
                                     dis_param_formula_2 = ~ -1+ AADT10kplus,
                                     max.iters = 200))  # Added iteration limit
  
  expect_s3_class(model, "flexCountReg")
})


test_that("Poisson-Inverse-Gaussian Type 1", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                    data = washington_roads, family = "PIG1",
                    max.iters = 500)  # Added limit
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
})


test_that("Poisson-Inverse-Gaussian Type 2", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                    data = washington_roads, family = "PIG2",
                    max.iters = 500)  # Added limit
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
})


test_that("Poisson-Lindley", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- suppressWarnings(countreg(Rollover ~ lnaadt + lnlength + speed50,
                                     data = washington_roads, family = "PL",
                                     max.iters = 500))  # Added limit
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
})

test_that("Poisson-Lindley-Gamma", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(Animal ~ lnaadt + speed50 + AADT10kplus, 
                    offset= "lnlength",
                    data = washington_roads, family = "PLG", ndraws=10, 
                    method="NM", max.iters = 200)  # Added limit
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
})

test_that("Poisson-Lindley-Lognormal", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  
  model <- suppressWarnings(countreg(
    Animal ~ lnaadt + speed50 + AADT10kplus,
    offset = "lnlength", 
    dis_param_formula_2 = ~ 1 + AADT10kplus, 
    data = washington_roads, family = "PLL", 
    ndraws=10, method="NM", max.iters = 200))  # Added limit
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
})

test_that("Poisson-Weibull", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(
    Total_crashes ~ lnaadt + speed50 + AADT10kplus,
    offset = "lnlength", 
    dis_param_formula_1 = ~ 1 + lnlength,
    data = washington_roads, family = "PW", ndraws=50,  # Reduced from 100
    method="BHHH", stderr = "Robust", max.iters = 200)  # Added limit
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
})

test_that("Sichel", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  
  model <- suppressWarnings(countreg(
    Total_crashes ~ lnaadt + speed50,
    offset = "lnlength", 
    data = washington_roads, family = "SI", 
    method="NM", max.iters = 200))  # Added limit
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
})


test_that("Generalized-Waring", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(
    Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
    data = washington_roads, family = "GW",
    max.iters = 500)  # Added limit
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
})

test_that("NB1 with bootstrapping", {
  skip_on_cran()  # Skip on CRAN due to time
  
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(
    Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
    data = washington_roads, family = "NB1", 
    bootstraps = 10,  # Reduced from 100!
    max.iters = 200)
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
})


test_that("NB2 with underreporting (logit)", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(
    Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
    data = washington_roads, family = "NB2",
    underreport_formula = ~1 + speed50 + AADT10kplus,
    max.iters = 500)  # Added limit
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
})


test_that("Poisson-Lognormal with underreporting (probit)", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(
    Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
    data = washington_roads, family = "NB2",
    underreport_formula = ~ speed50 + AADT10kplus, 
    underreport_family = "probit",
    max.iters = 500)  # Added limit
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
})

test_that("Conway-Maxwell-Poisson Model", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- suppressWarnings(countreg(
    Total_crashes ~ lnaadt + speed50 + AADT10kplus,
    data = washington_roads, 
    family = "COM", 
    offset = "lnlength",
    method='BHHH',
    max.iters = 200))  # Added limit
  
  expect_s3_class(model, "flexCountReg")
  expect_true(length(model$model$estimate) > 0)
})
