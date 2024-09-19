test_that("NB-2 model runs and returns correct output", {
  data("washington_roads")
  washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + ShouldWidth04 + AADTover10k,
               data = washington_roads, family = 'nb2', method = 'BHHH', max.iters = 3000)
  
  summary(model, confint_level=0.8, digits=4)
  
  expect_s3_class(model, "flexCountReg")  # Check the return class
  expect_true(length(model$model$estimate) > 0)  # Ensure estimates are returned
})

test_that("NB-1", {
  #data("washington_roads")
 # washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
  model1 <- countreg(Total_crashes ~ lnaadt + speed50,
                    data = washington_roads, family = "NB1",
                    offset = "lnlength",
                    stderr = "robust")
  
  expect_true(length(model1$model$estimate) > 0)  # Ensure estimates are returned
})


test_that("NB-p with sample weights", {
  set.seed(12345)
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  washington_roads$wgt <- runif(nrow(washington_roads))
  model2 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                    data = washington_roads, family = "NBP", 
                    weights = "wgt", max.iters=3000)
  
  summary(model2)
  
  expect_true(length(model2$model$estimate) > 0)  # Ensure estimates are returned
})

test_that("Poisson-Lognormal", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                    data = washington_roads, family = "PLN", ndraws=10)
  
  expect_s3_class(model, "flexCountReg")  # Check the return class
  expect_true(length(model$model$estimate) > 0)  # Ensure estimates are returned
})

test_that("Poisson Generalized-Exponential", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  
  model <- suppressWarnings(countreg(Animal ~ lnaadt,
                          data = washington_roads, family = "PGE", ndraws=5,
                          method="NM", 
                          offset = "lnlength",
                          dis_param_formula_2 = ~ -1+ AADT10kplus))
  
  expect_s3_class(model, "flexCountReg")  # Check the return class
})


test_that("Poisson-Inverse-Gaussian Type 1", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                    data = washington_roads, family = "PIG1")
  
  expect_s3_class(model, "flexCountReg")  # Check the return class
  expect_true(length(model$model$estimate) > 0)  # Ensure estimates are returned
})


test_that("Poisson-Inverse-Gaussian Type 2", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                    data = washington_roads, family = "PIG2")
  
  expect_s3_class(model, "flexCountReg")  # Check the return class
  expect_true(length(model$model$estimate) > 0)  # Ensure estimates are returned
})


test_that("Poisson-Lindley", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- suppressWarnings(countreg(Rollover ~ lnaadt + lnlength + speed50,
                    data = washington_roads, family = "PL"))
  
  expect_s3_class(model, "flexCountReg")  # Check the return class
  expect_true(length(model$model$estimate) > 0)  # Ensure estimates are returned
})

test_that("Poisson-Lindley-Gamma", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(Animal ~ lnaadt + speed50 + AADT10kplus, 
                     offset= "lnlength",
                    data = washington_roads, family = "PLG", ndraws=10, method="NM")
  
  expect_s3_class(model, "flexCountReg")  # Check the return class
  expect_true(length(model$model$estimate) > 0)  # Ensure estimates are returned
})

test_that("Poisson-Lindley-Lognormal", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(Animal ~ lnaadt + speed50 + AADT10kplus,
                     offset = "lnlength", 
                    dis_param_formula_2 = ~ 1 + AADT10kplus, 
                    data = washington_roads, family = "PLL", 
                    ndraws=10, method="NM")
  
  expect_s3_class(model, "flexCountReg")  # Check the return class
  expect_true(length(model$model$estimate) > 0)  # Ensure estimates are returned
})

test_that("Poisson-Weibull", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(Total_crashes ~ lnaadt + speed50 + AADT10kplus,
                     offset = "lnlength", 
                    dis_param_formula_1 = ~ 1 + lnlength,
                    data = washington_roads, family = "PW", ndraws=100, 
                    method="BHHH", stderr = "Robust")
  
  expect_s3_class(model, "flexCountReg")  # Check the return class
  expect_true(length(model$model$estimate) > 0)  # Ensure estimates are returned
})

test_that("Sichel", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(Total_crashes ~ lnaadt + speed50,
                    offset = "lnlength", 
                    data = washington_roads, family = "SI")
  
  expect_s3_class(model, "flexCountReg")  # Check the return class
  expect_true(length(model$model$estimate) > 0)  # Ensure estimates are returned
})


test_that("Generalized-Waring", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                    data = washington_roads, family = "GW")
  
  expect_s3_class(model, "flexCountReg")  # Check the return class
  expect_true(length(model$model$estimate) > 0)  # Ensure estimates are returned
})

test_that("NB1 with boostrapping", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                    data = washington_roads, family = "NB1", bootstraps = 100)
  
  expect_s3_class(model, "flexCountReg")  # Check the return class
  expect_true(length(model$model$estimate) > 0)  # Ensure estimates are returned
})


test_that("NB2 with underreporting (logit)", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                    data = washington_roads, family = "NB2",
                    underreport_formula = ~1 + speed50 + AADT10kplus)
  
  expect_s3_class(model, "flexCountReg")  # Check the return class
  expect_true(length(model$model$estimate) > 0)  # Ensure estimates are returned
})


test_that("Poisson-Lognormal with underreporting (probit)", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                    data = washington_roads, family = "NB2",
                    underreport_formula = ~ speed50 + AADT10kplus, underreport_family = "probit")
  
  expect_s3_class(model, "flexCountReg")  # Check the return class
  expect_true(length(model$model$estimate) > 0)  # Ensure estimates are returned
})

test_that("Conway-Maxwell-Poisson Model)", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- countreg(countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                             data = washington_roads, family = "COM"))
  
  expect_s3_class(model, "flexCountReg")  # Check the return class
  expect_true(length(model$model$estimate) > 0)  # Ensure estimates are returned
})

test_that("Poisson-Lindley RP", {
  data("washington_roads")
  washington_roads$AADT10kplus <- ifelse(washington_roads$AADT>10000,1,0)
  model <- rppLind(Animal ~ lnlength + lnaadt,
                   rpar_formula = ~ -1 + speed50,
                   data = washington_roads,
                   ndraws = 10,
                   correlated = FALSE,
                   rpardists = c(speed50="n"),
                   method = "nm",
                   print.level = 2)
  
  expect_s3_class(model, "flexCountReg")  # Check the return class
  expect_true(length(model$model$estimate) > 0)  # Ensure estimates are returned
})
