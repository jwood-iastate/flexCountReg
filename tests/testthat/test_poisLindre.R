test_that("Random Effects Poisson-Lindley model runs and returns correct output", {
  data("washington_roads")
  
  # Wrap the function call in suppressWarnings()
  model <- suppressWarnings(poisLind.re(Animal ~ lnaadt + lnlength ,
                                        data=washington_roads,
                                        group_var="ID",
                                        method="nm",
                                        max.iters = 1000))
  
  expect_s3_class(model, "flexCountReg")  # Check the return class
  expect_true(length(model$model$estimate) > 0)  # Ensure estimates are returned
  expect_named(model$model$estimate, c("(Intercept)", "lnaadt", "lnlength",  "ln(theta)"))  # Check names of estimates
})
