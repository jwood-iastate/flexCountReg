test_that("Sichel regression coefficients", {
  data("washington_roads")
  
  suppressMessages({
    sichel.mod <- sichel(Total_crashes ~ lnaadt + lnlength,
                         data = washington_roads,
                         method = "NM",
                         max.iters = 1000)
  })
  coefs <- coef(sichel.mod)
  test_value <- sum(abs(coefs - c(-9.21, 1.12, 0.74, 4.05, 2.5))) < 0.02
  
  testthat::expect_true(test_value)
})




