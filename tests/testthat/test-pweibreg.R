# Test the Poisson-Weibull model with alpha_formula
test_that("Poisson-Weibull model runs correctly with alpha_formula and sigma_formula", {
  data("washington_roads")
  model <- pwiebreg(Total_crashes ~ lnaadt + lnlength,
                    ndraws = 10,
                    data = washington_roads,
                    alpha_formula = ~ -1 + lnaadt,
                    sigma_formula = ~ lnaadt,
                    method = 'NM')
  
  expect_named(model$model$estimate, c("(Intercept)", "lnaadt", "lnlength", "ln(alpha):lnaadt", "ln(sigma):(Intercept)", "ln(sigma):lnaadt"))
})
