# Test the Poisson-Weibull model with alpha_formula
test_that("Poisson-Weibull model runs correctly with alpha_formula and sigma_formula", {
  data("washington_roads")
  washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
  model <- pwiebreg(Total_crashes ~ offset(lnaadt) + lnlength,
                    ndraws = 1500,
                    data = washington_roads,
                    alpha_formula = ~ -1 + lnaadt,
                    sigma_formula = ~ lnaadt,
                    method = 'NM')
  
  expect_named(model$model$estimate, c("(Intercept)", "lnlength", "ln(alpha):lnaadt", "ln(sigma):(Intercept)", "ln(sigma):lnaadt"))
})