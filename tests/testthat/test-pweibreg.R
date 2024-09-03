# Test the Poisson-Weibull model with alpha_formula
test_that("Poisson-Weibull model runs correctly with alpha_formula", {
  data("washington_roads")
  washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
  model <- pwiebreg(Total_crashes ~ lnaadt + lnlength,
                    ndraws = 100,
                    alpha_formula = ~ lnaadt,
                    data = washington_roads,
                    method = 'NM', 
                    max.iters = 3000)
  
  expect_named(model$model$estimate, c("lnlength", "lnaadt", "ln(alpha):(Intercept)", "ln(alpha):lnaadt", "ln(sigma)"))
})

# Test the Poisson-Weibull model with sigma_formula
test_that("Poisson-Weibull model runs correctly with sigma_formula", {
  data("washington_roads")
  washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
  model <- pwiebreg(Total_crashes ~ lnaadt + lnlength,
                    ndraws = 100,
                    sigma_formula  = ~ lnaadt,
                    data = washington_roads,
                    method = 'NM', 
                    max.iters = 3000)
  
  expect_named(model$model$estimate, c("(Intercept)","lnlength", "lnaadt", "ln(alpha)", "ln(sigma):(Intercept)", "ln(sigma):lnaadt"))
})
