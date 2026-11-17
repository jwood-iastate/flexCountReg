# With a fixed dispersion start of alpha = 1, quasi-Newton optimizers could
# step onto the flat region alpha -> 0 and stop there with a falsely
# "converged" Poisson fit (e.g. NB2 + BFGS on warpbreaks: logLik -228.48
# instead of -195.56). The moment-based start must reach the maximum.

test_that("NB2 with BFGS reaches the glm.nb maximum", {
  cases <- list(
    list(f = count ~ spray, d = datasets::InsectSprays),
    list(f = breaks ~ wool * tension, d = datasets::warpbreaks),
    list(f = Total_crashes ~ lnaadt + speed50 + offset(lnlength),
         d = flexCountReg::washington_roads))
  for (cs in cases) {
    ref <- MASS::glm.nb(cs$f, data = cs$d)
    for (m in c("BFGS", "NR")) {
      fit <- countreg(cs$f, data = cs$d, family = "NB2", method = m)
      expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)),
                   tolerance = 1e-6, label = paste(deparse(cs$f), m))
      expect_equal(coef(fit)[["ln(alpha)"]], log(1 / ref$theta),
                   tolerance = 1e-3, label = paste(deparse(cs$f), m))
    }
  }
})

test_that("the moment start uses ratios of sums and is clipped", {
  mstart <- flexCountReg:::.fcr_moment_disp_start
  y <- c(0, 5, 12, 20)   # overdispersed relative to mu
  mu <- c(1, 2, 4, 8)
  r <- (y - mu)^2 - y
  expect_equal(mstart("NB2", y, mu), log(sum(r) / sum(mu^2)))
  expect_equal(mstart("NB1", y, mu), log(sum(r) / sum(mu)))
  expect_equal(mstart("NB2", mu, mu), log(1e-3))   # underdispersed: floor
  expect_null(mstart("POISSON", y, mu))
})
