# Time-varying effects: tv(), temporal_effects(), temporal_tests()

wr <- function() {
  data("washington_roads", package = "flexCountReg", envir = environment())
  washington_roads
}

nb2 <- function(f, d, ...) {
  countreg(f, data = d, family = "NB2", method = "BFGS", ...)
}

test_that("tv() builds the step, level and hinge bases", {
  yr <- c(2016, 2017, 2018, 2017)
  x <- c(1, 2, 3, 4)
  s <- tv(x, yr)
  expect_equal(colnames(s), c("yr2016", "yr_step2017", "yr_step2018"))
  expect_equal(unclass(s)[, 1:3, drop = FALSE],
               cbind(x, x * c(0, 1, 1, 1), x * c(0, 0, 1, 0)),
               ignore_attr = TRUE)
  expect_equal(unclass(tv(x, yr, "level")), cbind(x * (yr == 2016),
                                                  x * (yr == 2017),
                                                  x * (yr == 2018)),
               ignore_attr = TRUE)
  # Hinge: max(t - k + 1, 0) with t the period's position
  h <- tv(x, yr, "hinge")
  expect_equal(unclass(h)[3, ], 3 * c(3, 2, 1), ignore_attr = TRUE)
  expect_equal(unclass(h)[4, ], 4 * c(2, 1, 0), ignore_attr = TRUE)
  expect_equal(attr(s, "levels"), c(2016, 2017, 2018))
  # Missing values give NA rows; unknown times are refused
  expect_true(all(is.na(tv(c(1, NA), c(2016, 2017), levels = 2016:2017)[2, ])))
  expect_error(tv(x, yr, levels = 2016:2017), "not in the estimation data")
  expect_equal(unclass(tv(1, 2020, levels = 2016:2018,
                          new_time = "nearest"))[1, ], c(1, 1, 1),
               ignore_attr = TRUE)
  expect_error(tv(factor(1:2), 1:2), "must be numeric")
  expect_error(tv(1:2, c(1, 1)), "fewer than two")
  expect_error(tv(1, 1:2, basis = "smooth"), "model formula")
})

test_that("the three bases give the same fit and the unconstrained effects", {
  w <- wr()
  f <- function(b) {
    nb2(stats::as.formula(paste0(
      "Total_crashes ~ tv(lnaadt, time = Year, basis = '", b,
      "') + speed50 + offset(lnlength)")), w)
  }
  ms <- f("step")
  ml <- f("level")
  mh <- f("hinge")
  expect_equal(names(coef(ms))[2:4],
               c("lnaadt:Year2016", "lnaadt:Year_step2017",
                 "lnaadt:Year_step2018"))
  expect_equal(as.numeric(logLik(ml)), as.numeric(logLik(ms)),
               tolerance = 1e-6)
  expect_equal(as.numeric(logLik(mh)), as.numeric(logLik(ms)),
               tolerance = 1e-6)
  w$fY <- factor(w$Year)
  mu <- nb2(Total_crashes ~ lnaadt:fY + speed50 + offset(lnlength), w)
  ref <- unname(coef(mu)[grep("lnaadt:fY", names(coef(mu)))])
  es <- temporal_effects(ms)
  expect_equal(es$period, c("2016", "2017", "2018"))
  expect_equal(es$part, rep("mean", 3))
  expect_equal(es$estimate, ref, tolerance = 1e-4)
  expect_equal(temporal_effects(mh)$estimate, ref, tolerance = 1e-4)
  expect_equal(temporal_effects(ml)$std.error, es$std.error,
               tolerance = 1e-3)
  # The tests are the same for every basis
  ts <- temporal_tests(ms)
  expect_equal(ts$test, c("all equal", "2017 vs 2016", "2018 vs 2017"))
  expect_equal(ts$df[1], 2)
  expect_equal(temporal_tests(mh)$p.value, ts$p.value, tolerance = 1e-3)
  # The step coefficients are the adjacent changes
  expect_equal(ts$estimate[2:3], unname(coef(ms)[3:4]), tolerance = 1e-8)
  # Prediction with new data reuses the training periods
  nd <- w[w$Year == 2018, ][1:3, ]
  expect_equal(predict(ms, newdata = nd), unname(fitted(ms)[w$Year == 2018][1:3]),
               tolerance = 1e-8, ignore_attr = TRUE)
  nd$Year <- 2025
  expect_error(predict(ms, newdata = nd), "not in the estimation data")
  # Marginal effects vary the variable, not the time
  me <- margEffTable(ms)
  expect_setequal(me$term, c("lnaadt", "speed50"))
})

test_that("tv(1, ...) and tv() in distribution-parameter formulas", {
  w <- wr()
  m <- nb2(Total_crashes ~ tv(1, time = Year) + lnaadt + offset(lnlength), w,
           dis_param_formula_1 = ~ tv(1, time = Year, basis = "level"))
  expect_true(all(c("(Intercept):Year_step2017",
                    "ln(alpha):(Intercept):Year2017") %in% names(coef(m))))
  e <- temporal_effects(m)
  expect_equal(unique(e$part), c("mean", "ln(alpha)"))
  cf <- coef(m)
  expect_equal(e$estimate[1:3],
               unname(cf["(Intercept)"] + c(0, cf[2], cf[2] + cf[3])),
               tolerance = 1e-10)
  w$fY <- factor(w$Year)
  ref <- nb2(Total_crashes ~ fY + lnaadt + offset(lnlength), w,
             dis_param_formula_1 = ~ fY)
  expect_equal(as.numeric(logLik(m)), as.numeric(logLik(ref)),
               tolerance = 1e-6)
  expect_equal(temporal_effects(m, "(Intercept)")$term[1], "(Intercept)")
  expect_error(temporal_effects(m, "speed50"), "No time-varying effect")
})

test_that("the smooth basis is a varying-coefficient smooth", {
  skip_if_not_installed("mgcv")
  w <- wr()
  m <- nb2(Total_crashes ~ tv(lnaadt, time = Year, basis = "smooth") +
             speed50 + offset(lnlength), w)
  expect_true(any(grepl("^s\\(Year\\):lnaadt", names(coef(m)))))
  e <- temporal_effects(m)
  expect_equal(e$basis, rep("smooth", 3))
  expect_equal(nrow(e), 3)
  expect_true(all(e$std.error > 0))
})

test_that("tv() works in the panel, random-parameter and hierarchical models", {
  skip_on_cran()
  w <- wr()
  f <- Total_crashes ~ tv(lnaadt, time = Year) + offset(lnlength)
  re <- suppressWarnings(renb(f, group_var = "ID", data = w,
                              method = "BFGS"))
  expect_equal(nrow(temporal_effects(re)), 3)
  h <- countreg_hier(update(f, . ~ . + (1 | ID)), data = w,
                     family = "POISSON")
  expect_equal(nrow(temporal_effects(h)), 3)
  expect_length(predict(h, newdata = w[1:3, ]), 3)
  r <- suppressWarnings(countreg_rp(
    Total_crashes ~ speed50 + offset(lnlength),
    rpar_formula = ~ -1 + tv(lnaadt, time = Year), data = w,
    family = "NB2", ndraws = 50, method = "BFGS", group_var = "ID"))
  e <- temporal_effects(r)
  expect_equal(unique(e$part), "random parameter mean")
  expect_equal(e$estimate[1], unname(coef(r)["lnaadt:Year2016:Mean"]))
})

test_that("models without tv() terms are refused", {
  m <- nb2(Total_crashes ~ lnaadt + offset(lnlength), wr())
  expect_error(temporal_effects(m), "no time-varying effects")
})
