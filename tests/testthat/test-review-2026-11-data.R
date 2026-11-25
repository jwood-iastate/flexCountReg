# Review corrections (2026-11): unused factor levels in the model data,
# temporal_effects() for multivariate fits and tv() smooths, and
# margEffTable() for multivariate fits.

wr <- function() {
  data("washington_roads", package = "flexCountReg", envir = environment())
  washington_roads
}

nb2 <- function(f, d, ...) {
  countreg(f, data = d, family = "NB2", method = "BFGS", ...)
}

test_that("unused factor levels are dropped, as in lm()", {
  w <- wr()
  w$sw <- factor(ifelse(w$ShouldWidth04 == 1, "wide", "narrow"),
                 levels = c("narrow", "wide", "unknown"))
  ref <- coef(glm(Total_crashes ~ lnaadt + sw + offset(lnlength),
                  data = w, family = poisson))
  # Declared level never used: main and auxiliary formulas
  m <- nb2(Total_crashes ~ lnaadt + sw + offset(lnlength), w,
           dis_param_formula_1 = ~ sw)
  expect_true("swwide" %in% names(coef(m)))
  expect_false(any(grepl("unknown", names(coef(m)))))
  bp <- m$model$blueprint
  expect_equal(bp$xlevels$sw, c("narrow", "wide"))
  expect_equal(bp$aux_xlevels$dis1$sw, c("narrow", "wide"))
  p <- countreg(Total_crashes ~ lnaadt + sw + offset(lnlength), data = w,
                family = "POISSON", method = "BFGS")
  expect_equal(unname(coef(p)), unname(ref), tolerance = 1e-4)
  # Prediction: training levels; a level not seen in training is refused
  expect_equal(predict(m, newdata = w[1:5, ]), fitted(m)[1:5],
               tolerance = 1e-8)
  nd <- w[1:2, ]
  nd$sw[1] <- "unknown"
  expect_error(predict(m, newdata = nd), "new level")
  # Level left only in a row removed for missing values
  w2 <- w
  w2$sw <- as.character(w2$sw)
  w2$sw[1] <- "unknown"
  w2$lnaadt[1] <- NA
  w2$sw <- factor(w2$sw)
  m2 <- suppressMessages(nb2(Total_crashes ~ lnaadt + sw + offset(lnlength),
                             w2))
  expect_equal(m2$model$blueprint$xlevels$sw, c("narrow", "wide"))
  # Level left unused by subsetting
  m3 <- nb2(Total_crashes ~ lnaadt + sw + offset(lnlength),
            w[w$lnaadt > 9, ])
  expect_true("swwide" %in% names(coef(m3)))
  # A factor made in the formula, and a factor by-variable of a smooth
  m4 <- nb2(Total_crashes ~ factor(speed50, levels = 0:2) +
              offset(lnlength), w)
  expect_equal(length(coef(m4)), 3)
  skip_if_not_installed("mgcv")
  w$sp <- factor(ifelse(w$speed50 == 1, "hi", "lo"),
                 levels = c("hi", "lo", "zz"))
  m5 <- nb2(Total_crashes ~ s(lnaadt, by = sp, k = 4) + sp +
              offset(lnlength), w)
  expect_false(any(grepl("zz", names(coef(m5)))))
  # Multivariate fits
  mv <- countreg_mv(cbind(Fatal_crashes, Injury_crashes) ~ lnaadt + sw +
                      offset(lnlength), data = w, family = "POISSON")
  expect_equal(sum(grepl("swwide", names(coef(mv)))), 2)
})

test_that("factor time variables of tv() with unused levels", {
  w <- wr()
  w$fY <- factor(w$Year, levels = 2015:2018)
  m <- nb2(Total_crashes ~ tv(lnaadt, time = fY) + speed50 +
             offset(lnlength), w)
  expect_equal(nrow(temporal_effects(m)), 3)
  expect_equal(predict(m, newdata = w[1:3, ]), fitted(m)[1:3],
               tolerance = 1e-8)
  nd <- w[1:2, ]
  nd$fY[1] <- "2015"
  expect_error(predict(m, newdata = nd), "not in the estimation data")
})

test_that("temporal_effects() and temporal_tests() for countreg_mv()", {
  w <- wr()
  mv <- countreg_mv(cbind(Fatal_crashes, Injury_crashes) ~
                      tv(lnaadt, time = Year) + offset(lnlength),
                    data = w, family = "POISSON")
  te <- temporal_effects(mv)
  expect_equal(names(te)[1], "outcome")
  expect_equal(te$outcome, rep(c("Fatal_crashes", "Injury_crashes"),
                               each = 3))
  expect_equal(te$period, rep(c("2016", "2017", "2018"), 2))
  # Step basis: the effects are cumulative sums of the coefficients
  cf <- coef(mv)
  b <- cf[c("Injury_crashes:lnaadt:Year2016",
            "Injury_crashes:lnaadt:Year_step2017",
            "Injury_crashes:lnaadt:Year_step2018")]
  expect_equal(te$estimate[4:6], unname(cumsum(b)), tolerance = 1e-10)
  tt <- temporal_tests(mv)
  expect_equal(names(tt)[1], "outcome")
  expect_equal(sum(tt$test == "all equal"), 2)
  expect_equal(nrow(temporal_effects(mv, term = "lnaadt")), 6)
})

test_that("tv(basis = \"smooth\") with any basis type", {
  skip_if_not_installed("mgcv")
  w <- wr()
  for (bs in c("cr", "ps", "tp")) {
    # (mgcv warns that a P-spline basis of 4 exceeds the 3 time values)
    m <- suppressWarnings(nb2(
      Total_crashes ~ tv(lnaadt, time = Year, basis = "smooth", bs = bs,
                         k = if (bs == "ps") 4 else 3) +
        speed50 + offset(lnlength), w))
    te <- temporal_effects(m)
    expect_equal(te$period, c("2016", "2017", "2018"), info = bs)
    expect_equal(te$term, rep("lnaadt", 3), info = bs)
    expect_true(all(te$std.error > 0), info = bs)
  }
  # tv() smooth in a distribution-parameter formula
  m2 <- suppressWarnings(nb2(
    Total_crashes ~ lnaadt + offset(lnlength), w,
    dis_param_formula_1 = ~ tv(1, time = Year, basis = "smooth", bs = "cr")))
  te2 <- temporal_effects(m2)
  expect_equal(unique(te2$part), "ln(alpha)")
  expect_equal(nrow(te2), 3)
  # x must be numeric
  w$sp <- factor(ifelse(w$speed50 == 1, "hi", "lo"))
  expect_error(nb2(Total_crashes ~ tv(sp, time = Year, basis = "smooth") +
                     offset(lnlength), w), "must be numeric")
})

test_that("smooths written with s() are not time-varying effects", {
  skip_if_not_installed("mgcv")
  w <- wr()
  set.seed(1)
  w$disc <- sample(1:8, nrow(w), TRUE)
  m <- nb2(Total_crashes ~ s(disc, k = 5) + lnaadt + offset(lnlength), w)
  expect_error(temporal_effects(m), "no time-varying effects")
  expect_error(temporal_tests(m), "no time-varying effects")
  w$sp <- factor(ifelse(w$speed50 == 1, "hi", "lo"))
  m2 <- nb2(Total_crashes ~ s(Year, k = 3, by = sp) + sp +
              tv(lnaadt, time = Year, basis = "smooth") + offset(lnlength),
            w)
  te <- temporal_effects(m2)
  expect_equal(te$term, rep("lnaadt", 3))
})

test_that("margEffTable() leaves out the time variable of a tv() smooth", {
  skip_if_not_installed("mgcv")
  w <- wr()
  m <- nb2(Total_crashes ~ tv(lnaadt, time = Year, basis = "smooth") +
             speed50 + offset(lnlength), w)
  me <- margEffTable(m)
  expect_setequal(me$term, c("lnaadt", "speed50"))
})

test_that("tv(1, time) in rpar_formula without a fixed intercept", {
  skip_on_cran()
  w <- wr()
  r <- suppressWarnings(countreg_rp(
    Total_crashes ~ -1 + lnaadt + offset(lnlength),
    rpar_formula = ~ tv(1, time = Year), data = w, family = "NB2",
    ndraws = 50, method = "BFGS"))
  te <- temporal_effects(r)
  expect_equal(te$term, rep("(Intercept)", 3))
  expect_equal(te$part, rep("random parameter mean", 3))
  cf <- coef(r)
  expect_equal(te$estimate[1], unname(cf["(Intercept):Mean"]))
  expect_equal(te$estimate[2], unname(cf["(Intercept):Mean"] +
                                        cf["(Intercept):Year_step2017:Mean"]))
})

test_that("tv() in an interaction is reported, not silently ignored", {
  w <- wr()
  m <- nb2(Total_crashes ~ tv(lnaadt, time = Year):speed50 + lnaadt +
             offset(lnlength), w)
  expect_error(temporal_effects(m), "interactions are not summarized")
  m2 <- nb2(Total_crashes ~ tv(lnaadt, time = Year):speed50 +
              tv(1, time = Year) + lnaadt + offset(lnlength), w)
  expect_warning(te <- temporal_effects(m2), "interactions")
  expect_equal(unique(te$term), "(Intercept)")
  expect_warning(temporal_tests(m2), "interactions")
})

test_that("margEffTable() for countreg_mv() gives effects per outcome", {
  w <- wr()
  mv <- countreg_mv(cbind(Fatal_crashes, Injury_crashes) ~ lnaadt +
                      speed50 + offset(lnlength), data = w, family = "NB2")
  me <- margEffTable(mv, measure = "ame", indicator = "discrete_change",
                     digits = 10)
  expect_equal(names(me)[1], "outcome")
  expect_equal(me$outcome, rep(c("Fatal_crashes", "Injury_crashes"),
                               each = 2))
  expect_equal(me$term, rep(c("lnaadt", "speed50"), 2))
  # Reference: finite differences of predict()
  h <- 1e-5
  p0 <- predict(mv, newdata = w)
  wp <- w
  wp$lnaadt <- wp$lnaadt + h
  ame <- colMeans((predict(mv, newdata = wp) - p0) / h)
  expect_equal(me$estimate[me$term == "lnaadt"], unname(ame),
               tolerance = 1e-4)
  w1 <- w
  w1$speed50 <- 1
  w0 <- w
  w0$speed50 <- 0
  dc <- colMeans(predict(mv, newdata = w1) - predict(mv, newdata = w0))
  expect_equal(me$estimate[me$term == "speed50"], unname(dc),
               tolerance = 1e-6)
  expect_true(all(me$std_error > 0))
  # Variables only in one outcome's equation
  mv2 <- countreg_mv(list(Fatal_crashes ~ lnaadt + offset(lnlength),
                          Injury_crashes ~ lnaadt + speed50 +
                            offset(lnlength)),
                     data = w, family = "POISSON")
  me2 <- margEffTable(mv2)
  expect_equal(me2$outcome, c("Fatal_crashes", "Injury_crashes",
                              "Injury_crashes"))
  expect_error(margEffTable(mv, se = "bootstrap"), "delta")
})
