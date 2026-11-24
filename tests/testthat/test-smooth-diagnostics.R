# Smooth types (by, re, fs, mrf) and the smooth diagnostics.

sd_data <- function() {
  data("washington_roads", package = "flexCountReg", envir = environment())
  w <- washington_roads
  set.seed(1)
  w$corr <- factor(sample(letters[1:6], nrow(w), TRUE))
  w$sp50 <- factor(w$speed50)
  w$reg <- factor(sample(paste0("r", 1:9), nrow(w), TRUE),
                  levels = paste0("r", 1:9))
  w
}

grid_nb <- function() {
  nb <- lapply(1:9, function(i) {
    r <- (i - 1) %/% 3
    c <- (i - 1) %% 3
    c(if (r > 0) i - 3, if (r < 2) i + 3, if (c > 0) i - 1, if (c < 2) i + 1)
  })
  stats::setNames(nb, paste0("r", 1:9))
}

test_that("by, re, fs and mrf smooths reproduce mgcv", {
  skip_if_not_installed("mgcv")
  w <- sd_data()
  nb <- grid_nb()
  cases <- list(
    list(Total_crashes ~ s(lnaadt, by = lnlength, k = 6) + offset(lnlength),
         1),
    list(Total_crashes ~ s(lnaadt, by = sp50, k = 6) + sp50 +
           offset(lnlength), c(1, 1)),
    list(Total_crashes ~ s(lnaadt, k = 6) + s(corr, bs = "re") +
           offset(lnlength), c(1, 1)),
    list(Total_crashes ~ s(lnaadt, corr, bs = "fs", k = 5) +
           offset(lnlength), c(1, 1, 1)),
    list(Total_crashes ~ s(lnaadt, k = 5) +
           s(reg, bs = "mrf", xt = list(nb = nb)) + offset(lnlength),
         c(1, 1)))
  rows <- c(5, 300, 1200)
  for (cs in cases) {
    f <- cs[[1]]
    environment(f) <- environment()
    g <- mgcv::gam(f, data = w, family = stats::poisson(), sp = cs[[2]])
    a <- countreg(f, data = w, family = "POISSON", sp = cs[[2]])
    b <- names(coef(a))
    expect_equal(coef(a), coef(g)[b], tolerance = 1e-6)
    expect_equal(predict(a, newdata = w[rows, ]), unname(fitted(g)[rows]),
                 tolerance = 1e-6)
  }
})

test_that("smooth_estimates and plots handle every smooth type", {
  skip_if_not_installed("mgcv")
  w <- sd_data()
  m <- countreg(Total_crashes ~ s(lnaadt, by = sp50, k = 6) + sp50 +
                  s(corr, bs = "re") + offset(lnlength), data = w,
                family = "POISSON", sp = c(1, 1, 1))
  e <- smooth_estimates(m, n = 15, nsim = 100, seed = 1)
  expect_setequal(unique(e$term),
                  c("s(lnaadt):sp500", "s(lnaadt):sp501", "s(corr)"))
  expect_identical(sum(e$term == "s(corr)"), 6L)
  expect_true(all(e$sp50[e$term == "s(lnaadt):sp501"] == "1"))
  # A numeric `by`: the varying coefficient, evaluated with by = 1
  mn <- countreg(Total_crashes ~ s(lnaadt, by = lnlength, k = 6) +
                   offset(lnlength), data = w, family = "POISSON", sp = 1)
  en <- smooth_estimates(mn, n = 10, nsim = 50, seed = 1)
  expect_true(all(en$lnlength == 1))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_silent(plot_smooths(m, nsim = 50))
  fs <- countreg(Total_crashes ~ s(lnaadt, corr, bs = "fs", k = 5) +
                   offset(lnlength), data = w, family = "POISSON",
                 sp = c(1, 1, 1))
  expect_silent(plot_smooths(fs, nsim = 50))
  te <- countreg(Total_crashes ~ te(lnaadt, lnlength, k = c(4, 4)),
                 data = w, family = "POISSON", sp = c(1, 1))
  g <- plot_surface(te)
  expect_true(all(c("lnaadt", "lnlength", "fit") %in% names(g)))
})

test_that("derivatives agree with differences of the smooth", {
  skip_if_not_installed("mgcv")
  w <- sd_data()
  m <- countreg(Total_crashes ~ s(lnaadt, k = 8) + offset(lnlength),
                data = w, family = "NB2")
  e <- smooth_estimates(m, n = 400, nsim = 50, seed = 1)
  d <- smooth_derivatives(m, n = 400, nsim = 50, seed = 1)
  num <- diff(e$fit) / diff(e$lnaadt)
  mid <- (d$fit[-1] + d$fit[-nrow(d)]) / 2
  expect_lt(max(abs(num - mid)), 1e-3 * max(abs(num)) + 1e-4)
  expect_true(all(d$supr - d$fit >= d$upr - d$fit))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_silent(plot_derivatives(m, nsim = 50))
})

test_that("smooth_check flags a basis that is too small", {
  skip_if_not_installed("mgcv")
  set.seed(2)
  d <- data.frame(x = stats::runif(1500, -2, 2))
  d$y <- stats::rpois(1500, exp(1 + sin(3 * d$x)))
  small <- countreg(y ~ s(x, k = 4), data = d, family = "POISSON")
  big <- countreg(y ~ s(x, k = 20), data = d, family = "POISSON")
  cs <- smooth_check(small, nperm = 200)
  cb <- smooth_check(big, nperm = 200)
  expect_identical(cs$residuals, "quantile")
  expect_lt(cs$k_index, 0.9)
  expect_lt(cs$p_value, 0.01)
  expect_gt(cb$k_index, 0.95)
  expect_gt(cb$p_value, 0.05)
})

test_that("concurvity and calibration tables", {
  skip_if_not_installed("mgcv")
  w <- sd_data()
  m <- countreg(Total_crashes ~ s(lnaadt, k = 6) + s(lnlength, k = 6) +
                  offset(lnlength), data = w, family = "NB2")
  cc <- concurvity_fcr(m)
  expect_identical(cc$term, c("s(lnaadt)", "s(lnlength)"))
  expect_true(all(cc$observed >= 0 & cc$observed <= 1))
  expect_true(all(cc$worst >= cc$observed - 1e-8 & cc$worst <= 1))
  cb <- calibration_by(m, "lnaadt", groups = 4)
  expect_identical(sum(cb$n), nrow(w))
  expect_equal(sum(cb$observed_mean * cb$n) / nrow(w),
               mean(w$Total_crashes))
  expect_true(all(cb$predicted_p0 > 0 & cb$predicted_p0 < 1))
  expect_error(calibration_by(m, "nope"), "numeric column")
})
