# Regression tests for renb() and poisLind.re(): shared data preparation,
# group handling, cluster bootstrap and the RENB likelihood.

test_that("renb log-likelihood is accurate for very large means", {
  x <- c(5, 1e5, 1e170)
  exact <- vapply(x, function(v) sum(log(v + 0:2)), numeric(1))
  expect_equal(.fcr_lgamma_diff(x, 3), exact)
  # A panel's probabilities sum to one
  p <- vapply(0:2000, function(y) exp(renb_ll(y, 50, 3, 2, 1)), numeric(1))
  expect_equal(sum(p), 1, tolerance = 1e-3)
})

test_that("renb with BFGS reaches the optimum (no runaway coefficients)", {
  m <- renb(Total_crashes ~ lnaadt + speed50, group_var = "ID",
            data = washington_roads, offset = "lnlength", method = "BFGS")
  expect_lt(abs(m$model$estimate["lnaadt"] - 1.12), 0.05)
  expect_gt(m$model$maximum, -1080)
  expect_lt(m$model$maximum, -1060)
})

test_that("group_var is used even when data has a panel_id column", {
  d <- washington_roads
  d$panel_id <- seq_len(nrow(d))
  for (fn in list(renb, poisLind.re)) {
    a <- fn(Total_crashes ~ lnaadt, group_var = "ID", data = d,
            method = "BFGS")
    b <- fn(Total_crashes ~ lnaadt, group_var = "ID",
            data = washington_roads, method = "BFGS")
    expect_equal(a$model$estimate, b$model$estimate)
    expect_equal(a$model$n_groups, length(unique(d$ID)))
  }
})

test_that("missing values, offsets and group_var are validated", {
  d <- washington_roads
  d$lnaadt[1:3] <- NA
  m <- suppressMessages(poisLind.re(Total_crashes ~ lnaadt, group_var = "ID",
                                    data = d, method = "BFGS"))
  expect_equal(m$model$nobs, nrow(d) - 3)
  expect_error(renb(Total_crashes ~ lnaadt + offset(lnlength),
                    group_var = "ID", data = washington_roads,
                    offset = "lnlength"), "not both")
  expect_error(renb(Total_crashes ~ lnaadt, data = washington_roads),
               "group_var")
  a <- renb(Total_crashes ~ lnaadt + offset(lnlength), group_var = "ID",
            data = washington_roads, method = "BFGS")
  b <- renb(Total_crashes ~ lnaadt, group_var = "ID", offset = "lnlength",
            data = washington_roads, method = "BFGS")
  expect_equal(a$model$estimate, b$model$estimate, tolerance = 1e-6)
})

test_that("bootstrap resamples whole clusters and SEs are named", {
  set.seed(8)
  prep <- prepare_model_data(Total_crashes ~ lnaadt, washington_roads,
                             group_var = "ID")
  sub <- .fcr_cluster_subset(prep, c(1, 1, 2))
  expect_equal(nlevels(sub$group), 3)
  expect_equal(length(sub$y), 2 * sum(prep$group == levels(prep$group)[1]) +
                 sum(prep$group == levels(prep$group)[2]))
  m <- poisLind.re(Total_crashes ~ lnaadt, group_var = "ID",
                   data = washington_roads, method = "BFGS", bootstraps = 10)
  expect_equal(m$model$vcov_type, "bootstrap")
  expect_identical(names(m$model$se), names(m$model$estimate))
  expect_equal(m$model$n_boot_ok, 10)
})
