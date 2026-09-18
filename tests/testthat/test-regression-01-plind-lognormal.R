# Error 01: validation, count support and cache safety in dplindlogn_cpp().
# Add this file; keep existing tests. Requires callr in DESCRIPTION/Suggests.

reg01_draws <- c(0.05, 0.2, 0.4, 0.6, 0.8, 0.95)
reg01_density <- function(x, mean = 2, theta = 3, sigma = 0.4, ...) {
  dplindLnorm(x, mean = mean, theta = theta, sigma = sigma,
             hdraws = reg01_draws, ...)
}

test_that("error 01: zero sigma is the Poisson-Lindley limit", {
  x <- 0:20
  for (mean in c(0.5, 2, 8)) {
    for (theta in c(0.5, 3)) {
      expect_equal(reg01_density(x, mean, theta, sigma = 0),
                   as.numeric(dplind(x, mean, theta)), tolerance = 1e-12)
    }
  }
})

test_that("error 01: missing counts propagate and unsupported counts have zero mass", {
  x <- c(NA_real_, NaN, -1, 0.5, Inf, -Inf, 0, 2)
  p <- reg01_density(x)
  lp <- reg01_density(x, log = TRUE)
  expect_true(all(is.na(p[1:2])))
  expect_true(all(is.na(lp[1:2])))
  expect_equal(p[3:6], rep(0, 4))
  expect_equal(lp[3:6], rep(-Inf, 4))
  expect_true(all(is.finite(p[7:8]) & p[7:8] > 0))
  expect_equal(exp(lp[7:8]), p[7:8], tolerance = 1e-12)
})

test_that("error 01: repeated and distinct cache keys match separate evaluations", {
  x <- c(0, 1, 5, 2, 8, 3, 1, 0)
  mu <- c(1, 2, 3, 4, 2, 1, 5, 2)
  theta <- c(1, 3, 2, 4, 1, 3, 2, 4)
  sigma <- c(0, 0.2, 0.7, 0.2, 0, 0.7, 0.4, 0.2)
  expected <- vapply(seq_along(x), function(i) {
    reg01_density(x[i], mu[i], theta[i], sigma[i])
  }, numeric(1))
  actual <- reg01_density(x, mu, theta, sigma)
  expect_equal(actual, expected, tolerance = 1e-12)
  order <- c(8, 3, 1, 7, 5, 2, 6, 4)
  expect_equal(reg01_density(x[order], mu[order], theta[order], sigma[order]),
               actual[order], tolerance = 1e-12)
})

test_that("error 01: scalar parameters and explicitly repeated parameters agree", {
  x <- 0:10
  expect_equal(reg01_density(x),
               reg01_density(x, rep(2, 11), rep(3, 11), rep(0.4, 11)),
               tolerance = 1e-12)
})

test_that("error 01: valid simulated probabilities normalize and log output agrees", {
  p <- reg01_density(0:500)
  expect_true(all(is.finite(p) & p >= 0 & p <= 1))
  expect_equal(sum(p), 1, tolerance = 1e-9)
  expect_equal(reg01_density(0:30, log = TRUE), log(p[1:31]),
               tolerance = 1e-12)
})

test_that("error 01: invalid parameter lengths are rejected", {
  expect_error(reg01_density(0:3, mean = c(1, 2)))
  expect_error(reg01_density(0:3, theta = c(1, 2)))
  expect_error(reg01_density(0:3, sigma = c(0.1, 0.2)))
})

test_that("error 01: invalid means and shapes are rejected before evaluation", {
  for (bad in c(0, -1, NA_real_, NaN, Inf, -Inf)) {
    expect_error(reg01_density(0:1, mean = bad))
    expect_error(reg01_density(0:1, theta = bad))
  }
})

test_that("error 01: invalid scalar standard deviations are rejected", {
  for (bad in c(-0.1, NA_real_, NaN, Inf, -Inf)) {
    expect_error(reg01_density(0:1, sigma = bad))
  }
})

test_that("error 01: empty and nonfinite normal draws are rejected", {
  kernel <- getFromNamespace("dplindlogn_cpp", "flexCountReg")
  expect_error(kernel(0:1, 2, 3, 0.4, numeric()))
  for (bad in c(NA_real_, NaN, Inf, -Inf)) {
    expect_error(kernel(0:1, 2, 3, 0.4, c(0, bad)))
  }
})

test_that("error 01: vectors over the old OpenMP threshold remain correct", {
  x <- rep(0:5, 201)  # 1,206 observations; old threshold was 1,000.
  expected <- rep(reg01_density(0:5), 201)
  expect_equal(reg01_density(x), expected, tolerance = 1e-12)
})

test_that("error 01: the former missing-sigma crash raises a recoverable R error", {
  # Load the CURRENT compiled DLL in a child, including during load_all/covr.
  # Loading an independently installed package here could test stale code.
  dll <- getLoadedDLLs()[["flexCountReg"]]
  expect_false(is.null(dll))
  if (is.null(dll)) return(invisible(NULL))

  result <- callr::r(function(dll_path) {
    loadNamespace("Rcpp")
    dll <- dyn.load(dll_path)
    native <- getNativeSymbolInfo("_flexCountReg_dplindlogn_cpp", PACKAGE = dll)
    caught <- tryCatch({
      .Call(native, c(0, 1), 2, 3, c(1, NA_real_), c(-1, 0, 1))
      list(error = FALSE, message = "No R error was raised")
    }, error = function(e) list(error = TRUE, message = conditionMessage(e)))
    # Also prove that a subsequent valid call works in the same process.
    valid <- .Call(native, c(0, 1), 2, 3, 0.4, c(-1, 0, 1))
    list(caught = caught, valid = valid)
  }, args = list(dll_path = dll[["path"]]), libpath = .libPaths(), timeout = 60)

  expect_true(result$caught$error)
  expect_match(result$caught$message, "sigma|standard deviation", ignore.case = TRUE)
  expect_length(result$valid, 2)
  expect_true(all(is.finite(result$valid) & result$valid > 0))
})
