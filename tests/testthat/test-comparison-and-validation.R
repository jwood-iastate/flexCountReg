# regCompTest() options and input validation of the panel estimators.

d <- washington_roads

test_that("regCompTest() base models of other families and their LR tests", {
  nb2 <- countreg(Total_crashes ~ lnaadt + speed50, d, offset = "lnlength",
                  method = "BFGS")
  # Intercept-only NB2 base model, fitted on the same rows and offset
  r <- regCompTest(nb2, basemodel = "NB2")
  b0 <- countreg(Total_crashes ~ 1, d, offset = "lnlength", method = "BFGS")
  expect_equal(r$LLbase, b0$model$maximum, tolerance = 1e-6)
  expect_equal(r$LRdof, 2)
  expect_equal(r$LR, 2 * (nb2$model$maximum - b0$model$maximum),
               tolerance = 1e-6)
  # Base model with the same variables: nothing to test (same parameters)
  rv <- regCompTest(nb2, basemodel = "NB2", variables = TRUE)
  expect_equal(rv$LLbase, nb2$model$maximum, tolerance = 1e-6)
  expect_true(is.na(rv$LR))
  expect_output(regCompTest(nb2, print = TRUE))
  expect_error(regCompTest(lm(Total_crashes ~ lnaadt, d)),
               "not a flexCountReg model")
})

test_that("regCompTest() warns for comparisons that are not nested", {
  a <- countreg(Total_crashes ~ lnaadt, d, family = "NB2", method = "BFGS")
  b <- countreg(Total_crashes ~ lnaadt, d, family = "NB1", method = "BFGS")
  expect_warning(r <- regCompTest(a, model2 = b), "not nested")
  expect_true(is.na(r$LR))
  # A badly converged larger model has a lower log-likelihood
  p <- countreg(Total_crashes ~ lnaadt, d, family = "Poisson",
                method = "BFGS")
  bad <- suppressWarnings(countreg(Total_crashes ~ lnaadt, d, family = "NB2",
                                   method = "NM", max.iters = 1,
                                   start.vals = c(-1, 0, 2)))
  expect_warning(r2 <- regCompTest(p, model2 = bad), "higher")
  expect_equal(r2$LR_pvalue, 1)
})

test_that("negmultinom() validates its inputs", {
  expect_error(negmultinom(Total_crashes ~ lnaadt, "ID", d,
                           family = c("NM1", "NM2")), "single string")
  expect_error(negmultinom(Total_crashes ~ lnaadt, "ID", d, method = "foo"),
               "Unsupported optimization method")
  m <- negmultinom(Total_crashes ~ lnaadt, "ID", d, method = "Nelder-Mead",
                   max.iters = 300)
  expect_equal(m$model$method, "NM")
  expect_error(negmultinom(Total_crashes ~ lnaadt, "ID", d,
                           start.vals = c(0, 1)), "must have length 3")
  expect_error(negmultinom(Total_crashes ~ lnaadt, "ID", d,
                           dis_param_formula_1 = y ~ x), "one-sided")
  expect_error(negmultinom(Total_crashes ~ lnaadt, "ID", as.matrix(d)),
               "data frame")
  expect_error(negmultinom(Total_crashes ~ lnaadt, "nope", d), "not found")
  expect_error(negmultinom(Total_crashes ~ lnaadt, "ID", d, offset = 1),
               "character")
  expect_error(negmultinom(Total_crashes ~ lnaadt, "ID", d, offset = "nope"),
               "not found")
  dn <- d
  dn$Total_crashes[1] <- -1
  expect_error(negmultinom(Total_crashes ~ lnaadt, "ID", dn), "nonnegative")
  dn$Total_crashes[1] <- 0.5
  # The starting Poisson fit also warns about the non-integer count
  ws <- character()
  withCallingHandlers(
    negmultinom(Total_crashes ~ lnaadt, "ID", dn, max.iters = 100),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("not integer", ws)))
  expect_error(predict(m, newdata = d, method = "Exact"), "Standard")
})

test_that("negmultinom() bootstrap with a dispersion formula", {
  set.seed(6)
  m <- negmultinom(Total_crashes ~ lnaadt, "ID", d,
                   dis_param_formula_1 = ~ speed50, stderr = "boot",
                   bootstraps = 5)
  expect_equal(dim(vcov(m)), rep(length(coef(m)), 2))
  expect_true(all(is.finite(sqrt(diag(vcov(m))))))
})

test_that("renb(), poisLind.re() and repois() validate their inputs", {
  for (fn in list(renb, poisLind.re, repois)) {
    expect_error(fn(Total_crashes ~ lnaadt, data = d), "group_var")
    expect_error(fn(Total_crashes ~ lnaadt, group_var = "ID", data = d,
                    method = "foo"), "`method` must be")
    expect_error(fn(Total_crashes ~ lnaadt, group_var = "ID", data = d,
                    bootstraps = 1), "at least 2")
  }
  dn <- d
  dn$Total_crashes[1] <- -1
  expect_error(repois(Total_crashes ~ lnaadt, group_var = "ID", data = dn),
               "nonnegative")
  # Every optimizer accepted by the check reaches the same maximum
  ll <- vapply(c("BFGS", "NR", "NM"), function(meth) {
    repois(Total_crashes ~ lnaadt, group_var = "ID", data = d,
           method = meth)$model$maximum
  }, numeric(1))
  expect_lt(diff(range(ll)), 1e-3)
  expect_s3_class(suppressWarnings(
    repois(Total_crashes ~ lnaadt, group_var = "ID", data = d,
           method = "BFGSR", max.iters = 20)), "repois")
})
