# Further prediction paths: probability matrices, objects without a stored
# blueprint, legacy fits, and Empirical Bayes with heterogeneity formulas.

test_that("type = 'prob' defaults to counts up to the largest training count", {
  m <- countreg(Total_crashes ~ lnaadt, washington_roads, family = "NB2",
                method = "BFGS")
  pm <- predict(m, type = "prob")
  expect_equal(ncol(pm), max(washington_roads$Total_crashes) + 1)
  expect_true(all(rowSums(pm) <= 1 + 1e-10))
  pm5 <- predict(m, type = "prob", max_count = 5)
  expect_equal(ncol(pm5), 6)
  expect_equal(pm5, pm[, 1:6], ignore_attr = TRUE)
})

test_that("types other than 'response' are for countreg() models only", {
  re <- renb(Total_crashes ~ lnaadt, washington_roads, group_var = "ID",
             method = "BFGS")
  expect_error(predict(re, type = "link"),
               "only available for countreg\\(\\) models")
})

test_that("the blueprint is rebuilt from the stored data when missing", {
  d <- washington_roads
  d$fac <- factor(ifelse(d$speed50 == 1, "fast", "slow"))
  m <- countreg(Total_crashes ~ lnaadt + fac, d, offset = "lnlength",
                method = "BFGS")
  ref <- predict(m, newdata = d[1:20, ])
  old <- m
  old$model$blueprint <- NULL
  expect_equal(predict(old, newdata = d[1:20, ]), ref)
  old$data <- NULL
  expect_error(predict(old, newdata = d[1:20, ]), "no stored blueprint")
})

test_that("predictions fail clearly when coefficients are missing", {
  m <- countreg(Total_crashes ~ lnaadt, washington_roads,
                underreport_formula = ~ speed50, method = "BFGS")
  bad <- m
  names(bad$model$estimate)[names(bad$model$estimate) ==
                              "Underreporting:speed50"] <- "renamed"
  expect_error(predict(bad, newdata = washington_roads[1:5, ]),
               "Coefficient\\(s\\) not found in the model: Underreporting:speed50")
})

test_that("probit underreporting fits from older versions keep their sign", {
  m <- countreg(Total_crashes ~ lnaadt, washington_roads,
                underreport_formula = ~ speed50,
                underreport_family = "probit", method = "BFGS")
  d <- washington_roads[1:10, ]
  p_new <- predict(m, newdata = d, type = "prob_reported")
  old <- m
  old$model$underreport_sign <- NULL
  # Older fits used P(reported) = 1 - pnorm(z'g)
  expect_equal(predict(old, newdata = d, type = "prob_reported"), 1 - p_new)
})

test_that("Empirical Bayes predictions work with heterogeneity formulas", {
  d <- washington_roads
  d$AADT10kplus <- as.numeric(d$AADT > 10000)
  rp <- countreg_rp(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
                    het_mean_formula = ~ AADT10kplus,
                    het_var_formula = ~ AADT10kplus,
                    data = d, family = "NB2", rpardists = c(speed50 = "n"),
                    ndraws = 50, offset = "lnlength", group_var = "ID",
                    method = "BFGS", max.iters = 200)
  eb <- predict(rp, newdata = d, method = "Individual")
  expect_length(eb, nrow(d))
  expect_true(all(is.finite(eb) & eb > 0))
  expect_true(all(attr(eb, "eb_sd")[, "speed50"] >= 0))
  # With no random variation left, EB equals the marginal prediction
  rp0 <- rp
  rp0$model$estimate["speed50:St.Dev"] <- 1e-8
  expect_equal(as.numeric(predict(rp0, newdata = d, method = "Individual")),
               predict(rp0, newdata = d, method = "Simulated"),
               tolerance = 1e-6)
  # A panel with all outcomes missing gets no EB prediction
  nd <- d[d$ID %in% unique(d$ID)[1:2], ]
  nd$Total_crashes <- NA
  p <- predict(rp, newdata = nd, method = "Individual")
  expect_true(all(is.na(p)))
  expect_true(all(is.na(attr(p, "eb_coefficients"))))
})
