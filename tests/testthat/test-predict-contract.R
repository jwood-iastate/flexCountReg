# Prediction contract: one prediction per row of newdata, NA for incomplete
# rows, training factor levels, offsets, prediction types, and Empirical Bayes
# predictions for random-parameter models.

make_data <- function() {
  d <- washington_roads
  d$fac <- factor(ifelse(d$speed50 == 1, "fast",
                         ifelse(d$AADT > 10000, "busy", "base")))
  d
}

test_that("predict keeps rows, returns NA for incomplete rows", {
  d <- make_data()
  m <- countreg(Total_crashes ~ lnaadt + fac, d, offset = "lnlength",
                dis_param_formula_1 = ~speed50, method = "BFGS")
  nd <- d[1:4, ]
  nd$lnaadt[2] <- NA
  p <- predict(m, newdata = nd)
  expect_length(p, 4)
  expect_true(is.na(p[2]))
  expect_equal(p[-2], predict(m, newdata = d[c(1, 3, 4), ]))

  nd2 <- d[1:4, ]
  nd2$speed50[3] <- NA   # only used by the dispersion formula
  expect_length(predict(m, newdata = nd2), 4)
})

test_that("predict uses the training factor levels", {
  d <- make_data()
  m <- countreg(Total_crashes ~ lnaadt + fac, d, method = "BFGS")
  i <- which(d$fac == "busy")[1]
  one <- d[i, , drop = FALSE]
  one$fac <- factor("busy")
  expect_equal(predict(m, newdata = one), predict(m, newdata = d)[i])
  # Same levels, different order
  re <- d[i, , drop = FALSE]
  re$fac <- factor("busy", levels = c("fast", "busy", "base"))
  expect_equal(predict(m, newdata = re), predict(m, newdata = d)[i])
})

test_that("offsets are required and used at prediction time", {
  d <- make_data()
  a <- countreg(Total_crashes ~ lnaadt + offset(lnlength), d, method = "BFGS")
  b <- countreg(Total_crashes ~ lnaadt, d, offset = "lnlength",
                method = "BFGS")
  expect_equal(predict(a, newdata = d), predict(b, newdata = d),
               tolerance = 1e-5)
  nd <- d[, setdiff(names(d), "lnlength")]
  expect_error(predict(a, newdata = nd), "lnlength")
  expect_error(predict(b, newdata = nd), "lnlength")
})

test_that("prediction types are consistent", {
  d <- make_data()
  m <- countreg(Total_crashes ~ lnaadt + lnlength, d,
                underreport_formula = ~ShouldWidth04, method = "BFGS")
  resp <- predict(m, newdata = d)
  tru <- predict(m, newdata = d, type = "true")
  pr <- predict(m, newdata = d, type = "prob_reported")
  lk <- predict(m, newdata = d, type = "link")
  expect_equal(resp, tru * pr)
  expect_equal(tru, exp(lk))
  expect_true(all(pr > 0 & pr < 1))

  m0 <- countreg(Total_crashes ~ lnaadt + lnlength, d, method = "BFGS")
  expect_equal(predict(m0, newdata = d, type = "true"),
               predict(m0, newdata = d))
  expect_true(all(predict(m0, newdata = d, type = "prob_reported") == 1))
})

test_that("PLL and PLN predictions include the lognormal mean factor", {
  d <- make_data()
  m <- countreg(Total_crashes ~ lnaadt + lnlength, d, family = "PLL",
                method = "BFGS", ndraws = 200)
  sig <- exp(m$model$estimate["ln(sigma)"])
  expect_equal(predict(m, newdata = d),
               exp(predict(m, newdata = d, type = "link")) * exp(sig^2 / 2))
})

test_that("the data argument is deprecated but still works", {
  d <- make_data()
  m <- countreg(Total_crashes ~ lnaadt, d, method = "BFGS")
  expect_warning(p <- predict(m, data = d[1:3, ]), "deprecated")
  expect_equal(p, predict(m, newdata = d[1:3, ]))
})

test_that("countreg.rp predictions use the offset", {
  d <- make_data()
  rp <- countreg.rp(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
                    data = d, family = "NB2", rpardists = c(speed50 = "n"),
                    ndraws = 100, offset = "lnlength", method = "BFGS")
  ps <- predict(rp, newdata = d, method = "Simulated")
  expect_equal(mean(ps), mean(d$Total_crashes), tolerance = 0.1)
})

test_that("Individual (Empirical Bayes) predictions follow Wood and Gayah", {
  d <- make_data()
  rp <- countreg.rp(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
                    data = d, family = "NB2", rpardists = c(speed50 = "n"),
                    ndraws = 100, offset = "lnlength", panel_id = "ID",
                    method = "BFGS")
  ps <- predict(rp, newdata = d, method = "Simulated")
  eb <- predict(rp, newdata = d, method = "Individual")
  expect_length(eb, nrow(d))
  expect_false(isTRUE(all.equal(as.numeric(eb), ps)))
  expect_equal(predict(rp, newdata = d, method = "EB"), eb)

  # Conditional coefficients are constant within a panel
  b <- attr(eb, "eb_coefficients")[, "speed50"]
  expect_true(all(tapply(b, d$ID, function(v) diff(range(v))) < 1e-12))
  expect_true(all(attr(eb, "eb_sd")[, "speed50"] >= 0))

  # eb_type = "mean" differs from the coefficient form
  ebm <- predict(rp, newdata = d, method = "Individual", eb_type = "mean")
  expect_length(ebm, nrow(d))

  # With no heterogeneity (sigma -> 0), EB equals the marginal prediction
  rp0 <- rp
  i <- grep("St.Dev", names(rp0$model$estimate))
  rp0$model$estimate[i] <- 1e-8
  expect_equal(as.numeric(predict(rp0, newdata = d, method = "Individual")),
               predict(rp0, newdata = d, method = "Simulated"),
               tolerance = 1e-8)

  # Requires the outcome and the panel id; NA rows stay NA
  nd <- d[d$ID %in% unique(d$ID)[1:3], ]
  expect_error(predict(rp, newdata = nd[, setdiff(names(nd), "Total_crashes")],
                       method = "Individual"), "requires the outcome")
  expect_error(predict(rp, newdata = nd[, setdiff(names(nd), "ID")],
                       method = "Individual"), "panel id")
  nd$lnaadt[2] <- NA
  p <- predict(rp, newdata = nd, method = "Individual")
  expect_true(is.na(p[2]))
  expect_false(anyNA(p[-2]))
  expect_error(predict(rp, newdata = nd, method = "Bogus"))
})
