# countreg_rp(): panel ids, correlated random parameters, distribution
# parameter formulas, standard-error options and argument validation.

d <- washington_roads
fit_rp <- function(...) {
  countreg_rp(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
              data = d, ndraws = 50, method = "BFGS", ...)
}

test_that("a panel-id vector is equivalent to a panel-id column", {
  a <- fit_rp(group_var = "ID")
  b <- fit_rp(group_var = d$ID)
  expect_equal(coef(a), coef(b))
  expect_equal(b$model$n_groups, length(unique(d$ID)))
  expect_error(fit_rp(group_var = d$ID[1:10]), "one value per row")
  # Empirical Bayes predictions match panels by the stored id column
  # (.fcr_panel_id), never by row position
  eb <- predict(b, method = "Individual")
  expect_length(eb, nrow(d))
  expect_equal(as.vector(predict(a, newdata = d, method = "Individual")),
               as.vector(eb), tolerance = 1e-6)
  rev_rows <- rev(seq_len(nrow(d)))
  expect_equal(as.vector(predict(b, newdata = b$data[rev_rows, ],
                                 method = "Individual")),
               as.vector(eb)[rev_rows], tolerance = 1e-12)
  expect_error(predict(b, newdata = d[1:10, ], method = "Individual"),
               "group_var")
})

test_that("correlated random parameters: Cholesky terms and Exact = Simulated", {
  expect_warning(
    m <- countreg_rp(Total_crashes ~ lnlength,
                     rpar_formula = ~ -1 + lnaadt + speed50, data = d,
                     rpardists = c(lnaadt = "t", speed50 = "n"),
                     correlated = TRUE, ndraws = 50, method = "BFGS",
                     max.iters = 200),
    "only normal")
  expect_equal(unname(m$model$rpardists), c("n", "n"))
  expect_equal(sum(grepl("^Chol:", names(coef(m)))), 3)
  m$model$ndraws <- 20000
  ps <- predict(m, newdata = d[1:50, ], method = "Simulated")
  pe <- predict(m, newdata = d[1:50, ], method = "Exact")
  expect_lt(max(abs(ps / pe - 1)), 0.01)
})

test_that("distribution-parameter formulas name their coefficients", {
  m1 <- fit_rp(dis_param_formula_1 = ~ ShouldWidth04)
  expect_true(all(c("ln(alpha):(Intercept)", "ln(alpha):ShouldWidth04") %in%
                    names(coef(m1))))
  m2 <- countreg_rp(Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50,
                    data = d, family = "NBP", dis_param_formula_2 = ~ speed50,
                    ndraws = 30, method = "BFGS", max.iters = 100)
  expect_true("ln(p):speed50" %in% names(coef(m2)))
  expect_length(predict(m2, newdata = d[1:5, ]), 5)
})

test_that("invalid specifications give clear errors or warnings", {
  expect_error(fit_rp(family = "NB3"), "Unknown `family`")
  expect_error(fit_rp(dis_param_formula_2 = ~ speed50),
               "no second distribution parameter")
  expect_error(countreg_rp(Total_crashes ~ lnaadt, rpar_formula = ~ 1,
                           data = d, ndraws = 20),
               "no random parameters")
  dw <- d
  dw$w <- rep(1:2, length.out = nrow(dw))
  expect_error(countreg_rp(Total_crashes ~ lnaadt,
                           rpar_formula = ~ -1 + speed50, data = dw,
                           weights = "w", group_var = "ID", ndraws = 20,
                           method = "BFGS", max.iters = 20),
               "vary within panels")
})

test_that("rpardists shorthands", {
  rp2 <- ~ -1 + speed50 + ShouldWidth04
  a <- countreg_rp(Total_crashes ~ lnaadt, rpar_formula = rp2, data = d,
                   rpardists = "u", ndraws = 20, method = "BFGS",
                   max.iters = 30)
  expect_equal(a$model$rpardists, c(speed50 = "u", ShouldWidth04 = "u"))
  b <- countreg_rp(Total_crashes ~ lnaadt, rpar_formula = rp2, data = d,
                   rpardists = c("n", "t"), ndraws = 20, method = "BFGS",
                   max.iters = 30)
  expect_equal(b$model$rpardists, c(speed50 = "n", ShouldWidth04 = "t"))
  expect_error(countreg_rp(Total_crashes ~ lnaadt, rpar_formula = rp2,
                           data = d, rpardists = c("n", "t", "u")),
               "one entry per random")
})

test_that("robust and bootstrap standard errors", {
  r <- fit_rp(group_var = "ID", stderr = "robust")
  expect_equal(r$model$vcov_type, "robust")
  expect_true(all(is.finite(sqrt(diag(vcov(r))))))
  set.seed(4)
  b <- fit_rp(group_var = "ID", stderr = "boot", bootstraps = 4,
              max.iters = 50)
  expect_equal(b$model$vcov_type, "bootstrap")
  expect_equal(b$model$n_boot_ok, 4)
  expect_identical(rownames(vcov(b)), names(coef(b)))
})

test_that("legacy heterogeneity-in-means fits predict with their own form", {
  m <- countreg_rp(Total_crashes ~ lnlength, rpar_formula = ~ -1 + speed50,
                   het_mean_formula = ~ ShouldWidth04, data = d,
                   ndraws = 50, method = "BFGS")
  old <- m
  old$model$draw_version <- NULL
  old$model$rpar <- NULL
  cf <- old$model$estimate
  names(cf)[names(cf) == "HetMean:speed50:ShouldWidth04"] <-
    "HetMean:ShouldWidth04"
  old$model$estimate <- cf
  nd <- d[1:6, ]
  # Old form: mu_i = mu * exp(z delta) for every random parameter
  mu_i <- cf["speed50:Mean"] * exp(cf["HetMean:ShouldWidth04"] *
                                     nd$ShouldWidth04)
  manual <- exp(cf["(Intercept)"] + cf["lnlength"] * nd$lnlength +
                  nd$speed50 * mu_i +
                  0.5 * nd$speed50^2 * cf["speed50:St.Dev"]^2)
  expect_equal(predict(old, newdata = nd, method = "Exact"), unname(manual))
})
