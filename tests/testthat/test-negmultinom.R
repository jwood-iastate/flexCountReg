# Tests for negmultinom() (NM1, NM2, NM-P)

sim_nm <- function(G = 300, Tn = 4, beta = c(0.2, 0.4, 0.3),
                   lnalpha = log(0.5), p = 2, seed = 42) {
  set.seed(seed)
  id <- rep(seq_len(G), each = Tn)
  x1 <- stats::rnorm(G * Tn)
  z <- rep(stats::rbinom(G, 1, 0.5), each = Tn)
  lnL <- log(stats::runif(G * Tn, 0.5, 2))
  mu <- exp(beta[1] + beta[2] * x1 + beta[3] * z + lnL)
  M <- tapply(mu, id, sum)
  r <- M^(2 - p) / exp(lnalpha)
  U <- stats::rgamma(G, shape = r, rate = r)
  data.frame(id = id, x1 = x1, z = z, lnL = lnL,
             y = stats::rpois(G * Tn, mu * U[id]))
}

nm_data <- sim_nm()
ns <- asNamespace("flexCountReg")

test_that("analytic gradient matches numerical gradient for all families", {
  for (fam in c("NM2", "NM1", "NMP")) {
    prep <- ns$.nm_prepare(y ~ x1 + z, "id", nm_data, "lnL",
                           ~ z, if (fam == "NMP") ~ z else NULL, fam)
    par <- c(0.1, 0.3, 0.2, -0.5, 0.2, if (fam == "NMP") c(1.4, 0.1))
    ga <- colSums(ns$.nm_grad(par, prep, fam))
    gn <- vapply(seq_along(par), function(j) {
      e <- replace(numeric(length(par)), j, 1e-6)
      (sum(ns$.nm_loglik(par + e, prep, fam)) -
         sum(ns$.nm_loglik(par - e, prep, fam))) / 2e-6
    }, numeric(1))
    expect_equal(unname(ga), gn, tolerance = 1e-5)
  }
})

test_that("NM2 log-likelihood equals the repois log-likelihood", {
  mod <- negmultinom(y ~ x1 + z, group_var = "id", data = nm_data,
                     offset = "lnL")
  pr <- ns$.repois_prepare_inputs(y ~ x1 + z, "id", nm_data, "lnL")
  ll_re <- ns$.repois_loglik(mod$model$estimate, pr$y, pr$X,
                             pr$offset_vec, pr$cluster_list)
  expect_equal(mod$model$maximum, ll_re, tolerance = 1e-8,
               ignore_attr = TRUE)
})

test_that("NM-P likelihood reduces to NM2 at p = 2 and NM1 at p = 1", {
  prep_p <- ns$.nm_prepare(y ~ x1 + z, "id", nm_data, NULL, NULL, NULL, "NMP")
  prep_2 <- ns$.nm_prepare(y ~ x1 + z, "id", nm_data, NULL, NULL, NULL, "NM2")
  par <- c(0.1, 0.3, 0.2, -0.5)
  expect_equal(sum(ns$.nm_loglik(c(par, 2), prep_p, "NMP")),
               sum(ns$.nm_loglik(par, prep_2, "NM2")))
  expect_equal(sum(ns$.nm_loglik(c(par, 1), prep_p, "NMP")),
               sum(ns$.nm_loglik(par, prep_2, "NM1")))
})

test_that("NM2 fit returns a well-formed flexCountReg object", {
  mod <- negmultinom(y ~ x1 + z, group_var = "id", data = nm_data,
                     offset = "lnL")
  expect_s3_class(mod, "flexCountReg")
  expect_equal(mod$model$modelType, "NegMultinom")
  expect_equal(mod$model$family, "NM2")
  expect_equal(names(mod$model$estimate),
               c("(Intercept)", "x1", "z", "ln(alpha)"))
  expect_true(all(is.finite(mod$model$se)))
  expect_equal(mod$model$n_groups, 300)
  expect_equal(length(mod$model$eb.fitted), nrow(nm_data))
  out <- capture.output(s <- summary(mod))
  expect_equal(nrow(s), 5) # 4 parameters + offset row
})

test_that("family names are cleaned and invalid names rejected", {
  mod <- negmultinom(y ~ x1, group_var = "id", data = nm_data,
                     family = "nm-p", max.iters = 200)
  expect_equal(mod$model$family, "NMP")
  expect_true("p" %in% names(mod$model$estimate))
  expect_error(negmultinom(y ~ x1, "id", nm_data, family = "NB2"),
               "must be one of")
})

test_that("ln(alpha) and p can be functions of covariates", {
  # NM-P needs a reasonably large panel to separate alpha from p
  big <- sim_nm(G = 2000, Tn = 5, p = 1.5, seed = 7)
  mod <- negmultinom(y ~ x1 + z, group_var = "id", data = big,
                     family = "NMP", dis_param_formula_1 = ~ z,
                     dis_param_formula_2 = ~ z)
  expect_equal(tail(names(mod$model$estimate), 4),
               c("ln(alpha):(Intercept)", "ln(alpha):z",
                 "p:(Intercept)", "p:z"))
  expect_equal(length(mod$model$alpha), 2000)
  expect_true(all(is.finite(mod$model$se)))
  expect_equal(length(unique(round(mod$model$p, 8))), 2)
})

test_that("within-cluster variation in dispersion covariates warns", {
  expect_warning(
    negmultinom(y ~ x1, group_var = "id", data = nm_data,
                dis_param_formula_1 = ~ x1),
    "vary within clusters")
})

test_that("dis_param_formula_2 is ignored (with warning) for NM1/NM2", {
  expect_warning(
    mod <- negmultinom(y ~ x1, group_var = "id", data = nm_data,
                       family = "NM1", dis_param_formula_2 = ~ z),
    "only used when")
  expect_false(any(grepl("^p", names(mod$model$estimate))))
})

test_that("EB predictions match the closed-form NM2 expression", {
  mod <- negmultinom(y ~ x1 + z, group_var = "id", data = nm_data,
                     offset = "lnL")
  mu <- predict(mod)
  expect_equal(mu, as.numeric(mod$model$fitted.values))
  alpha <- exp(mod$model$estimate[["ln(alpha)"]])
  yR <- ave(nm_data$y, nm_data$id, FUN = sum)
  muR <- ave(mu, nm_data$id, FUN = sum)
  eb_manual <- unname(mu * (1 + alpha * yR) / (1 + alpha * muR))
  expect_equal(predict(mod, method = "EB"), eb_manual)
  expect_equal(as.numeric(mod$model$eb.fitted), eb_manual)

  # EB on newdata uses the clusters present in newdata
  nd <- nm_data[nm_data$id %in% 1:10, ]
  expect_equal(length(predict(mod, newdata = nd, method = "EB")), nrow(nd))
  nd$y <- NULL
  expect_error(predict(mod, newdata = nd, method = "EB"), "require")
  expect_equal(length(predict(mod, newdata = nd)), nrow(nd))
})

test_that("offset() in the formula is equivalent to the offset argument", {
  m1 <- negmultinom(y ~ x1 + z, "id", nm_data, offset = "lnL")
  m2 <- negmultinom(y ~ x1 + z + offset(lnL), "id", nm_data)
  expect_equal(m1$model$maximum, m2$model$maximum, tolerance = 1e-6)
  expect_equal(predict(m1), predict(m2), tolerance = 1e-5)
})

test_that("robust and bootstrap standard errors are produced", {
  mr <- negmultinom(y ~ x1, "id", nm_data, stderr = "robust")
  expect_true(all(is.finite(mr$model$se)))
  set.seed(1)
  mb <- negmultinom(y ~ x1, "id", nm_data, stderr = "boot", bootstraps = 5)
  expect_true(all(is.finite(mb$model$se)))
  expect_error(negmultinom(y ~ x1, "id", nm_data, stderr = "boot"),
               "bootstraps")
})

test_that("missing values are dropped consistently", {
  d <- nm_data
  d$x1[c(1, 50)] <- NA
  d$id[7] <- NA
  mod <- negmultinom(y ~ x1, "id", d)
  expect_equal(mod$model$nobs, nrow(d) - 3)
  eb <- predict(mod, method = "EB")
  expect_true(all(is.na(eb[c(1, 7, 50)])))
})

test_that("multiple grouping variables define clusters by interaction", {
  d <- nm_data
  d$id_a <- (d$id - 1) %/% 10
  d$id_b <- (d$id - 1) %% 10
  m1 <- negmultinom(y ~ x1, "id", d)
  m2 <- negmultinom(y ~ x1, c("id_a", "id_b"), d)
  expect_equal(m1$model$maximum, m2$model$maximum, tolerance = 1e-6)
})

test_that("margEffTable works with negmultinom models", {
  mod <- negmultinom(y ~ x1 + z, "id", nm_data, offset = "lnL")
  me <- margEffTable(mod, vars = c("x1", "z"))
  expect_true(all(is.finite(me$std_error)))
})
