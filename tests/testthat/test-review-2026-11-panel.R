# Regression tests for the November 2026 review of the panel estimators:
# analytic scores and accurate Hessians for repois() and poislind_re(),
# robust covariances of countreg_fe() fits, mean_fe predictions, and the
# matching of panel scores to rows.

data("washington_roads", package = "flexCountReg")

re_f <- Total_crashes ~ lnaadt + speed50 + offset(lnlength)

# Weights constant within each segment, and the same data with each segment
# repeated w times as separate panels
re_weighted_data <- function() {
  d <- washington_roads
  set.seed(3)
  wid <- sample(1:3, nlevels(d$ID), TRUE)
  d$w <- wid[as.integer(d$ID)]
  dd <- d[rep(seq_len(nrow(d)), d$w), ]
  dd$rep <- stats::ave(seq_len(nrow(dd)), dd$ID, dd$Year, FUN = seq_along)
  dd$IDr <- paste(dd$ID, dd$rep)
  list(d = d, dd = dd)
}

test_that("repois() Hessian standard errors use the analytic scores", {
  skip_on_cran()
  d <- washington_roads
  rp <- repois(re_f, group_var = "ID", data = d, method = "BFGS")
  nm <- negmultinom(re_f, group_var = "ID", data = d)
  # Same model: same estimates and standard errors (previously ln(alpha)
  # had 0.2125 instead of 0.2003)
  expect_equal(unname(coef(rp)), unname(coef(nm)), tolerance = 1e-4)
  expect_equal(unname(rp$model$se), unname(nm$model$se), tolerance = 1e-3)
  expect_equal(unname(sqrt(diag(vcov(rp, type = "robust")))),
               unname(sqrt(diag(vcov(nm, type = "robust")))),
               tolerance = 1e-3)
  # The stored per-panel scores are the analytic ones
  prep <- prepare_model_data(re_f, d, group_var = "ID")
  cl <- split(seq_along(prep$y), prep$group, drop = TRUE)
  G <- .repois_grad(coef(rp), as.numeric(prep$y), prep$X, prep$offset, cl)
  expect_equal(unname(as.matrix(estfun(rp))), G, tolerance = 1e-8)
  # The default Nelder-Mead fit gets the same standard errors
  rpn <- repois(re_f, group_var = "ID", data = d)
  expect_equal(unname(rpn$model$se), unname(nm$model$se), tolerance = 5e-3)
})

test_that("repois() with panel weights equals the expanded data", {
  skip_on_cran()
  wd <- re_weighted_data()
  rpw <- repois(re_f, group_var = "ID", data = wd$d, method = "BFGS",
                weights = "w")
  rpe <- repois(re_f, group_var = "IDr", data = wd$dd, method = "BFGS")
  expect_equal(as.numeric(logLik(rpw)), as.numeric(logLik(rpe)),
               tolerance = 1e-6)
  expect_equal(unname(coef(rpw)), unname(coef(rpe)), tolerance = 1e-4)
  # Previously 0.291 (weighted) and 0.978 (expanded) for the intercept
  expect_equal(unname(rpw$model$se), unname(rpe$model$se), tolerance = 1e-4)
  expect_equal(unname(rpw$model$se[1]), 0.3477, tolerance = 1e-3)
  # Previously 0.283 and 3.35
  expect_equal(unname(sqrt(diag(vcov(rpw, type = "robust")))),
               unname(sqrt(diag(vcov(rpe, type = "robust")))),
               tolerance = 1e-4)
  expect_equal(unname(sqrt(diag(vcov(rpw, type = "robust")))[1]), 0.3995,
               tolerance = 1e-3)
})

test_that("repois() bootstrap refits still work with the analytic score", {
  skip_on_cran()
  d <- washington_roads[washington_roads$ID %in%
                          levels(washington_roads$ID)[1:60], ]
  set.seed(1)
  rb <- suppressWarnings(repois(Total_crashes ~ lnaadt + offset(lnlength),
                                group_var = "ID", data = d,
                                method = "BFGS", stderr = "boot",
                                bootstraps = 5))
  expect_true(all(is.finite(rb$model$se)))
  expect_gt(rb$model$n_boot_ok, 0)
})

test_that("poislind_re() has an analytic score and an accurate Hessian", {
  skip_on_cran()
  d <- washington_roads
  pl <- poislind_re(re_f, group_var = "ID", data = d, method = "BFGS")
  prep <- prepare_model_data(re_f, d, group_var = "ID")
  y <- as.numeric(prep$y)
  ll <- function(b) .fcr_poislind_re_ll(b, y, prep$X, prep$group,
                                        prep$offset)
  # Per-panel scores against numerical derivatives, at and off the estimates
  for (b in list(coef(pl), coef(pl) + c(0.05, -0.01, 0.02, 0.5))) {
    expect_equal(.fcr_poislind_re_grad(b, y, prep$X, prep$group,
                                       prep$offset),
                 unname(numDeriv::jacobian(ll, b)), tolerance = 1e-6)
  }
  # Hessian standard errors equal those of an accurate numerical Hessian
  # (previously the intercept had 0.700 instead of 0.497)
  H <- numDeriv::hessian(function(b) sum(ll(b)), coef(pl))
  expect_equal(unname(pl$model$se), sqrt(diag(solve(-H))), tolerance = 1e-4)
  expect_equal(unname(as.matrix(estfun(pl))),
               .fcr_poislind_re_grad(coef(pl), y, prep$X, prep$group,
                                     prep$offset), tolerance = 1e-8)
})

test_that("poislind_re() with panel weights equals the expanded data", {
  skip_on_cran()
  wd <- re_weighted_data()
  f <- Animal ~ lnaadt + speed50 + offset(lnlength)
  pw <- poislind_re(f, group_var = "ID", data = wd$d, method = "BFGS",
                    weights = "w")
  pe <- poislind_re(f, group_var = "IDr", data = wd$dd, method = "BFGS")
  expect_equal(unname(coef(pw)), unname(coef(pe)), tolerance = 1e-4)
  expect_equal(unname(pw$model$se), unname(pe$model$se), tolerance = 1e-3)
})

fe_data <- function() {
  d <- washington_roads
  set.seed(1)
  d$x1 <- stats::rnorm(nrow(d))
  d$w <- sample(1:3, nrow(d), TRUE)
  d$y <- stats::rpois(nrow(d), exp(0.5 + 0.3 * d$x1 + d$lnlength))
  d
}

test_that("vcov(type = 'robust') works for countreg_fe() fits", {
  skip_on_cran()
  d <- fe_data()
  f1 <- suppressMessages(countreg_fe(y ~ x1 + lnaadt + offset(lnlength),
                                     fixef = ~ ID, data = d))
  V <- vcov(f1, type = "robust")
  expect_true(all(is.finite(V)))
  # Groups are the units: clustered by group without G / (G - 1)
  G <- f1$model$n_clusters[[1]]
  expect_equal(V, vcov(f1, type = "cluster", cluster = ~ ID) * (G - 1) / G,
               ignore_attr = TRUE, tolerance = 1e-10)
  expect_equal(sqrt(diag(V)), sqrt(diag(sandwich::sandwich(f1))),
               ignore_attr = TRUE, tolerance = 1e-10)

  # Weights that vary within groups: the weighted fit equals the expanded
  # data (rows repeated w times within their group)
  d2 <- d[rep(seq_len(nrow(d)), d$w), ]
  for (fam in c("POISSON", "NB2")) {
    fw <- suppressMessages(countreg_fe(y ~ x1 + lnaadt, fixef = ~ ID,
                                       data = d, weights = "w",
                                       family = fam))
    fe <- suppressMessages(countreg_fe(y ~ x1 + lnaadt, fixef = ~ ID,
                                       data = d2, family = fam))
    k <- 1:2
    expect_equal(sqrt(diag(vcov(fw, type = "robust")))[k],
                 sqrt(diag(vcov(fe, type = "robust")))[k],
                 tolerance = if (fam == "POISSON") 1e-6 else 1e-3)
  }

  # Demeaned (several dimensions): rows are the units, weights frequency
  f3 <- suppressMessages(countreg_fe(y ~ x1 + lnlength,
                                     fixef = ~ ID + Year, data = d,
                                     weights = "w"))
  f3e <- suppressMessages(countreg_fe(y ~ x1 + lnlength,
                                      fixef = ~ ID + Year, data = d2))
  expect_equal(sqrt(diag(vcov(f3, type = "robust"))),
               sqrt(diag(vcov(f3e, type = "robust"))), tolerance = 1e-6)
})

test_that("countreg_fe() robust SEs equal fixest's", {
  skip_on_cran()
  skip_if_not_installed("fixest")
  d <- fe_data()
  ss <- fixest::ssc(adj = FALSE, cluster.adj = FALSE)
  f1 <- suppressMessages(countreg_fe(y ~ x1 + lnaadt + offset(lnlength),
                                     fixef = ~ ID, data = d))
  g1 <- fixest::fepois(y ~ x1 + lnaadt + offset(lnlength) | ID, data = d,
                       vcov = ~ ID, ssc = ss, notes = FALSE)
  expect_equal(sqrt(diag(vcov(f1, type = "robust"))), fixest::se(g1),
               ignore_attr = TRUE, tolerance = 1e-6)
  f3 <- suppressMessages(countreg_fe(y ~ x1 + lnlength,
                                     fixef = ~ ID + Year, data = d))
  g3 <- fixest::fepois(y ~ x1 + lnlength | ID + Year, data = d,
                       vcov = "hetero", ssc = ss, notes = FALSE)
  expect_equal(sqrt(diag(vcov(f3, type = "robust"))), fixest::se(g3),
               ignore_attr = TRUE, tolerance = 1e-6)
})

test_that("mean_fe predictions do not depend on the other rows", {
  skip_on_cran()
  d <- fe_data()
  f1 <- suppressMessages(countreg_fe(y ~ x1, fixef = ~ ID, data = d))
  nd <- d[c(2, 3, 5), ]
  nd$ID <- factor(c("2", "3", "zzz"))
  p_all <- predict(f1, newdata = nd, newgroup = "mean_fe")[3]
  p_one <- predict(f1, newdata = nd[3, ], newgroup = "mean_fe")
  nd2 <- nd
  nd2$ID <- factor(c("2", "2", "zzz"))
  expect_equal(p_all, p_one)
  expect_equal(predict(f1, newdata = nd2, newgroup = "mean_fe")[3], p_one)
  # The mean of the effects, weighted by the groups' estimation sizes
  a <- f1$model$fixef$fe1
  n_g <- table(factor(d$ID[f1$model$rows_used], levels = names(a)))
  expect_equal(p_one, exp(sum(coef(f1) * nd$x1[3]) +
                            stats::weighted.mean(a, as.vector(n_g))),
               ignore_attr = TRUE)
  # Fits without the stored sizes use the plain mean
  f0 <- f1
  f0$model$fixef_n <- NULL
  expect_equal(predict(f0, newdata = nd[3, ], newgroup = "mean_fe"),
               exp(sum(coef(f1) * nd$x1[3]) + mean(a)), ignore_attr = TRUE)
})

test_that("panel scores are matched to rows by panel id with one row each", {
  skip_on_cran()
  set.seed(5)
  n <- 200
  dd <- data.frame(x = stats::rnorm(n), id = sample(sprintf("p%03d", 1:n)),
                   cl = rep(1:20, each = 10))
  dd$y <- stats::rnbinom(n, mu = exp(0.5 + 0.3 * dd$x), size = 2)
  m <- repois(y ~ x, group_var = "id", data = dd, method = "BFGS")
  nb <- countreg(y ~ x, data = dd, family = "NB2", method = "BFGS")
  # With one row per panel, the RE Poisson is the NB2 model: the clustered
  # standard errors agree only when the panel scores reach the right rows
  expect_equal(sqrt(diag(vcov(m, type = "cluster", cluster = ~ cl)))[1:2],
               sqrt(diag(vcov(nb, type = "cluster", cluster = ~ cl)))[1:2],
               ignore_attr = TRUE, tolerance = 1e-4)
  expect_equal(.fcr_unit_index(m, n), m$model$panel_index)
})
