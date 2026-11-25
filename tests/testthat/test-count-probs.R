# Count probabilities (chosen counts) and intervals of every estimator.

q_ <- function(x) suppressMessages(suppressWarnings(x))
cp_data <- function() {
  ids <- unique(washington_roads$ID)[1:150]
  d <- washington_roads[washington_roads$ID %in% ids, ]
  d$ID <- as.integer(as.character(d$ID))
  set.seed(4)
  d$c1 <- stats::rbinom(nrow(d), d$Total_crashes, 0.3)
  d$c2 <- d$Total_crashes - d$c1
  d
}
cp_fits <- function(d) {
  f <- Total_crashes ~ lnaadt
  g <- cbind(c1, c2) ~ lnaadt
  list(
    countreg = countreg(f, d, method = "BFGS"),
    zip = countreg(f, d, family = "POISSON", zero_formula = ~ 1,
                   method = "BFGS"),
    hurdle = countreg(f, d, family = "POISSON", zero_formula = ~ 1,
                      zero_type = "hurdle", method = "BFGS"),
    rp = countreg_rp(f, ~ -1 + lnlength, data = d, ndraws = 50),
    rp_uniform = countreg_rp(f, ~ -1 + lnlength, data = d, ndraws = 50,
                             rpardists = c(lnlength = "u")),
    mix = countreg_mix(f, d, n_starts = 2, seed = 1),
    fe = countreg_fe(f, ~ ID, d),
    hier = countreg_hier(Total_crashes ~ lnaadt + (1 | ID), d,
                         family = "POISSON"),
    hier_gamma = countreg_hier(Total_crashes ~ lnaadt + (1 | ID), d,
                               family = "POISSON", re_dist = "gamma",
                               re_draws = 100),
    alloc = countreg_alloc(g, d),
    alloc_shared = countreg_alloc(g, d, family = "POISSON",
                                  shared_formula = ~ 1, shared_nodes = 9),
    mv = countreg_mv(g, d),
    mv_rp = countreg_mv(g, d, family = "POISSON", rpar_formula = ~ 1,
                        correlation = "none", group_var = "ID",
                        rpar_draws = 30),
    renb = renb(f, "ID", d, method = "BFGS"),
    repois = repois(f, "ID", d, method = "BFGS"),
    poislind = poislind_re(f, "ID", d, method = "BFGS"),
    nm1 = negmultinom(f, "ID", d, family = "NM1"),
    nm2 = negmultinom(f, "ID", d))
}
cp_d <- cp_data()
cp_m <- q_(cp_fits(cp_d))

for (nm in names(cp_m)) {
  test_that(paste("count probabilities:", nm), {
    m <- cp_m[[nm]]
    nd <- if (nm == "fe") cp_d[m$model$rows_used[1:5], ] else cp_d[1:5, ]
    ty <- if (grepl("alloc", nm)) "count_prob" else "prob"
    extra <- if (grepl("hier", nm)) list(marginal = TRUE) else list()
    if (grepl("hier", nm)) nd$ID[1] <- max(cp_d$ID) + 1   # a new group
    pr <- function(...) {
      P <- do.call(predict, c(list(m, newdata = nd, type = ty, ...), extra))
      if (is.list(P)) P else list(P)
    }
    k <- 0:300
    P <- pr(counts = k)
    mu <- as.matrix(do.call(predict, c(list(m, newdata = nd), extra)))
    for (j in seq_along(P)) {
      expect_equal(rowSums(P[[j]]), rep(1, 5), tolerance = 1e-6)
      expect_equal(as.vector(P[[j]] %*% k), as.vector(mu[, j]),
                   tolerance = if (nm == "rp") 0.01 else
                     if (nm == "rp_uniform") 0.01 else
                     if (nm == "hier_gamma") 0.01 else 1e-5)
    }
    # A chosen range equals those columns; tails complete the rows
    a <- pr(counts = c(4, 2, 3, 3), overflow = TRUE, underflow = TRUE)[[1]]
    expect_equal(colnames(a), c("<2", "2", "3", "4", ">4"))
    expect_equal(unname(a[, 2:4]), unname(P[[1]][, 3:5]), tolerance = 1e-8)
    expect_equal(rowSums(a), rep(1, 5), tolerance = 1e-8)
    expect_equal(pr(max_count = 3), pr(counts = 0:3))
    expect_equal(dim(pr(counts = 7)[[1]]), c(5L, 1L))
    expect_error(pr(counts = 0:2, max_count = 2), "not both")
    expect_error(pr(counts = -1), "whole numbers")
    expect_error(pr(counts = 1.5), "whole numbers")
    # A row with a missing covariate has no probabilities
    nd$lnaadt[2] <- NA
    z <- pr(counts = 0:2)[[1]]
    expect_true(all(is.na(z[2, ])))
    expect_false(anyNA(z[-2, ]))
  })

  test_that(paste("intervals:", nm), {
    m <- cp_m[[nm]]
    nd <- cp_d[1:4, ]
    if (nm == "fe") {
      expect_error(predict(m, nd, interval = "confidence"), "group effects")
      return(invisible())
    }
    iv <- function(...) {
      v <- q_(predict(m, newdata = nd, ndraws = 25, ...))
      if (is.data.frame(v)) list(v) else v
    }
    mu <- as.matrix(predict(m, newdata = nd))
    ci <- iv(interval = "confidence", level = 0.9)
    pi <- iv(interval = "prediction", level = 0.9)
    expect_length(ci, ncol(mu))
    for (j in seq_along(ci)) {
      expect_equal(ci[[j]]$fit, as.vector(mu[, j]))
      expect_true(all(ci[[j]]$lwr < ci[[j]]$upr))
      expect_true(all(pi[[j]]$lwr == round(pi[[j]]$lwr)))
      expect_true(all(pi[[j]]$upr >= ci[[j]]$fit - 1))
    }
    expect_equal(iv(interval = "confidence", level = 0.9), ci)
    nd$lnaadt[2] <- NA
    expect_true(is.na(iv(interval = "confidence")[[1]]$upr[2]))
    expect_true(is.na(iv(interval = "prediction")[[1]]$upr[2]))
  })
}

test_that("confidence intervals agree with the delta method", {
  m <- cp_m$repois
  nd <- cp_d[1:3, ]
  v <- predict(m, nd, interval = "confidence", ndraws = 4000)
  X <- cbind(1, nd$lnaadt)
  se <- sqrt(rowSums((X %*% vcov(m)[1:2, 1:2]) * X))
  expect_equal(v$lwr, v$fit * exp(-1.96 * se), tolerance = 0.02)
  expect_equal(v$upr, v$fit * exp(1.96 * se), tolerance = 0.02)
})

test_that("known marginal distributions", {
  nd <- cp_d[1:4, ]
  m <- cp_m$repois
  expect_equal(unname(predict(m, nd, type = "prob", counts = 0:3)),
               outer(predict(m, nd), 0:3, function(u, k)
                 dnbinom(k, size = exp(-coef(m)[["lnalpha"]]), mu = u)))
  # Allocation: independent thinning of a Poisson total is Poisson
  a <- q_(countreg_alloc(cbind(c1, c2) ~ lnaadt, cp_d, family = "POISSON"))
  P <- predict(a, nd, type = "count_prob", counts = 0:4)
  mu <- predict(a, nd)
  expect_equal(unname(P$c1), outer(unname(mu[, "c1"]), 0:4, function(u, k)
    dpois(k, u)), tolerance = 1e-7)
  P <- predict(a, nd, type = "total_prob", counts = 0:200, conditional = TRUE)
  expect_equal(unname(P[, 1]), rep(0, 4))
  expect_equal(rowSums(P), rep(1, 4))
  expect_equal(as.vector(P %*% 0:200),
               as.vector(predict(a, nd, type = "total_positive")))
  fr <- q_(countreg_alloc(cbind(c1, c2) ~ lnaadt, cp_d,
                          likelihood = "fractional"))
  expect_error(predict(fr, nd, type = "count_prob"), "fractional")
  expect_equal(rowSums(predict(fr, nd, type = "total_prob", counts = 0:200)),
               rep(1, 4), tolerance = 1e-8)
})

test_that("hierarchical prediction intervals include new-group variation", {
  m <- cp_m$hier
  nd <- cp_d[1:3, ]
  nd$ID <- max(cp_d$ID) + 1:3
  wide <- q_(predict(m, nd, interval = "prediction", ndraws = 40))
  narrow <- q_(predict(m, nd, interval = "prediction", ndraws = 40,
                       marginal = FALSE))
  expect_true(all(wide$upr >= narrow$upr))
})

test_that("parameters at a boundary are held fixed", {
  g <- NULL
  m <- cp_m$repois
  m$model$vcov["lnalpha", "lnalpha"] <- 1e6
  expect_message(g <- .fcr_general_draws(m, 50, NULL), "lnalpha")
  expect_true(all(g$draws[, "lnalpha"] == coef(m)[["lnalpha"]]))
})

test_that("rootogram() works for every estimator with one outcome", {
  for (nm in setdiff(names(cp_m), c("alloc", "alloc_shared", "mv", "mv_rp",
                                    "fe"))) {
    r <- rootogram(cp_m[[nm]], plot = FALSE, overflow = TRUE)
    expect_equal(sum(r$observed), nrow(cp_d))
    expect_equal(sum(r$expected), nrow(cp_d), tolerance = 1e-5)
  }
})

test_that("margEffTable() gives effects on each allocation category", {
  m <- cp_m$alloc
  e <- as.data.frame(margEffTable(m, digits = 8))
  expect_equal(e$outcome, c("c1", "c2", "total"))
  expect_equal(sum(e$estimate[1:2]), e$estimate[3], tolerance = 1e-5)
  expect_true(all(e$std_error > 0))
  # The total's mean is exp(x'b): AME = b * mean(mu)
  expect_equal(e$estimate[3],
               coef(m)[[2]] * mean(predict(m, type = "total")),
               tolerance = 1e-4)
  expect_error(margEffTable(m, se = "bootstrap"), "several outcomes")
})
