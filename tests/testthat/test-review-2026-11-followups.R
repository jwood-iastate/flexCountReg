# Follow-ups of the November 2026 review

test_that("countreg_fe(): frequency weights count as observations of a group", {
  set.seed(4)
  G <- 60
  d <- data.frame(id = rep(seq_len(G), each = 3), x = stats::rnorm(3 * G))
  d$y <- stats::rnbinom(nrow(d), size = 2,
                        mu = exp(0.3 * d$x + stats::rnorm(G)[d$id]))
  # Some groups are a single row with weight 3 (three observations)
  d <- d[!(d$id %in% 1:10 & duplicated(d$id)), ]
  d$w <- ifelse(d$id %in% 1:10, 3, 1)
  ex <- d[rep(seq_len(nrow(d)), d$w), ]
  fw <- suppressMessages(countreg_fe(y ~ x, fixef = ~ id, data = d,
                                     family = "NB2", method = "concentrated",
                                     weights = "w"))
  fe <- suppressMessages(countreg_fe(y ~ x, fixef = ~ id, data = ex,
                                     family = "NB2",
                                     method = "concentrated"))
  expect_equal(as.numeric(logLik(fw)), as.numeric(logLik(fe)),
               tolerance = 1e-6)
  expect_equal(coef(fw), coef(fe), tolerance = 1e-4)
})

test_that("countreg_hier(): random-slope factors with unused levels", {
  set.seed(5)
  d <- data.frame(g = rep(1:30, each = 5), x = stats::rnorm(150))
  d$f <- factor(sample(c("a", "b"), 150, replace = TRUE),
                levels = c("a", "b", "c"))
  d$y <- stats::rpois(150, exp(0.2 + 0.3 * d$x + stats::rnorm(30, 0, 0.4)[d$g]))
  m <- countreg_hier(y ~ x + (1 + f | g), data = d, family = "POISSON")
  expect_true(all(is.finite(coef(m))))
})

test_that("countreg_alloc(): marginal predictions use enough draws", {
  set.seed(2)
  n <- 400
  d <- data.frame(x = stats::rnorm(n), site = rep(seq_len(n / 4), each = 4))
  tot <- stats::rpois(n, exp(-0.2 + 0.3 * d$x + 0.6 * stats::rnorm(n / 4)[d$site]))
  d$a <- stats::rbinom(n, tot, 0.6)
  d$b <- tot - d$a
  m <- countreg_alloc(cbind(a, b) ~ x, data = d, family = "POISSON",
                      alloc_formula = ~ x, group_var = "site",
                      rpar_total = ~ 1, rpar_draws = 50, seed = 1)
  cf <- coef(m)
  s <- abs(cf["sd:Total:(Intercept)"])
  exact <- exp(cf["(Intercept)"] + cf["x"] * d$x[1:3] + s^2 / 2)
  # Gauss-Hermite integration: exact up to rounding
  expect_equal(unname(predict(m, newdata = d[1:3, ], type = "total")),
               unname(exact), tolerance = 1e-8)
})

test_that("margEffTable() stops clearly when the expected count is infinite", {
  set.seed(6)
  d <- data.frame(x = stats::runif(300, 0.2, 1.5))
  d$y <- stats::rpois(300, exp(0.2 + 0.4 * d$x))
  m <- suppressWarnings(countreg_mix(
    y ~ x, data = d, components = 2, rpar_formula = ~ x, rpardists = "ln",
    rpar_draws = 20, n_starts = 1, seed = 1))
  expect_error(suppressWarnings(margEffTable(m)), "infinite")
})
