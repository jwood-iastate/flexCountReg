# Cluster-robust covariances, post-fit vcov() types, and the sandwich
# package's estfun() / bread() methods.

cl_data <- function() {
  data("washington_roads", package = "flexCountReg", envir = environment())
  w <- washington_roads
  set.seed(1)
  ids <- unique(w$ID)
  # Corridors contain whole sites; years vary within sites
  w$corr <- sample(1:40, length(ids), TRUE)[match(w$ID, ids)]
  w$yr <- sample(1:6, nrow(w), TRUE)
  w
}

test_that("cluster covariances match sandwich::vcovCL()", {
  w <- cl_data()
  m <- countreg(Total_crashes ~ lnaadt + lnlength, data = w, family = "NB2",
                method = "BFGS")
  v1 <- vcov(m, type = "cluster", cluster = ~ ID)
  expect_equal(v1, sandwich::vcovCL(m, cluster = w$ID, type = "HC0",
                                    cadjust = TRUE),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_identical(attr(v1, "n_clusters"), c(ID = length(unique(w$ID))))
  # Two-way (Cameron, Gelbach and Miller 2011)
  expect_equal(vcov(m, type = "cluster", cluster = ~ ID + yr),
               sandwich::vcovCL(m, cluster = w[, c("ID", "yr")],
                                type = "HC0", cadjust = TRUE),
               tolerance = 1e-10, ignore_attr = TRUE)
  # Other ways to give the clusters
  expect_equal(vcov(m, type = "cluster", cluster = "ID"), v1)
  expect_equal(vcov(m, type = "cluster", cluster = w$ID), v1,
               ignore_attr = TRUE)
  # Clusters of single rows: the robust sandwich times n / (n - 1)
  n <- nrow(w)
  expect_equal(vcov(m, type = "cluster", cluster = seq_len(n)),
               vcov(m, type = "robust") * n / (n - 1),
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(sandwich::sandwich(m), vcov(m, type = "robust"),
               tolerance = 1e-12)
  expect_equal(vcov(m, type = "hessian"), vcov(m))
  expect_identical(vcov(m, type = "fit"), vcov(m))
})

test_that("stderr = \"cluster\" at fit time", {
  w <- cl_data()
  m <- countreg(Total_crashes ~ lnaadt + lnlength, data = w, family = "NB2",
                method = "BFGS", cluster = ~ ID)
  expect_identical(m$model$vcov_type, "cluster")
  expect_identical(m$model$stderr, "cluster")
  expect_identical(m$model$cluster, "ID")
  expect_equal(vcov(m), vcov(m, type = "cluster", cluster = ~ ID),
               ignore_attr = TRUE)
  expect_null(attr(vcov(m), "n_clusters"))
  expect_output(print(summary(m)), "cluster-robust \\(507 clusters of ID\\)")
  expect_error(countreg(Total_crashes ~ lnaadt, data = w, stderr = "cluster",
                        method = "BFGS"), "needs `cluster`")
  expect_error(countreg(Total_crashes ~ lnaadt, data = w, stderr = "robust",
                        cluster = ~ ID, method = "BFGS"), "is used with")
  expect_error(countreg(Total_crashes ~ lnaadt, data = w, cluster = ~ nope,
                        method = "BFGS"), "not found")
  w$bad <- w$ID
  w$bad[3] <- NA
  expect_error(countreg(Total_crashes ~ lnaadt, data = w, cluster = ~ bad,
                        method = "BFGS"), "missing values")
  # Estimators without fit-time clustering reject stderr = "cluster"
  expect_error(renb(Total_crashes ~ lnaadt, data = w, group_var = "ID",
                    stderr = "cluster"), "should be one of")
})

test_that("panel models cluster their panel scores", {
  skip_on_cran()
  w <- cl_data()
  rp <- suppressWarnings(countreg_rp(
    Total_crashes ~ lnaadt, rpar_formula = ~ -1 + speed50, data = w,
    group_var = "ID", ndraws = 30, cluster = ~ corr))
  expect_identical(rp$model$vcov_type, "cluster")
  expect_equal(vcov(rp), vcov(rp, type = "cluster", cluster = ~ corr),
               ignore_attr = TRUE)
  # One cluster id per panel, in the order of the panel ids
  cu <- w$corr[!duplicated(w$ID)][order(unique(w$ID))]
  expect_equal(vcov(rp), sandwich::vcovCL(rp, cluster = cu, type = "HC0",
                                          cadjust = TRUE),
               tolerance = 1e-10, ignore_attr = TRUE)
  G <- length(unique(w$ID))
  expect_equal(vcov(rp, type = "cluster", cluster = ~ ID),
               vcov(rp, type = "robust") * G / (G - 1),
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_error(vcov(rp, type = "cluster", cluster = ~ yr),
               "varies within a panel")
  re <- suppressWarnings(renb(Total_crashes ~ lnaadt, data = w,
                              group_var = "ID"))
  v <- vcov(re, type = "cluster", cluster = ~ corr)
  expect_true(all(is.finite(sqrt(diag(v)))))
  expect_identical(attr(v, "n_clusters"), c(corr = 40L))
  rpo <- repois(Total_crashes ~ lnaadt, data = w, group_var = "ID")
  expect_equal(nrow(sandwich::estfun(rpo)), G)
  expect_true(all(is.finite(diag(vcov(rpo, type = "cluster",
                                      cluster = ~ corr)))))
})

test_that("weights: frequency copies share their row's cluster", {
  w <- cl_data()[1:400, ]
  set.seed(3)
  w$f <- sample(1:3, nrow(w), TRUE)
  m <- countreg(Total_crashes ~ lnaadt, data = w, family = "NB2",
                method = "BFGS", weights = "f")
  # The expanded data, clustered by site, gives the same covariance
  e <- w[rep(seq_len(nrow(w)), w$f), ]
  me <- countreg(Total_crashes ~ lnaadt, data = e, family = "NB2",
                 method = "BFGS", start.vals = coef(m))
  expect_equal(vcov(m, type = "cluster", cluster = ~ ID),
               vcov(me, type = "cluster", cluster = ~ ID),
               tolerance = 1e-4, ignore_attr = TRUE)
  # Sampling weights: invariant to rescaling the weights
  s1 <- countreg(Total_crashes ~ lnaadt, data = w, family = "NB2",
                 method = "BFGS", weights = "f", weight_type = "sampling",
                 cluster = ~ ID)
  w$f100 <- 100 * w$f
  s2 <- countreg(Total_crashes ~ lnaadt, data = w, family = "NB2",
                 method = "BFGS", weights = "f100", weight_type = "sampling",
                 cluster = ~ ID, start.vals = coef(s1))
  expect_equal(vcov(s1), vcov(s2), tolerance = 1e-4)
})

test_that("estfun() and bread() follow the sandwich package's conventions", {
  w <- cl_data()
  m <- countreg(Total_crashes ~ lnaadt, data = w, family = "POISSON",
                method = "BFGS")
  U <- sandwich::estfun(m)
  expect_identical(dim(U), c(nrow(w), 2L))
  expect_identical(colnames(U), names(coef(m)))
  expect_equal(colSums(U), c(0, 0), tolerance = 1e-3, ignore_attr = TRUE)
  expect_equal(sandwich::bread(m), nrow(w) * vcov(m, type = "hessian"))
  # Poisson scores are (y - mu) x, as for glm()
  g <- stats::glm(Total_crashes ~ lnaadt, data = w, family = poisson())
  expect_equal(sandwich::vcovCL(m, cluster = w$ID, type = "HC0"),
               sandwich::vcovCL(g, cluster = w$ID, type = "HC0"),
               tolerance = 1e-4, ignore_attr = TRUE)
  skip_if_not_installed("lmtest")
  ct <- lmtest::coeftest(m, vcov. = vcov(m, type = "cluster",
                                         cluster = ~ ID))
  expect_equal(unname(ct[, 2]),
               unname(sqrt(diag(vcov(m, type = "cluster", cluster = ~ ID)))))
})

test_that("a model without stored data needs cluster ids as a vector", {
  w <- cl_data()
  m <- countreg(Total_crashes ~ lnaadt, data = w, family = "NB2",
                method = "BFGS", keep_data = FALSE)
  expect_error(vcov(m, type = "cluster", cluster = ~ ID), "keep_data")
  expect_true(all(is.finite(diag(vcov(m, type = "cluster",
                                      cluster = w$ID)))))
  expect_error(vcov(m, type = "cluster"), "needs `cluster`")
  expect_error(vcov(m, type = "robust", cluster = ~ ID), "is used with")
})
