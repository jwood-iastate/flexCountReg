# tidymodels integration: the count_reg() model and its flexCountReg engine.

test_that("count_family() is a dials parameter over count distributions", {
  skip_if_not_installed("dials")
  p <- count_family()
  expect_s3_class(p, "qual_param")
  expect_true("NB2" %in% p$values)
})

test_that("parsnip engine fits and predicts with the tidymodels contract", {
  skip_if_not_installed("parsnip")
  spec <- count_reg(family = "NB2") |>
    parsnip::set_engine("flexCountReg", method = "BFGS")
  f <- parsnip::fit(spec, Total_crashes ~ lnaadt + lnlength + speed50,
                    data = washington_roads)
  expect_s3_class(f, "model_fit")
  p <- predict(f, new_data = washington_roads[1:5, ])
  expect_s3_class(p, "tbl_df")
  expect_named(p, ".pred")
  expect_equal(nrow(p), 5)
  m <- MASS::glm.nb(Total_crashes ~ lnaadt + lnlength + speed50,
                    washington_roads)
  expect_equal(p$.pred, unname(fitted(m)[1:5]), tolerance = 1e-3)
  
  nd <- washington_roads[1:5, ]
  nd$lnaadt[2] <- NA
  p2 <- predict(f, new_data = nd)
  expect_equal(nrow(p2), 5)
  expect_true(is.na(p2$.pred[2]))
  expect_false(anyNA(p2$.pred[-2]))
  
  # The family is a main argument
  f2 <- parsnip::fit(count_reg(family = "Poisson") |>
                       parsnip::set_engine("flexCountReg", method = "BFGS"),
                     Total_crashes ~ lnaadt, data = washington_roads)
  expect_equal(f2$fit$model$family, "POISSON")
})

test_that("count_reg_fit passes case weights as frequency weights", {
  d <- washington_roads
  w <- rep(1:3, length.out = nrow(d))
  a <- count_reg_fit(Total_crashes ~ lnaadt, d, family = "NB2", weights = w,
                     method = "BFGS")
  d$w <- w
  b <- countreg(Total_crashes ~ lnaadt, d, weights = "w", method = "BFGS")
  expect_equal(coef(a), coef(b))
})
