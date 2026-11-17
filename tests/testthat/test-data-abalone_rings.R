test_that("abalone_rings has the documented structure", {
  e <- new.env()
  data("abalone_rings", package = "flexCountReg", envir = e)
  d <- e$abalone_rings
  expect_s3_class(d, "data.frame")
  expect_identical(dim(d), c(4177L, 9L))
  expect_identical(names(d), c("sex", "length", "diameter", "height",
                               "whole_weight", "shucked_weight",
                               "viscera_weight", "shell_weight", "rings"))
  expect_identical(levels(d$sex), c("F", "I", "M"))
  expect_identical(as.vector(table(d$sex)), c(1307L, 1342L, 1528L))
  expect_type(d$rings, "integer")
  expect_identical(range(d$rings), c(1L, 29L))
  expect_false(anyNA(d))
  expect_true(all(vapply(d[2:8], is.double, logical(1))))
  expect_equal(sum(d$height == 0), 2)
  expect_equal(mean(d$rings), 9.933684, tolerance = 1e-6)
})

test_that("a Poisson fit to abalone_rings matches glm()", {
  e <- new.env()
  data("abalone_rings", package = "flexCountReg", envir = e)
  d <- e$abalone_rings
  f <- rings ~ sex + length + diameter + shell_weight
  fit <- countreg(f, data = d, family = "POISSON", method = "BFGS")
  ref <- stats::glm(f, data = d, family = stats::poisson())
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)),
               tolerance = 1e-6)
  expect_equal(as.numeric(predict(fit, type = "response")),
               as.numeric(fitted(ref)), tolerance = 1e-5)
})
