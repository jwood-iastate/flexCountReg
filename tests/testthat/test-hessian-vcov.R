test_that("flat likelihood directions get NA variances, not zero", {
  # Hessian of a 3-parameter fit that is flat in the third parameter
  H <- -diag(c(4, 2, 0))
  dimnames(H) <- list(c("a", "b", "c"), c("a", "b", "c"))
  fit <- list(hessian = H, estimate = c(a = 1, b = 2, c = 3))
  expect_warning(V <- flexCountReg:::.fcr_hessian_vcov(fit),
                 "flat in a direction involving c")
  expect_equal(diag(V)[1:2], c(1 / 4, 1 / 2))
  expect_true(all(is.na(V[3, ])) && all(is.na(V[, 3])))
  # A flat direction shared by two parameters (a + b not identified)
  H2 <- -matrix(c(1, 1, 0, 1, 1, 0, 0, 0, 3), 3)
  fit2 <- list(hessian = H2, estimate = c(a = 0, b = 0, c = 0))
  expect_warning(V2 <- flexCountReg:::.fcr_hessian_vcov(fit2), "a, b")
  expect_true(is.na(V2[1, 1]) && is.na(V2[2, 2]))
  expect_equal(V2[3, 3], 1 / 3)
  # Nonsingular Hessians are inverted as before, without a warning
  expect_silent(V3 <- flexCountReg:::.fcr_hessian_vcov(
    list(hessian = -diag(2))))
  expect_equal(V3, diag(2))
})
