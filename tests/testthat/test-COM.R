# PDF ----

test_that("COM PDF", {
  
  val <- dcom(3, lambda=2, nu=3)
  
  expect_equal(val, 0.01046772, tolerance = 0.000001)
})

test_that("COM PDF using mu", {
  
  val <- dcom(3, mu=1.15, nu=3)
  
  expect_equal(val, 0.03134405, tolerance = 0.000001)
})

test_that("COM PDF using vectors", {
  
  val <- dcom(c(0,1,2,3,4,5), mu=c(.1,.2,.3,.4,.5,.6), nu=c(0.1,0.3, 0.2, 0.1, 0.7, 1))
  expected <- c(0.9085289259, 0.1475683696, 0.0403047769, 0.0157318300, 0.0028058045, 0.0003556299)
  expect_equal(val, expected, tolerance = 0.000001)
  
  
  val <- dcom(c(0,1,2,3,4,5), mu=0.5, nu=c(0.1,0.3, 0.2, 0.1, 0.7, 1))
  expected <- c(0.6587951083, 0.2486166660, 0.0770023700, 0.0238262615, 0.0028058045, 0.0001579507)
  expect_equal(val, expected, tolerance = 0.000001)
  
  val <- dcom(c(0,1,2,3,4,5), mu=c(.1,.2,.3,.4,.5,.6), nu=1.1)
  expected <- c(0.9045414667, 0.1656661440, 0.0322230252, 0.0064112827, 0.0012901245, 0.0002613611)
  expect_equal(val, expected, tolerance = 0.000001)
  
})


test_that("COM PDF using vectors (lambda)", {
  
  val <- dcom(c(0,1,2,3,4,5), lambda=c(.1,.2,.3,.4,.5,.6), nu=c(0.1,0.3, 0.2, 0.1, 0.7, 1))
  expected <- c(0.9007014296, 0.1615683242, 0.0558568176, 0.0328045144, 0.0039801598, 0.0003556299)
  expect_equal(val, expected, tolerance = 0.000001)
  
  
  val <- dcom(c(0,1,2,3,4,5), lambda=0.5, nu=c(0.1,0.3, 0.2, 0.1, 0.7, 1))
  expected <- c(00.5218466915, 0.2761512603, 0.1172372363, 0.0545302919, 0.0039801598, 0.0001579507)
  expect_equal(val, expected, tolerance = 0.000001)
  
  val <- dcom(c(0,1,2,3,4,5), lambda=c(.1,.2,.3,.4,.5,.6), nu=1.1)
  expected <- c(0.9051349952, 0.1639578793, 0.0311933627, 0.0060071028, 0.0011583352, 0.0002227429)
  expect_equal(val, expected, tolerance = 0.000001)
  
})

test_that("COM PDF (error 1)", {
  
  ## Negative outcome
  expect_error(dcom(-3, lambda=2, nu=3))
})

test_that("COM PDF (error 2)", {
  
  ## Negative lambda
  expect_error(dcom(3, lambda=-2, nu=3))
})

test_that("COM PDF (error 3)", {
  
  ## Negative nu
  expect_error(dcom(3, lambda=2, nu=-3))
})


test_that("COM PDF using vectors - errors", {
  
  # different number of mu
  expect_error(dcom(c(0,1,2,3,4,5), mu=c(.1,.2,.3,.4,.5), nu=c(0.1,0.3, 0.2, 0.1, 0.7, 1)))
  
  # different number of nu
  expect_error(dcom(c(0,1,2,3,4,5), mu=c(.1,.2,.3,.4,.5, 0.2), nu=c(0.1,0.3, 0.2, 0.1, 0.7, 1,0.5)))
  
  # different number of mu and nu
  expect_error(dcom(c(0,1,2,3,4,5), mu=c(.1,.2,.3,.4,.5), nu=c(0.1,0.3, 1)))
  
  # negative value in x
  expect_error(dcom(c(0,1,2,-3,4,5), mu=c(.1,.2,.3,.4,.5, 0.6), nu=c(0.1,0.3, 0.2, 0.1, 0.7, 1)))
  
  # negative value in mu
  expect_error(dcom(c(0,1,2,3,4,5), mu=c(.1,.2,-.3,.4,.5, 0.6), nu=c(0.1,0.3, 0.2, 0.1, 0.7, 1)))
  
  # negative value in nu
  expect_error(dcom(c(0,1,2,3,4,5), mu=c(.1,.2,.3,.4,.5, 0.6), nu=c(0.1,0.3, -0.2, 0.1, 0.7, 1)))
  
  # negative value for mu
  expect_error(dcom(c(0,1,2,3,4,5), mu=-2, nu=c(0.1,0.3, 0.2, 0.1, 0.7, 1)))
  
  # negative value for nu
  expect_error(dcom(c(0,1,2,3,4,5), mu=c(.1,.2,.3,.4,.5, 0.6), nu=-2))
  
})

test_that("COM CDF", {
  
  pdf_vals <- pcom(2, mu=0.9, nu=0.85)
  pdf_expected <- 0.9306465
  
  expect_true(isTRUE(all.equal(pdf_vals, pdf_expected, tolerance = 0.000001)))
  
  pdf_vals <- pcom(c(0, 1, 3, 9, 10), mu=0.9, nu=0.85)
  pdf_expected <- c(0.4205842, 0.7697908, 0.9831408, 0.9999998, 1.0000000)
  
  expect_true(isTRUE(all.equal(pdf_vals, pdf_expected, tolerance = 0.000001)))
  
  pdf_vals <- pcom(c(0, 1, 3, 9, 10), mu=c(0.9,0.3,1,2,0.19), nu=0.85)
  pdf_expected <- c(0.4205842, 0.9606628, 0.9766657, 0.9998778, 1.0000000)
  
  expect_true(isTRUE(all.equal(pdf_vals, pdf_expected, tolerance = 0.000001)))
  
  pdf_vals <- pcom(c(0, 1, 3, 9, 10), mu=c(0.9,0.3,1,2,0.19), nu=c(0.9,0.3,1,2,0.19))
  pdf_expected <- c(0.4158138, 0.9515524, 0.9810118, 0.9999999, 1.0000000)
  
  expect_true(isTRUE(all.equal(pdf_vals, pdf_expected, tolerance = 0.000001)))
  
  pdf_vals <- pcom(c(0, 1, 3, 9, 10), mu=0.95, nu=c(0.9,0.3,1,2,0.19))
  pdf_expected <- c(0.3964257, 0.7524078, 0.9839256, 1.0000000, 0.9999473)
  
  expect_true(isTRUE(all.equal(pdf_vals, pdf_expected, tolerance = 0.000001)))
  
  
  
  pdf_vals <- pcom(2, lambda=0.9, nu=0.85)
  pdf_expected <- 0.9159908
  
  expect_true(isTRUE(all.equal(pdf_vals, pdf_expected, tolerance = 0.000001)))
  
  pdf_vals <- pcom(c(0, 1, 3, 9, 10), lambda=0.9, nu=0.85)
  pdf_expected <- c(0.3898869, 0.7407852, 0.9779688, 0.9999996, 1.0000000)
  
  expect_true(isTRUE(all.equal(pdf_vals, pdf_expected, tolerance = 0.000001)))
  
  pdf_vals <- pcom(c(0, 1, 3, 9, 10), lambda=c(0.9,0.3,1,2,0.19), nu=0.85)
  pdf_expected <- c(0.3898869, 0.9584149, 0.9689737, 0.9995927, 1.0000000)
  
  expect_true(isTRUE(all.equal(pdf_vals, pdf_expected, tolerance = 0.000001)))
  
  pdf_vals <- pcom(c(0, 1, 3, 9, 10), lambda=c(0.9,0.3,1,2,0.19), nu=c(0.9,0.3,1,2,0.19))
  pdf_expected <- c(0.3957888, 0.9334406, 0.9810118, 1.0000000, 1.0000000)
  
  expect_true(isTRUE(all.equal(pdf_vals, pdf_expected, tolerance = 0.000001)))
  
  pdf_vals <- pcom(c(0, 1, 3, 9, 10), lambda=0.95, nu=c(0.9,0.3,1,2,0.19))
  pdf_expected <- c(0.3753800, 0.5011549, 0.9839256, 1.0000000, 0.9895617)
  
  expect_true(isTRUE(all.equal(pdf_vals, pdf_expected, tolerance = 0.000001)))
})



test_that("COM CDF errors", {
  
  # Negative values
  expect_error(pcom(-2, mu=0.9, nu=0.85))
  expect_error(pcom(2, mu=-0.9, nu=0.85))
  expect_error(pcom(2, mu=0.9, nu=-0.85))
  expect_error(pcom(-2, lambda=0.9, nu=0.85))
  expect_error(pcom(2, lambda=-0.9, nu=0.85))
  expect_error(pcom(2, lambda=0.9, nu=-0.85))
  
  # Negative values with vectors
  expect_error(pcom(c(0, 1, 3, -9, 10), mu=0.9, nu=0.85))
  expect_error(pcom(c(0, 1, 3, 9, 10), mu=c(-0.9,0.3,1,2,0.19), nu=0.85))
  expect_error(pcom(c(0, 1, 3, 9, 10), mu=c(0.9,0.3,1,2,0.19), nu=c(0.9,-0.3,1,2,0.19)))
  expect_error(pcom(c(0, 1, -3, 9, 10), mu=c(-0.9,0.3,1,2,0.19), nu=c(0.9,0.3,1,2,-0.19)))
  expect_error(pcom(c(0, 1, 3, 9, 10), mu=0.9, nu=-0.85))
  expect_error(pcom(c(0, 1, 3, 9, 10), mu=-0.9, nu=0.85))
  expect_error(pcom(c(0, 1, 3, -9, 10), lambda=0.9, nu=0.85))
  expect_error(pcom(c(0, 1, 3, 9, 10), lambda=c(-0.9,0.3,1,2,0.19), nu=0.85))
  expect_error(pcom(c(0, 1, 3, 9, 10), lambda=c(0.9,0.3,1,2,0.19), nu=c(0.9,-0.3,1,2,0.19)))
  expect_error(pcom(c(0, 1, -3, 9, 10), lambda=c(-0.9,0.3,1,2,0.19), nu=c(0.9,0.3,1,2,-0.19)))
  expect_error(pcom(c(0, 1, 3, 9, 10), lambda=0.9, nu=-0.85))
  expect_error(pcom(c(0, 1, 3, 9, 10), lambda=-0.9, nu=0.85))
  
  # Incorrect number of values in vectors
  expect_error(pcom(c(0, 1, 3, 9, 10), mu=c(0.9,1), nu=0.85))
  expect_error(pcom(c(0, 1, 3, 9, 10), mu=0.9, nu=c(0.85,1,1,3)))
  expect_error(pcom(c(0, 1, 3, 9, 10), lambda=c(0.9,1), nu=0.85))
  expect_error(pcom(c(0, 1, 3, 9, 10), lambda=0.9, nu=c(0.85,1,1,3)))
})



test_that("COM CDF (errors)", {
  
  ## negative outcome
  expect_error(pcom(c(0, -1, 3, 9, 10), mu=0.9, nu=0.85))
  
  ## negative mean
  expect_error(ppcom(c(0, 1, 3, 9, 10), mu=-0.9, nu=0.85))
  
  ## negative lambda
  expect_error(pcom(c(0, 1, 3, 9, 10), lambda=-0.9, nu=0.85))
  
  ## negative mu
  expect_error(pcom(c(0, 1, 3, 9, 10), mu=0.9, nu=-0.85))
})

# Quantiles ----

test_that("COM Quantiles", {
  
  quant <- qcom(seq(0.1,0.9,0.1), mu=0.5, nu=1.5)
  quant_expected <- c(0, 0, 0, 0, 0, 1, 1, 1, 1)
  out <- all(quant_expected == quant)
  
  expect_true(out)
})


test_that("COM Quantiles with lambda", {
  
  quant <- qcom(seq(0.1,0.9,0.1), lambda=1.5, nu=1.5)
  quant_expected <- c(0, 0, 1, 1, 1, 1, 2, 2, 2)
  out <- all(quant_expected == quant)
  
  expect_true(out)
})

test_that("COM Quantiles with lambda", {
  
  quant <- qcom(seq(0.1,0.9,0.1), mu=0.5, nu=1.5)
  quant_expected <- c(0, 0, 0, 0, 0, 1, 1, 1, 1)
  out <- all(quant_expected == quant)
  
  expect_true(out)
})


test_that("COM Quantiles (errors)", {
  
  ## QUantile larger than 1
  expect_error(qcom(seq(0.1,1.1,0.1), lambda=1.5, nu=1.5))
  
  ## lambda is negative
  expect_error(qcom(seq(0.1,0.9,0.1), lambda=-1.5, nu=1.5))
  
  ## mu is negative
  expect_error(qcom(seq(0.1,0.9,0.1), mu=-1.5, nu=1.5))
  
  ## nu is negative
  expect_error(qcom(seq(0.1,0.9,0.1), mu=1.5, nu=-1.5))
})

# Random samples ----

test_that("COM Samples", {
  
  set.seed(666)
  comsamples <- rcom(1e6, mu = 0.72, nu = 0.7)
  sample_mean <- mean(comsamples)
  diff <- abs(sample_mean - 0.72) 
  
  testthat::expect_true(diff < 0.04)
})


test_that("COM Samples (errors)", {
  
  ## too many n values
  expect_error(rcom(c(1,5), mu = 0.72, nu = 0.7))
  
  ## negative mu
  expect_error(rcom(5, mu = -3, nu = 0.7))
  
  ## negative lambda
  expect_error(rcom(5, lambda = -3, nu = 0.7))
  
  ## negative nu
  expect_error(rcom(5, mu = 0.72, nu = -0.7))
})

