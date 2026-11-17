test_that("stream_salamanders has the documented structure", {
  e <- new.env()
  data("stream_salamanders", package = "flexCountReg", envir = e)
  d <- e$stream_salamanders
  expect_s3_class(d, "data.frame")
  expect_identical(dim(d), c(644L, 9L))
  expect_identical(names(d), c("site", "mined", "cover", "sample", "DOP",
                               "Wtemp", "DOY", "spp", "count"))
  expect_identical(nlevels(d$site), 23L)
  expect_true(all(table(d$site) == 28L))
  expect_identical(levels(d$mined), c("no", "yes"))
  expect_identical(levels(d$spp),
                   c("GP", "PR", "DM", "EC-A", "EC-L", "DES-L", "DF"))
  expect_identical(sort(unique(d$sample)), 1:4)
  expect_false(anyDuplicated(d[c("site", "sample", "spp")]) > 0)
  expect_type(d$count, "integer")
  expect_true(all(d$count >= 0))
  expect_false(anyNA(d))
  expect_equal(sum(d$count), 852L)
  # mined is constant within site, sample covariates within site x sample
  expect_true(all(tapply(as.integer(d$mined), d$site,
                         function(v) length(unique(v))) == 1))
  expect_true(all(tapply(d$DOY, paste(d$site, d$sample),
                         function(v) length(unique(v))) == 1))
})
