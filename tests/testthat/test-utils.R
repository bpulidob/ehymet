test_that("'funspline' works for single curve, unidimensional data", {
  t <- seq(0, 1, length.out = 100)
  curves <- matrix(sin(2*pi*t), nrow = 1)
  result <- funspline(curves)

  expect_equal(dim(result$smooth), c(1, length(t)))
  expect_equal(dim(result$deriv),  c(1, length(t)))
  expect_equal(dim(result$deriv2), c(1, length(t)))
})

test_that("'funspline' works for multiple curves, unidimensional data", {
  t <- seq(0, 1, length.out = 100)
  curves <- matrix(c(sin(2*pi*t), cos(2*pi*t)), nrow = 2)
  result <- funspline(curves)

  expect_equal(dim(result$smooth), c(2, length(t)))
  expect_equal(dim(result$deriv),  c(2, length(t)))
  expect_equal(dim(result$deriv2), c(2, length(t)))
})

test_that("'funspline' works for multiple curves, multidimensional data", {
  t <- seq(0, 1, length.out = 100)
  curves <- array(c(sin(2*pi*t), cos(2*pi*t), exp(t)), dim = c(2, length(t), 3))
  result <- funspline(curves)

  expect_equal(dim(result$smooth), c(2, length(t), 3))
  expect_equal(dim(result$deriv),  c(2, length(t), 3))
  expect_equal(dim(result$deriv2), c(2, length(t), 3))
})

test_that("'funspline' works with non-uniform time sequence", {
  # should it be working?
  t <- c(seq(0, 0.5, length.out = 50), seq(0.6, 1, length.out = 50))
  curves <- matrix(sin(2*pi*t), nrow = 1)
  result <- funspline(curves)

  expect_equal(dim(result$smooth), c(1, length(t)))
  expect_equal(dim(result$deriv),  c(1, length(t)))
  expect_equal(dim(result$deriv2), c(1, length(t)))
})

test_that("'funspline' throws an error with invalid input", {
  t <- seq(0, 1, length.out = 100)
  curves <- "invalid"

  expect_error(funspline(curves))
})

test_that("the 'check_n_cores' function works", {
  # the proper behavior of this function can only be tested on windows
  skip_on_os(c("mac", "linux", "solaris"))

  n_cores <- 542
  check_n_cores(n_cores)

  expect_equal(n_cores, 1)
})
