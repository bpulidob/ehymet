test_that("'EI' works for a matrix", {
  curves <- array(c(1,0,1,0,0,1,0,1), dim=c(2,4))
  expect_equal(EI(curves), 1-c(0.5,0.5))
})

test_that("'EI' works for a 3-dimensional array", {
  curves <- array(c(1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1), dim=c(2,4,2))
  expect_equal(EI(curves), 1-c(0.5,0.5))
})

test_that("'EI' does not work for a 4-dimensional array and (for example) for the print function", {
  curves <- array(c(1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1), dim=c(2,2,2,2))
  expect_error(EI(curves))
  expect_error(EI(print))
})

test_that("'HI' works for a matrix", {
  curves <- array(c(1,0,1,0,0,1,0,1), dim=c(2,4))
  expect_equal(HI(curves), c(0.5,0.5))
})

test_that("'HI' works for a 3-dimensional array", {
  curves <- array(c(1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1), dim=c(2,4,2))
  expect_equal(HI(curves), c(0.5,0.5))
})

test_that("'HI' does not work for a string", {
  expect_error(HI("I will throw an error"))
})

test_that("'MEI' works for a matrix", {
  curves <- array(c(1,0,1,0,0,1,0,1), dim=c(2,4))
  expect_equal(MEI(curves), 1-c(0.75,0.75))
})

test_that("'MEI' works for a 3-dimensional array", {
  curves <- array(c(1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1), dim=c(2,4,2))
  expect_equal(MEI(curves), 1-c(0.75,0.75))
})

test_that("'MEI' does not work for a complex number", {
  expect_error(MEI(1 + 2i))
})

test_that("'MHI' works for a matrix", {
  curves <- array(c(1,0,1,0,0,1,0,1), dim=c(2,4))
  expect_equal(MHI(curves), c(0.75,0.75))
})

test_that("'MHI' works for a 3-dimensional array", {
  curves <- array(c(1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1), dim=c(2,4,2))
  expect_equal(MHI(curves), c(0.75,0.75))
})

test_that("'MEI' does not work for a boolean", {
  expect_error(MEI(FALSE))
})

test_that("'generate_indices' works for a matrix", {
  curves <- array(c(1,0,1,0,0,1,0,1), dim=c(2,4))

  res <- generate_indices(curves, indices = c("EI", "MHI"))

  expect_equal(dim(res), c(2,6))
  expect_equal(res[2,1], 0.5)
})

test_that("'generate_indices' works for a 3-dimnensional array", {
  curves <- array(c(1,2,3, 3,2,1, 5,2,3, 9,8,7, -1,-5,-6, 2,3,0, -1,0,2, -1,-2,0),
                  dim = c(3,4,2))

  res <- suppressWarnings(generate_indices(curves, indices = c("EI", "MHI")))

  expect_equal(dim(res), c(3,6))
  expect_equal(res[2,1], 0.666666666666666740681)
})
