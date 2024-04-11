test_that("'EI' works for a matrix", {
  curves <- array(c(1,0,1,0,0,1,0,1), dim=c(2, 4))
  expect_equal(EI(curves), 1-c(0.5,0.5))
})

test_that("'EI' works for a 3-dimensional array", {
  curves <- array(c(1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1), dim=c(2, 4, 2))
  expect_equal(EI(curves), 1-c(0.5,0.5))
})

test_that("'EI' does not work for a 4-dimensional array and (for example) for the print function", {
  curves <- array(c(1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1), dim=c(2, 2, 2, 2))
  expect_error(EI(curves))
  expect_error(EI(print))
})

test_that("'HI' works for a matrix", {
  curves <- array(c(1,0,1,0,0,1,0,1), dim=c(2, 4))
  expect_equal(HI(curves), c(0.5,0.5))
})

test_that("'HI' works for a 3-dimensional array", {
  curves <- array(c(1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1), dim=c(2, 4, 2))
  expect_equal(HI(curves), c(0.5,0.5))
})

test_that("'HI' does not work for a string", {
  expect_error(HI("I will throw an error"))
})

test_that("MEI works", {
  curves <- array(c(1,0,1,0,0,1,0,1), dim=c(2,4))
  expect_equal(MEI(curves), 1-c(0.75,0.75))
})

test_that("MHI works", {
  curves <- array(c(1,0,1,0,0,1,0,1), dim=c(2,4))
  expect_equal(MHI(curves), c(0.75,0.75))
})

test_that("mulMEI works", {
  curves <- array(c(1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1), dim=c(2,4,2))
  expect_equal(mulMEI(curves), 1-c(0.75,0.75))
})

test_that("mulMHI works", {
  curves <- array(c(1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1), dim=c(2,4,2))
  expect_equal(mulMHI(curves), c(0.75,0.75))
})
