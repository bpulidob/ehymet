test_that("the parameter checking is working as expected", {
  data <- ehyclus_example_data()
  curves <- data$curves
  vars_combinations <- data$vars_combinations

  # Repeated element in 'indices'
  expect_error(
    EHyClus(curves, vars_combinations, indices = c("EI", "EI", "HI", "MEI", "MHI"))
  )

  # Empty list in 'l_method_hierarch'
  expect_error(
    EHyClus(curves, vars_combinations, l_method_hierarch = c())
  )

  # Non-valid argument in 'l_dist_hierarch'
  expect_error(
    EHyClus(curves, vars_combinations, l_dist_hierarch = c("euclidean", "i_do_not_exist"))
  )

  # 'vars_combinations' not being a list (TIENE QUE SER LIST !!!!!)
  expect_error(
    EHyClus(curves, unlist(vars_combinations))
  )
})

test_that("the checking related to 'vars_combinations' is doing its work", {
  data <- ehyclus_example_data()
  curves <- data$curves

  vars_wrong_name <- list(c("dtaEI", "ddtaET", "d2dtaEI"))
  vars_singular   <- list(c("dtaEI", "ddtaEI", "d2dtaEI"))
  vars_wrong_and_correct <- list(c("dtaEI", "ddtaET", "d2dtaEI"), c("dtaEI", "dtaMEI"))

  # The case with only one variable should throw a warning but still works
  # expect_warning(EHyClus(curves = curves, vars_combinations = vars_but_only_one))
  # TODO Error: unable to find an inherited method for function 'affinMult' for signature 'kernel = "rbfkernel", x = "numeric"'
  # same happening sometimes in "the 'n_clusters' parameter is working as expected"


  # Both the case with wrong names and with a singular (or almost singular)
  # matrix should throw an error
  expect_error(EHyClus(curves = curves, vars_combinations = vars_wrong_name))
  expect_error(EHyClus(curves = curves, vars_combinations = vars_singular))

  # Also, an empty list should throw an error
  expect_error(EHyClus(curves = curves, vars_combinations = list()))

  # As one of the vars combinations is correct but there is also one incorrect,
  # this should work but a warning is expected
  expect_warning(EHyClus(curves = curves, vars_combinations = vars_wrong_and_correct))

})

# test_that("the 'n_clusters' parameter is working as expected", {
#   data <- ehyclus_example_data()
#   curves <- data$curves
#   vars_combinations <- data$vars_combinations
#
#   res <- EHyClus(curves, vars_combinations, n_clusters = 3)
#   expect_equal(
#     max(res$cluster[[1]][[1]]$cluster),
#     3
#   )
# })

