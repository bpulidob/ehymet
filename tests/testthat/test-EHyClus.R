test_that("the parameter checking is working as expected", {
  set.seed(42)

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
  set.seed(42)

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
  expect_error(
    suppressWarnings(EHyClus(curves = curves, k = 10, vars_combinations = vars_wrong_name))
  )

  expect_error(
    suppressWarnings(EHyClus(curves = curves, k = 10, vars_combinations = vars_singular))
  )

  # Also, an empty list should throw an error
  expect_error(
    suppressWarnings(EHyClus(curves = curves, k = 10, vars_combinations = list()))
  )

  # As one of the vars combinations is correct but there is also one incorrect,
  # this should work but a warning is expected
  expect_warning(EHyClus(curves = curves, k = 10, vars_combinations = vars_wrong_and_correct))

})

test_that("the 'n_clusters' parameter is working as expected", {
  set.seed(42)

  data <- ehyclus_example_data(n = 10)
  curves <- data$curves
  vars_combinations <- data$vars_combinations

  res <- EHyClus(curves, vars_combinations = vars_combinations, n_clusters = 3)
  expect_equal(
    max(res$cluster[[1]][[1]]$cluster),
    3
  )
})

test_that("the 'n_cores' parameter is not breaking anything", {
  set.seed(42)

  data <- ehyclus_example_data(n = 10)
  curves <- data$curves
  vars_combinations <- data$vars_combinations

  expect_no_error(EHyClus(curves, vars_combinations = vars_combinations, n_cores = 2))
})

test_that("metrics are correctly created when 'true_labels' is given to 'EHyClus'", {
  n <- 100
  labels <- rep(c(1,2), each = n)

  vars1  <- c("dtaMEI","ddtaMEI")
  vars2  <- c("dtaMEI","d2dtaMEI")

  curves <- sim_model_ex2(n = n)
  res <- EHyClus(curves, vars_combinations = list(vars1, vars2), true_labels = labels)

  expect_equal(names(res), c("cluster", "metrics"))
  expect_equal(dim(res$metrics), c(32, 4))
})


test_that("the 'get_best_vars_combinations' is giving the expected results", {
  set.seed(42)

  curves  <- sim_model_ex1()
  grid <- c(1, 2)
  k  <- 30
  indices <- c("EI", "HI", "MEI", "MHI")

  expected_best_combinations <- list(
    c("dtaHI", "dtaMEI"),
    c("dtaHI", "dtaMHI")
  )

  ind_curves <- generate_indices(curves, grid = grid, k = k, indices = indices)

  expect_error(get_best_vars_combinations(ind_curves = ind_curves, top_n = -2))
  expect_error(get_best_vars_combinations(ind_curves = ind_curves, top_n = 1.432534534))

  top_n <- 2
  best_combinatons <- get_best_vars_combinations(ind_curves = ind_curves, top_n = top_n)

  expect_length(best_combinatons, top_n)

  expect_equal(best_combinatons, expected_best_combinations)
})

test_that("giving an integer number to 'vars_combinations' works", {
  set.seed(33)

  curves <- sim_model_ex1()

  # By default, only one combination should be used
  res <- EHyClus(curves, k = 30)
  expect_length(attr(res, "vars_combinations"), 1)

  # Changing the 'vars_combinations' parameter should increase the number of combinations
  res <- EHyClus(curves, vars_combinations = 2, k = 30)
  expect_length(attr(res, "vars_combinations"), 2)

  indices <- c("EI", "MHI")
  humongous_vars_combinations_number <- 9999999999
  max_vars_combinations_number <- 2^6 - 6 - 1
  res <- suppressWarnings(
    EHyClus(curves,
            vars_combinations = humongous_vars_combinations_number,
            indices = indices
    )
  )

  expect_true(length(attr(res, "vars_combinations")) <= max_vars_combinations_number)
})
