test_that("the parameter checking is working as expected", {
  set.seed(42)

  data <- ehyclus_example_data()
  curves <- data$curves
  vars_combinations <- data$vars_combinations

  # Empty list in 'l_method_hierarch'
  expect_error(
    EHyClus(curves, vars_combinations, l_method_hierarch = c())
  )

  # Non-valid argument in 'l_dist_hierarch'
  expect_error(
    EHyClus(curves, vars_combinations, l_dist_hierarch = c("euclidean", "i_do_not_exist"))
  )

  # 'vars_combinations' not being a list
  expect_error(
    EHyClus(curves, unlist(vars_combinations))
  )
})

test_that("the checking related to 'vars_combinations' is doing its work", {
  set.seed(42)

  data <- ehyclus_example_data()
  curves <- data$curves

  vars_wrong_name <- list(c("dtaEI", "ddtaET", "d2dtaEI"))
  vars_singular <- list(c("dtaEI", "ddtaEI", "d2dtaEI"))
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
  labels <- rep(c(1, 2), each = n)

  vars1 <- c("dtaMEI", "ddtaMEI")
  vars2 <- c("dtaMEI", "d2dtaMEI")

  curves <- sim_model_ex2(n = n)
  res <- EHyClus(curves, vars_combinations = list(vars1, vars2), true_labels = labels)

  expect_equal(names(res), c("cluster", "metrics"))
  expect_equal(dim(res$metrics), c(32, 4))
})

test_that("the 'generic_vars_combinations' is returning an object of the expected lenght", {
  curves_unidimensional <- sim_model_ex2(i_sim = 1)
  vars_combinations_unidimensional <-
    generic_vars_combinations(length(dim(curves_unidimensional)) == 3)
  expect_length(vars_combinations_unidimensional, 18)

  curves_multidimensional <- sim_model_ex2(i_sim = 3)
  vars_combinations_multidimensional <-
    generic_vars_combinations(length(dim(curves_multidimensional)) == 3)
  expect_length(vars_combinations_multidimensional, 15)
})

test_that("the 'EHyClus' function works without providing 'vars_combinations'", {
  curves <- sim_model_ex2(n = 5, i_sim = 3)
  res <- EHyClus(curves)

  vars_combinations <- attr(res, "vars_combinations")

  expect_length(vars_combinations, 15)
})

test_that("the 'only_best' parameter works", {
  set.seed(32)

  n <- 5
  labels <- rep(c(1, 2), each = n)

  vars1 <- c("dtaMEI", "ddtaMEI")
  vars2 <- c("dtaMEI", "d2dtaMEI")

  curves <- sim_model_ex2(n = n)
  res <- EHyClus(curves, vars_combinations = list(vars1, vars2), true_labels = labels,
                 only_best = TRUE)

  expect_equal(dim(res$metrics), c(1, 4))
  expect_equal(length(res$cluster), 1)
})
