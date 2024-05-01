test_that("the parameter checking is working as expected", {
  data <- ehyclus_example_data()
  curves <- data$curves
  t <- data$t
  vars_combinations <- data$vars_combinations

  # Repeated element in 'indices'
  expect_error(
    EHyClus(curves, t, vars_combinations, indices = c("EI", "EI", "HI", "MEI", "MHI"))
  )

  # Empty list in 'l_method_hierarch'
  expect_error(
    EHyClus(curves, t, vars_combinations, l_method_hierarch = c())
  )

  # Non-valid argument in 'l_dist_hierarch'
  expect_error(
    EHyClus(curves, t, vars_combinations, l_dist_hierarch = c("euclidean", "i_do_not_exist"))
  )

  # 'vars_combinations' not being a list (TIENE QUE SER LIST !!!!!)
  expect_error(
    EHyClus(curves, t, unlist(vars_combinations))
  )
})

# test_that("the 'n_clusters' parameter is working as expected", {
#   data <- ehyclus_example_data()
#   curves <- data$curves
#   t <- data$t
#   vars_combinations <- data$vars_combinations
#
#   res <- EHyClus(curves, t, vars_combinations, n_clusters = 3)
#   expect_equal(
#     max(res$cluster[[1]]$cluster),
#     3
#   )
# })

