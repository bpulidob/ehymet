test_that("the parameter checking is working as expected", {
  data <- ehyclus_example_data()
  curves <- data$curves
  t <- data$t
  vars_list <- data$vars_list

  # Repeated element in 'indices'
  expect_error(
    EHyClus(curves, t, vars_list, indices = c("EI", "EI", "HI", "MEI", "MHI"))
  )

  # Empty list in 'l_method_hierarch'
  expect_error(
    EHyClus(curves, t, vars_list, l_method_hierarch = c())
  )

  # Non-valid argument in 'l_dist_hierarch'
  expect_error(
    EHyClus(curves, t, vars_list, l_dist_hierarch = c("euclidean", "i_do_not_exist"))
  )

  # 'vars_list' not being a list (TIENE QUE SER LIST !!!!!)
  expect_error(
    EHyClus(curves, t, unlist(vars_list))
  )
})

test_that("the 'n_clusters' parameter is working as expected", {
  data <- ehyclus_example_data()
  curves <- data$curves
  t <- data$t
  vars_list <- data$vars_list

  res <- EHyClus(curves, t, vars_list, n_clusters = 4)
  expect_equal(
    max(res$cluster$hierarch_single_euclidean_dtaEIdtaMEI$cluster),
    4
  )
})

