test_that("the 'get_result_names' function is working as expexted for 'clustInd_hierarch'", {
  data <- sim_model_ex1()
  data_ind <- generate_indices(data, grid_ll = 0, grid_ul = 1)
  res <- clustInd_hierarch(
    ind_data = data_ind,
    vars_combinations = list(c("dtaEI", "dtaMEI"), c("dtaMEI", "ddtaMEI", "d2dtaMEI")),
    method_list = c("single", "complete")
  )
  expected <- c(
    "hierarch_single_euclidean_dtaEIdtaMEI",
    "hierarch_single_euclidean_dtaMEIddtaMEId2dtaMEI",
    "hierarch_complete_euclidean_dtaEIdtaMEI",
    "hierarch_complete_euclidean_dtaMEIddtaMEId2dtaMEI",
    "hierarch_single_manhattan_dtaEIdtaMEI",
    "hierarch_single_manhattan_dtaMEIddtaMEId2dtaMEI",
    "hierarch_complete_manhattan_dtaEIdtaMEI",
    "hierarch_complete_manhattan_dtaMEIddtaMEId2dtaMEI"
  )

  expect_equal(names(res), expected)
})


test_that("the 'get_result_names' function is working as expexted for 'clustInd_kmeans'", {
  data <- sim_model_ex1()
  data_ind <- generate_indices(data, grid_ll = 0, grid_ul = 1)
  res <- clustInd_kmeans(
    ind_data = data_ind,
    vars_combinations = list(c("dtaEI", "dtaMEI"), c("dtaMEI", "ddtaMEI", "d2dtaMEI")),
    dist_vector = c("euclidean", "mahalanobis")
  )
  expected <- c(
    "kmeans_euclidean_dtaEIdtaMEI",
    "kmeans_euclidean_dtaMEIddtaMEId2dtaMEI",
    "kmeans_mahalanobis_dtaEIdtaMEI",
    "kmeans_mahalanobis_dtaMEIddtaMEId2dtaMEI"
  )

  expect_equal(names(res), expected)
})
