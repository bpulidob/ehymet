test_that("the 'get_result_names' function is working as expexted for 'clustInd_hierarch'", {
  data <- sim_model_ex1()
  data_ind <- ind(data, t=seq(0, 1, length = 30))
  res <- clustInd_hierarch(ind_data = data_ind,
                           vars_list = list(c("dtaEI", "dtaMEI"), c("dtaHI")),
                           method_list = c("single", "complete"))
  expected <- c(
    "hierarch_single_euclidean_dtaEIdtaMEI",
    "hierarch_single_euclidean_dtaHI",
    "hierarch_complete_euclidean_dtaEIdtaMEI",
    "hierarch_complete_euclidean_dtaHI",
    "hierarch_single_manhattan_dtaEIdtaMEI",
    "hierarch_single_manhattan_dtaHI",
    "hierarch_complete_manhattan_dtaEIdtaMEI",
    "hierarch_complete_manhattan_dtaHI"
  )

  expect_equal(names(res), expected)
})


test_that("the 'get_result_names' function is working as expexted for 'clustInd_kmeans'", {
  data <- sim_model_ex1()
  data_ind <- ind(data, t=seq(0, 1, length = 30))
  res <- clustInd_kmeans(ind_data = data_ind,
                         vars_list = list(c("dtaMEI"), c("dtaHI")),
                         dist_list = c("euclidean", "mahalanobis"))
  expected <- c(
    "kmeans_euclidean_dtaMEI",
    "kmeans_euclidean_dtaHI",
    "kmeans_mahalanobis_dtaMEI",
    "kmeans_mahalanobis_dtaHI"
  )

  expect_equal(names(res), expected)
})

