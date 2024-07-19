test_that("the external validation metrics are computted", {
  set.seed(33)
  vars <- list(c("dtaEI", "dtaMEI"))
  data <- sim_model_ex1()
  true_labels <- c(rep(1, 50), rep(2, 50))
  data_ind <- generate_indices(data)
  clus_kmeans <- clustInd_kmeans(data_ind, vars)
  cluskmeans_mahalanobis_dtaEIdtaMEI <- clus_kmeans$kmeans_mahalanobis_dtaEIdtaMEI$cluster
  validation_results <- clustering_validation(cluskmeans_mahalanobis_dtaEIdtaMEI, true_labels)

  expect_equal(names(validation_results), c("Purity", "Fmeasure", "RI", "ARI"))
  expect_length(validation_results, 4)
})

test_that("the internal validation metrics are computed", {
  curves <- sim_model_ex1(n = 10)
  vars_combinations <- list(c("dtaEI", "dtaMEI"), c("dtaHI", "dtaMHI"))
  res <- EHyClus(curves, vars_combinations = vars_combinations)

  expect_named(
    res$cluster$hierarch$hierarch_single_euclidean_dtaEIdtaMEI$internal_metrics,
    c("davies_bouldin", "dunn", "silhouette", "infomax")
  )
})
