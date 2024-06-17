#' Create a table containing three validation metrics for clustering: Purity, F-measure
#' and Rand Index (RI)
#'
#' @param true_labels Atomic vector with the true labels of the data.
#' @param clusters The clusters predicted by the clustering method.
#'
#' @return A \code{table} containing values for Purity, F-measure and RI.
#'
#' @examples
#' set.seed(1221)
#' vars1 <- c("dtaEI", "dtaMEI")
#' data <- ehymet::sim_model_ex1()
#' true_labels <- c(rep(1, 50), rep(2, 50))
#' data_ind <- generate_indices(data)
#' clus_kmeans <- ehymet::clustInd_kmeans(data_ind, list(vars1))
#' cluskmeans_mahalanobis_dtaEIdtaMEI <- clus_kmeans$kmeans_mahalanobis_dtaEIdtaMEI$cluster
#' clustering_validation(true_labels, cluskmeans_mahalanobis_dtaEIdtaMEI)
#'
#' @export
clustering_validation <- function(true_labels, clusters) {
  if (is.integer(true_labels)) {
    true_labels <- as.numeric(true_labels)
  }

  if (is.integer(clusters)) {
    clusters <- as.numeric(clusters)
  }

  if (!is.vector(true_labels) || !is.numeric(true_labels)) {
    stop("'true_labels' should be a numeric vector", call. = FALSE)
  }

  if (!is.vector(clusters) || !is.numeric(clusters)) {
    stop("'clusters' should be a numeric vector", call. = FALSE)
  }

  if (length(true_labels) != length(clusters)) {
    stop("the length of the 'true_labels' vector should be equal to the length of the clusters vector", call. = FALSE)
  }

  tbl <- table(clusters, true_labels) # contingency table
  conv_df <- as.data.frame.matrix(tbl)

  # Purity
  res_purity <- sum(apply(conv_df, 1, max)) / length(true_labels)

  # True positives(tp), false positives(fp), true negatives(tn),
  # and false negatives(fn)

  # (needed for calculating RI)
  tp_plus_fp <- sum(choose(rowSums(conv_df), 2))
  tp_plus_fn <- sum(choose(colSums(conv_df), 2))
  tp <- sum(choose(as.vector(as.matrix(conv_df)), 2))
  fp <- tp_plus_fp - tp
  fn <- tp_plus_fn - tp
  tn <- choose(sum(as.vector(as.matrix(conv_df))), 2) - tp - fp - fn

  # Precision and recall (needed for calculating Fmeasure)
  prec <- tp / (tp + fp) # Precision
  rec <- tp / (tp + fn) # Recall

  # (Purity, Fmeasure, RI)
  res <- c(
    round(res_purity, 4),
    round(2 * ((prec * rec) / (prec + rec)), 4),
    round((tp + tn) / (tp + fp + fn + tn), 4)
  )
  names(res) <- c("Purity", "Fmeasure", "RI")

  # Adjusted Rand Index (ARI)
  n <- length(true_labels)
  total_pairs <- choose(n, 2)

  a_i <- rowSums(conv_df)
  b_j <- colSums(conv_df)

  index_sum <- sum(choose(a_i, 2)) + sum(choose(b_j, 2))
  product_sum <- sum(choose(as.vector(as.matrix(conv_df)), 2))

  expected_ri <- (sum(choose(a_i, 2)) * sum(choose(b_j, 2))) / total_pairs
  max_ri <- 0.5 * (sum(choose(a_i, 2)) + sum(choose(b_j, 2)))

  ari <- (product_sum - expected_ri) / (max_ri - expected_ri)

  res <- c(res, round(ari, 4))
  names(res) <- c("Purity", "Fmeasure", "RI", "ARI")

  as.table(res)
}

