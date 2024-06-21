#' Create a table containing three validation metrics for clustering: Purity, F-measure,
#' Rand Index (RI) and Adjusted Rand Index (ARI)
#'
#' @param true_labels Atomic vector with the true labels of the data.
#' @param clusters The clusters predicted by the clustering method.
#' @param digits Number of digits for rounding.
#'
#' @return A \code{list} containing values for Purity, F-measure, RI and ARI.
#'
#' @examples
#' set.seed(1221)
#' vars <- list(c("dtaEI", "dtaMEI"))
#' data <- sim_model_ex1()
#' true_labels <- c(rep(1, 50), rep(2, 50))
#' data_ind <- generate_indices(data)
#' clus_kmeans <- clustInd_kmeans(data_ind, vars)
#' cluskmeans_mahalanobis_dtaEIdtaMEI <- clus_kmeans$kmeans_mahalanobis_dtaEIdtaMEI$cluster
#' clustering_validation(cluskmeans_mahalanobis_dtaEIdtaMEI, true_labels)
#'
#' @export
clustering_validation <- function(clusters, true_labels, digits = 4) {
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

  # Purity
  res_purity <- sum(apply(tbl, 1, max)) / length(true_labels)

  # True positives(tp), false positives(fp), true negatives(tn),
  # and false negatives(fn)

  # (needed for calculating RI)
  tp_plus_fp <- sum(choose(rowSums(tbl), 2))
  tp_plus_fn <- sum(choose(colSums(tbl), 2))
  tp <- sum(choose(as.vector(as.matrix(tbl)), 2))
  fp <- tp_plus_fp - tp
  fn <- tp_plus_fn - tp
  tn <- choose(sum(as.vector(as.matrix(tbl))), 2) - tp - fp - fn

  # Precision and recall (needed for calculating Fmeasure)
  prec <- tp / (tp + fp) # Precision
  rec <- tp / (tp + fn) # Recall

  # F1 score
  res_fmeasure <- 2 * ((prec * rec) / (prec + rec))

  # RI
  res_ri <- (tp + tn) / (tp + fp + fn + tn)

  # Adjusted Rand Index (ARI)
  n_pairs <- choose(length(true_labels),2)
  nij_sum <- sum(choose(tbl, 2))
  ai_sum <- sum(choose(rowSums(tbl),2))
  bj_sum <- sum(choose(colSums(tbl),2))
  aibj_term <- (ai_sum*bj_sum)/n_pairs

  res_ari <- (nij_sum-aibj_term)/(0.5*(ai_sum+bj_sum)-aibj_term)

  list(
    Purity = round(res_purity, digits),
    Fmeasure = round(res_fmeasure, digits),
    RI = round(res_ri, digits),
    ARI = round(res_ari, digits)
  )

}

