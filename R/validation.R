#' Create a table containing three validation metrics: Purity, F-measure
#' and Rand Index (RI)
#'
#' @param true_labels Vector of true labels for validation
#' @param clusters Number of clusters to create
#'
#' @return A table containing values for Purity, F-measure and RI
#'
#' @examples
#' set.seed(1221)
#' vars1 <- c("dtaEI", "dtaMEI")
#' data <- ehymet::sim_model_ex1()
#' true_labels <- c(rep(1, 50), rep(2, 50))
#' data_ind <- ehymet::ind(data)
#' clus_kmeans <- ehymet::clustInd_kmeans(data_ind, list(vars1))
#' cluskmeans_mahalanobis_dtaEIdtaMEI <- clus_kmeans$kmeans_mahalanobis_dtaEIdtaMEI$cluster
#' valid(true_labels, cluskmeans_mahalanobis_dtaEIdtaMEI)
#'
#' @export
valid <- function(true_labels, clusters) {
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
    stop("clusters should be a numeric vector", call. = FALSE)
  }

  if (length(true_labels) != length(clusters)) {
    stop("The length of the true_labels vector should be equal to the length of the clusters vector", call. = FALSE)
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
  res <- c(round(res_purity, 4), round(2 * ((prec * rec) / (prec + rec)), 4),
           round((tp + tn) / (tp + fp + fn + tn), 4))
  names(res) <- c("Purity", "Fmeasure", "RI")

  as.table(res)
}
