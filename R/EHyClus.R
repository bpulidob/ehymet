#' EHyClus method for clustering. It creates a multivariate dataset containing
#' the epigraph, hypograph and/or modified version on the curves and derivatives
#' and perform hierarchical clustering, kmeans, kernel kmeans, support vector
#' clustering and spectral clustering
#'
#' @param curves Dataset containing the curves to apply a clustering algorithm.
#' The functional dataset can be one dimensional (\eqn{n \times p}) where n is the number of
#' curves and p the number of time points, or multidimensional (\eqn{n \times p \times k}) where k
#' represents the number of dimensions in the data
#' @param grid_ll lower limit of the grid.
#' @param grid_ul upper limit of the grid.
#' @param vars_list \code{list} containing one or more combinations of indexes in
#' \code{ind_data}. If it is non-named, the names of the variables are set to
#' vars1, ..., varsk, where k is the number of elements in \code{vars_list}.
#' @param clustering_methods character vector specifying at least one of the following
#' clustering methods to be computed: "hierarch", "kmeans", "kkmeans", "svc", "spc".
#' @param nbasis Number of basis for the B-splines.
#' @param norder Order of the B-splines.
#' @param indices Names of the indices that need to be generated. They should be
#' one or more between 'EI', 'HI', 'MEI' and 'MHI'. Depending on the dimension on the data
#' they are calculated for one or multiple dimension
#' @param l_method_hierarch \code{list} of clustering methods for hierarchical
#' clustering
#' @param l_dist_hierarch \code{list} of distances for hierarchical clustering
#' @param l_dist_kmeans \code{list} of distances for kmeans clustering
#' @param l_kernel \code{list} of kernels
#' @param l_method_svc \code{list} of clustering methods for support vector clustering
#' @param n_clusters Number of clusters to create
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
#' @param colapse It is a boolean. If it is true a dataframe with metrics values
#' is generated. If \code{true_labels} is True the dataframe contains Purity,
#' F-measure, RI and Time, and if it is False, only Time.
#' @param verbose If \code{TRUE}, the function will print logs for about the execution of
#' some clustering methods. Defaults to \code{FALSE}.
#' @param num_cores Number of cores to do parallel computation. 1 by default,
#' which mean no parallel execution.
#' @param ... Ignored.
#'
#' @return A list containing the clustering partition for each method and indexes
#' combination and a data frame containing the time elapsed for obtaining a
#' clustering partition of the indexes dataset for each methodology
#'
#' @examples
#' vars1 <- c("dtaEI", "dtaMEI"); vars2 <- c("dtaHI", "dtaMHI")
#' varsl <- list(vars1, vars2)
#' data <- ehymet::sim_model_ex1()
#' EHyClus(data, varsl, grid_ll = 0, grid_ul = 1)
#'
#' @export
EHyClus <- function(curves, vars_list, nbasis = 30,  n_clusters = 2, norder = 4,
                    clustering_methods = c("hierarch", "kmeans", "kkmeans", "svc", "spc"),
                    indices            = c("EI", "HI", "MEI", "MHI"),
                    l_method_hierarch  = c("single", "complete", "average", "centroid", "ward.D2"),
                    l_dist_hierarch    = c("euclidean", "manhattan"),
                    l_dist_kmeans      = c("euclidean", "mahalanobis"),
                    l_kernel           = c("rbfdot", "polydot"),
                    l_method_svc       = c("kmeans", "kernkmeans"),
                    grid_ll = 0, grid_ul = 1,
                    true_labels = NULL, colapse = FALSE, verbose = FALSE, num_cores = 1, ...) {

  # vars_list TIENE QUE SER LIST !!!!!
  if (!is.list(vars_list)) {
    stop("input 'vars_list' must be a list.", call. = FALSE)
  }

  # list that maps each clustering method to its corresponding function
  default_clustering_methods <- list(
    "hierarch" = clustInd_hierarch,
    "kmeans"   = clustInd_kmeans,
    "kkmeans"  = clustInd_kkmeans,
    "svc"      = clustInd_svc,
    "spc"      = clustInd_spc
  )

  # Constants definition
  INDICES            <- c("EI", "HI", "MEI", "MHI")
  METHOD_HIERARCH    <- c("single", "complete", "average", "centroid", "ward.D2")
  DIST_HIERARCH      <- c("euclidean", "manhattan")
  DIST_KMEANS        <- c("euclidean", "mahalanobis")
  KERNEL             <- c("rbfdot", "polydot")
  METHOD_SVC         <- c("kmeans", "kernkmeans")
  CLUSTERING_METHODS <- names(default_clustering_methods)

  check_list_parameter(clustering_methods, CLUSTERING_METHODS, "clustering_method")
  check_list_parameter(indices, INDICES, "indices")
  check_list_parameter(l_method_hierarch, METHOD_HIERARCH, "l_method_hierarch")
  check_list_parameter(l_dist_hierarch, DIST_HIERARCH, "l_dist_hierarch")
  check_list_parameter(l_dist_kmeans, DIST_KMEANS, "l_dist_kmeans")
  check_list_parameter(l_kernel, KERNEL, "l_kernel")
  check_list_parameter(l_method_svc, METHOD_SVC, "l_method_svc")

  # Generate the dataset with the indexes
  ind_curves <- ind(curves, grid_ll = grid_ll, grid_ul = grid_ul, nbasis, norder, indices)

  # common arguments for all the clustering methods that are implemented
  # in the package
  common_clustering_arguments <- list(
    "ind_data"    = ind_curves,
    "vars_list"   = vars_list,
    "n_cluster"   = n_clusters,
    "true_labels" = true_labels,
    "colapse"     = colapse,
    "num_cores"   = num_cores
  )

  cluster <- list()
  for (method in clustering_methods) {
    method_args <- switch(method,
      "hierarch" = append(common_clustering_arguments, list(method_list = l_method_hierarch, dist_list = l_dist_hierarch)),
      "kmeans"   = append(common_clustering_arguments, list(dist_list   = l_dist_kmeans)),
      "kkmeans"  = append(common_clustering_arguments, list(kernel_list = l_kernel)),
      "svc"      = append(common_clustering_arguments, list(method_list = l_method_svc)),
      "spc"      = append(common_clustering_arguments, list(kernel_list = l_kernel))
    )

    cluster[[method]] <- if (verbose) {
      do.call(default_clustering_methods[[method]], method_args)
    } else {
      suppressMessages(quiet(do.call(default_clustering_methods[[method]], method_args)))
    }
  }

  if (colapse) {
    metrics <- do.call(rbind, sapply(cluster, "[[", "metrics"))
    result  <- list("cluster" = cluster, "metrics" = metrics)
  } else {
    result <- list("cluster" = cluster)
  }

  class(result) <- c("EHyClus", class(result))
  attr(result, "n_clusters") <- n_clusters

  result
}

print.EHyClus <- function(x, ...) {
  cat("Clustering methods used:", paste(names(x$cluster), collapse = ", "), "\n")
  cat("Number of clusters:", attr(x, "n_clusters"))
  cat("More and more and more things.........\n")
  cat("......................................")

  invisible(x)
}

