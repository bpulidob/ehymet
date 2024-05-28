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
#' @param vars_combinations \code{integer} or \code{list}.
#' If \code{integer}, the method will automatically determine the best combinations
#' of variables. As many combinations will be selected as the value of the variable.
#' If \code{list},  each element of the list should be an atomic \code{vector} of strings with the
#' names of the variables. Combinations with non-valid variable names will be discarded.
#' If the list is non-named, the names of the variables are set to
#' vars1, ..., varsk, where k is the number of elements in \code{vars_combinations}.
#' Default to an \code{integer} with value \code{1}, i.e. it only uses the theoretically
#' best combination.
#' @param clustering_methods character vector specifying at least one of the following
#' clustering methods to be computed: "hierarch", "kmeans", "kkmeans", "spc".
#' @param nbasis Number of basis for the B-splines. If , it will
#' be automatically set.
#' @param bs A two letter chatacter string indicating the (penalized) smoothing
#' basis to use. See \code{\link{smooth.terms}}.
#' @param indices Names of the indices that need to be generated. They should be
#' one or more between 'EI', 'HI', 'MEI' and 'MHI'. Depending on the dimension on the data
#' they are calculated for one or multiple dimension
#' @param l_method_hierarch \code{list} of clustering methods for hierarchical
#' clustering
#' @param l_dist_hierarch \code{list} of distances for hierarchical clustering
#' @param l_dist_kmeans \code{list} of distances for kmeans clustering
#' @param l_kernel \code{list} of kernels
#' @param n_clusters Number of clusters to create
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
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
#' data <- sim_model_ex1()
#' EHyClus(data, varsl, grid_ll = 0, grid_ul = 1)
#'
#' @export
EHyClus <- function(curves, vars_combinations = 1, nbasis = 30, n_clusters = 2, bs = "cr",
                    clustering_methods = c("hierarch", "kmeans", "kkmeans", "spc"),
                    indices            = c("EI", "HI", "MEI", "MHI"),
                    l_method_hierarch  = c("single", "complete", "average", "centroid", "ward.D2"),
                    l_dist_hierarch    = c("euclidean", "manhattan"),
                    l_dist_kmeans      = c("euclidean", "mahalanobis"),
                    l_kernel           = c("rbfdot", "polydot"),
                    grid_ll = 0, grid_ul = 1,
                    true_labels = NULL, verbose = FALSE, num_cores = 1, ...) {

  if (!is.list(vars_combinations) && !is.numeric(vars_combinations)) {
    stop("input 'vars_combinations' must be a list or an integer number", call. = FALSE)
  }

  if (is.list(vars_combinations) && !length(vars_combinations)) {
    stop("input 'vars_combinations' is an empty list", call. = FALSE)
  }

  if (!is.null(true_labels) && length(true_labels) != dim(curves)[1]) {
    stop("'true labels' should have the same length as the number of curves", call. = FALSE)
  }

  if (!is.numeric(nbasis) || nbasis %% 1 != 0) {
    stop("'nbasis' should be an integer number", call. = FALSE)
  }

  # list that maps each clustering method to its corresponding function
  default_clustering_methods <- list(
    "hierarch" = clustInd_hierarch,
    "kmeans"   = clustInd_kmeans,
    "kkmeans"  = clustInd_kkmeans,
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

  # Generate the dataset with the indexes
  if (nbasis) {
    ind_curves <- generate_indices(curves, nbasis, grid_ll = grid_ll, grid_ul = grid_ul, bs = bs, indices = indices)
  } else {
    ind_curves <- generate_indices(curves, grid_ll = grid_ll, grid_ul = grid_ul, bs = bs, indices = indices)
  }

  if (!is.list(vars_combinations)) {
    max_n <- 2^length(ind_curves) - length(ind_curves) - 1 # power set - 1-variable combinations - empty set
    if (vars_combinations > max_n) {
      warning(paste0("The maximum number for 'vars_combinations' in this setting is ", max_n))
      vars_combinations <- max_n
    }

    vars_combinations <- get_best_vars_combinations(ind_curves, vars_combinations)
  }

  # Check for correct vars combinations
  vars_combinations_to_remove <- check_vars_combinations(vars_combinations, ind_curves)

  if (length(vars_combinations_to_remove)) {
    vars_combinations <- vars_combinations[-vars_combinations_to_remove]
  }


  # common arguments for all the clustering methods that are implemented
  # in the package
  common_clustering_arguments <- list(
    "ind_data"          = ind_curves,
    "vars_combinations" = vars_combinations,
    "n_cluster"         = n_clusters,
    "true_labels"       = true_labels,
    "num_cores"         = num_cores
  )

  cluster <- list()
  for (method in clustering_methods) {
    method_args <- switch(method,
      "hierarch" = append(common_clustering_arguments, list(method_list = l_method_hierarch, dist_list = l_dist_hierarch)),
      "kmeans"   = append(common_clustering_arguments, list(dist_list   = l_dist_kmeans)),
      "kkmeans"  = append(common_clustering_arguments, list(kernel_list = l_kernel)),
      "spc"      = append(common_clustering_arguments, list(kernel_list = l_kernel))
    )

    cluster[[method]] <- if (verbose) {
      do.call(default_clustering_methods[[method]], method_args)
    } else {
      suppressMessages(quiet(do.call(default_clustering_methods[[method]], method_args)))
    }
  }

  if (!is.null(true_labels)) {
    methods <- c()
    metrics <- data.frame(Purity = numeric(0), Fmeasure = numeric(0), RI = numeric(0), Time = numeric(0))
    for (clustering_method in names(cluster)) {
      for (method in names(cluster[[clustering_method]])) {
        methods <- c(methods, method)
        metrics <- rbind(metrics,
                         c(cluster[[clustering_method]][[method]][["valid"]], cluster[[clustering_method]][[method]][["time"]]))
      }
    }
    names(metrics) <- c("Purity", "Fmeasure", "RI", "Time")
    rownames(metrics) <- methods

    result <- list("cluster" = cluster, "metrics" = metrics)
  } else {
    result <- list("cluster" = cluster)
  }

  class(result) <- c("EHyClus", class(result))

  attr(result, "n_clusters")        <- n_clusters
  attr(result, "vars_combinations") <- vars_combinations

  result
}

# print.EHyClus <- function(x, ...) {
#   cat("Clustering methods used:", paste(names(x$cluster), collapse = ", "), "\n")
#   cat("Number of clusters:", attr(x, "n_clusters"))
#   cat("More and more and more things.........\n")
#   cat("......................................")
#
#   invisible(x)
# }


#' Search for the best combinations of variables
#'
#' @param ind_curves Dataset with indexes from a functional dataset in one or multiple
#' dimensions.
#' @param top_n Number of desired variable combinations.
#'
#' @return \code{top_n} combinations of variables
#'
#' @noRd
get_best_vars_combinations <- function(ind_curves, top_n) {
  if (top_n %% 1 != 0 || top_n < 1) {
    stop("'vars_combinations' must be an integer greater than 1", call. = FALSE)
  }

  vars <- names(ind_curves)
  all_vars_combinations <- do.call(c, lapply(2:length(vars), utils::combn, x = vars, simplify = FALSE))
  dets <- lapply(all_vars_combinations, function(combination) det(stats::cov(ind_curves[, combination])))

  best_n <- sort(unlist(dets), index.return = TRUE, decreasing = TRUE)$ix[1:top_n]

  all_vars_combinations[best_n]
}


#' Check all combinations of variables and found the non-valid ones
#'
#' @param vars_combinations \code{list} containing one or more combination of variables.
#' @param ind_curves dataset with indices from a functional dataset in one or multiple
#' dimensions.
#'
#' @return Atomic vector with the index of the non-valid combinations of variables.
#'
#' @noRd
check_vars_combinations <- function(vars_combinations, ind_curves) {
  vars_combinations_to_remove <- c()

  vars_empty           <- c()
  vars_invalid_name    <- c()
  vars_almost_singular <- c()


  for (i in seq_along(vars_combinations)) {
    if (length(vars_combinations[[i]]) == 0) {
      vars_combinations_to_remove <- c(vars_combinations_to_remove, i)
      vars_empty <- c(vars_empty, i)

      next
    }

    if (length(vars_combinations[[i]]) == 1) {
      warning(paste0("Combination of varaibles '", vars_combinations[[i]],
                     "' with index ", i, " is only one variable, which ",
                     "does not have much sense in this context...")
      )
    }

    if (!all(vars_combinations[[i]] %in% names(ind_curves))) {
      vars_combinations_to_remove <- c(vars_combinations_to_remove, i)
      vars_invalid_name <- c(vars_invalid_name, i)

      next
    }

    if (det(stats::var(ind_curves[,vars_combinations[[i]]])) == 0) {
      vars_combinations_to_remove <- c(vars_combinations_to_remove, i)
      vars_almost_singular <- c(vars_almost_singular, i)
    }
  }

  if (length(vars_empty)) {
    warning(paste("Index/indices ", paste0(vars_empty, collapse = ", "), "of 'vars_combinations' is/are empty.",
                  "Removing them..."))
  }

  if (length(vars_invalid_name)) {
    warning(paste("Invalid variable name in 'vars_combinations' for index/indices ",
                  paste0(vars_invalid_name, collapse = ", "),
                  ". Removing them..."))
  }

  if (length(vars_almost_singular)) {
    warning(paste("Combination/s of variables with index/indices", paste0(vars_almost_singular, collapse = ", "),
                  "is/are singular or almost singular. Removing them..."))
  }

  if (length(vars_combinations_to_remove)) {
    warning(paste("Combination/s of variable/s with index", paste0(vars_combinations_to_remove, collapse = ", "),
                  "are not valid. Excluding them from any computation..."))
  }

  if (length(vars_combinations_to_remove) == length(vars_combinations)) {
    stop("none of the combinations provided in 'vars_combinations' is valid.", call. = FALSE)
  }

  vars_combinations_to_remove
}
