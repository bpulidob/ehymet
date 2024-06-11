#' Clustering using Epigraph and Hypograph indices
#'
#' It creates a multivariate dataset containing
#' the epigraph, hypograph and/or its modified versions on the curves and derivatives
#' and then perform hierarchical clustering, kmeans, kernel kmeans, and spectral clustering
#'
#' @param curves Dataset containing the curves to apply a clustering algorithm.
#' The functional dataset can be one dimensional (\eqn{n \times p}) where \eqn{n} is the number of
#' curves and \code{p} the number of time points, or multidimensional (\eqn{n \times p \times q}) where \eqn{q}
#' represents the number of dimensions in the data
#' @param vars_combinations If \code{list},  each element of the list should be an
#' atomic \code{vector} of strings with the
#' names of the variables. Combinations with non-valid variable names will be discarded.
#' If the list is non-named, the names of the variables are set to
#' vars1, ..., varsk, where k is the number of elements in \code{vars_combinations}.
#' If not provided, generic combinations of variables will be used. They will not be
#' the same for uni-dimensional and multi-dimensional problems.
#' @param clustering_methods character vector specifying at least one of the following
#' clustering methods to be computed: "hierarch", "kmeans", "kkmeans" or "spc".
#' @param k Number of basis functions for the B-splines. If equals to \code{0}, the number
#' of basis functions will be automatically selected.
#' @param bs A two letter character string indicating the (penalized) smoothing
#' basis to use. See \code{\link{smooth.terms}}.
#' @param l_method_hierarch \code{list} of clustering methods for hierarchical
#' clustering.
#' @param l_dist_hierarch \code{list} of distances for hierarchical clustering.
#' @param l_dist_kmeans \code{list} of distances for kmeans clustering.
#' @param l_kernel \code{list} of kernels for kkmeans or spc.
#' @param grid Atomic vector of type numeric with two elements: the lower limit and the upper
#' limit of the evaluation grid. If not provided, it will be selected automatically.
#' @param n_clusters Number of clusters to generate.
#' @param true_labels Numeric vector of true labels for validation. If provided,
#' more metrics are computed in the final result.
#' @param only_best \code{logical} value. If \code{TRUE} and \code{true_labels}
#' is provided, the function will return only the result for the best clustering
#' method based on the Rand Index. Defaults to \code{FALSE}.
#' @param verbose If \code{TRUE}, the function will print logs for about the execution of
#' some clustering methods. Defaults to \code{FALSE}.
#' @param n_cores Number of cores to do parallel computation. 1 by default,
#' which mean no parallel execution. Must be an integer number greater than 1.
#'
#' @return A \code{list} containing the clustering partition for each method and indices
#' combination and, if \code{true_labels} is provided a data frame containing the time elapsed for obtaining a
#' clustering partition of the indices dataset for each methodology. Also, the number of
#' generated clusters and the combinations of variables used can be seen as attributes
#' of this object.
#'
#' @examples
#' vars_combinations <- list(c("dtaEI", "dtaMEI"), c("dtaHI", "dtaMHI"))
#' curves <- sim_model_ex1(n = 10)
#' EHyClus(curves, vars_combinations = vars_combinations)
#'
#' @export
EHyClus <- function(curves, vars_combinations, k = 30, n_clusters = 2, bs = "cr",
                    clustering_methods = c("hierarch", "kmeans", "kkmeans", "spc"),
                    l_method_hierarch = c("single", "complete", "average", "centroid", "ward.D2"),
                    l_dist_hierarch = c("euclidean", "manhattan"),
                    l_dist_kmeans = c("euclidean", "mahalanobis"),
                    l_kernel = c("rbfdot", "polydot"),
                    grid,
                    true_labels = NULL, only_best = FALSE, verbose = FALSE, n_cores = 1) {
  if (length(dim(curves)) > 3 || length(dim(curves)) < 2) {
    stop("'curves' should be 2-dimensional or 3-dimensional", call. = FALSE)
  }

  if (!missing(vars_combinations) && !is.list(vars_combinations) && !is.numeric(vars_combinations)) {
    stop("input 'vars_combinations' must be a list", call. = FALSE)
  }

  if (!missing(vars_combinations) && is.list(vars_combinations) && !length(vars_combinations)) {
    stop("input 'vars_combinations' is an empty list", call. = FALSE)
  }

  if (!is.null(true_labels) && length(true_labels) != dim(curves)[1]) {
    stop("'true labels' should have the same length as the number of curves", call. = FALSE)
  }

  if (!is.numeric(k) || k %% 1 != 0) {
    stop("'k' should be an integer number", call. = FALSE)
  }

  if (!missing(grid)) {
    if (length(grid) != 2 && !is.numeric(grid)) {
      stop("'grid' should be a numeric atomic vector with two elements", call. = FALSE)
    } else if (diff(grid) <= 0) {
      stop("the second element of 'grid' should be greater than the first one", call. = FALSE)
    } else if (grid[1] < 1) {
      stop("the first element of 'grid' should be equal or greater than one", call. = FALSE)
    }
  }

  if (missing(vars_combinations)) {
    vars_combinations <-
      generic_vars_combinations(length(dim(curves)) == 3)
  }

  # list that maps each clustering method to its corresponding function
  default_clustering_methods <- list(
    "hierarch" = clustInd_hierarch,
    "kmeans"   = clustInd_kmeans,
    "kkmeans"  = clustInd_kkmeans,
    "spc"      = clustInd_spc
  )

  # Constants definition
  METHOD_HIERARCH <- c("single", "complete", "average", "centroid", "ward.D2")
  DIST_HIERARCH <- c("euclidean", "manhattan")
  DIST_KMEANS <- c("euclidean", "mahalanobis")
  KERNEL <- c("rbfdot", "polydot")
  METHOD_SVC <- c("kmeans", "kernkmeans")
  CLUSTERING_METHODS <- names(default_clustering_methods)

  check_list_parameter(clustering_methods, CLUSTERING_METHODS, "clustering_method")
  check_list_parameter(l_method_hierarch, METHOD_HIERARCH, "l_method_hierarch")
  check_list_parameter(l_dist_hierarch, DIST_HIERARCH, "l_dist_hierarch")
  check_list_parameter(l_dist_kmeans, DIST_KMEANS, "l_dist_kmeans")
  check_list_parameter(l_kernel, KERNEL, "l_kernel")

  n_cores <- check_n_cores(n_cores)

  # Generate the dataset with the indices

  generate_indices_parameters <- list(
    curves  = curves,
    k       = k,
    bs      = bs,
    indices = c("EI", "HI", "MEI", "MHI")
  )

  if (k) {
    generate_indices_parameters[["k"]] <- k
  }

  if (!missing(grid)) {
    generate_indices_parameters[["grid"]] <- grid
  }

  ind_curves <- do.call(generate_indices, generate_indices_parameters)

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
    "n_cores"           = n_cores
  )

  cluster <- list()
  for (method in clustering_methods) {
    method_args <- switch(method,
      "hierarch" = append(common_clustering_arguments, list(method_list = l_method_hierarch, dist_vector = l_dist_hierarch)),
      "kmeans"   = append(common_clustering_arguments, list(dist_vector = l_dist_kmeans)),
      "kkmeans"  = append(common_clustering_arguments, list(kernel_list = l_kernel)),
      "spc"      = append(common_clustering_arguments, list(kernel_list = l_kernel))
    )

    cluster_res <- tryCatch({
      if (verbose) {
        do.call(default_clustering_methods[[method]], method_args)
      } else {
        suppressMessages(quiet(do.call(default_clustering_methods[[method]], method_args)))
      }
    },
    error = function(x) NA
    )

    if (!all(is.na(cluster_res))) {
      cluster[[method]] <- cluster_res
    }

  }

  if (!is.null(true_labels)) {
    methods <- c()
    metrics <- data.frame(Purity = numeric(0), Fmeasure = numeric(0), RI = numeric(0), Time = numeric(0))
    for (clustering_method in names(cluster)) {
      for (method in names(cluster[[clustering_method]])) {
        methods <- c(methods, method)
        metrics <- rbind(
          metrics,
          c(cluster[[clustering_method]][[method]][["valid"]], cluster[[clustering_method]][[method]][["time"]])
        )
      }
    }
    names(metrics) <- c("Purity", "Fmeasure", "RI", "Time")
    rownames(metrics) <- methods

    metrics <- metrics[order(metrics$RI, decreasing = TRUE),]

    if (only_best) {
      metrics <- metrics[1,]
      best_clustering <- NA

      for (clustering_method in cluster) {
        if (rownames(metrics)[[1]] %in% names(clustering_method)) {
          best_clustering <- clustering_method[rownames(metrics)[[1]]]
          break
        }
      }

      cluster <- best_clustering
    }

    result <- list("cluster" = cluster, "metrics" = metrics)
  } else {
    result <- list("cluster" = cluster)
  }

  class(result) <- c("EHyClus", class(result))

  attr(result, "n_clusters") <- n_clusters
  attr(result, "vars_combinations") <- vars_combinations

  result
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

  vars_empty <- c()
  vars_invalid_name <- c()
  vars_almost_singular <- c()


  for (i in seq_along(vars_combinations)) {
    if (length(vars_combinations[[i]]) == 0) {
      vars_combinations_to_remove <- c(vars_combinations_to_remove, i)
      vars_empty <- c(vars_empty, i)

      next
    }

    if (length(vars_combinations[[i]]) == 1) {
      warning(paste0(
        "Combination of varaibles '", vars_combinations[[i]],
        "' with index ", i, " is only one variable, which ",
        "does not have much sense in this context..."
      ))
    }

    if (!all(vars_combinations[[i]] %in% names(ind_curves))) {
      vars_combinations_to_remove <- c(vars_combinations_to_remove, i)
      vars_invalid_name <- c(vars_invalid_name, i)

      next
    }

    if (det(stats::var(ind_curves[, vars_combinations[[i]]])) == 0) {
      vars_combinations_to_remove <- c(vars_combinations_to_remove, i)
      vars_almost_singular <- c(vars_almost_singular, i)
    }
  }

  if (length(vars_empty)) {
    warning(paste(
      "Index/indices ", paste0(vars_empty, collapse = ", "), "of 'vars_combinations' is/are empty.",
      "Removing them..."
    ))
  }

  if (length(vars_invalid_name)) {
    warning(paste(
      "Invalid variable name in 'vars_combinations' for index/indices ",
      paste0(vars_invalid_name, collapse = ", "),
      ". Removing them..."
    ))
  }

  if (length(vars_almost_singular)) {
    warning(paste(
      "Combination/s of variables with index/indices", paste0(vars_almost_singular, collapse = ", "),
      "is/are singular or almost singular. Removing them..."
    ))
  }

  if (length(vars_combinations_to_remove)) {
    warning(paste(
      "Combination/s of variable/s with index", paste0(vars_combinations_to_remove, collapse = ", "),
      "are not valid. Excluding them from any computation..."
    ))
  }

  if (length(vars_combinations_to_remove) == length(vars_combinations)) {
    stop("none of the combinations provided in 'vars_combinations' is valid.", call. = FALSE)
  }

  vars_combinations_to_remove
}


#' Return the default combinations of variables
#'
#' @param multidimensional \code{logical} determining if the vars_combinations
#' are for a uni-dimensional or multi-dimensional dataset.
#'
#' @return \code{list} with combinations of variables.
#'
#' @noRd
generic_vars_combinations <- function(multidimensional = TRUE) {
  if (multidimensional) {
    return(
      list(
        c("dtaMEI","dtaMHI"),
        c("ddtaMEI","ddtaMHI"),
        c("d2dtaMEI","d2dtaMHI"),
        c("dtaMEI", "dtaMHI", "ddtaMEI", "ddtaMHI"),
        c("dtaMEI", "dtaMHI", "d2dtaMEI", "d2dtaMHI"),
        c("ddtaMEI", "ddtaMHI", "d2dtaMEI", "d2dtaMHI"),
        c("dtaMEI", "dtaMHI", "ddtaMEI", "ddtaMHI", "d2dtaMEI","d2dtaMHI"),
        c("dtaMEI","ddtaMEI"),
        c("dtaMEI","d2dtaMEI"),
        c("ddtaMEI","d2dtaMEI"),
        c("dtaMEI","ddtaMEI","d2dtaMEI"),
        c("dtaMHI","ddtaMHI"),
        c("dtaMHI","d2dtaMHI"),
        c("ddtaMHI","d2dtaMHI"),
        c("dtaMHI","ddtaMHI","d2dtaMHI")
      )
    )
  } else {
    return(
      list(
        c("dtaEI","dtaHI"),
        c("ddtaEI","ddtaHI"),
        c("d2dtaEI","d2dtaHI"),
        c(c("dtaEI","dtaHI"),c("ddtaEI","ddtaHI")),
        c(c("dtaEI","dtaHI"),c("d2dtaEI","d2dtaHI")),
        c(c("ddtaEI","ddtaHI"),c("d2dtaEI","d2dtaHI")),
        c(c(c("dtaEI","dtaHI"),c("ddtaEI","ddtaHI")),c("d2dtaEI","d2dtaHI")),
        c("dtaMEI","ddtaMEI"),
        c("dtaMEI","d2dtaMEI"),
        c("ddtaMEI","d2dtaMEI"),
        c(c("dtaMEI","ddtaMEI"),"d2dtaMEI"),
        c(c("dtaEI","dtaHI"),"dtaMEI"),
        c(c("ddtaEI","ddtaHI"),"ddtaMEI"),
        c(c("d2dtaEI","d2dtaHI"),"d2dtaMEI"),
        c(c(c("dtaEI","dtaHI"),c("ddtaEI","ddtaHI")),c("dtaMEI","ddtaMEI")),
        c(c(c("dtaEI","dtaHI"),c("d2dtaEI","d2dtaHI")),c("dtaMEI","d2dtaMEI")),
        c(c(c("ddtaEI","ddtaHI"),c("d2dtaEI","d2dtaHI")),c("ddtaMEI","d2dtaMEI")),
        c(c(c(c("dtaEI","dtaHI"),c("ddtaEI","ddtaHI")),c("d2dtaEI","d2dtaHI")), c(c("dtaMEI","ddtaMEI"),"d2dtaMEI"))
      )
    )
  }
}
