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
#' clustering methods to be computed: "hierarch", "kmeans", "kkmeans" or "spc".
#' @param k Number of basis functions for the B-splines. If equals to \code{0}, the number
#' of basis functions will be automatically selected.
#' @param bs A two letter character string indicating the (penalized) smoothing
#' basis to use. See \code{\link{smooth.terms}}.
#' @param indices Names of the indices that need to be generated. They should be
#' one or more between "EI", "HI", "MEI" and "MHI". Depending on the dimension on the data,
#' its one-dimensional version or multi-dimensional version is computed.
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
#' @param verbose If \code{TRUE}, the function will print logs for about the execution of
#' some clustering methods. Defaults to \code{FALSE}.
#' @param n_cores Number of cores to do parallel computation. 1 by default,
#' which mean no parallel execution. Must be an integer number greater than 1.
#'
#' @return A \code{list} containing the clustering partition for each method and indexes
#' combination and, if \code{true_labels} is provided a data frame containing the time elapsed for obtaining a
#' clustering partition of the indexes dataset for each methodology.
#'
#' @examples
#' vars_combinations <- list(c("dtaEI", "dtaMEI"), c("dtaHI", "dtaMHI"))
#' curves <- sim_model_ex1(n = 10)
#' EHyClus(curves, vars_combinations = vars_combinations)
#'
#' @export
EHyClus <- function(curves, vars_combinations = 1, k = 30, n_clusters = 2, bs = "cr",
                    clustering_methods = c("hierarch", "kmeans", "kkmeans", "spc"),
                    indices = c("EI", "HI", "MEI", "MHI"),
                    l_method_hierarch = c("single", "complete", "average", "centroid", "ward.D2"),
                    l_dist_hierarch = c("euclidean", "manhattan"),
                    l_dist_kmeans = c("euclidean", "mahalanobis"),
                    l_kernel = c("rbfdot", "polydot"),
                    grid,
                    true_labels = NULL, verbose = FALSE, n_cores = 1) {
  if (!is.list(vars_combinations) && !is.numeric(vars_combinations)) {
    stop("input 'vars_combinations' must be a list or an integer number", call. = FALSE)
  }

  if (is.list(vars_combinations) && !length(vars_combinations)) {
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

  # list that maps each clustering method to its corresponding function
  default_clustering_methods <- list(
    "hierarch" = clustInd_hierarch,
    "kmeans"   = clustInd_kmeans,
    "kkmeans"  = clustInd_kkmeans,
    "spc"      = clustInd_spc
  )

  # Constants definition
  INDICES <- c("EI", "HI", "MEI", "MHI")
  METHOD_HIERARCH <- c("single", "complete", "average", "centroid", "ward.D2")
  DIST_HIERARCH <- c("euclidean", "manhattan")
  DIST_KMEANS <- c("euclidean", "mahalanobis")
  KERNEL <- c("rbfdot", "polydot")
  METHOD_SVC <- c("kmeans", "kernkmeans")
  CLUSTERING_METHODS <- names(default_clustering_methods)

  check_list_parameter(clustering_methods, CLUSTERING_METHODS, "clustering_method")
  check_list_parameter(indices, INDICES, "indices")
  check_list_parameter(l_method_hierarch, METHOD_HIERARCH, "l_method_hierarch")
  check_list_parameter(l_dist_hierarch, DIST_HIERARCH, "l_dist_hierarch")
  check_list_parameter(l_dist_kmeans, DIST_KMEANS, "l_dist_kmeans")
  check_list_parameter(l_kernel, KERNEL, "l_kernel")

  n_cores <- check_n_cores(n_cores)

  # Generate the dataset with the indexes

  generate_indices_parameters <- list(
    curves  = curves,
    k       = k,
    bs      = bs,
    indices = indices
  )

  if (k) {
    generate_indices_parameters[["k"]] <- k
  }

  if (!missing(grid)) {
    generate_indices_parameters[["grid"]] <- grid
  }

  ind_curves <- do.call(generate_indices, generate_indices_parameters)

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
    "n_cores"           = n_cores
  )

  cluster <- list()
  for (method in clustering_methods) {
    method_args <- switch(method,
      "hierarch" = append(common_clustering_arguments, list(method_list = l_method_hierarch, dist_list = l_dist_hierarch)),
      "kmeans"   = append(common_clustering_arguments, list(dist_list = l_dist_kmeans)),
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
        metrics <- rbind(
          metrics,
          c(cluster[[clustering_method]][[method]][["valid"]], cluster[[clustering_method]][[method]][["time"]])
        )
      }
    }
    names(metrics) <- c("Purity", "Fmeasure", "RI", "Time")
    rownames(metrics) <- methods

    result <- list("cluster" = cluster, "metrics" = metrics)
  } else {
    result <- list("cluster" = cluster)
  }

  class(result) <- c("EHyClus", class(result))

  attr(result, "n_clusters") <- n_clusters
  attr(result, "vars_combinations") <- vars_combinations

  result
}

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
