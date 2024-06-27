#' Perform hierarchical clustering for a given combination of indices,
#' method and distance
#'
#' @param ind_data Dataframe containing indices applied to the original data and
#' its first and second derivatives.
#' @param vars vector with a combinations of indices in \code{ind_data}.
#' @param method The aglomerative method to be used in hierarchical clustering.
#' @param dist The distance method to be used.
#' @param n_cluster Number of clusters to create.
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
#'
#' @return A \code{list} containing clustering results and execution time.
#'
#' @noRd
clustInd_hierarch_aux <- function(ind_data, vars, method = "single",
                                  dist = "euclidean", n_cluster = 2,
                                  true_labels = NULL) {
  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.")
  }

  # Check if variables exist in the data frame
  if (!all(vars %in% names(ind_data))) {
    stop("Invalid variable name.")
  }

  if (!all(method %in% c("single", "complete", "average", "centroid", "ward.D2"))) {
    stop("Invalid method name.")
  }

  if (!all(dist %in% c("euclidean", "manhattan"))) {
    stop("Invalid distance name.")
  }

  t0 <- Sys.time()

  # Perform hierarchical clustering
  d <- stats::dist(ind_data[, vars], method = dist)
  met <- stats::hclust(d, method = method)
  clus <- stats::cutree(met, k = n_cluster)

  t1 <- Sys.time()

  # Calculate execution time
  et <- data.frame(difftime(t1, t0, "secs"))

  if (is.null(true_labels)) {
    res <- list("cluster" = clus, "time" = as.numeric(et))
  } else {
    valid <- clustering_validation(clus, true_labels)
    res <- list("cluster" = clus, "valid" = valid, "time" = as.numeric(et))
  }

  return(res)
}

#' Perform hierarchical clustering for a different combinations of indices,
#' method and distance
#'
#' @param ind_data Dataframe containing indices applied to the original data and
#' its first and second derivatives. See \link{generate_indices}.
#' @param vars_combinations \code{list} containing one or more combinations of indices in
#' \code{ind_data}. If it is non-named, the names of the variables are set to
#' vars1, ..., varsk, where k is the number of elements in \code{vars_combinations}.
#' @param method_list \code{list} of clustering methods.
#' @param dist_vector \code{list} of distance metrics.
#' @param n_cluster number of clusters to generate.
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
#' @param n_cores Number of cores to do parallel computation. 1 by default,
#' which mean no parallel execution.
#'
#' @return A \code{list} containing hierarchical clustering results
#' for each configuration.
#'
#' @examples
#' vars1 <- c("dtaEI", "dtaMEI")
#' vars2 <- c("dtaHI", "dtaMHI")
#' data <- ehymet::sim_model_ex1()
#' data_ind <- generate_indices(data)
#' clustInd_hierarch(data_ind, list(vars1, vars2))
#'
#' @export
clustInd_hierarch <- function(ind_data, vars_combinations,
                              method_list = c(
                                "single", "complete", "average",
                                "centroid", "ward.D2"
                              ),
                              dist_vector = c("euclidean", "manhattan"),
                              n_cluster = 2, true_labels = NULL,
                              n_cores = 1) {
  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("input 'ind_data' must be a data frame.", call. = FALSE)
  }

  if (!is.list(vars_combinations)) {
    stop("input 'vars_combinations' must be a list.", call. = FALSE)
  }

  n_cores <- check_n_cores(n_cores)

  # Check for correct vars combinations
  vars_combinations_to_remove <- check_vars_combinations(vars_combinations, ind_data)

  if (length(vars_combinations_to_remove)) {
    vars_combinations <- vars_combinations[-vars_combinations_to_remove]
  }

  # Check if indices, methods and distances lists are provided
  if (!is.character(method_list) ||
    !is.character(dist_vector) || length(vars_combinations) == 0 ||
    length(method_list) == 0 || length(dist_vector) == 0) {
    stop("invalid 'method_list' or 'dist_vector'. Both must be non-empty character vectors.", call. = FALSE)
  }

  if (is.null(names(vars_combinations))) {
    names(vars_combinations) <- paste0("vars", seq_along(vars_combinations))
  }

  # Generate all the possible combinations of indices, methods and distances
  parameter_combinations <- expand.grid(
    vars = names(vars_combinations),
    method = method_list, distance = dist_vector
  )

  tl_null <- is.null(true_labels)

  n_comb <- nrow(parameter_combinations)

  result <- parallel::mclapply(1:n_comb, function(i) {
    vars <- vars_combinations[[parameter_combinations$vars[i]]]
    met <- parameter_combinations$method[i]
    dist <- parameter_combinations$distance[i]

    clustInd_hierarch_aux(ind_data, vars, met, dist, n_cluster, true_labels)
  }, mc.cores = n_cores)

  result_name <- get_result_names("hierarch", parameter_combinations, vars_combinations)
  names(result) <- result_name

  return(result)
}


init_kplus2 <- function(ind_data, n_cluster){
  init_center <- rep(0, n_cluster)
  # Choose the first center at random from data points
  init_center[1] <- sample(1:nrow(ind_data), 1)

  for(i in 2:n_cluster) {
    # For each data point x not choosen yet, compute D(x), the distance between
    # x and the nearest center that has already been chosen
    dist_center_min <- apply(ind_data,1,function(row){
      dist_centers = sapply(init_center[1:(i-1)],
                            function(center) {
                              sqrt(sum((row - ind_data[center,])^2))
                            })
      min(dist_centers)
    })

    # Choose one new data point at random as a new center, using a weighted
    # probability distribution where a point x is chosen with probability
    # proportional to D(x)^2
    probabilities <- dist_center_min^2/sum(dist_center_min)
    zero_prob = which(dist_center_min == 0)
    init_center[i] = sample(seq_along(probabilities)[-zero_prob], 1,
                            prob=probabilities[-zero_prob])
  }
  return(ind_data[init_center,])
}

#' k-means clustering with Mahalanobis distance
#'
#' @param ind_data Dataframe containing indices applied to the original data and
#' its first and second derivatives.
#' @param n_cluster Number of clusters to create.
#' @param init Initialization strategy.
#'
#' @return k-means clustering with Mahalanobis distance output
#'
#' @noRd
kmeans_mahal <- function(ind_data, n_cluster, init = "random") {
  # Check if input is numeric matrix or array
  if (!(is.numeric(ind_data) || is.matrix(ind_data) || is.array(ind_data) ||
    is.data.frame(ind_data))) {
    stop("Input must be a numeric matrix, array or data frame.")
  }

  # Check if the initialization method given can be used
  if (!init %in% c("random", "kmeans++")) {
    stop("Invalid initialization method.")
  }

  # Convert data to matrix
  if (!is.matrix(ind_data)) {
    data_matrix <- as.matrix(ind_data)
  } else {
    data_matrix <- ind_data
  }

  # Cholesky transformation
  c <- chol(stats::var(data_matrix))
  data_transform <- data_matrix %*% solve(c)

  # vector containing the clustering partition
  centers_init <- if (init == "random") n_cluster else
    init_kplus2(data_transform, n_cluster)
  km <- stats::kmeans(data_transform,
    centers = centers_init, iter.max = 1000,
    nstart = 100
  )$cluster
  return(km)
}


#' Perform kmeans clustering for a given combination of indices and distance
#'
#' @param ind_data Dataframe containing indices applied to the original data and
#' its first and second derivatives. See \link{generate_indices}.
#' @param vars vector with a combinations of indices in \code{ind_data}
#' @param dist The distance method to be used.
#' @param n_cluster Number of clusters to create.
#' @param true_labels Vector of true labels for validation.
#' (if it is not known true_labels is set to NULL)
#'
#' @return A list containing clustering results and execution time.
#'
#' @noRd
clustInd_kmeans_aux <- function(ind_data, vars, dist = "euclidean",
                                n_cluster = 2, init = "random",
                                true_labels = NULL) {
  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.")
  }

  # Check if variables exist in the data frame
  if (!all(vars %in% names(ind_data))) {
    stop("Invalid variable names.")
  }

  # Check if the distance given can be used
  if (!dist %in% c("euclidean", "mahalanobis")) {
    stop("Invalid distance.")
  }

  # Check if the initialization method given can be used
  if (!init %in% c("random", "kmeans++")) {
    stop("Invalid initialization method.")
  }

  t0 <- Sys.time()

  if (dist == "euclidean") {
    centers_init <- if (init == "random") n_cluster else
      init_kplus2(ind_data[, vars], n_cluster)
    clus <- stats::kmeans(ind_data[, vars],
      centers = centers_init,
      iter.max = 1000, nstart = 100
    )$cluster
  } else {
    clus <- kmeans_mahal(ind_data[, vars], n_cluster, init)
  }
  t1 <- Sys.time()
  t <- data.frame(difftime(t1, t0, "secs"))

  if (is.null(true_labels)) {
    res <- list("cluster" = clus, "time" = as.numeric(t))
  } else {
    valid <- clustering_validation(clus, true_labels)
    res <- list("cluster" = clus, "valid" = valid, "time" = as.numeric(t))
  }
  return(res)
}

#' K-means clustering with indices
#'
#' Perform k-means clustering for a different combinations of indices and
#' distances.
#'
#' @param ind_data Dataframe containing indices applied to the original data and
#' its first and second derivatives. See \link{generate_indices}.
#' @param vars_combinations \code{list} containing one or more combinations of indices in
#' \code{ind_data}. If it is non-named, the names of the variables are set to
#' vars1, ..., varsk, where k is the number of elements in \code{vars_combinations}.
#' @param dist_vector Atomic vector of distance metrics. The possible values are,
#' "euclidean", "mahalanobis" or both.
#' @param n_cluster Number of clusters to create.
#' @param true_labels Vector of true labels for validation.
#' (if it is not known true_labels is set to NULL)
#' @param n_cores Number of cores to do parallel computation. 1 by default,
#' which mean no parallel execution.
#' @return A list containing hierarchical clustering results
#' for each configuration
#'
#' @return A list containing kmeans clustering results for each configuration
#'
#' @examples
#' vars1 <- c("dtaEI", "dtaMEI")
#' vars2 <- c("dtaHI", "dtaMHI")
#' data <- ehymet::sim_model_ex1()
#' data_ind <- generate_indices(data)
#' clustInd_kmeans(data_ind, list(vars1, vars2))
#'
#' @export
clustInd_kmeans <- function(ind_data, vars_combinations,
                            dist_vector = c("euclidean", "mahalanobis"),
                            n_cluster = 2, init = "random", true_labels = NULL,
                            n_cores = 1) {
  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.", call. = FALSE)
  }

  if (!is.list(vars_combinations)) {
    stop("Input 'vars_combinations' must be a list.", call. = FALSE)
  }

  n_cores <- check_n_cores(n_cores)

  DIST_VECTOR <- c("euclidean", "mahalanobis")
  check_list_parameter(dist_vector, DIST_VECTOR, "dist")

  # Check for correct vars combinations
  vars_combinations_to_remove <- check_vars_combinations(vars_combinations, ind_data)

  if (length(vars_combinations_to_remove)) {
    vars_combinations <- vars_combinations[-vars_combinations_to_remove]
  }

  # Check if indices, methods and distances lists are provided
  if (length(vars_combinations) == 0 || length(dist_vector) == 0) {
    stop("Invalid 'vars_combinations' or 'dist_vector'. Both must be non-empty
         character vectors.", call. = FALSE)
  }

  if (is.null(names(vars_combinations))) {
    names(vars_combinations) <- paste0("vars", seq_along(vars_combinations))
  }

  # Generate all the possible combinations of indices, methods and distances
  parameter_combinations <- expand.grid(
    vars = names(vars_combinations),
    distance = dist_vector
  )

  tl_null <- is.null(true_labels)

  n_comb <- nrow(parameter_combinations)

  result <- parallel::mclapply(1:n_comb, function(i) {
    vars <- vars_combinations[[parameter_combinations$vars[i]]]
    dist <- parameter_combinations$distance[i]

    clustInd_kmeans_aux(
      ind_data = ind_data, vars = vars, dist = dist,
      n_cluster = n_cluster, init = init, true_labels = true_labels
    )
  }, mc.cores = n_cores)

  result_name <- get_result_names("kmeans", parameter_combinations, vars_combinations)
  names(result) <- result_name

  return(result)
}

#' Perform kernel kmeans clustering for a given combination of indices
#' and distance.
#'
#' @param ind_data Dataframe containing indices applied to the original data and
#' its first and second derivatives. See \link{generate_indices}.
#' @param vars vector with a combinations of indices in \code{ind_data}
#' @param kernel The kernel method to be used
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
#'
#' @return A list containing clustering results and execution time.
#'
#' @noRd
clustInd_kkmeans_aux <- function(ind_data, vars, kernel = "rbfdot",
                                 n_cluster = 2, true_labels = NULL, ...) {
  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.")
  }

  # Check if variables exist in the data frame
  if (!all(vars %in% names(ind_data))) {
    stop("Invalid variable names.")
  }

  # Check if the kernel given can be used
  if (!kernel %in% c("rbfdot", "polydot")) {
    stop("Invalid kernel.")
  }


  t0 <- Sys.time()
  kkmeans_out <- kernlab::kkmeans(as.matrix(ind_data[, vars]),
    centers = n_cluster, kernel = kernel, ...
  )
  clus <- kkmeans_out@.Data
  t1 <- Sys.time()
  t <- data.frame(difftime(t1, t0, "secs"))

  if (is.null(true_labels)) {
    res <- list("cluster" = clus, "time" = as.numeric(t))
  } else {
    valid <- clustering_validation(clus, true_labels)
    res <- list("cluster" = clus, "valid" = valid, "time" = as.numeric(t))
  }
  return(res)
}

#' Kernel k-means clustering using indices
#'
#' Perform kernel kmeans clustering for a different combinations of indices
#' and kernel
#'
#' @param ind_data Dataframe containing indices applied to the original data and
#' its first and second derivatives. See \link{generate_indices}.
#' @param vars_combinations \code{list} containing one or more combinations of indices in
#' \code{ind_data}. If it is non-named, the names of the variables are set to
#' vars1, ..., varsk, where k is the number of elements in \code{vars_combinations}.
#' @param kernel_list List of kernels
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
#' @param n_cores Number of cores to do parallel computation. 1 by default,
#' which mean no parallel execution.
#'
#' @return A \code{list} containing kernel-kmeans clustering results for each configuration.
#'
#' @examples
#' vars1 <- c("dtaEI", "dtaMEI")
#' vars2 <- c("dtaHI", "dtaMHI")
#' data <- ehymet::sim_model_ex1()
#' data_ind <- generate_indices(data)
#' clustInd_kkmeans(data_ind, list(vars1, vars2))
#'
#' @export
clustInd_kkmeans <- function(ind_data, vars_combinations,
                             kernel_list = c("rbfdot", "polydot"),
                             n_cluster = 2, true_labels = NULL,
                             n_cores = 1) {
  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.", call. = FALSE)
  }

  if (!is.list(vars_combinations)) {
    stop("Input 'vars_combinations' must be a list.", call. = FALSE)
  }

  n_cores <- check_n_cores(n_cores)

  # Check for correct vars combinations
  vars_combinations_to_remove <- check_vars_combinations(vars_combinations, ind_data)

  if (length(vars_combinations_to_remove)) {
    vars_combinations <- vars_combinations[-vars_combinations_to_remove]
  }

  # Check if indices, and kernel lists are provided
  if (!is.character(kernel_list) || length(vars_combinations) == 0 || length(kernel_list) == 0) {
    stop("Invalid 'kernel_list' or 'vars_combinations'. Both must be non-empty
         character vectors.", call. = FALSE)
  }

  if (is.null(names(vars_combinations))) {
    names(vars_combinations) <- paste0("vars", seq_along(vars_combinations))
  }

  # Generate all the possible combinations of indices, methods and distances
  parameter_combinations <- expand.grid(
    vars = names(vars_combinations),
    kernel = kernel_list
  )

  tl_null <- is.null(true_labels)

  n_comb <- nrow(parameter_combinations)

  result <- parallel::mclapply(1:n_comb, function(i) {
    vars <- vars_combinations[[parameter_combinations$vars[i]]]
    kern <- as.character(parameter_combinations$kernel[i])

    clustInd_kkmeans_aux(ind_data, vars, kern, n_cluster, true_labels)
  }, mc.cores = n_cores)

  result_name <- get_result_names("kkmeans", parameter_combinations, vars_combinations)
  names(result) <- result_name

  return(result)
}

#' Perform spectral clustering for a given combination of indices and kernel
#'
#' @param ind_data Dataframe containing indices applied to the original data and
#' its first and second derivatives. See \link{generate_indices}.
#' @param vars vector with a combinations of indices in \code{ind_data}
#' @param kernel The kernel method to be used
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
#'
#' @return A list containing clustering results and execution time.
#'
#' @noRd
clustInd_spc_aux <- function(ind_data, vars, kernel = "rbfdot", n_cluster = 2,
                             true_labels = NULL, ...) {
  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.", call. = FALSE)
  }

  # Check if variables exist in the data frame
  if (!all(vars %in% names(ind_data))) {
    stop("Invalid variable names.", call. = FALSE)
  }

  t0 <- Sys.time()

  # create target variable
  nrow_ind_data <- nrow(ind_data)
  specc_output <- kernlab::specc(as.matrix(ind_data[, vars]),
    centers = n_cluster, kernel = kernel, ...
  )
  clus <- specc_output@.Data

  t1 <- Sys.time()
  t <- data.frame(difftime(t1, t0, "secs"))

  if (is.null(true_labels)) {
    res <- list("cluster" = clus, "time" = as.numeric(t))
  } else {
    valid <- clustering_validation(clus, true_labels)
    res <- list("cluster" = clus, "valid" = valid, "time" = as.numeric(t))
  }
  return(res)
}

#' Spectral clustering using indices
#'
#' Perform spectral clustering for a different combinations of indices
#' and kernels
#'
#' @param ind_data Dataframe containing indices applied to the original data and
#' its first and second derivatives. See \link{generate_indices}.
#' @param vars_combinations \code{list} containing one or more combinations of indices in
#' \code{ind_data}. If it is non-named, the names of the variables are set to
#' vars1, ..., varsk, where k is the number of elements in \code{vars_combinations}.
#' @param kernel_list List of kernels
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
#' @param n_cores Number of cores to do parallel computation. 1 by default,
#' which mean no parallel execution.
#'
#' @return A list containing kkmeans clustering results for each configuration
#'
#' @examples
#' vars1 <- c("dtaEI", "dtaMEI")
#' vars2 <- c("dtaHI", "dtaMHI")
#' data <- ehymet::sim_model_ex1()
#' data_ind <- generate_indices(data)
#' clustInd_spc(data_ind, list(vars1, vars2))
#'
#' @export
clustInd_spc <- function(ind_data, vars_combinations,
                         kernel_list = c("rbfdot", "polydot"),
                         n_cluster = 2, true_labels = NULL,
                         n_cores = 1) {
  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.", call. = FALSE)
  }

  if (!is.list(vars_combinations)) {
    stop("Input 'vars_combinations' must be a data frame.", call. = FALSE)
  }

  n_cores <- check_n_cores(n_cores)

  # Check for correct vars combinations
  vars_combinations_to_remove <- check_vars_combinations(vars_combinations, ind_data)

  if (length(vars_combinations_to_remove)) {
    vars_combinations <- vars_combinations[-vars_combinations_to_remove]
  }

  # Check if indices, methods and distances lists are provided
  if (!is.character(kernel_list) || length(vars_combinations) == 0 ||
    length(kernel_list) == 0) {
    stop("Invalid 'kernel_list' or 'vars_combinations'. Both must be non-empty
         character vectors.", call. = FALSE)
  }

  if (is.null(names(vars_combinations))) {
    names(vars_combinations) <- paste0("vars", seq_along(vars_combinations))
  }

  # Generate all the possible combinations of indices, methods and distances
  parameter_combinations <- expand.grid(
    vars = names(vars_combinations),
    kernel = kernel_list
  )

  tl_null <- is.null(true_labels)

  n_comb <- nrow(parameter_combinations)

  result <- parallel::mclapply(1:n_comb, function(i) {
    vars <- vars_combinations[[parameter_combinations$vars[i]]]
    kern <- as.character(parameter_combinations$kernel[i])

    clustInd_spc_aux(ind_data, vars, kern, n_cluster, true_labels)
  }, mc.cores = n_cores)

  result_name <- get_result_names("spc", parameter_combinations, vars_combinations)
  names(result) <- result_name

  return(result)
}
