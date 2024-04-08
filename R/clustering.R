#' Perform hierarchical clustering for a given combination of indexes,
#' method and distance
#'
#' @param ind_data Dataframe containing indexes applied to the original data and
#' its first and second derivatives
#' @param vars vector with a combinations of indexes in \code{ind_data}
#' @param method The agglomerative method to be used in hierarchical clustering
#' @param dist The distance method to be used
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
#'
#' @return A list containing clustering results and execution time.
#' @noRd
#'
clustInd_hierarch_aux <- function(ind_data, vars, method = "single",
                                  dist = "euclidean", n_cluster = 2,
                                  true_labels=NULL){

  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.")
  }

  # Check if variables exist in the data frame
  if (!all(vars %in% names(ind_data))) {
    stop("Invalid variable name.")
  }

  if(!all(method %in% c("single","complete","average", "centroid","ward.D2"))) {
    stop("Invalid method name.")
  }

  if(!all(dist %in% c("euclidean", "manhattan"))) {
    stop("Invalid distance name.")
  }

  t0 <- Sys.time()

  # Perform hierarchical clustering
  d <- stats::dist(ind_data[,vars], method = dist)
  met <- stats::hclust(d, method = method)
  clus <- stats::cutree(met, k = n_cluster)

  t1 <- Sys.time()

  # Calculate execution time
  et <- data.frame(difftime(t1,t0,'secs'))

  if (is.null(true_labels)) {
    res <- list("cluster" = clus, "time" = as.numeric(et))
  } else {
    valid <- valid(true_labels, clus)
    res <- list("cluster" = clus, "valid" = valid, "time" = as.numeric(et))
  }

  return(res)
}

#' Perform hierarchical clustering for a different combinations of indexes,
#' method and distance
#'
#' @param ind_data Dataframe containing indexes applied to the original data and
#' its first and second derivatives
#' @param vars_list List containing one or more combinations of indexes in
#' \code{ind_data}
#' @param name_vars A vector with names for \code{vars_list}. NULL by default
#' in which case names are set to vars1, ..., varsk, where k is the number of
#' elements in \code{vars_list}.
#' @param method_list List of clustering methods
#' @param dist_list List of distance metrics
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
#' @param colapse It is a boolean. If it is true a dataframe with metrics values
#' is generated. If \code{true_labels} is True the dataframe contains Purity,
#' F-measure, RI and Time, and if it is False, only Time.
#' @param num_cores Number of cores to do parallel computation. 1 by default,
#' which mean no parallel execution.
#' @return A list containing hierarchical clustering results
#' for each configuration
#' @export
#' @examples
#' vars1 <- c("dtaEI", "dtaMEI")
#' vars2 <- c("dtaHI", "dtaMHI")
#' data <- ehymet::sim_model_ex1()
#' data_ind <- ehymet::ind(data, t=seq(0, 1, length = 30))
#' clustInd_hierarch(data_ind, list(vars1, vars2))
clustInd_hierarch <- function(ind_data, vars_list, name_vars = NULL,
                              method_list = c("single","complete","average",
                                              "centroid","ward.D2"),
                              dist_list = c("euclidean", "manhattan"),
                              n_cluster=2, true_labels = NULL, colapse = FALSE,
                              num_cores=1) {

  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("input 'ind_data' must be a data frame.", call. = FALSE)
  }

  if(!is.list(vars_list)) {
    stop("input 'vars_list' must be a list.", call. = FALSE)
  }

  if(!is.null(name_vars) && !length(vars_list) == length(name_vars)){
    stop("'name_vars' and 'vars_list' should have the same length.", call. = FALSE)
  }

  # Check if indices, methods and distances lists are provided
  if (!is.character(method_list) ||
      !is.character(dist_list) || length(vars_list) == 0 ||
      length(method_list) == 0 || length(dist_list) == 0) {
    stop("invalid 'method_list' or 'dist_list'. Both must be non-empty character vectors.", call. = FALSE)
  }

  if(is.null(name_vars)) name_vars <- paste0("vars", seq_along(vars_list))
  names(vars_list) <- name_vars

  # Generate all the possible combinations of indices, methods and distances
  parameter_combinations <- expand.grid(vars = names(vars_list),
                                        method = method_list, distance = dist_list)

  tl_null <- is.null(true_labels)

  n_comb <- nrow(parameter_combinations)

  result <- parallel::mclapply(1:n_comb, function(i) {
    vars <- vars_list[[parameter_combinations$vars[i]]]
    met  <- parameter_combinations$method[i]
    dist <- parameter_combinations$distance[i]

    clustInd_hierarch_aux(ind_data, vars, met, dist, n_cluster, true_labels)
}, mc.cores = num_cores)

  # result_name <- apply(parameter_combinations, 1, function(row) {
  #   paste("hierarch", row["method"], row["distance"],
  #         paste(get(as.character(row["vars"])), collapse = ""),
  #         sep = "_")
  # })
  #
  # names(result) <- result_name

  if (colapse) {
    result <- list("list" = result,
                   "metrics" = result_to_table(result, colapse))
  }


  return(result)
}


#' k-means clustering with Mahalanobis distance
#'
#' @param ind_data Dataframe containing indexes applied to the original data and
#' its first and second derivatives
#' @param n_cluster Number of clusters to create
#'
#' @return k-means clustering with Mahalanobis distance output
#' @noRd
kmeans_mahal <- function(ind_data, n_cluster){

  # Check if input is numeric matrix or array
  if (!(is.numeric(ind_data) || is.matrix(ind_data) || is.array(ind_data) ||
        is.data.frame(ind_data))) {
    stop("Input must be a numeric matrix, array or data frame.")
  }

  # Convert data to matrix
  if (! is.matrix(ind_data)) {
    data_matrix <- as.matrix(ind_data)
  } else {
    data_matrix <- ind_data
  }

  # Cholesky transformation
  c <- chol(stats::var(data_matrix))
  data_transform <- data_matrix %*% solve(c)

  #vector containing the clustering partition
  km <- stats::kmeans(data_transform, centers=n_cluster, iter.max=1000,
                      nstart=100)$cluster
  return(km)
}

#' Perform kmeans clustering for a given combination of indexes and distance
#'
#' @param ind_data Dataframe containing indexes applied to the original data and
#' its first and second derivatives
#' @param vars vector with a combinations of indexes in \code{ind_data}
#' @param dist The distance method to be used
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
#'
#' @return A list containing clustering results and execution time.
#' @noRd
clustInd_kmeans_aux <- function(ind_data, vars, dist = "euclidean",
                                n_cluster = 2, true_labels = NULL){

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

  t0 <- Sys.time()

  if(dist=="euclidean"){
    clus <- stats::kmeans(ind_data[,vars], centers = n_cluster,
                          iter.max = 1000, nstart = 100)$cluster
  } else{
    clus <- kmeans_mahal(ind_data[,vars], n_cluster)
  }
  t1 <- Sys.time()
  t <- data.frame(difftime(t1,t0,'secs'))

  if (is.null(true_labels)) {
    res <- list("cluster" = clus, "time" = as.numeric(t))
  } else {
    valid <- valid(true_labels, clus)
    res <- list("cluster" = clus, "valid" = valid, "time" = as.numeric(t))
  }
  return(res)
}

#' Perform hierarchical clustering for a different combinations of indexes,
#' method and distance
#'
#' @param ind_data Dataframe containing indexes applied to the original data and
#' its first and second derivatives
#' @param vars_list List containing one or more combinations of indexes in
#' \code{ind_data}
#' @param name_vars A vector with names for \code{vars_list}. NULL by default
#' in which case names are set to vars1, ..., varsk, where k is the number of
#' elements in \code{vars_list}.
#' @param dist_list List of distance metrics
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
#' @param colapse It is a boolean. If it is true a dataframe with metrics values
#' is generated. If \code{true_labels} is True the dataframe contains Purity,
#' F-measure, RI and Time, and if it is False, only Time.
#' @param num_cores Number of cores to do parallel computation. 1 by default,
#' which mean no parallel execution.
#' @return A list containing hierarchical clustering results
#' for each configuration
#'
#' @return A list containing kmeans clustering results for each configuration
#' @export
#' @examples
#' vars1 <- c("dtaEI", "dtaMEI")
#' vars2 <- c("dtaHI", "dtaMHI")
#' data <- ehymet::sim_model_ex1()
#' data_ind <- ehymet::ind(data, t=seq(0, 1, length = 30))
#' clustInd_kmeans(data_ind, list(vars1, vars2))
clustInd_kmeans <- function(ind_data, vars_list, name_vars = NULL,
                            dist_list = c("euclidean", "mahalanobis"),
                            n_cluster = 2, true_labels = NULL, colapse = FALSE,
                            num_cores = 1) {

  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.", call. = FALSE)
  }

  if(!is.list(vars_list)) {
    stop("Input 'vars_list' must be a list.", call. = FALSE)
  }

  if(!is.null(name_vars) && !length(vars_list) == length(name_vars)){
    stop("'name_vars' and 'vars_list' should have the same length.", call. = FALSE)
  }

  # Check if indices, methods and distances lists are provided
  if (length(vars_list) == 0 || length(dist_list) == 0) {
    stop("Invalid 'vars_list' or 'dist_list'. Both must be non-empty
         character vectors.", call. = FALSE)
  }

  if(is.null(name_vars)) name_vars <- paste0("vars", seq_along(vars_list))

  names(vars_list) <- name_vars

  # Generate all the possible combinations of indices, methods and distances
  parameter_combinations <- expand.grid(vars = names(vars_list),
                                        distance = dist_list)

  tl_null <- is.null(true_labels)

  n_comb <- nrow(parameter_combinations)

  result <- parallel::mclapply(1:n_comb, function(i) {
    vars <- vars_list[[parameter_combinations$vars[i]]]
    dist <- parameter_combinations$distance[i]

    clustInd_kmeans_aux(ind_data = ind_data, vars =vars, dist = dist,
                           n_cluster = n_cluster, true_labels = true_labels)
  }, mc.cores = num_cores)

  # result_name <- apply(parameter_combinations, 1, function(row) {
  #   paste("kmeans", row["distance"], paste(get(as.character(row["vars"])),
  #                                          collapse = ""), sep = "_")
  # })
  #
  # names(result) <- result_name

  if (colapse) {
    result <- list("list" = result,
                   "metrics" = result_to_table(result, colapse))
  }

  return(result)
}

#' Perform kernel kmeans clustering for a given combination of indexes
#' and distance
#'
#' @param ind_data Dataframe containing indexes applied to the original data and
#' its first and second derivatives
#' @param vars vector with a combinations of indexes in \code{ind_data}
#' @param kernel The kernel method to be used
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
#'
#' @return A list containing clustering results and execution time.
#' @noRd
clustInd_kkmeans_aux <- function(ind_data, vars, kernel = "rbfdot",
                                 n_cluster = 2, true_labels = NULL, ...){

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
  kkmeans_out <- kernlab::kkmeans(as.matrix(ind_data[,vars]),
                                  centers= n_cluster, kernel = kernel, ...)
  clus <- kkmeans_out@.Data
  t1 <- Sys.time()
  t <- data.frame(difftime(t1,t0,'secs'))

  if (is.null(true_labels)) {
    res <- list("cluster" = clus, "time" = as.numeric(t))
  } else {
    valid <- valid(true_labels, clus)
    res <- list("cluster" = clus, "valid" = valid, "time" = as.numeric(t))
  }
  return(res)
}

#' Perform kernel kmeans clustering for a different combinations of indexes
#' and kernel
#'
#' @param ind_data Dataframe containing indexes applied to the original data and
#' its first and second derivatives
#' @param vars_list List containing one or more combinations of indexes in
#' \code{ind_data}
#' @param name_vars A vector with names for \code{vars_list}. NULL by default
#' in which case names are set to vars1, ..., varsk, where k is the number of
#' elements in \code{vars_list}.
#' @param kernel_list List of kernels
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
#' @param colapse It is a boolean. If it is true a dataframe with metrics values
#' is generated. If \code{true_labels} is True the dataframe contains Purity,
#' F-measure, RI and Time, and if it is False, only Time.
#' @param num_cores Number of cores to do parallel computation. 1 by default,
#' which mean no parallel execution.
#' @param ... Additional arguments (unused)
#'
#' @return A list containing kkmeans clustering results for each configuration
#' @export
#' @examples
#' vars1 <- c("dtaEI", "dtaMEI")
#' vars2 <- c("dtaHI", "dtaMHI")
#' data <- ehymet::sim_model_ex1()
#' data_ind <- ehymet::ind(data, t=seq(0, 1, length = 30))
#' clustInd_kkmeans(data_ind, list(vars1, vars2))
clustInd_kkmeans <- function(ind_data, vars_list, name_vars = NULL,
                             kernel_list = c("rbfdot", "polydot"),
                             n_cluster = 2, true_labels = NULL, colapse = FALSE,
                             num_cores = 1, ...) {

  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.", call. = FALSE)
  }

  if(!is.list(vars_list)) {
    stop("Input 'vars_list' must be a list.", call. = FALSE)
  }

  if(!is.null(name_vars) && !(length(vars_list) == length(name_vars))){
    stop("'name_vars' and 'vars_list' should have the same length.", call. = FALSE)
  }

  # Check if indices, and kernel lists are provided
  if (!is.character(kernel_list) || length(vars_list) == 0 || length(kernel_list) == 0) {
    stop("Invalid 'kernel_list' or 'vars_list'. Both must be non-empty
         character vectors.", call. = FALSE)
  }

  if(is.null(name_vars)) name_vars <- paste0("vars", seq_along(vars_list))
  names(vars_list) <- name_vars

  # Generate all the possible combinations of indices, methods and distances
  parameter_combinations <- expand.grid(vars = names(vars_list),
                                        kernel = kernel_list)

  tl_null <- is.null(true_labels)

  n_comb <- nrow(parameter_combinations)

  result <- parallel::mclapply(1:n_comb, function(i) {
    vars <- vars_list[[parameter_combinations$vars[i]]]
    kern <- as.character(parameter_combinations$kernel[i])

    clustInd_kkmeans_aux(ind_data, vars, kern, n_cluster, true_labels)
  }, mc.cores = num_cores)

  # result_name <- apply(parameter_combinations, 1, function(row) {
  #   paste("kkmeans", row["kernel"],
  #         paste(get(as.character(row["vars"])), collapse = ""),
  #         sep = "_")
  # })
  #
  # names(result) <- result_name

  if (colapse) {
    result <- list("list" = result,
                   "metrics" = result_to_table(result, colapse))
  }

  return(result)
}

#' Perform support vector clustering for a given combination of indexes
#' and kernel
#'
#' @param ind_data Dataframe containing indexes applied to the original data and
#' its first and second derivatives
#' @param vars vector with a combinations of indexes in \code{ind_data}
#' @param method The initialization method to be used
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
#'
#' @return A list containing clustering results and execution time.
#' @noRd
clustInd_svc_aux <- function(ind_data, vars, method = "kmeans", n_cluster = 2,
                             true_labels = NULL, ...){

  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.")
  }

  # Check if variables exist in the data frame
  if (!all(vars %in% names(ind_data))) {
    stop("Invalid variable names.")
  }

  t0 <- Sys.time()

  # create target variable
  nrow_ind_data <- nrow(ind_data)
  r1 <- round(nrow_ind_data/n_cluster)
  r2 <- nrow_ind_data - (n_cluster-1) * r1
  y <-  c(rep(1:(n_cluster-1), each =r1), rep(n_cluster, r2))
  clusterSVM_output <- SwarmSVM::clusterSVM(ind_data[,vars], y, n_cluster,
                                            cluster.method=method, ...)
  clus <- clusterSVM_output$label

  t1 <- Sys.time()
  t <- data.frame(difftime(t1,t0,'secs'))

  if (is.null(true_labels)) {
    res <- list("cluster" = clus, "time" = as.numeric(t))
  } else {
    valid <- valid(true_labels, clus)
    res <- list("cluster" = clus, "valid" = valid, "time" = as.numeric(t))
  }
  return(res)
}

#' Perform support vector clustering for a different combinations of indexes
#' and distance
#'
#' @param ind_data Dataframe containing indexes applied to the original data and
#' its first and second derivatives
#' @param vars_list List containing one or more combinations of indexes in
#' \code{ind_data}
#' @param name_vars A vector with names for \code{vars_list}. NULL by default
#' in which case names are set to vars1, ..., varsk, where k is the number of
#' elements in \code{vars_list}.
#' @param method_list List of methods
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
#' @param colapse It is a boolean. If it is true a dataframe with metrics values
#' is generated. If \code{true_labels} is True the dataframe contains Purity,
#' F-measure, RI and Time, and if it is False, only Time.
#' @param num_cores Number of cores to do parallel computation. 1 by default,
#' which mean no parallel execution.
#' @param ... Additional arguments (unused)
#'
#' @return A list containing kkmeans clustering results for each configuration
#' @export
#' @examples
#' vars1 <- c("dtaEI", "dtaMEI")
#' vars2 <- c("dtaHI", "dtaMHI")
#' data <- ehymet::sim_model_ex1()
#' data_ind <- ehymet::ind(data, t=seq(0, 1, length = 30))
#' clustInd_svc(data_ind, list(vars1, vars2))
clustInd_svc <- function(ind_data, vars_list, name_vars = NULL,
                         method_list = c("kmeans", "kernkmeans"),
                         n_cluster = 2, true_labels = NULL, colapse = FALSE,
                         num_cores = 1, ...) {

  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.", call. = FALSE)
  }

  if (!is.list(vars_list)) {
    stop("Input 'vars_list' must be a data frame.", call. = FALSE)
  }

  if (!is.null(name_vars) & !length(vars_list) == length(name_vars)){
    stop("'name_vars' and 'vars_list' should have the same length.", call. = FALSE)
  }

  # Check if indices, methods and distances lists are provided
  if ( !is.character(method_list) ||  length(vars_list) == 0 ||
       length(method_list) == 0) {
    stop("Invalid 'method_list' or 'vars_list'. Both must be non-empty
         character vectors.", call. = FALSE)
  }

  if(is.null(name_vars)) name_vars <- paste0("vars", seq_along(vars_list))
  names(vars_list) <- name_vars

  # Generate all the possible combinations of indices, methods and distances
  parameter_combinations <- expand.grid(vars = names(vars_list),
                                        method = method_list)

  tl_null <- is.null(true_labels)

  n_comb <- nrow(parameter_combinations)

  result <- parallel::mclapply(1:n_comb, function(i) {
    vars <- vars_list[[parameter_combinations$vars[i]]]
    met <- as.character(parameter_combinations$method[i])

    clustInd_svc_aux(ind_data, vars, met, n_cluster, true_labels)
  }, mc.cores = num_cores)

  # result_name <- apply(parameter_combinations, 1, function(row) {
  #   paste("svc", row["method"],
  #         paste(get(as.character(row["vars"])), collapse = ""),
  #         sep = "_")
  # })
  #
  # names(result) <- result_name

  if(colapse) result <- list("list" = result,
                             "metrics" = result_to_table(result, colapse))

  return(result)
}

#' Perform spectral clustering for a given combination of indexes and kernel
#'
#' @param ind_data Dataframe containing indexes applied to the original data and
#' its first and second derivatives
#' @param vars vector with a combinations of indexes in \code{ind_data}
#' @param kernel The kernel method to be used
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
#'
#' @return A list containing clustering results and execution time.
#' @noRd
clustInd_spc_aux <- function(ind_data, vars, kernel = "rbfdot", n_cluster = 2,
                             true_labels=NULL, ...){

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
  specc_output <- kernlab::specc(as.matrix(ind_data[,vars]),
                                 centers = n_cluster,  kernel = kernel, ...)
  clus <- specc_output@.Data

  t1 <- Sys.time()
  t <- data.frame(difftime(t1,t0,'secs'))

  if (is.null(true_labels)) {
    res <- list("cluster" = clus, "time" = as.numeric(t))
  } else {
    valid <- valid(true_labels, clus)
    res <- list("cluster" = clus, "valid" = valid, "time" = as.numeric(t))
  }
  return(res)
}

#' Perform spectral clustering for a different combinations of indexes
#' and kernels
#'
#' @param ind_data Dataframe containing indexes applied to the original data and
#' its first and second derivatives
#' @param vars_list List containing one or more combinations of indexes in
#' \code{ind_data}
#' @param name_vars A vector with names for \code{vars_list}. NULL by default
#' in which case names are set to vars1, ..., varsk, where k is the number of
#' elements in \code{vars_list}.
#' @param kernel_list List of kernels
#' @param n_cluster Number of clusters to create
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
#' @param colapse It is a boolean. If it is true a dataframe with metrics values
#' is generated. If \code{true_labels} is True the dataframe contains Purity,
#' F-measure, RI and Time, and if it is False, only Time.
#' @param num_cores Number of cores to do parallel computation. 1 by default,
#' which mean no parallel execution.
#' @param ... Additional arguments (unused)
#'
#' @return A list containing kkmeans clustering results for each configuration
#' @export
#' @examples
#' vars1 <- c("dtaEI", "dtaMEI")
#' vars2 <- c("dtaHI", "dtaMHI")
#' data <- ehymet::sim_model_ex1()
#' data_ind <- ehymet::ind(data, t=seq(0, 1, length = 30))
#' clustInd_spc(data_ind, list(vars1, vars2))
clustInd_spc <- function(ind_data, vars_list, name_vars = NULL,
                         kernel_list = c("rbfdot", "polydot"),
                         n_cluster = 2, true_labels = NULL, colapse = FALSE,
                         num_cores = 1, ...) {

  # Check if input is a data frame
  if (!is.data.frame(ind_data)) {
    stop("Input 'ind_data' must be a data frame.", call. = FALSE)
  }

  if(!is.list(vars_list)) {
    stop("Input 'vars_list' must be a data frame.", call. = FALSE)
  }

  if(!is.null(name_vars) & !length(vars_list) == length(name_vars)){
    stop("'name_vars' and 'vars_list' should have the same length.", call. = FALSE)
  }

  # Check if indices, methods and distances lists are provided
  if (!is.character(kernel_list) || length(vars_list) == 0 ||
      length(kernel_list) == 0) {
    stop("Invalid 'kernel_list' or 'vars_list'. Both must be non-empty
         character vectors.", call. = FALSE)
  }

  if(is.null(name_vars)) name_vars <- paste0("vars", 1:length(vars_list))
  names(vars_list) <- name_vars

  # Generate all the possible combinations of indices, methods and distances
  parameter_combinations <- expand.grid(vars = names(vars_list),
                                        kernel = kernel_list)

  tl_null <- is.null(true_labels)

  n_comb <- nrow(parameter_combinations)

  result <- parallel::mclapply(1:n_comb, function(i) {
    vars <- vars_list[[parameter_combinations$vars[i]]]
    kern <- as.character(parameter_combinations$kernel[i])

    clustInd_spc_aux(ind_data, vars, kern, n_cluster, true_labels)
  }, mc.cores = num_cores)

  # result_name <- apply(parameter_combinations, 1, function(row) {
  #   paste("spc", row["kernel"],
  #         paste(get(as.character(row["vars"])), collapse = ""),
  #         sep = "_")
  # })
  #
  # names(result) <- result_name

  if(colapse) result <- list("list" = result,
                             "metrics" = result_to_table(result, colapse))

  return(result)
}
