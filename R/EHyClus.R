#' EHyClus method for clustering. It creates a multivariate dataset containing
#' the epigraph, hypograph and/or modified version on the curves and derivatives
#' and perform hierarchical clustering, kmeans, kernel kmeans, support vector
#' clustering and spectral clustering
#'
#' @param curves Dataset containing the curves to apply a clustering algorithm.
#' The functional dataset can be one dimensional (nxp) where n is the number of
#' curves and p the number of time points, or multidimensional (nxpxk) where k
#' represents the number of dimensions in the data
#' @param t Grid
#' @param vars_list List containing one or more combinations of indexes in
#' \code{ind_data}
#' @param name_vars A vector with names for \code{vars_list}. NULL by default
#' in which case names are set to vars1, ..., varsk, where k is the number of
#' elements in \code{vars_list}.
#' @param nbasis Number of basis for the B-splines
#' @param norder Order of the B-splines
#' @param indices Names ofthe indices that need to be generated. They should be
#' one or more between EI, HI, MEI, MHI. Depending on the dimension on the data
#' they are calculated for one or multiple dimension
#' @param l_method_hierarch List of clustering methods for hierarchical
#' clustering
#' @param l_dist_hierarch List of distances for hierarchical clustering
#' @param l_dist_kmeans List of distances for kmeans clustering
#' @param l_kernel List of kernels
#' @param l_method_svc List of clustering methods for support vector clustering
#' @param n_clusters Number of clusters to create
#' @param true_labels Vector of true labels for validation
#' (if it is not known true_labels is set to NULL)
#' @param colapse It is a boolean. If it is true a dataframe with metrics values
#' is generated. If \code{true_labels} is True the dataframe contains Purity,
#' F-measure, RI and Time, and if it is False, only Time.
#' @param num_cores Number of cores to do parallel computation. 1 by default,
#' which mean no parallel execution.
#' @param ... Additional arguments (unused)
#'
#' @return A list containing the clustering partition for each method and indexes
#' combination and a data frame containing the time elapsed for obtaining a
#' clustering partition of the indexes dataset for each methodology
#' @export
#' @examples
#' vars1 <- c("dtaEI", "dtaMEI"); vars2 <- c("dtaHI", "dtaMHI")
#' varsl <- list(vars1, vars2)
#' data <- ehymet::sim_model_ex1()
#' t <- seq(0, 1, length = 30)
#' EHyClus(data, t, varsl)
EHyClus <- function(curves, t, vars_list, name_vars = NULL, nbasis = 30,
                    norder = 4, indices = c("EI", "HI", "MEI", "MHI"),
                    l_method_hierarch = c("single","complete","average",
                                          "centroid","ward.D2"),
                    l_dist_hierarch = c("euclidean", "manhattan"),
                    l_dist_kmeans =  c("euclidean", "mahalanobis"),
                    l_kernel = c("rbfdot", "polydot"),
                    l_method_svc = c("kmeans", "kernkmeans"),
                    n_clusters = 2, true_labels = NULL, colapse = FALSE,
                    num_cores = 1, ...){

  # Constants definition
  INDICES         <- c("EI", "HI", "MEI", "MHI")
  METHOD_HIERARCH <- c("single", "complete", "average", "centroid", "ward.D2")
  DIST_HIERARCH   <- c("euclidean", "manhattan")
  DIST_KMEANS     <- c("euclidean", "mahalanobis")
  KERNEL          <- c("rbfdot", "polydot")
  METHOD_SVC      <- c("kmeans", "kernkmeans")

  if(!is.list(vars_list)) {
    stop("Input 'vars_list' must be a list.", call. = FALSE)
  }

  check_ehyclus_parameter(indices, INDICES, "indices")
  check_ehyclus_parameter(l_method_hierarch, METHOD_HIERARCH, "l_method_hierarch")
  check_ehyclus_parameter(l_dist_hierarch, DIST_HIERARCH, "l_dist_hierarch")
  check_ehyclus_parameter(l_dist_kmeans, DIST_KMEANS, "l_dist_kmeans")
  check_ehyclus_parameter(l_kernel, KERNEL, "l_kernel")
  check_ehyclus_parameter(l_method_svc, METHOD_SVC, "l_method_svc")

  # Check indices names
  if(!all(indices %in% c("EI", "HI", "MEI", "MHI"))){
    stop("Indices should be one or more of the following: EI, HI, MEI, MHI")
  }

  # vars_list TIENE QUE SER LIST !!!!!

  # Generate the dataset with the indexes
  ind_curves <- ind(curves, t, nbasis, norder, indices)

  # Hierarchical clustering
  cl_hierarch <- clustInd_hierarch(ind_data = ind_curves, vars_list = vars_list,
                                   name_vars = name_vars,
                                   method_list = l_method_hierarch,
                                   dist_list = l_dist_hierarch,
                                   n_cluster = n_clusters, true_labels = true_labels,
                                   colapse = colapse, num_cores = num_cores)

  # kmeans
  cl_kmeans <- clustInd_kmeans(ind_data = ind_curves, vars_list = vars_list,
                               dist_list = l_dist_kmeans, n_cluster = n_clusters,
                               true_labels = true_labels, colapse = colapse,
                               num_cores = num_cores)

  # kernel kmeans
  cl_kkmeans <- clustInd_kkmeans(ind_data = ind_curves, vars_list = vars_list,
                                 kernel_list = l_kernel, n_cluster = n_clusters,
                                 true_labels = true_labels, colapse = colapse,
                                 num_cores = num_cores)

  # support vector clustering
  cl_svc <- clustInd_svc(ind_data = ind_curves, vars_list = vars_list,
                         method_list = l_method_svc, n_cluster = n_clusters,
                         true_labels = true_labels, colapse = colapse,
                         num_cores = num_cores)

  # spectral clustering
  cl_spc <- clustInd_spc(ind_data = ind_curves, vars_list = vars_list,
                         kernel_list = l_kernel, n_cluster = n_clusters,
                         true_labels = true_labels, colapse = colapse,
                         num_cores = num_cores)

  cluster <- c(cl_hierarch, cl_kmeans, cl_kkmeans, cl_svc, cl_spc)

  if (colapse) {
    metrics <- rbind(cl_hierarch$metrics, cl_kmeans$metrics, cl_kkmeans$metrics,
                     cl_svc$metrics, cl_spc$metrics)
    result <- list("cluster" = cluster, "metrics" = metrics)

  } else {
    result <- list("cluster" = cluster)
  }


  return(result)
}


#' @noRd
check_ehyclus_parameter <- function(argument, parameter_values, parameter_name) {
  if (length(argument) == 0) {
    stop("Parameter '", parameter_name, "' should have at least one element.", call. = FALSE)
  }

  if (any(duplicated(argument))) {
    stop("Duplicated argument in '", parameter_name,"'.", call. = FALSE)
  }

  indices <- pmatch(argument, parameter_values)
  if (any(is.na(indices))) {
    stop("Invalid argument in '", parameter_name, "': ", paste(argument[is.na(indices)], collapse = ", "), ".",
         call. = FALSE)
  }
}
