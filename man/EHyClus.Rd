% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{EHyClus}
\alias{EHyClus}
\title{EHyClus method for clustering. It creates a multivariate dataset containing
the epigraph, hypograph and/or modified version on the curves and derivatives
and perform hierarchical clustering, kmeans, kernel kmeans, support vector
clustering and spectral clustering}
\usage{
EHyClus(
  curves,
  t,
  vars_list,
  name_vars = NULL,
  nbasis = 30,
  norder = 4,
  indices = c("EI", "HI", "MEI", "MHI"),
  l_method_hierarch = c("single", "complete", "average", "centroid", "ward.D2"),
  l_dist_hierarch = c("euclidean", "manhattan"),
  l_dist_kmeans = c("euclidean", "mahalanobis"),
  l_kernel = c("rbfdot", "polydot"),
  l_method_svc = c("kmeans", "kernkmeans"),
  n_clusters = 2,
  true_labels = NULL,
  colapse = FALSE,
  num_cores = 1,
  ...
)
}
\arguments{
\item{curves}{Dataset containing the curves to apply a clustering algorithm.
The functional dataset can be one dimensional (nxp) where n is the number of
curves and p the number of time points, or multidimensional (nxpxk) where k
represents the number of dimensions in the data}

\item{t}{Grid}

\item{vars_list}{List containing one or more combinations of indexes in
\code{ind_data}}

\item{name_vars}{A vector with names for \code{vars_list}. NULL by default
in which case names are set to vars1, ..., varsk, where k is the number of
elements in \code{vars_list}.}

\item{nbasis}{Number of basis for the B-splines}

\item{norder}{Order of the B-splines}

\item{indices}{Names ofthe indices that need to be generated. They should be
one or more between EI, HI, MEI, MHI. Depending on the dimension on the data
they are calculated for one or multiple dimension}

\item{l_method_hierarch}{List of clustering methods for hierarchical
clustering}

\item{l_dist_hierarch}{List of distances for hierarchical clustering}

\item{l_dist_kmeans}{List of distances for kmeans clustering}

\item{l_kernel}{List of kernels}

\item{l_method_svc}{List of clustering methods for support vector clustering}

\item{n_clusters}{Number of clusters to create}

\item{true_labels}{Vector of true labels for validation
(if it is not known true_labels is set to NULL)}

\item{colapse}{It is a boolean. If it is true a dataframe with metrics values
is generated. If \code{true_labels} is True the dataframe contains Purity,
F-measure, RI and Time, and if it is False, only Time.}

\item{num_cores}{Number of cores to do parallel computation. 1 by default,
which mean no parallel execution.}

\item{...}{Additional arguments (unused)}
}
\value{
A list containing the clustering partition for each method and indexes
combination and a data frame containing the time elapsed for obtaining a
clustering partition of the indexes dataset for each methodology
}
\description{
EHyClus method for clustering. It creates a multivariate dataset containing
the epigraph, hypograph and/or modified version on the curves and derivatives
and perform hierarchical clustering, kmeans, kernel kmeans, support vector
clustering and spectral clustering
}
\examples{
vars1 <- c("dtaEI", "dtaMEI"); vars2 <- c("dtaHI", "dtaMHI")
varsl <- list(vars1, vars2)
data <- ehymet::sim_model_ex1()
t <- seq(0, 1, length = 30)
EHyClus(data, t, varsl)
}
