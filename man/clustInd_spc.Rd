% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{clustInd_spc}
\alias{clustInd_spc}
\title{Perform spectral clustering for a different combinations of indexes
and kernels}
\usage{
clustInd_spc(
  ind_data,
  vars_list,
  name_vars = NULL,
  kernel_list = c("rbfdot", "polydot"),
  n_cluster = 2,
  true_labels = NULL,
  colapse = FALSE,
  num_cores = 1,
  ...
)
}
\arguments{
\item{ind_data}{Dataframe containing indexes applied to the original data and
its first and second derivatives}

\item{vars_list}{List containing one or more combinations of indexes in
\code{ind_data}}

\item{name_vars}{A vector with names for \code{vars_list}. NULL by default
in which case names are set to vars1, ..., varsk, where k is the number of
elements in \code{vars_list}.}

\item{kernel_list}{List of kernels}

\item{n_cluster}{Number of clusters to create}

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
A list containing kkmeans clustering results for each configuration
}
\description{
Perform spectral clustering for a different combinations of indexes
and kernels
}
\examples{
vars1 <- c("dtaEI", "dtaMEI")
vars2 <- c("dtaHI", "dtaMHI")
data <- ehymet::sim_model_ex1()
data_ind <- ehymet::ind(data, t=seq(0, 1, length = 30))
clustInd_spc(data_ind, list(vars1, vars2))
}
