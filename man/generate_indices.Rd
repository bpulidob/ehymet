% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indices.R
\name{generate_indices}
\alias{generate_indices}
\title{Create a dataset with indexes from a functional dataset in one or multiple
dimensions}
\usage{
generate_indices(
  curves,
  k,
  bs = "cr",
  grid_ll = 0,
  grid_ul = 1,
  indices = c("EI", "HI", "MEI", "MHI"),
  ...
)
}
\arguments{
\item{curves}{\code{matrix} with dimension \eqn{n \times p} in the case of a
one-dimensional functional dataset, or \code{array} of dimension
\eqn{n \times p \times k} in the case of a multivariate functional dataset.
\eqn{n} represents the number of curves, \eqn{p} the number of values along
the curve, and in the second case, \eqn{k} is the number of dimensions.}

\item{k}{Number of basis functions for the B-splines. If not provided, it will
be automatically set.}

\item{bs}{A two letter chatacter string indicating the (penalized) smoothing
basis to use. See \code{\link{smooth.terms}}.}

\item{grid_ll}{Lower limit of the grid.}

\item{grid_ul}{Upper limit of the grid.}

\item{indices}{Set of indices to be applied to the dataset. They should be
any between EI, HI, MEI and MHI}

\item{...}{Additional arguments (unused)}
}
\value{
A dataframe containing the indexes provided in \code{indices} for
original data, first and second derivatives
}
\description{
Create a dataset with indexes from a functional dataset in one or multiple
dimensions
}
\examples{
x1 <- array(c(1,2,3, 3,2,1, 5,2,3, 9,8,7, -1,-5,-6, 2,3,0, -1,0,2, -1,-2,0),
dim = c(3,4,2))
generate_indices(x1, k = 4)

x2 <- matrix(c(1,2,3,3,2,1,5,2,3,9,8,7), nrow = 3, ncol  = 4)
generate_indices(x2, k = 4)

}