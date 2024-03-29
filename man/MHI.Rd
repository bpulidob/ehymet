% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indices.R
\name{MHI}
\alias{MHI}
\title{Modified Hypograph Index (MHI) for a univariate functional dataset.}
\usage{
MHI(curves)
}
\arguments{
\item{curves}{A matrix where each row represents a curve, and each column
represents values along the curve.}
}
\value{
Return a numeric vector containing the MHI for each curve
}
\description{
The Modified Hypograph Index of a curve x is the proportion of "time"
the curves in the sample are below x.
}
\examples{
x <- matrix(c(1,2,3,3,2,1,5,2,3,9,8,7), ncol = 3, nrow = 4)
MHI(x)
}
