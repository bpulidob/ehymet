#' Epigraph Index (EI) for a univariate functional dataset.
#'
#' The Epigraph Index of a curve x is one minus the proportion of curves
#' in the sample that are above x.
#'
#' @param curves A \code{matrix} where each row represents a curve, and each column
#' represents values along the curve or an \code{array} with dimension \eqn{n \times p \times k} with
#' \eqn{n} curves, \eqn{p} values along the curve, and \eqn{k} dimensions.
#'
#' @return Return a numeric vector containing the EI for each curve.
#'
#' @examples
#' x <- matrix(c(1,2,3,3,2,1,5,2,3,9,8,7), ncol = 3, nrow = 4)
#' EI(x)
#'
#' y <- array(c(1,2,3, 3,2,1, 5,2,3, 9,8,7, -1,-5,-6, 2,3,0, -1,0,2, -1,-2,0),
#' dim = c(3,4,2))
#' EI(y)
#'
#' @export
EI <- function(curves) {
  UseMethod("EI")
}

#' @export
EI.matrix <- function(curves){
  n_curves <- dim(curves)[1]
  l_curves <- dim(curves)[2]

  index <- apply(curves,1, function(y)
    sum(apply(curves,1,function(x)
      sum(x>=y)==l_curves)))/n_curves

  return(1 - index)
}

#' @export
EI.default <- function(curves) {
  if (length(dim(curves)) != 3 || is.null(dim(curves))) {
    stop("'curves' should be a matrix or a 3-dimensional array", call. = FALSE)
  }

  n_curves <- dim(curves)[1]
  l_curves <- dim(curves)[2]
  d_curves <- dim(curves)[3]

  index <- colSums(Reduce('*',lapply(1:d_curves, function(k)
    sapply(1:n_curves, function(j)
      colSums(sapply(1:n_curves, function(i)
        curves[i,,k] >= curves[j,,k]))==l_curves))))

  return(1 - index/n_curves)
}

#' Hypograph Index (HI) for a univariate functional dataset.
#'
#' The Hypograph Index of a curve x is the proportion of curves in the sample
#' that are below x.
#'
#' @param curves A matrix where each row represents a curve, and each column
#' represents values along the curve.
#'
#' @return Return a numeric vector containing the HI for each curve
#' @export
#'
#' @examples
#' x <- matrix(c(1,2,3,3,2,1,5,2,3,9,8,7), ncol = 3, nrow = 4)
#' HI(x)
HI <- function(curves){
  n_curves <- dim(curves)[1]
  l_curves <- dim(curves)[2]
  index <- apply(curves,1, function(y)
    sum(apply(curves,1,function(x)
      sum(x<=y)==l_curves)))/n_curves
  return (index)
}

#' Modified Epigraph Index (MEI) for a univariate functional dataset.
#'
#' The Modified Epigraph Index of a curve x is one minus the proportion of
#' "time" the curves in the sample are above x
#'
#' @param curves A matrix where each row represents a curve, and each column
#' represents values along the curve.
#'
#' @return Return a numeric vector containing the MEI for each curve
#' @export
#'
#' @examples
#' x <- matrix(c(1,2,3,3,2,1,5,2,3,9,8,7),ncol = 3, nrow = 4)
#' MEI(x)
MEI <- function(curves){
  n_curves <- dim(curves)[1]
  l_curves <- dim(curves)[2]
  rankm <- apply(curves, 2, function(y) (rank(y, ties.method = "min")))
  n_a <- n_curves-rankm+1
  index <- rowSums(n_a)/(n_curves*l_curves)
  return(1-index)
}

#' Modified Hypograph Index (MHI) for a univariate functional dataset.
#'
#' The Modified Hypograph Index of a curve x is the proportion of "time"
#' the curves in the sample are below x.
#'
#' @param curves A matrix where each row represents a curve, and each column
#' represents values along the curve.
#'
#' @return Return a numeric vector containing the MHI for each curve
#' @export
#'
#' @examples
#' x <- matrix(c(1,2,3,3,2,1,5,2,3,9,8,7), ncol = 3, nrow = 4)
#' MHI(x)
MHI <- function(curves){
  n_curves <- dim(curves)[1]
  l_curves <- dim(curves)[2]
  rankm <- apply(curves, 2, function(y) (rank(y, ties.method = "max")))
  index <- rowSums(rankm)/(n_curves*l_curves)
  return(index)
}

#' Hypograph Index (HI) for a multivariate functional dataset.
#'
#' @param curves An array with dimension \eqn{n \times p \times k} with
#' \eqn{n} curves, \eqn{p} values along the curve, and \eqn{k} dimensions.
#'
#' @return A numeric vector containing the HI for each multivariate curve
#' @export
#'
#' @examples
#' x <- array(c(1,2,3, 3,2,1, 5,2,3, 9,8,7, -1,-5,-6, 2,3,0, -1,0,2, -1,-2,0),
#' dim = c(3,4,2))
#' mulHI(x)
mulHI <- function(curves) {
  n_curves <- dim(curves)[1]
  l_curves <- dim(curves)[2]
  d_curves <- dim(curves)[3]

  index <- colSums(Reduce('*',lapply(1:d_curves, function(k)
    sapply(1:n_curves, function(j)
      colSums(sapply(1:n_curves, function(i)
        curves[i,,k] <= curves[j,,k]))==l_curves))))

  return (index/n_curves)
}

#' Modified epigraph Index (MEI) for a multivariate functional dataset.
#'
#' @param curves An array with dimension \eqn{n \times p \times k} with
#' \eqn{n} curves, \eqn{p} values along the curve, and \eqn{k} dimensions.
#'
#' @return A numeric vector containing the MEI for each multivariate curve
#' @export
#'
#' @examples
#' x <- array(c(1,2,3, 3,2,1, 5,2,3, 9,8,7, -1,-5,-6, 2,3,0, -1,0,2, -1,-2,0),
#' dim = c(3,4,2))
#' mulMEI(x)
mulMEI <- function(curves) {
  n_curves <- dim(curves)[1]
  l_curves <- dim(curves)[2]
  d_curves <- dim(curves)[3]

  index <- sapply(1:n_curves, function(j)
    sum(Reduce('*', lapply(1:d_curves, function(k)
      sapply(1:n_curves, function(i)
        curves[i,,k] >= curves[j,,k])))))

  return (1 - index/(n_curves*l_curves))
}

#' Modified hypograph Index (MHI) for a multivariate functional dataset.
#'
#' @param curves An array with dimension \eqn{n \times p \times k} with
#' \eqn{n} curves, \eqn{p} values along the curve, and \eqn{k} dimensions.
#'
#' @return A numeric vector containing the MHI for each multivariate curve
#' @export
#'
#' @examples
#' x <- array(c(1,2,3, 3,2,1, 5,2,3, 9,8,7, -1,-5,-6, 2,3,0, -1,0,2, -1,-2,0),
#' dim = c(3,4,2))
#' mulMHI(x)
mulMHI <- function(curves) {
  n_curves <- dim(curves)[1]
  l_curves <- dim(curves)[2]
  d_curves <- dim(curves)[3]

  index <- sapply(1:n_curves, function(j)
    sum(Reduce('*', lapply(1:d_curves, function(k)
      sapply(1:n_curves, function(i)
        curves[i,,k] <= curves[j,,k])))))

  return (index/(n_curves*l_curves))
}

#' Create a dataset with indexes from a functional dataset in one or multiple
#' dimensions
#'
#' @param curves A matrix with dimension \eqn{n \times p} in the case of a
#' one-dimensional functional dataset, or a matrix of dimension
#' \eqn{n \times p \times k} in the case of a multivariate functional dataset.
#' \eqn{n} represents the number of curves, \eqn{p} the number of values along
#' the curve, and in the second case, \eqn{k} is the number of dimensions.
#' @param t Grid
#' @param nbasis Number of basis for the B-splines
#' @param norder Order of the B-splines
#' @param indices Set of indices to be applied to the dataset. They should be
#' any between EI, HI, MEI and MHI
#' @param ... Additional arguments (unused)
#'
#' @return A dataframe containing the indexes provided in \code{indices} for
#' original data, first and second derivatives
#' @export
#'
#' @examples
#' x1 <- array(c(1,2,3, 3,2,1, 5,2,3, 9,8,7, -1,-5,-6, 2,3,0, -1,0,2, -1,-2,0),
#' dim = c(3,4,2))
#' ind(x1, t = seq(0,1,length=4), nbasis = 4)
#'
#' x2 <- matrix(c(1,2,3,3,2,1,5,2,3,9,8,7), nrow = 3, ncol  = 4)
#' ind(x2, t = seq(0,1,length=4), nbasis = 4)
ind <- function(curves, t, nbasis=25, norder=4,
                     indices = c("EI", "HI", "MEI", "MHI"), ...){

  # Check the curves dimension
  curves_dim <- length(dim(curves))

  # Check indices names
  if(!all(indices %in% c("EI", "HI", "MEI", "MHI"))){
    stop("Indices should be one or more of the following: EI, HI, MEI, MHI")
  }

  if(curves_dim == 2) indices2 <- indices
  else if(curves_dim == 3) indices2 <- paste0("mul", indices)
  else  stop("Invalid number of dimensions")

  # Smoothed data and derivatives
  fun_data <- funspline(curves = curves, t = t, nbasis = nbasis,
                        norder = norder, ...)

  # Initialize an empty data frame to store the results
  ind_data <- as.data.frame(matrix(NA, nrow = nrow(fun_data$smooth), ncol = 0))

  # Loop through the list of functions and apply them to the smoothed and
  # its first and second derivatives
  for (i in 1:length(indices)) {
    smooth_col <- paste0("dta", indices[i])
    deriv_col <- paste0("ddta", indices[i])
    deriv2_col <- paste0("d2dta", indices[i])

    smooth_result <- get(indices2[i])(fun_data$smooth)
    deriv_result <- get(indices2[i])(fun_data$deriv)
    deriv2_result <- get(indices2[i])(fun_data$deriv2)

    ind_data <- cbind(ind_data,
                      stats::setNames(data.frame(smooth_result, deriv_result,
                                                 deriv2_result),
                               c(smooth_col, deriv_col, deriv2_col)))
  }

  return(ind_data)
}



