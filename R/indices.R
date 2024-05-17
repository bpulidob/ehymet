#' Epigraph Index (EI) for a univariate functional dataset.
#'
#' The Epigraph Index of a curve x is one minus the proportion of curves
#' in the sample that are above x.
#'
#' @param curves \code{matrix} where each row represents a curve, and each column
#' represents values along the curve or \code{array} with dimension
#' \eqn{n \times p \times k} with \eqn{n} curves, \eqn{p} values along the curve, and
#' \eqn{k} dimensions.
#'
#' @param ... Ignored.
#'
#' @return numeric \code{vector} containing the EI for each curve.
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
EI <- function(curves, ...) {
  UseMethod("EI")
}

#' @export
EI.matrix <- function(curves, ...) {
  n_curves <- dim(curves)[1]
  l_curves <- dim(curves)[2]

  index <- apply(curves,1, function(y)
    sum(apply(curves,1,function(x)
      sum(x>=y)==l_curves)))/n_curves

  return(1 - index)
}

#' @export
EI.array <- function(curves, ...) {
  if (length(dim(curves)) != 3) {
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

#' @export
EI.default <- function(curves, ...) {
  stop("'curves' should be a matrix or a 3-dimensional array", call. = FALSE)
}

#' Hypograph Index (HI) for a univariate functional dataset.
#'
#' The Hypograph Index of a curve x is the proportion of curves in the sample
#' that are below x.
#'
#' @param curves \code{matrix} where each row represents a curve, and each column
#' represents values along the curve or \code{array} with dimension
#' \eqn{n \times p \times k} with \eqn{n} curves, \eqn{p} values along the curve, and
#' \eqn{k} dimensions.
#'
#' @param ... Ignored.
#'
#' @return \code{numeric vector} containing the HI for each curve.
#'
#' @examples
#' x <- matrix(c(1,2,3,3,2,1,5,2,3,9,8,7), ncol = 3, nrow = 4)
#' HI(x)
#'
#' y <- array(c(1,2,3, 3,2,1, 5,2,3, 9,8,7, -1,-5,-6, 2,3,0, -1,0,2, -1,-2,0),
#' dim = c(3,4,2))
#' HI(y)
#'
#' @export
HI <- function(curves, ...) {
  UseMethod("HI")
}

#' @export
HI.matrix <- function(curves, ...) {
  n_curves <- dim(curves)[1]
  l_curves <- dim(curves)[2]

  index <- apply(curves,1, function(y)
    sum(apply(curves,1,function(x)
      sum(x<=y)==l_curves)))/n_curves

  return(index)
}

#' @export
HI.array <- function(curves, ...) {
  if (length(dim(curves)) != 3) {
    stop("'curves' should be a matrix or a 3-dimensional array", call. = FALSE)
  }

  n_curves <- dim(curves)[1]
  l_curves <- dim(curves)[2]
  d_curves <- dim(curves)[3]

  index <- colSums(Reduce('*',lapply(1:d_curves, function(k)
    sapply(1:n_curves, function(j)
      colSums(sapply(1:n_curves, function(i)
        curves[i,,k] <= curves[j,,k]))==l_curves))))

  return(index/n_curves)
}

#' @export
HI.default <- function(curves, ...) {
  stop("'curves' should be a matrix or a 3-dimensional array", call. = FALSE)
}

#' Modified Epigraph Index (MEI) for a univariate functional dataset.
#'
#' The Modified Epigraph Index of a curve x is one minus the proportion of
#' "time" the curves in the sample are above x
#'
#' @param curves \code{matrix} where each row represents a curve, and each column
#' represents values along the curve or an \code{array} with dimension
#' \eqn{n \times p \times k} with \eqn{n} curves, \eqn{p} values along the curve, and
#' \eqn{k} dimensions.
#'
#' @param ... Ignored.
#'
#' @return \code{numeric vector} containing the MEI for each curve.
#'
#' @examples
#' x <- matrix(c(1,2,3,3,2,1,5,2,3,9,8,7),ncol = 3, nrow = 4)
#' MEI(x)
#' y <- array(c(1,2,3, 3,2,1, 5,2,3, 9,8,7, -1,-5,-6, 2,3,0, -1,0,2, -1,-2,0),
#' dim = c(3,4,2))
#' MEI(y)
#'
#' @export
MEI <- function(curves, ...) {
  UseMethod("MEI")
}

#' @export
MEI.matrix <- function(curves, ...) {
  n_curves <- dim(curves)[1]
  l_curves <- dim(curves)[2]
  rankm <- apply(curves, 2, function(y) (rank(y, ties.method = "min")))
  n_a <- n_curves-rankm+1
  index <- rowSums(n_a)/(n_curves*l_curves)
  return(1 - index)
}

#' @export
MEI.array <- function(curves, ...) {
  if (length(dim(curves)) != 3) {
    stop("'curves' should be a matrix or a 3-dimensional array", call. = FALSE)
  }

  n_curves <- dim(curves)[1]
  l_curves <- dim(curves)[2]
  d_curves <- dim(curves)[3]

  index <- sapply(1:n_curves, function(j)
    sum(Reduce('*', lapply(1:d_curves, function(k)
      sapply(1:n_curves, function(i)
        curves[i,,k] >= curves[j,,k])))))

  return (1 - index / (n_curves * l_curves))
}

#' @export
MEI.default <- function(curves, ...) {
  stop("'curves' should be a matrix or a 3-dimensional array", call. = FALSE)
}

#' Modified Hypograph Index (MHI) for a univariate functional dataset.
#'
#' The Modified Hypograph Index of a curve x is the proportion of "time"
#' the curves in the sample are below x.
#'
#' @param curves \code{matrix} where each row represents a curve, and each column
#' represents values along the curve or an \code{array} with dimension
#' \eqn{n \times p \times k} with \eqn{n} curves, \eqn{p} values along the curve, and
#' \eqn{k} dimensions.
#'
#' @param ... Ignored.
#'
#' @return \code{numeric vector} containing the MHI for each curve.
#'
#' @examples
#' x <- matrix(c(1,2,3,3,2,1,5,2,3,9,8,7), ncol = 3, nrow = 4)
#' MHI(x)
#' y <- array(c(1,2,3, 3,2,1, 5,2,3, 9,8,7, -1,-5,-6, 2,3,0, -1,0,2, -1,-2,0),
#' dim = c(3,4,2))
#' MHI(y)
#'
#' @export
MHI <- function(curves, ...) {
  UseMethod("MHI")
}

#' @export
MHI.matrix <- function(curves, ...) {
  n_curves <- dim(curves)[1]
  l_curves <- dim(curves)[2]
  rankm <- apply(curves, 2, function(y) (rank(y, ties.method = "max")))
  index <- rowSums(rankm) / (n_curves * l_curves)
  return(index)
}

#' @export
MHI.array <- function(curves, ...) {
  if (length(dim(curves)) != 3) {
    stop("'curves' should be a matrix or a 3-dimensional array", call. = FALSE)
  }

  n_curves <- dim(curves)[1]
  l_curves <- dim(curves)[2]
  d_curves <- dim(curves)[3]

  index <- sapply(1:n_curves, function(j)
    sum(Reduce('*', lapply(1:d_curves, function(k)
      sapply(1:n_curves, function(i)
        curves[i,,k] <= curves[j,,k])))))

  return (index / (n_curves * l_curves))
}

#' @export
MHI.default <- function(curves, ...) {
  stop("'curves' should be a matrix or a 3-dimensional array", call. = FALSE)
}

#' Create a dataset with indexes from a functional dataset in one or multiple
#' dimensions
#'
#' @param curves \code{matrix} with dimension \eqn{n \times p} in the case of a
#' one-dimensional functional dataset, or \code{array} of dimension
#' \eqn{n \times p \times k} in the case of a multivariate functional dataset.
#' \eqn{n} represents the number of curves, \eqn{p} the number of values along
#' the curve, and in the second case, \eqn{k} is the number of dimensions.
#' @param nbasis Number of basis for the B-splines.
#' @param norder Order of the B-splines.
#' @param grid_ll Lower limit of the grid.
#' @param grid_ul Upper limit of the grid.
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
#' ind(x1, nbasis = 4)
#'
#' x2 <- matrix(c(1,2,3,3,2,1,5,2,3,9,8,7), nrow = 3, ncol  = 4)
#' ind(x2, nbasis = 4)
#'
#' @export
ind <- function(curves, nbasis = 30, norder = 4, grid_ll = 0, grid_ul = 1,
                     indices = c("EI", "HI", "MEI", "MHI"), ...) {
  # define indices constant
  INDICES <- c("EI", "HI", "MEI", "MHI")
  curves_dim <- length(dim(curves))

  # stop conditions
  if (!(curves_dim %in% c(2, 3)) || is.null(curves_dim)) {
    stop("'curves' should be a matrix or a 3-dimensional array", call. = FALSE)
  }

  check_list_parameter(indices, INDICES, "indices")

  # Smoothed data and derivatives
  fun_data <- funspline(curves = curves, nbasis = nbasis,
                        norder = norder, grid_ll = grid_ll,
                        grid_ul = grid_ul, ...)

  # Initialize an empty data frame to store the results
  ind_data <- as.data.frame(matrix(NA, nrow = nrow(fun_data$smooth), ncol = 0))

  # Loop through the list of functions and apply them to the smoothed and
  # its first and second derivatives
  for (index in indices) {
    smooth_col <- paste0("dta", index)
    deriv_col  <- paste0("ddta", index)
    deriv2_col <- paste0("d2dta", index)

    smooth_result <- get(index)(fun_data$smooth)
    deriv_result  <- get(index)(fun_data$deriv)
    deriv2_result <- get(index)(fun_data$deriv2)

    ind_data <- cbind(ind_data,
                      stats::setNames(data.frame(smooth_result, deriv_result,
                                                 deriv2_result),
                               c(smooth_col, deriv_col, deriv2_col)))
  }

  ind_data
}

