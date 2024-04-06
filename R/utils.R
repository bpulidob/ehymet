#' Smooth data and calculate first and second derivatives
#'
#' @param curves A matrix where each row represents a curve, and each column
#' represents values along the curve.
#' @param t Grid
#' @param nbasis Number of basis for the B-splines
#' @param norder Order of the B-splines
#' @param ... Additional arguments (unused)
#'
#' @return A list containing smoothed data, first and second derivatives
#' @noRd
#'
funspline <- function(curves, t, nbasis, norder, ...){
  #Create B-spline basis
  basisobj <- fda::create.bspline.basis(rangeval = c(min(t), max(t)),
                                        nbasis = nbasis, norder = norder, ...)

  curves_dim <- length(dim(curves))

  if(curves_dim == 2){

    # Smooth data using B-spline basis
    ys <-  fda::smooth.basis(argvals = t, y = t(curves), fdParobj = basisobj)

    # Evaluate smoothed data and derivatives
    smooth <- t(fda::eval.fd(t,ys$fd,0)) # smoothed data
    deriv <- t(fda::eval.fd(t,ys$fd,1)) # first derivatives
    deriv2 <- t(fda::eval.fd(t,ys$fd,2)) # second derivatives
  } else if(curves_dim == 3){
    n_curves <- dim(curves)[1]
    l_curves <- dim(curves)[2]
    d_curves <- dim(curves)[3]

    # Initialize empty dataframes to store the results
    smooth <- array(rep(NaN,n_curves*l_curves),dim=c(n_curves,l_curves,d_curves))
    deriv <- array(rep(NaN,n_curves*l_curves),dim=c(n_curves,l_curves,d_curves))
    deriv2 <- array(rep(NaN,n_curves*l_curves),dim=c(n_curves,l_curves,d_curves))

    for(d in 1:dim(curves)[3]){
      # Smooth data using B-spline basis
      ys <-  fda::smooth.basis(argvals = t, y = t(curves[,,d]),
                               fdParobj = basisobj)

      # Evaluate smoothed data and derivatives
      smooth[,,d] <- t(fda::eval.fd(t,ys$fd,0)) # smoothed data
      deriv[,,d] <- t(fda::eval.fd(t,ys$fd,1)) # first derivatives
      deriv2[,,d] <- t(fda::eval.fd(t,ys$fd,2)) # second derivatives
    }
  } else{
    stop("Invalid number of dimensions")
  }

  # Return a list containing the data and derivatives
  res <- list(
    "smooth" = smooth,
    "deriv" = deriv,
    "deriv2" = deriv2
    )

  return(res)
}

#' Transform metrics results from clustering functions in cluster.R to a dataframe
#'
#' @param res list containing clustering partition and metric for different
#' combinations
#' @param tl_null a bool to indicate weather metrics other than time ane or not
#' available
#'
#' @return Dataframe
#' @noRd
result_to_table <- function(res, tl_null){
  name_res <- names(res)
  len_res <- length(name_res)
  if(tl_null){
      metrics_df <-
        data.frame(Time = sapply(1:len_res, function (i) res[[i]][[2]]))
      row.names(metrics_df) <- name_res
  } else{
    metrics_df <-
      data.frame(t(sapply(1:len_res, function (i) c(res[[i]][[2]],
                                                    "Time" = res[[i]][[3]]))))
    row.names(metrics_df) <- name_res
    }
  return(metrics_df)
}

#' Checks for list function arguments
#'
#' Checks that a list given as argument to a function is not empty,
#' has no repeated values and all its elements are within a bounded set.
#'
#' @param argument Input argument of the function.
#' @param parameter_values Possible values that may appear in the list.
#' @param parameter_name Name of the parameter.
#'
#' @noRd
check_list_parameter <- function(argument, parameter_values, parameter_name) {
  if (length(argument) == 0) {
    stop("parameter '", parameter_name, "' should have at least one element.", call. = FALSE)
  }

  if (any(duplicated(argument))) {
    stop("duplicated argument in '", parameter_name,"'.", call. = FALSE)
  }

  indices <- pmatch(argument, parameter_values)
  if (any(is.na(indices))) {
    stop("invalid argument in '", parameter_name, "': ", paste(argument[is.na(indices)], collapse = ", "), ".",
         call. = FALSE)
  }
}
