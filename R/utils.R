#' Smooth data and calculate first and second derivatives
#'
#' @param curves \code{matrix} where each row represents a curve, and each column
#' represents values along the curve or \code{array} with dimension
#' \eqn{n \times p \times k} with \eqn{n} curves, \eqn{p} values along the curve, and
#' \eqn{k} dimensions.
#' @param nbasis Number of basis for the B-splines
#' @param norder Order of the B-splines
#' @param grid_ll lower limit of the grid.
#' @param grid_ul upper limit of the grid.
#' @param ... Ignored.
#'
#' @return A list containing smoothed data, first and second derivatives
#' @noRd
funspline <- function(curves, nbasis, norder, grid_ll = 0, grid_ul = 1, ...) {
  #Create B-spline basis
  basisobj <- fda::create.bspline.basis(rangeval = c(grid_ll, grid_ul),
                                        nbasis = nbasis, norder = norder, ...)

  curves_dim <- length(dim(curves))
  grid <- seq(from = grid_ll, to = grid_ul, length.out = dim(curves)[2])

  if(curves_dim == 2){

    # Smooth data using B-spline basis
    ys <-  fda::smooth.basis(argvals = grid, y = t(curves), fdParobj = basisobj)

    # Evaluate smoothed data and derivatives
    smooth <- t(fda::eval.fd(grid ,ys$fd,0)) # smoothed data
    deriv <- t(fda::eval.fd(grid ,ys$fd,1)) # first derivatives
    deriv2 <- t(fda::eval.fd(grid ,ys$fd,2)) # second derivatives
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
      ys <-  fda::smooth.basis(argvals = grid, y = t(curves[,,d]),
                               fdParobj = basisobj)

      # Evaluate smoothed data and derivatives
      smooth[,,d] <- t(fda::eval.fd(grid ,ys$fd,0)) # smoothed data
      deriv[,,d] <- t(fda::eval.fd(grid ,ys$fd,1)) # first derivatives
      deriv2[,,d] <- t(fda::eval.fd(grid ,ys$fd,2)) # second derivatives
    }
  } else {
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

#' Generate the name for the results of the clustering methods
#' @noRd
get_result_names <- function(method_name, parameter_combinations, vars_combinations) {
  args <- list(method_name)
  for (combination in parameter_combinations[-1]) {
    args <- append(args, list(combination))
  }

  args <- append(args, list(rep(sapply(vars_combinations, function(x) paste0(x, collapse = "")),
                                times = nrow(parameter_combinations) / length(vars_combinations))
  ))
  args[["sep"]] <- "_"
  do.call(paste, args)
}

#' Suppress outputs from cat (by Hadley Wickham)
#' @noRd
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}


#' Check all combinations of variables and found the non-valid ones
#'
#' @param vars_combinations \code{list} containing one or more combination of variables.
#' @param ind_curves dataset with indices from a functional dataset in one or multiple
#' dimensions.
#'
#' @return Atomic vector with the index of the non-valid combinations of variables.
#'
#' @noRd
check_vars_combinations <- function(vars_combinations, ind_curves) {
  vars_combinations_to_remove <- c()
  for (i in seq_along(vars_combinations)) {
    if (length(vars_combinations[[i]]) == 0) {
      vars_combinations_to_remove <- c(vars_combinations_to_remove, i)
      warning(paste0("Index '", i, "' of 'vars_combinations' is empty.",
                     "Removing it..."))
      next
    }

    if (length(vars_combinations[[i]]) == 1) {
      warning(paste0("Combination of varaibles '", vars_combinations[[i]],
                     "' with index ", i, " is only one variable, which ",
                     "does not have much sense in this context...")
      )
    }

    if (!all(vars_combinations[[i]] %in% names(ind_curves))) {
      vars_combinations_to_remove <- c(vars_combinations_to_remove, i)
      warning(paste0("Invalid variable name in 'vars_combinations' for index ", i,
                     ". Removing combination..."))

      next
    }

    if (det(stats::var(ind_curves[,vars_combinations[[i]]])) == 0) {
      vars_combinations_to_remove <- c(vars_combinations_to_remove, i)

      warning(paste0("Combination of variables '",
                     paste0(vars_combinations[i], collapse = ", "),
                     "' with index ", i, " is singular or almost singular.\n",
                     "Excluding it from any computation...")
      )
    }
  }

  if (length(vars_combinations_to_remove) == length(vars_combinations)) {
    stop("none of the combinations provided in 'vars_combinations' is valid.", call. = FALSE)
  }

  vars_combinations_to_remove
}
