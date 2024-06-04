#' Smooth the data and computes the first and second derivatives
#'
#' @param curves \code{matrix} with dimension \eqn{n \times p} in the case of a
#' one-dimensional functional dataset, or \code{array} of dimension
#' \eqn{n \times p \times k} in the case of a multivariate functional dataset.
#' \eqn{n} represents the number of curves, \eqn{p} the number of values along
#' the curve, and in the second case, \eqn{k} is the number of dimensions.
#' @param k Number of basis functions for the B-splines. If equals to 0, the number
#' of basis functions will be automatically selected.
#' @param bs A two letter character string indicating the (penalized) smoothing
#' basis to use. See \code{\link{smooth.terms}}.
#' @param grid Atomic vector of type numeric with two elements: the lower limit and the upper
#' limit of the evaluation grid. If not provided, it will be selected automatically.
#'
#' @return A list containing smoothed data, first and second derivatives
#' @noRd
funspline <- function(curves, k, bs = "cr", grid) {
  curves_dim <- length(dim(curves))

  tfb_params <- list(bs = bs)

  if (!missing(k)) {
    tfb_params[["k"]] <- k
  }

  if (!missing(grid)) {
    grid <- seq(grid[1], grid[2], length.out = dim(curves)[2])
    tfb_params[["arg"]] <- grid
  }

  if (curves_dim == 2) {
    tfb_params[["data"]] <- curves

    ys <- suppressMessages(do.call(tf::tfb, tfb_params))

    # Evaluate smoothed data and derivatives
    smooth <- as.matrix(ys) # smoothed data

    if (!missing(grid)) {
      deriv  <- as.matrix(tf::tf_derive(ys, arg = grid, order = 1))
      deriv2 <- as.matrix(tf::tf_derive(ys, arg = grid, order = 2))
    } else {
      deriv  <- as.matrix(tf::tf_derive(ys, order = 1))
      deriv2 <- as.matrix(tf::tf_derive(ys, order = 2))
    }

  } else {
    n_curves <- dim(curves)[1]
    l_curves <- dim(curves)[2]
    d_curves <- dim(curves)[3]

    # Initialize empty dataframes to store the results
    smooth <- array(rep(NaN, n_curves * l_curves),   dim = c(n_curves, l_curves, d_curves))
    deriv  <- array(rep(NaN, n_curves * (l_curves)), dim = c(n_curves, l_curves, d_curves))
    deriv2 <- array(rep(NaN, n_curves * (l_curves)), dim = c(n_curves, l_curves, d_curves))

    for (d in seq_len(dim(curves)[3])) {
      tfb_params[["data"]] <- curves[,,d]

      # Smooth data using B-spline basis
      ys <- suppressMessages(do.call(tf::tfb, tfb_params))

      # Evaluate smoothed data and derivatives
      smooth[,,d] <- as.matrix(ys) # smoothed data
      if (!missing(grid)) {
        deriv[,,d]  <- as.matrix(tf::tf_derive(ys, arg = grid, order = 1))
        deriv2[,,d] <- as.matrix(tf::tf_derive(ys, arg = grid, order = 2))
      } else {
        deriv[,,d]  <- as.matrix(tf::tf_derive(ys, order = 1))
        deriv2[,,d] <- as.matrix(tf::tf_derive(ys, order = 2))
      }
    }
  }

  list(
    "smooth" = smooth,
    "deriv"  = deriv,
    "deriv2" = deriv2
  )
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

#' Maps the index name to its corresponding function
#'
#' Used to ensure that we call our index functions inside \code{generate_indices}
#'
#' @noRd
map_index_name_to_function <- function(index) {
  switch(index,
    "EI"  = EI,
    "HI"  = HI,
    "MEI" = MEI,
    "MHI" = MHI
  )
}

#' Disable multicore on windows without killing the execution
#' @noRd
check_n_cores <- function(n_cores) {
  if (n_cores < 1 || n_cores %% 1 != 0) {
    stop("'n_cores' must be an integer number greater than '1'", call. = FALSE)
  }

  if (.Platform$OS.type != "unix" && n_cores > 1) {
    warning("Running this function using multiples cores is only supported on unix systems.",
            "Setting 'n_cores' parameter to '1'")
    return(1)
  }

  return(n_cores)
}
