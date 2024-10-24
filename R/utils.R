#' Smooth the data and computes the first and second derivatives
#'
#' @param curves \code{matrix} with dimension \eqn{n \times p} in the case of a
#' one-dimensional functional dataset, or \code{array} of dimension
#' \eqn{n \times p \times q} in the case of a multivariate functional dataset.
#' \eqn{n} represents the number of curves, \eqn{p} the number of values along
#' the curve, and in the second case, \eqn{q} is the number of dimensions.
#' @param k Number of basis functions for the B-splines. If equals to 0, the number
#' of basis functions will be automatically selected.
#' @param bs A two letter character string indicating the (penalized) smoothing
#' basis to use. See \code{\link{smooth.terms}}.
#'
#' @return A list containing smoothed data, first and second derivatives
#'
#' @noRd
funspline <- function(curves, k, bs = "cr") {
  curves_dim <- length(dim(curves))

  tfb_params <- list(bs = bs)

  if (!missing(k)) {
    tfb_params[["k"]] <- k
  }

  if (curves_dim == 2) {
    tfb_params[["data"]] <- curves

    ys <- suppressMessages(do.call(tf::tfb, tfb_params))

    # Evaluate smoothed data and derivatives
    smooth <- as.matrix(ys) # smoothed data

    deriv <- as.matrix(tf::tf_derive(ys, order = 1))
    deriv2 <- as.matrix(tf::tf_derive(ys, order = 2))

  } else {
    n_curves <- dim(curves)[1]
    l_curves <- dim(curves)[2]
    d_curves <- dim(curves)[3]

    # Initialize empty dataframes to store the results
    smooth <- array(rep(NaN, n_curves * l_curves), dim = c(n_curves, l_curves, d_curves))
    deriv <- array(rep(NaN, n_curves * (l_curves)), dim = c(n_curves, l_curves, d_curves))
    deriv2 <- array(rep(NaN, n_curves * (l_curves)), dim = c(n_curves, l_curves, d_curves))

    for (d in seq_len(dim(curves)[3])) {
      tfb_params[["data"]] <- curves[, , d]

      # Smooth data using B-spline basis
      ys <- suppressMessages(do.call(tf::tfb, tfb_params))

      # Evaluate smoothed data and derivatives
      smooth[, , d] <- as.matrix(ys) # smoothed data
      deriv[, , d]  <- as.matrix(tf::tf_derive(ys, order = 1))
      deriv2[, , d] <- as.matrix(tf::tf_derive(ys, order = 2))

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
    stop("duplicated argument in '", parameter_name, "'.", call. = FALSE)
  }

  indices <- pmatch(argument, parameter_values)
  if (any(is.na(indices))) {
    stop("invalid argument in '", parameter_name, "': ", paste(argument[is.na(indices)], collapse = ", "), ".",
      call. = FALSE
    )
  }
}

#' Generate the name for the results of the clustering methods
#'
#' @noRd
get_result_names <- function(method_name, parameter_combinations, vars_combinations) {
  args <- list(method_name)
  for (combination in parameter_combinations[-1]) {
    args <- append(args, list(combination))
  }

  args <- append(args, list(rep(sapply(vars_combinations, function(x) paste0(x, collapse = "")),
    times = nrow(parameter_combinations) / length(vars_combinations)
  )))
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
#'
#' @noRd
check_n_cores <- function(n_cores) {
  if (n_cores < 1 || n_cores %% 1 != 0) {
    stop("'n_cores' must be an integer number greater than '1'", call. = FALSE)
  }

  if (.Platform$OS.type != "unix" && n_cores > 1) {
    warning(
      "Running this function using multiples cores is only supported on unix systems.",
      "Setting 'n_cores' parameter to '1'"
    )
    return(1)
  }

  n_cores
}

#' Replacement for caret::findCorrelation
#'
#' @noRd
findCorrelation <- function(cor_matrix, cutoff = 0.75) {
  # Get the average absolute correlation for each variable
  avg_cor <- colMeans(abs(cor_matrix - diag(rep(1, ncol(cor_matrix)))))
  var_names <- colnames(cor_matrix)
  remove_vars <- integer()

  # While there are still correlations above the cutoff
  while(TRUE) {
    # Get absolute correlations excluding diagonal
    cor_tri <- abs(cor_matrix)
    diag(cor_tri) <- 0

    # Find maximum correlation
    max_cor <- max(cor_tri, na.rm = TRUE)
    if(max_cor < cutoff) break

    # Find variable pairs with maximum correlation
    max_indices <- which(cor_tri == max_cor, arr.ind = TRUE)

    # If no more high correlations found, break
    if(nrow(max_indices) == 0) break

    # For each pair, compare mean correlations
    for(i in 1:nrow(max_indices)) {
      rows <- max_indices[i, 1]
      cols <- max_indices[i, 2]

      # Remove variable with higher average correlation
      if(avg_cor[rows] > avg_cor[cols]) {
        if(!(rows %in% remove_vars)) {
          remove_vars <- c(remove_vars, rows)
          # Zero out correlations for removed variable
          cor_matrix[rows, ] <- cor_matrix[, rows] <- 0
        }
      } else {
        if(!(cols %in% remove_vars)) {
          remove_vars <- c(remove_vars, cols)
          # Zero out correlations for removed variable
          cor_matrix[cols, ] <- cor_matrix[, cols] <- 0
        }
      }
    }
  }

  sort(remove_vars)
}


#' Select an unique combination of indices
#'
#' @noRd
select_var_ind <- function(data_ind) {
  # Calculate ratio of distinct values
  rat <- sapply(data_ind, \(x) length(unique(x)) / length(x))
  rat_var <- names(rat)[rat > 0.5]

  data_ind_f_rat <- data_ind[, rat_var, drop = FALSE]

  cor_data <- stats::cor(data_ind_f_rat)

  # Get positions of highly correlated variables
  cor_var_pos <- findCorrelation(cor_data)

  # Get the names of variables to remove
  if (length(cor_var_pos) > 0) {
    cor_var <- colnames(data_ind_f_rat)[cor_var_pos]
  } else {
    cor_var <- character(0)
  }

  # Get final selected variable names
  ind_f_rat_cor_var <- setdiff(colnames(data_ind_f_rat), cor_var)

  # Return the final indices
  data_ind_f_rat[, ind_f_rat_cor_var, drop = FALSE]
}

