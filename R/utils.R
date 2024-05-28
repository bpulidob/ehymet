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
