ehyclus_example_data <- function() {
  vars1 <- c("dtaEI", "dtaMEI")
  vars2 <- c("dtaHI", "dtaMHI")
  vars_combinations <- list(vars1, vars2)
  curves <- sim_model_ex1(seed = 42)

  list(
    curves = curves,
    vars_combinations = vars_combinations
  )
}
