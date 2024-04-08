ehyclus_example_data <- function() {
  vars1 <- c("dtaEI", "dtaMEI")
  vars2 <- c("dtaHI", "dtaMHI")
  vars_list <- list(vars1, vars2)
  curves <- sim_model_ex1()
  t <- seq(0, 1, length = 30)

  list(
    curves = curves,
    t = t,
    vars_list = vars_list
  )
}
