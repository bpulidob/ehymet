#' Function for generating functional data in one dimension
#'
#' Each dataset has 2 groups with \code{n} curves each, defined in the interval
#' \eqn{t=[0, 1]} with \code{p} equidistant points. The first \code{n} curves are
#' generated fron the following model
#' \eqn{X_1(t)=E_1(t)+e(t)} where \eqn{E_1(t)=E_1(X(t))=30t^{ \frac{3}{2}}(1-t)}
#' is the mean function and \eqn{e(t)} is a centered Gaussian process with
#' covariance matrix \eqn{Cov(e(t_i),e(t_j))=0.3 \exp(-\frac{\lvert t_i-t_j \rvert}{0.3})}
#' The remaining 50 functions are generated from model \code{i_sim} with
#' \code{i_sim} \eqn{\in \{1, \ldots, 8\}.}
#' The first three models contain changes in the mean, while the covariance
#' matrix does not change. Model 4 and 5 are obtained by multiplying the
#' covariance matrix by a constant. Model 6 is obtained from adding to
#' \eqn{E_1(t)} a centered Gaussian process \eqn{h(t)} whose covariance matrix
#' is given by \eqn{Cov(e(t_i),e(t_j))=0.5 \exp (-\frac{\lvert t_i-t_j\rvert}{0.2})}.
#' Model 7 and 8 are obtained by a different mean function.
#' \describe{
#'    \item{Model 1.}{\eqn{X_1(t)=30t^{\frac{3}{2}}(1-t)+0.5+e(t).}}
#'    \item{Model 2.}{\eqn{X_2(t)=30t^{\frac{3}{2}}(1-t)+0.75+e(t).}}
#'    \item{Model 3.}{\eqn{X_3(t)=30t^{\frac{3}{2}}(1-t)+1+e(t).}}
#'    \item{Model 4.}{\eqn{X_4(t)=30t^{\frac{3}{2}}(1-t)+2 e(t).}}
#'    \item{Model 5.}{\eqn{X_5(t)=30t^{\frac{3}{2}}(1-t)+0.25 e(t).}}
#'    \item{Model 6.}{\eqn{X_6(t)=30t^{\frac{3}{2}}(1-t)+ h(t).}}
#'    \item{Model 7.}{\eqn{X_7(t)=30t{(1-t)}^2+ h(t).}}
#'    \item{Model 8.}{\eqn{X_8(t)=30t{(1-t)}^2+ e(t).}}
#' }
#'
#' @param n Number of curves to generate for each of the two groups.
#' Set to 50 by default
#' @param p Number of grid points of the curves.
#' Curves are generated over the interval \eqn{[0, 1]}.
#' Set to 30 grid point by default
#' @param i_sim Integer set to \eqn{1, \ldots, 8}
#' @param seed A seed to set for reproducibility. NULL by default in which case a seed is not set.
#' @return data matrix of size \eqn{2n \times p}
#' @export
#'
#' @examples
#' sm1 <- sim_model_ex1()
#' dim(sm1)

sim_model_ex1 <- function(n = 50, p = 30, i_sim = 1, seed = NULL){

  if (!(i_sim %in% 1:8))
    stop("Argument \"i_sim\" should be one between 1 and 8.")

  # set seed
  if(!is.null(seed)) set.seed(seed)

  t_interval <- seq(0, 1, length = p)

  S <- function(t) exp(-1/0.3 * abs(outer(t,t,"-")))
  sigma <- 0.3 * (t_interval |> S())
  # centered Gaussian process
  egauss <- MASS::mvrnorm(n ,rep(0, p), sigma)

  S2 <- function(t) exp(-1/0.2 * abs(outer(t,t,"-")))
  sigma2 <- 0.5 * (t_interval |> S2())
  # centered Gaussian process
  hgauss <- MASS::mvrnorm(n ,rep(0, p), sigma2)

  f1 <- (30 * t_interval^(3/2) * (1 - t_interval))
  matr <- matrix(rep(f1, each = n), nrow = n)
  dat1 <- matr + egauss

  f2 <- (30 * t_interval * (1 - t_interval))
  matr2 <- matrix(rep(f2, each = n), nrow = n)

  if(i_sim == 1){
    dat2 <- dat1 + 0.5
  } else if(i_sim == 2){
    dat2 <- dat1 + 0.75
  } else if(i_sim == 3){
    dat2 <- dat1 + 1
  } else if(i_sim == 4){
    dat2 <- matr + 2*egauss
  } else if(i_sim == 5){
    dat2 <- matr + 0.25*egauss
  } else if(i_sim == 6){
    dat2 <- matr + hgauss
  } else if(i_sim == 7){
    dat2 <- matr2 + hgauss
  } else{
    dat2 <- matr2 + egauss
  }

  dat <- rbind(dat1, dat2)

  return(dat)
}


#' Function for generating functional data in one or multiple dimension
#'
#' @param n Number of curves to generate for each of the two groups.
#' Set to 50 by default
#' @param p Number of grid points of the curves.
#' Curves are generated over the interval \eqn{[0, 1]}.
#' Set to 150 grid point by default
#' @param i_sim Integer set to \eqn{1, \ldots, 4}
#' @param seed A seed to set for reproducibility.
#' NULL by default in which case a seed is not set.
#'
#' @return data matrix of size \eqn{2n \times p} if \eqn{\code{i_sim} \in {1,2}}
#' or an array of dimensions
#' \eqn{2n \times p \times 2} if \eqn{\code{i_sim} \in {3, 4}}
#' @export
#'
#' @examples
#' sm1 <- sim_model_ex2()
#' dim(sm1)
#'
#' sm4 <- sim_model_ex2(i_sim=4)
#' dim(sm4)
#'
#' @importFrom dplyr "%>%"
sim_model_ex2 <- function(n = 50, p = 150, i_sim = 1, seed = NULL){

  if (!(i_sim %in% 1:4))
    stop("Argument \"i_sim\" should be between 1 and 4.")

  # set seed
  if(!is.null(seed)) set.seed(seed)

  t_interval <- seq(0, 1, length = p)

  K <- 100
  if(i_sim %in% c(1,2)){
    E1 <- t_interval * (1 - t_interval)
  } else{
    E1 <- rbind(t_interval * (1 - t_interval), 4*t_interval^2 * (1 - t_interval))
  }


  rho <- ifelse(1:K < 4, 1 / (1:K + 1), 1 / (1:K + 1)^2)

  theta <- sapply(1:K, function(k) {
    if (k %% 2 == 0) sqrt(2) * sin(k * pi * t_interval)
    else if (k %% 2 != 0 && k != 1) sqrt(2) * cos((k - 1) * pi * t_interval)
    else rep(1, p)
  }) %>% t()

  m <- sqrt(rho)*theta

  if(i_sim ==1){
    E2 <- E1 + colSums(m[1:4,])
  } else if(i_sim == 2){
    E2 <- E1 + colSums(m[5:K,])
  } else if(i_sim == 3){
    E2 <- (t(E1) + colSums(m[1:4,])) |> t()
  } else{
    E2 <- (t(E1) + colSums(m[5:K,])) |> t()
  }

  if(i_sim %in% c(1,2)){
    z1 <- (matrix(stats::rnorm(n*K), K, n) * sqrt(rho)) |> t()
    z2 <- (matrix(stats::rnorm(n*K), K, n) * sqrt(rho)) |> t()

    X1 <- (t(z1 %*% theta) + E1) |> t()
    X2 <- (t(z2 %*% theta) + E2) |> t()
    X <- rbind(X1,X2)
  } else{
    z1 <- (array(stats::rnorm(n*K*2), dim =c(K, n, 2)) * sqrt(rho)) |>
      aperm(perm = c(2,1,3))
    z2 <- (array(stats::rnorm(n*K*2), dim =c(K, n, 2)) * sqrt(rho)) |>
      aperm(perm = c(2,1,3))

    X1 <- unlist(lapply(1:2, function(k) t(z1[ , , k] %*% theta) + E1[k, ])) |>
      array(dim = c(p, n, 2)) |> aperm(perm = c(2,1,3))
    X2 <- unlist(lapply(1:2, function(k) t(z2[ , , k] %*% theta) + E2[k, ])) |>
      array(dim = c(p, n, 2)) |> aperm(perm = c(2,1,3))
    X <- unlist(lapply(1:2, function(k) rbind(X1[ , , k], X2[ , , k]))) |>
      array(dim = c(2*n, p, 2))
  }

  return(X)
}


#' Function for generating plots of one dimensional functional datasets
#'
#' @param data data matrix of size \eqn{n \times p}
#' @param true_labels array containing the true groups in which the data should
#' be classified
#'
#' @return A plot
#' @export
#' @importFrom dplyr "%>%"
#'
#' @examples
#' dat1 <- sim_model_ex1()
#' true_labels <- c(rep(1,50), rep(2,50))
#' plt_fun(dat1, true_labels)
plt_fun <- function(data, true_labels){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("package 'ggplot2' is required for this functionality", call. = FALSE)
  }

  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("package 'tidyr' is required for this functionality", call. = FALSE)
  }

  if (!(length(dim(data)) == 2)) {
    stop("This function can be only used with 2-dimensional datasets.", call. = FALSE)
  }

  df <-  dplyr::as_tibble(data)
  t_interval <- seq(0, 1, length = ncol(data))
  names(df) <- as.character(t_interval)
  df$id <- 1:nrow(data)
  df$Order <- true_labels
  df_long<- df %>% tidyr::pivot_longer(-c(id, Order), names_to="variable", values_to="values") %>%
    dplyr::mutate(variable=as.numeric(variable))
  pa <- df_long %>% ggplot2::ggplot(ggplot2::aes(x=variable, y=values,group=id, color=factor(Order)))

  plt <- pa +
    ggplot2::geom_line(linewidth=0.1)+
    ggplot2::scale_color_brewer(palette = "Set1")+
    # scale_color_manual(values=c("#CC6600","#3399FF")) +
    # ggtitle("MEI. First dimension")+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))+
    ggplot2::ylab("") + ggplot2::xlab("") +
    ggplot2::theme(legend.position = "none")
  return(plt)
}

