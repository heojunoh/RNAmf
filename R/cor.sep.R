#' cor.sep
#'
#' calculating separable matern kernel
#'
#' @param X vector or matrix of input.
#' @param x vector or matrix of new input. Default is NULL
#' @param theta lengthscale parameter. It should have the length of ncol(X).
#' @param nu numerical value of smoothness hyperparameter. It should be 0.5, 1.5, 2.5, 3.5, or 4.5.
#' @param derivative logical indicating for its first derivative(derivative=1)
#' @noRd
#' @keywords internal
#' @return A covariance matrix of matern kernel.

cor.sep <- function(X, x=NULL, theta, nu, derivative=0){
d <- NCOL(X)
n <- NROW(X)
nu <- rep(nu, d)
if(is.null(x)){
  K <- matrix(1, n, n)
  for(i in 1:d){
    R <- sqrt(distance(X[,i]/theta[i]))
    K <- K * matern.kernel(R, nu=nu[i], derivative=derivative)
  }
}else{
  n.new <- NROW(x)
  K <- matrix(1, n, n.new)
  for(i in 1:d){
    R <- sqrt(distance(X[,i]/theta[i], x[,i]/theta[i]))
    K <- K * matern.kernel(R, nu=nu[i], derivative=derivative)
  }
}
return(K)
}
