#' closed
#'
#' Fitting the model with two fidelity levels
#'
#' @param X1 vector or matrix of input locations for the low fidelity level.
#' @param y1 vector of response values for the low fidelity level.
#' @param X2 vector or matrix of input locations for the high fidelity level.
#' @param y2 vector of response values for the high fidelity level.
#' @param kernel chracter specifying kernel type to be used, to be chosen between "sqex"(squared exponential), "matern1.5", or "matern2.5"
#' @param constant logical indicating for constant mean (constant=TRUE) or zero mean (constant=FALSE). Default is FALSE.
#' @return A list containing fitted models f1 and f2, constant, and kernel structure:
#' \itemize{
#'   \item \code{fit1}: list of fitted model for (X1, y1).
#'   \item \code{fit2}: list of fitted model for ((X2, f1(X2)), y2).
#'   \item \code{constant}: copy of constant.
#'   \item \code{kernel}: copy of kernel.
#' }
#' @export
#'

closed <- function(X1, y1, X2, y2, kernel, constant=FALSE){
  if(kernel=="sqex"){
    if(constant){
      fit1 <- GP(X1, y1, constant=TRUE)
      fit2 <- GP(cbind(X2, pred.GP(fit1, X2)$mu), y2, constant=TRUE)
    }else{
      fit1 <- GP(X1, y1)
      fit2 <- GP(cbind(X2, pred.GP(fit1, X2)$mu), y2)
    }
  }else if(kernel=="matern1.5"){
    if(constant){
      fit1 <- matGP(X1, y1, nu=1.5, constant=TRUE)
      fit2 <- matGP(cbind(X2, pred.matGP(fit1, X2)$mu), y2, nu=1.5, constant=TRUE)
    }else{
      fit1 <- matGP(X1, y1, nu=1.5)
      fit2 <- matGP(cbind(X2, pred.matGP(fit1, X2)$mu), y2, nu=1.5)
    }
  }else if(kernel=="matern2.5"){
    if(constant){
      fit1 <- matGP(X1, y1, nu=2.5, constant=TRUE)
      fit2 <- matGP(cbind(X2, pred.matGP(fit1, X2)$mu), y2, nu=2.5, constant=TRUE)
    }else{
      fit1 <- matGP(X1, y1, nu=2.5)
      fit2 <- matGP(cbind(X2, pred.matGP(fit1, X2)$mu), y2, nu=2.5)
    }
  }

  return(list(fit1=fit1, fit2=fit2, constant=constant, kernel=kernel))
}
