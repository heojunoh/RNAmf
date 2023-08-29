#' RNAmf2
#'
#' Fitting the model with three fidelity levels
#'
#' @param X1 vector or matrix of input locations for the low fidelity level.
#' @param y1 vector of response values for the low fidelity level.
#' @param X2 vector or matrix of input locations for the medium fidelity level.
#' @param y2 vector of response values for the medium fidelity level.
#' @param X3 vector or matrix of input locations for the high fidelity level.
#' @param y3 vector of response values for the high fidelity level.
#' @param kernel chracter specifying kernel type to be used, to be chosen between "sqex"(squared exponential), "matern1.5", or "matern2.5"
#' @param constant logical indicating for constant mean (constant=TRUE) or zero mean (constant=FALSE). Default is FALSE.
#' @return A list containing fitted models f1, f2 and f3, constant, and kernel structure:
#' \itemize{
#'   \item \code{fit.RNAmf1}: list of fitted model for low fidelity level (X1, y1) and medium fidelity level ((X2, f1(X2)), y2).
#'   \item \code{fit3}: list of fitted model for ((X2, f2(X3, f1(X3))), y3).
#'   \item \code{constant}: copy of constant.
#'   \item \code{kernel}: copy of kernel.
#' }
#' @export
#'

RNAmf2 <- function(X1, y1, X2, y2, X3, y3, kernel, constant=FALSE){
  # if(all(X2 %in% X1) == FALSE){
  #   stop("X2 is not nested by X1")
  # }
  #
  # if(all(X3 %in% X2) == FALSE){
  #   stop("X2 is not nested by X1")
  # }

  if(kernel=="sqex"){
    if(constant){
      fit.RNAmf1 <- RNAmf(X1, y1, X2, y2, kernel="sqex", constant=TRUE)

      fit1 <- fit.RNAmf1$fit1
      fit2 <- fit.RNAmf1$fit2
      fit3 <- GP(cbind(X3, pred.GP(fit2, cbind(X3, pred.GP(fit1, X3)$mu))$mu), y3, constant=TRUE)
    }else{
      fit.RNAmf1 <- RNAmf(X1, y1, X2, y2, kernel="sqex")

      fit1 <- fit.RNAmf1$fit1
      fit2 <- fit.RNAmf1$fit2
      fit3 <- GP(cbind(X3, pred.GP(fit2, cbind(X3, pred.GP(fit1, X3)$mu))$mu), y3)
    }
  }else if(kernel=="matern1.5"){
    if(constant){
      fit.RNAmf1 <- RNAmf(X1, y1, X2, y2, kernel="matern1.5", constant=TRUE)

      fit1 <- fit.RNAmf1$fit1
      fit2 <- fit.RNAmf1$fit2
      fit3 <- matGP(cbind(X3, pred.matGP(fit2, cbind(X3, pred.matGP(fit1, X3)$mu))$mu), y3, nu=1.5, constant=TRUE)
    }else{
      fit.RNAmf1 <- RNAmf(X1, y1, X2, y2, kernel="matern1.5")

      fit1 <- fit.RNAmf1$fit1
      fit2 <- fit.RNAmf1$fit2
      fit3 <- matGP(cbind(X3, pred.matGP(fit2, cbind(X3, pred.matGP(fit1, X3)$mu))$mu), y3, nu=1.5)
    }
  }else if(kernel=="matern2.5"){
    if(constant){
      fit.RNAmf1 <- RNAmf(X1, y1, X2, y2, kernel="matern2.5", constant=TRUE)

      fit1 <- fit.RNAmf1$fit1
      fit2 <- fit.RNAmf1$fit2
      fit3 <- matGP(cbind(X3, pred.matGP(fit2, cbind(X3, pred.matGP(fit1, X3)$mu))$mu), y3, nu=2.5, constant=TRUE)
    }else{
      fit.RNAmf1 <- RNAmf(X1, y1, X2, y2, kernel="matern2.5")

      fit1 <- fit.RNAmf1$fit1
      fit2 <- fit.RNAmf1$fit2
      fit3 <- matGP(cbind(X3, pred.matGP(fit2, cbind(X3, pred.matGP(fit1, X3)$mu))$mu), y3, nu=2.5)
    }
  }

  return(list(fit.RNAmf1=fit.RNAmf1, fit3=fit3, constant=constant, kernel=kernel))
}

