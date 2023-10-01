#' Fitting the Recursive non-additive model with two fidelity levels.
#'
#' @description The function fits RNA models with designs of two fidelity levels,
#' including estimation of hyperparameters. The estimation methods are based on MLE.
#' Possible kernel choices are squared exponential, Matern 1.5 and Matern 2.5.
#' The function returns fitted model at level 1 and 2, whether constant mean or not, and kernel choice.
#'
#' @seealso \code{\link{predRNAmf_two_level}} for prediction.
#'
#' @details Consider the model
#' \eqn{\begin{cases}
#' & f_1(\mathbf{x}) = W_1(\mathbf{x}),\\
#' & f_2(\mathbf{x}) = W_2(\mathbf{x}, f_1(\mathbf{x})),
#' \end{cases}}
#' where \eqn{W(\mathbf{x})} is GP model.
#' Hyperparameters \eqn{(\alpha, \tau^2, \mathbf{\theta})} are estimated by
#' maximizing the log-likelihood via an optimization algorithm "L-BFGS-B".
#' For details, see Heo and Sung (2023+, <arXiv:2309.11772>).
#'
#' @param X1 vector or matrix of input locations for the low fidelity level.
#' @param y1 vector of response values for the low fidelity level.
#' @param X2 vector or matrix of input locations for the high fidelity level.
#' @param y2 vector of response values for the high fidelity level.
#' @param kernel chracter specifying kernel type to be used, to be chosen between "sqex"(squared exponential), "matern1.5", or "matern2.5". Default is "sqex".
#' @param constant logical indicating for constant mean (constant=TRUE) or zero mean (constant=FALSE). Default is FALSE.
#' @return A list containing fitted models f1 and f2, constant, and kernel structure:
#' \itemize{
#'   \item \code{fit1}: list of fitted model for (X1, y1).
#'   \item \code{fit2}: list of fitted model for ((X2, f1(X2)), y2).
#'   \item \code{constant}: copy of constant.
#'   \item \code{kernel}: copy of kernel.
#' }
#' @export
#' @examples
#' library(MuFiCokriging)
#' library(lhs)
#'
#' ### Perdikaris function ###
#' f1 <- function(x)
#' {
#'   sin(8*pi*x)
#' }
#'
#' f2 <- function(x)
#' {
#'   (x-sqrt(2))*(sin(8*pi*x))^2
#' }
#'
#' ### training data ###
#' n1 <- 13; n2 <- 8
#'
#' X1 <- maximinLHS(n1, 1)
#' X2 <- maximinLHS(n2, 1)
#'
#' NestDesign <- NestedDesignBuild(design = list(X1,X2))
#'
#' X1 <- NestDesign$PX
#' X2 <- ExtractNestDesign(NestDesign,2)
#'
#' y1 <- f1(X1)
#' y2 <- f2(X2)
#'
#' ### test data ###
#' x <- seq(0,1,length.out=1000)
#'
#' ### fitting ###
#' fit.RNAmf <- RNAmf_two_level(X1, y1, X2, y2, kernel="sqex", constant=TRUE)
#' pred.RNAmf <- predRNAmf_two_level(fit.RNAmf, x)
#' predy <- pred.RNAmf$mu
#' predsig2 <- pred.RNAmf$sig2
#'
#' ### RMSE ###
#' sqrt(mean((predy-f2(x))^2))
#'

RNAmf_two_level <- function(X1, y1, X2, y2, kernel="sqex", constant=FALSE){
  if(checknested(X1, X2) == FALSE){
    stop("X2 is not nested by X1")
  }

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
