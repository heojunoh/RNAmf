#' Fitting the model with three fidelity levels
#'
#' @description The function fits RNA models with designs of three fidelity levels,
#' including estimation of hyperparameters. The estimation methods are based on MLE.
#' Possible kernel choices are squared exponential, Matern 1.5 and Matern 2.5.
#' The function returns fitted model by \code{\link{RNAmf_two_level}}, fitted model at level 3, whether constant mean or not, and kernel choice.
#'
#' @seealso \code{\link{predRNAmf_three_level}} for prediction.
#'
#' @details Consider the model
#' \eqn{\begin{cases}
#' & f_1(\mathbf{x}) = W_1(\mathbf{x}),\\
#' & f_l(\mathbf{x}) = W_l(\mathbf{x}, f_{l-1}(\mathbf{x})) \quad\text{for}\quad l=2,3,
#' \end{cases}}
#' where \eqn{W(\mathbf{x})} is GP model.
#' Hyperparameters \eqn{(\alpha, \tau^2, \mathbf{\theta})} are estimated by
#' maximizing the log-likelihood via an optimization algorithm "L-BFGS-B".
#' For details, see Heo and Sung (2023+, <arXiv:2309.11772>).
#'
#' @param X1 vector or matrix of input locations for the low fidelity level.
#' @param y1 vector of response values for the low fidelity level.
#' @param X2 vector or matrix of input locations for the medium fidelity level.
#' @param y2 vector of response values for the medium fidelity level.
#' @param X3 vector or matrix of input locations for the high fidelity level.
#' @param y3 vector of response values for the high fidelity level.
#' @param kernel chracter specifying kernel type to be used, to be chosen between "sqex"(squared exponential), "matern1.5", or "matern2.5". Default is "sqex".
#' @param constant logical indicating for constant mean (constant=TRUE) or zero mean (constant=FALSE). Default is FALSE.
#' @return A list containing fitted models f1, f2 and f3, constant, and kernel structure:
#' \itemize{
#'   \item \code{fit.RNAmf_two_level}: list of fitted model for low fidelity level (X1, y1) and medium fidelity level ((X2, f1(X2)), y2).
#'   \item \code{fit3}: list of fitted model for ((X2, f2(X3, f1(X3))), y3).
#'   \item \code{constant}: copy of constant.
#'   \item \code{kernel}: copy of kernel.
#' }
#' @export
#' @examples
#' library(MuFiCokriging)
#' library(lhs)
#'
#' ### Branin function ###
#' branin <- function(xx){
#'   x1 <- xx[1]
#'   x2 <- xx[2]
#'
#'   (-1.275*x1^2/pi^2+5*x1/pi+x2-6)^2 + (10-5/(4*pi))*cos(x1)+ 10
#' }
#'
#' braninm <- function(xx)
#' {
#'   x1 <- xx[1]
#'   x2 <- xx[2]
#'
#'   10*sqrt((-1.275*(x1+2)^2/pi^2+5*(x1+2)/pi+(x2+2)-6)^2 + (10-5/(4*pi))*cos((x1+2))+ 10) +
#'   2*(x1-0.5) - 3*(3*x2-1) - 1
#' }
#'
#' braninl <- function(xx)
#' { x1 <- xx[1]
#' x2 <- xx[2]
#'
#' 10*sqrt((-1.275*(1.2*x1+0.4)^2/pi^2+5*(1.2*x1+0.4)/pi+(1.2*x2+0.4)-6)^2 +
#' (10-5/(4*pi))*cos((1.2*x1+0.4))+ 10) + 2*(1.2*x1+1.9) - 3*(3*(1.2*x2+2.4)-1) - 1 - 3*x2 + 1
#' }
#'
#' output.branin <- function(x){
#'   factor_range <- list("x1" = c(-5, 10), "x2" = c(0, 15))
#'
#'   for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
#'   branin(x[1:2])
#' }
#'
#' output.braninl <- function(x){
#'   factor_range <- list("x1" = c(-5, 10), "x2" = c(0, 15))
#'
#'   for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
#'   braninl(x[1:2])
#' }
#'
#' output.braninm <- function(x){
#'   factor_range <- list("x1" = c(-5, 10), "x2" = c(0, 15))
#'
#'   for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
#'   braninm(x[1:2])
#' }
#'
#'
#' ### training data ###
#' n1 <- 20; n2 <- 15; n3 <- 10
#'
#' X1 <- maximinLHS(n1, 2)
#' X2 <- maximinLHS(n2, 2)
#' X3 <- maximinLHS(n3, 2)
#'
#' NestDesign <- NestedDesignBuild(design = list(X1,X2,X3))
#'
#' X1 <- NestDesign$PX
#' X2 <- ExtractNestDesign(NestDesign,2)
#' X3 <- ExtractNestDesign(NestDesign,3)
#'
#' y1 <- apply(X1,1,output.braninl)
#' y2 <- apply(X2,1,output.braninm)
#' y3 <- apply(X3,1,output.branin)
#'
#' ### test data ###
#' x <- maximinLHS(1000, 2)
#'
#' ### fitting ###
#' fit.RNAmf <- RNAmf_three_level(X1, y1, X2, y2, X3, y3, kernel="sqex", constant=TRUE)
#' pred.RNAmf <- predRNAmf_three_level(fit.RNAmf, x)
#' predy <- pred.RNAmf$mu
#' predsig2 <- pred.RNAmf$sig2
#'
#' ### RMSE ###
#' sqrt(mean((predy-apply(x,1,output.branin))^2))
#'

RNAmf_three_level <- function(X1, y1, X2, y2, X3, y3, kernel="sqex", constant=FALSE){
  if(checknested(X1, X2) == FALSE){
    stop("X2 is not nested by X1")
  }
  if(checknested(X2, X3) == FALSE){
    stop("X3 is not nested by X2")
  }

  if(kernel=="sqex"){
    if(constant){
      fit.RNAmf_two_level <- RNAmf_two_level(X1, y1, X2, y2, kernel="sqex", constant=TRUE)

      fit1 <- fit.RNAmf_two_level$fit1
      fit2 <- fit.RNAmf_two_level$fit2
      fit3 <- GP(cbind(X3, pred.GP(fit2, cbind(X3, pred.GP(fit1, X3)$mu))$mu), y3, constant=TRUE)
    }else{
      fit.RNAmf_two_level <- RNAmf_two_level(X1, y1, X2, y2, kernel="sqex")

      fit1 <- fit.RNAmf_two_level$fit1
      fit2 <- fit.RNAmf_two_level$fit2
      fit3 <- GP(cbind(X3, pred.GP(fit2, cbind(X3, pred.GP(fit1, X3)$mu))$mu), y3)
    }
  }else if(kernel=="matern1.5"){
    if(constant){
      fit.RNAmf_two_level <- RNAmf_two_level(X1, y1, X2, y2, kernel="matern1.5", constant=TRUE)

      fit1 <- fit.RNAmf_two_level$fit1
      fit2 <- fit.RNAmf_two_level$fit2
      fit3 <- matGP(cbind(X3, pred.matGP(fit2, cbind(X3, pred.matGP(fit1, X3)$mu))$mu), y3, nu=1.5, constant=TRUE)
    }else{
      fit.RNAmf_two_level <- RNAmf_two_level(X1, y1, X2, y2, kernel="matern1.5")

      fit1 <- fit.RNAmf_two_level$fit1
      fit2 <- fit.RNAmf_two_level$fit2
      fit3 <- matGP(cbind(X3, pred.matGP(fit2, cbind(X3, pred.matGP(fit1, X3)$mu))$mu), y3, nu=1.5)
    }
  }else if(kernel=="matern2.5"){
    if(constant){
      fit.RNAmf_two_level <- RNAmf_two_level(X1, y1, X2, y2, kernel="matern2.5", constant=TRUE)

      fit1 <- fit.RNAmf_two_level$fit1
      fit2 <- fit.RNAmf_two_level$fit2
      fit3 <- matGP(cbind(X3, pred.matGP(fit2, cbind(X3, pred.matGP(fit1, X3)$mu))$mu), y3, nu=2.5, constant=TRUE)
    }else{
      fit.RNAmf_two_level <- RNAmf_two_level(X1, y1, X2, y2, kernel="matern2.5")

      fit1 <- fit.RNAmf_two_level$fit1
      fit2 <- fit.RNAmf_two_level$fit2
      fit3 <- matGP(cbind(X3, pred.matGP(fit2, cbind(X3, pred.matGP(fit1, X3)$mu))$mu), y3, nu=2.5)
    }
  }

  return(list(fit.RNAmf_two_level=fit.RNAmf_two_level, fit3=fit3, constant=constant, kernel=kernel))
}

