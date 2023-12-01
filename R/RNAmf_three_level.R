#' Fitting the model with three fidelity levels
#'
#' @description The function fits RNA models with designs of three fidelity levels.
#' The estimation method is based on MLE.
#' Possible kernel choices are squared exponential, Matern kernel with smoothness parameter 1.5 and 2.5.
#' The function returns fitted model by \code{\link{RNAmf_two_level}}, fitted model at level 3, whether constant mean or not, and kernel choice.
#'
#' @seealso \code{\link{predict.RNAmf}} for prediction.
#'
#' @details Consider the model
#' \eqn{\begin{cases}
#' & f_1(\bm{x}) = W_1(\bm{x}),\\
#' & f_l(\bm{x}) = W_l(\bm{x}, f_{l-1}(\bm{x})) \quad\text{for}\quad l=2,3,
#' \end{cases}}
#' where \eqn{f_l} is the simulation code at fidelity level \eqn{l}, and
#' \eqn{W_l(\bm{x}) \sim GP(\alpha_l, \tau_l^2 K_l(\bm{x}, \bm{x}'))} is GP model.
#' Hyperparameters \eqn{(\alpha_l, \tau_l^2, \bm{\theta_l})} are estimated by
#' maximizing the log-likelihood via an optimization algorithm "L-BFGS-B".
#' For \code{constant=FALSE}, \eqn{\alpha_l=0}.
#'
#' Covariance kernel is defined as:
#' \eqn{K_l(\bm{x}, \bm{x}')=\prod^d_{j=1}\phi(x_j,x'_j;\theta_{lj})} with
#' \eqn{\phi(x, x';\theta) = \exp \left( -\frac{ \left( x - x' \right)^2}{\theta}  \right)}
#' for squared exponential kernel; \code{kernel="sqex"},
#' \eqn{\phi(x,x';\theta) =\left( 1+\frac{\sqrt{3}|x- x'|}{\theta} \right) \exp \left( -\frac{\sqrt{3}|x- x'|}{\theta} \right)}
#' for Matern kernel with the smoothness parameter of 1.5; \code{kernel="matern1.5"} and
#' \eqn{\phi(x, x';\theta) = \left( 1+\frac{\sqrt{5}|x-x'|}{\theta} +\frac{5(x-x')^2}{3\theta^2} \right) \exp \left( -\frac{\sqrt{5}|x-x'|}{\theta} \right)}
#' for Matern kernel with the smoothness parameter of 2.5; \code{kernel="matern2.5"}.
#'
#' For details, see Heo and Sung (2023+, <arXiv:2309.11772>).
#'
#' @param X1 vector or matrix of input locations for the low fidelity level.
#' @param y1 vector of response values for the low fidelity level.
#' @param X2 vector or matrix of input locations for the medium fidelity level.
#' @param y2 vector of response values for the medium fidelity level.
#' @param X3 vector or matrix of input locations for the high fidelity level.
#' @param y3 vector of response values for the high fidelity level.
#' @param kernel character specifying kernel type to be used, to be chosen between \code{"sqex"}(squared exponential), \code{"matern1.5"}, or \code{"matern2.5"}. Default is \code{"sqex"}.
#' @param constant logical indicating for constant mean of GP (\code{constant=TRUE}) or zero mean (\code{constant=FALSE}). Default is \code{TRUE}.
#' @param ... for compatibility with \code{optim}.
#' @return
#' \itemize{
#'   \item \code{fit.RNAmf_two_level}: a class \code{RNAmf} object fitted by \code{RNAmf_two_level}. It contains a list of \eqn{\begin{cases} & \text{\code{fit1} for } (X_1, y_1),\\ & \text{\code{fit2} for } ((X_2, f_1(X_2)), y_2), \end{cases}}. See \code{\link{RNAmf_two_level}}.
#'   \item \code{fit3}: list of fitted model for \eqn{((X_2, f_2(X_3, f_1(X_3))), y_3)}.
#'   \item \code{constant}: copy of \code{constant}.
#'   \item \code{kernel}: copy of \code{kernel}.
#'   \item \code{level}: a level of the fidelity. It returns 3.
#'   \item \code{time}: a scalar of the time for the computation.
#' }
#' @usage RNAmf_three_level(X1, y1, X2, y2, X3, y3, kernel = "sqex", constant = TRUE, ...)
#' @export
#'

RNAmf_three_level <- function(X1, y1, X2, y2, X3, y3, kernel = "sqex", constant = TRUE, ...) {
  t1 <- proc.time()[3]
  if (checknested(X1, X2) == FALSE) {
    stop("X2 is not nested by X1")
  }
  if (checknested(X2, X3) == FALSE) {
    stop("X3 is not nested by X2")
  }

  if (kernel == "sqex") {
    if (constant) {
      fit.RNAmf_two_level <- RNAmf_two_level(X1, y1, X2, y2, kernel = "sqex", constant = TRUE)

      fit1 <- fit.RNAmf_two_level$fit1
      fit2 <- fit.RNAmf_two_level$fit2
      fit3 <- GP(cbind(X3, pred.GP(fit2, cbind(X3, pred.GP(fit1, X3)$mu))$mu), y3, constant = TRUE)
    } else {
      fit.RNAmf_two_level <- RNAmf_two_level(X1, y1, X2, y2, kernel = "sqex")

      fit1 <- fit.RNAmf_two_level$fit1
      fit2 <- fit.RNAmf_two_level$fit2
      fit3 <- GP(cbind(X3, pred.GP(fit2, cbind(X3, pred.GP(fit1, X3)$mu))$mu), y3)
    }
  } else if (kernel == "matern1.5") {
    if (constant) {
      fit.RNAmf_two_level <- RNAmf_two_level(X1, y1, X2, y2, kernel = "matern1.5", constant = TRUE)

      fit1 <- fit.RNAmf_two_level$fit1
      fit2 <- fit.RNAmf_two_level$fit2
      fit3 <- matGP(cbind(X3, pred.matGP(fit2, cbind(X3, pred.matGP(fit1, X3)$mu))$mu), y3, nu = 1.5, constant = TRUE)
    } else {
      fit.RNAmf_two_level <- RNAmf_two_level(X1, y1, X2, y2, kernel = "matern1.5")

      fit1 <- fit.RNAmf_two_level$fit1
      fit2 <- fit.RNAmf_two_level$fit2
      fit3 <- matGP(cbind(X3, pred.matGP(fit2, cbind(X3, pred.matGP(fit1, X3)$mu))$mu), y3, nu = 1.5)
    }
  } else if (kernel == "matern2.5") {
    if (constant) {
      fit.RNAmf_two_level <- RNAmf_two_level(X1, y1, X2, y2, kernel = "matern2.5", constant = TRUE)

      fit1 <- fit.RNAmf_two_level$fit1
      fit2 <- fit.RNAmf_two_level$fit2
      fit3 <- matGP(cbind(X3, pred.matGP(fit2, cbind(X3, pred.matGP(fit1, X3)$mu))$mu), y3, nu = 2.5, constant = TRUE)
    } else {
      fit.RNAmf_two_level <- RNAmf_two_level(X1, y1, X2, y2, kernel = "matern2.5")

      fit1 <- fit.RNAmf_two_level$fit1
      fit2 <- fit.RNAmf_two_level$fit2
      fit3 <- matGP(cbind(X3, pred.matGP(fit2, cbind(X3, pred.matGP(fit1, X3)$mu))$mu), y3, nu = 2.5)
    }
  }

  model <- list()

  model$fit.RNAmf_two_level <- fit.RNAmf_two_level
  # model$fit1 <- fit.RNAmf_two_level$fit1
  # model$fit2 <- fit.RNAmf_two_level$fit2
  model$fit3 <- fit3
  model$constant <- constant
  model$kernel <- kernel
  model$level <- 3
  model$time <- proc.time()[3] - t1

  class(model) <- "RNAmf"

  return(model)
}
