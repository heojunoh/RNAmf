#' pred.GP
#'
#' predictive posterior mean and variance with squared exponential kernel.
#'
#' @param fit an object of class GP.
#' @param xnew vector or matrix of new input locations to predict.
#'
#' @return A list predictive posterior mean and variance:
#' \itemize{
#'   \item \code{mu}: vector of predictive posterior mean.
#'   \item \code{sig2}: vector of predictive posterior variance.
#' }
#'
#' @importFrom plgp covar.sep
#' @noRd
#' @keywords internal
#' @examples
#' \dontrun{
#' library(lhs)
#' ### synthetic function ###
#' f1 <- function(x)
#' {
#'   sin(8*pi*x)
#' }
#'
#' ### training data ###
#' n1 <- 15
#'
#' X1 <- maximinLHS(n1, 1)
#' y1 <- f1(X1)
#'
#' fit1 <- GP(X1, y1)
#'
#' ### test data ###
#' x <- seq(0,1,0.01)
#' pred.GP(fit1, x)
#' }

pred.GP <- function(fit, xnew){
  constant <- fit$constant

  if(constant){
    xnew <- as.matrix(xnew)

    Xscale <- fit$Xscale
    Ki <- fit$Ki
    theta <- fit$theta
    g <- fit$g
    X <- fit$X
    y <- fit$y
    tau2hat <- fit$tau2hat
    mu.hat <- fit$mu.hat

    if(Xscale) xnew <- t((t(xnew)-attr(X,"scaled:center"))/attr(X,"scaled:scale"))

    KXX <- covar.sep(xnew, d=theta, g=g)
    KX <- covar.sep(xnew, X, d=theta, g=0)

    mup2 <- mu.hat + KX %*% Ki %*% (y - mu.hat)
    Sigmap2 <- pmax(0, diag(tau2hat*(KXX - KX %*% Ki %*% t(KX))))

    return(list(mu=mup2, sig2=Sigmap2))
  }else{
    xnew <- as.matrix(xnew)

    Xscale <- fit$Xscale
    Yscale <- fit$Yscale
    Ki <- fit$Ki
    theta <- fit$theta
    g <- fit$g
    X <- fit$X
    y <- fit$y
    tau2hat <- fit$tau2hat

    if(Xscale) xnew <- t((t(xnew)-attr(X,"scaled:center"))/attr(X,"scaled:scale"))

    KXX <- covar.sep(xnew, d=theta, g=g)
    KX <- covar.sep(xnew, X, d=theta, g=0)

    if(Yscale) mup2 <- KX %*% Ki %*% (y + attr(y, "scaled:center")) else mup2 <- KX %*% Ki %*% y
    Sigmap2 <- pmax(0, diag(tau2hat*(KXX - KX %*% Ki %*% t(KX))))

    return(list(mu=mup2, sig2=Sigmap2))
  }
}
