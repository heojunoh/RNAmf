#' pred.matGP
#'
#' predictive posterior mean and variance with matern kernel.
#'
#' @param fit an object of class matGP.
#' @param xnew vector or matrix of new input locations to predict.
#'
#' @return A list predictive posterior mean and variance:
#' \itemize{
#'   \item \code{mu}: vector of predictive posterior mean.
#'   \item \code{sig2}: vector of predictive posterior variance.
#' }
#'
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
#' fit1 <- matGP(X1, y1, nu=2.5)
#'
#' ### test data ###
#' x <- seq(0,1,0.01)
#' pred.matGP(fit1, x)
#' }

pred.matGP <- function(fit, xnew){
  constant <- fit$constant

  if(constant){
    xnew <- as.matrix(xnew)

    Xscale <- fit$Xscale
    Ki <- fit$Ki
    theta <- fit$theta
    nu <- fit$nu
    g <- fit$g
    X <- fit$X
    y <- fit$y
    tau2hat <- fit$tau2hat
    mu.hat <- fit$mu.hat

    if(Xscale) xnew <- t((t(xnew)-attr(X,"scaled:center"))/attr(X,"scaled:scale"))

    KXX <- cor.sep(xnew, theta=theta, nu=nu)
    KX <- t(cor.sep(X, xnew, theta=theta, nu=nu))

    mup2 <- mu.hat + KX %*% Ki %*% (y - mu.hat)
    Sigmap2 <- pmax(0, diag(tau2hat*(KXX + diag(g,nrow(xnew)) - KX %*% Ki %*% t(KX))))

    return(list(mu=mup2, sig2=Sigmap2))
  }else{
    xnew <- as.matrix(xnew)

    Xscale <- fit$Xscale
    Yscale <- fit$Yscale
    Ki <- fit$Ki
    theta <- fit$theta
    nu <- fit$nu
    g <- fit$g
    X <- fit$X
    y <- fit$y
    tau2hat <- fit$tau2hat

    if(Xscale) xnew <- t((t(xnew)-attr(X,"scaled:center"))/attr(X,"scaled:scale"))

    KXX <- cor.sep(xnew, theta=theta, nu=nu)
    KX <- t(cor.sep(X, xnew, theta=theta, nu=nu))

    if(Yscale) mup2 <- KX %*% Ki %*% (y + attr(y, "scaled:center")) else mup2 <- KX %*% Ki %*% y
    Sigmap2 <- pmax(0, diag(tau2hat*(KXX + diag(g,nrow(xnew)) - KX %*% Ki %*% t(KX))))

    return(list(mu=mup2, sig2=Sigmap2))
  }
}
