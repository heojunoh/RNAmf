#' matern.kernel
#'
#' calculating matern kernel with corresponding smoothness parameter
#'
#'
#' @param r vector or matrix of input.
#' @param nu numerical value of smoothness hyperparameter. It should be 0.5, 1.5, 2.5, 3.5, or 4.5.
#' @param derivative logical indicating for its first derivative(derivative=1)
#' @noRd
#' @keywords internal
#' @return A value from matern kernel.

matern.kernel <- function(r,nu,derivative=0){
  if(nu==1/2){ #nu=0.5
    if(derivative == 0) out <- exp(-r)
    if(derivative == 1) out <- -exp(-r)
  }else if(nu==3/2){ #nu=1.5
    if(derivative == 0) out <- (1+r*sqrt(3)) * exp(-r*sqrt(3))
    if(derivative == 1) out <- -3*r*exp(-sqrt(3)*r)
  }else if(nu==5/2){ #nu=2.5
    if(derivative == 0) out <- (1+r*sqrt(5)+5*r^2/3) * exp(-r*sqrt(5))
    if(derivative == 1) out <- -(r*(5^(3/2)*r+5)*exp(-sqrt(5)*r))/3
  }else if(nu==7/2){ #nu=3.5
    if(derivative == 0) out <- (1+r*sqrt(7)+2.8*r^2+7/15*sqrt(7)*r^3) * exp(-r*sqrt(7))
    if(derivative == 1) out <- -(r*(49*r^2+3*7^(3/2)*r+21)*exp(-sqrt(7)*r))/15
  }else if(nu==9/2){ #nu=4.5
    if(derivative == 0) out <- (1+r*sqrt(9)+27*r^2/7+18/7*r^3+27/35*r^4) * exp(-r*3)
    if(derivative == 1) out <- -((81*r^4+162*r^3+135*r^2+45*r)*exp(-3*r))/35
  }

  return(out)
}

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

#' matGP
#'
#' fitting the model with matern kernel.
#'
#' @param X vector or matrix of input locations.
#' @param y vector of response values.
#' @param g nugget parameter. Default is 1.490116e-08.
#' @param nu numerical value of smoothness hyperparameter. It should be 0.5, 1.5, 2.5, 3.5, or 4.5.
#' @param lower lower bound of theta.
#' @param upper upper bound of theta.
#' @param Xscale logical indicating whether to scale X or not. Default is TRUE.
#' @param Yscale logical indicating whether to scale y or not. Only used if constant=FALSE. Default is TRUE.
#' @param constant logical indicating for constant mean (constant=TRUE) or zero mean (constant=FALSE). Default is FALSE.
#'
#' @return A list containing hyperparameters, covariance inverse matrix, X, y and logical inputs:
#' \itemize{
#'   \item \code{theta}: vector of lengthscale hyperparameter.
#'   \item \code{nu}: copy of nu.
#'   \item \code{g}: copy of g.
#'   \item \code{Ki}: matrix of covariance inverse.
#'   \item \code{mu.hat}: optimized constant mean. If constant=FALSE, 0.
#'   \item \code{X}: copy of X. If Xscale=TRUE, scaled X.
#'   \item \code{y}: copy of y. If Yscale=TRUE, scaled y.
#'   \item \code{tau2hat}: estimated scale hyperparameter.
#'   \item \code{Xscale}: copy of Xscale.
#'   \item \code{Yscale}: copy of Yscale.
#'   \item \code{constant}: copy of constant.
#' }
#'
#' @importFrom stats optim quantile
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
#' matGP(X1, y1, nu=2.5)
#' }

matGP <- function(X, y, nu=2.5, g=sqrt(.Machine$double.eps),
                  lower=rep(0.1, ncol(X)), upper=rep(100,ncol(X)),
                  Xscale=TRUE, Yscale=TRUE, constant=FALSE){
  if(constant){

    if(is.null(dim(X))) X <- matrix(X, ncol = 1)

    parbound <- function(XX){
      XX <- matrix(XX, ncol = 1)

      lower = - quantile(distance(XX)[lower.tri(distance(XX))], 0.05)/log(0.01) * (apply(XX, 2, range)[2,] - apply(XX, 2, range)[1,])^2
      upper = - quantile(distance(XX)[lower.tri(distance(XX))], 0.95)/log(0.5) * (apply(XX, 2, range)[2,] - apply(XX, 2, range)[1,])^2

      return(c(lower, upper))
    }

    # # darg way
    # init <- rep(sort(distance(X))[sort(distance(X))!=0][0.1*length(sort(distance(X))[sort(distance(X))!=0])], ncol(X))
    # lower <- min(sort(distance(X))[sort(distance(X))!=0])
    # upper <- max(sort(distance(X))[sort(distance(X))!=0])

    # hetGP way * 10^(d*2/3)
    XX <- (X - matrix(apply(X, 2, range)[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) %*% diag(1/(apply(X, 2, range)[2,] - apply(X, 2, range)[1,]), ncol(X))
    lower <- max(10^(ncol(XX)*2/3)*apply(XX, 2, parbound)[1,], 0.2)
    upper <- 10^(ncol(XX)*2/3)*apply(XX, 2, parbound)[2,]
    init <- sqrt(lower*upper)

    # # hetGP way
    # XX <- (X - matrix(apply(X, 2, range)[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) %*% diag(1/(apply(X, 2, range)[2,] - apply(X, 2, range)[1,]), ncol(X))
    # lower <- - quantile(distance(XX)[lower.tri(distance(XX))], 0.05)/log(0.01) * (apply(XX, 2, range)[2,] - apply(XX, 2, range)[1,])^2
    # upper <- - quantile(distance(XX)[lower.tri(distance(XX))], 0.95)/log(0.5) * (apply(XX, 2, range)[2,] - apply(XX, 2, range)[1,])^2
    # init <- sqrt(lower*upper)

    if(Xscale){
      X <- scale(X, center = TRUE, scale = TRUE)
    }else{
      attr(X,"scaled:center") <- rep(0, ncol(X))
      attr(X,"scaled:scale") <- rep(1, ncol(X))
    }

    n <- length(y)

    nlsep <- function(par, X, Y)
    {
      theta <- par # lengthscale
      K <- cor.sep(X, theta=theta, nu=nu)
      Ki <- solve(K+diag(g,n))
      ldetK <- determinant(K, logarithm=TRUE)$modulus

      one.vec <- matrix(1,ncol=1,nrow=n)
      mu.hat <- drop((t(one.vec)%*%Ki%*%Y)/(t(one.vec)%*%Ki%*%one.vec))

      tau2hat <- drop(t(Y-mu.hat)%*%Ki%*%(Y-mu.hat)/n)
      ll <- - (n/2)*log(tau2hat) - (1/2)*ldetK
      return(-ll)
    }

    # gradnlsep <- function(par, X, Y)
    # {
    #   theta <- par
    #   n <- length(Y)
    #   K <- cor.sep(X, theta=theta, nu=nu)
    #   Ki <- solve(K+diag(g,n))
    #   KiY <- Ki %*% Y
    #
    #   ## loop over theta components
    #   dlltheta <- rep(NA, length(theta))
    #   # dK <- matern.kernel(R, nu=nu, derivative = 1)
    #   dK <- cor.sep(X, theta=theta, nu=nu, derivative = 1)
    #   for(j in 1:length(dlltheta)) {
    #     dotK <- dK * (-distance(X[,j])/(theta[j]^3)/R)
    #     diag(dotK) <- rep(0, n)
    #     dlltheta[j] <- (n/2) * t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) -
    #       (1/2)*sum(diag(Ki %*% dotK))
    #   }
    #
    #   return(-c(dlltheta))
    # }

    outg <- optim(init, nlsep, #gradnlsep,
                  method="L-BFGS-B", lower=lower, upper=upper, X=X, Y=y)

    theta <- outg$par

    # R <- sqrt(distance(t(t(X)/theta)))
    # K <- matern.kernel(R, nu=nu)
    K <- cor.sep(X, theta=theta, nu=nu)
    Ki <- solve(K+diag(g,n))
    one.vec <- matrix(1,ncol=1,nrow=n)
    mu.hat <- drop((t(one.vec)%*%Ki%*%y)/(t(one.vec)%*%Ki%*%one.vec))
    tau2hat <- drop(t(y-mu.hat) %*% Ki %*% (y-mu.hat) / nrow(X))

    return(list(theta = theta, nu=nu, g=g, Ki=Ki, mu.hat=mu.hat, X = X, y = y, tau2hat=tau2hat, Xscale=Xscale, constant=constant))
  }else{

    if(is.null(dim(X))) X <- matrix(X, ncol = 1)

    parbound <- function(XX){
      XX <- matrix(XX, ncol = 1)

      lower = - quantile(distance(XX)[lower.tri(distance(XX))], 0.05)/log(0.01) * (apply(XX, 2, range)[2,] - apply(XX, 2, range)[1,])^2
      upper = - quantile(distance(XX)[lower.tri(distance(XX))], 0.95)/log(0.5) * (apply(XX, 2, range)[2,] - apply(XX, 2, range)[1,])^2

      return(c(lower, upper))
    }

    # # darg way
    # init <- rep(sort(distance(X))[sort(distance(X))!=0][0.1*length(sort(distance(X))[sort(distance(X))!=0])], ncol(X))
    # lower <- min(sort(distance(X))[sort(distance(X))!=0])
    # upper <- max(sort(distance(X))[sort(distance(X))!=0])

    # hetGP way * 10^(d*2/3)
    XX <- (X - matrix(apply(X, 2, range)[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) %*% diag(1/(apply(X, 2, range)[2,] - apply(X, 2, range)[1,]), ncol(X))
    lower <- max(10^(ncol(XX)*2/3)*apply(XX, 2, parbound)[1,], 0.2)
    upper <- 10^(ncol(XX)*2/3)*apply(XX, 2, parbound)[2,]
    init <- sqrt(lower*upper)

    # # hetGP way
    # XX <- (X - matrix(apply(X, 2, range)[1,], nrow = nrow(X), ncol = ncol(X), byrow = TRUE)) %*% diag(1/(apply(X, 2, range)[2,] - apply(X, 2, range)[1,]), ncol(X))
    # lower <- - quantile(distance(XX)[lower.tri(distance(XX))], 0.05)/log(0.01) * (apply(XX, 2, range)[2,] - apply(XX, 2, range)[1,])^2
    # upper <- - quantile(distance(XX)[lower.tri(distance(XX))], 0.95)/log(0.5) * (apply(XX, 2, range)[2,] - apply(XX, 2, range)[1,])^2
    # init <- sqrt(lower*upper)

    if(Xscale){
      X <- scale(X, center = TRUE, scale = TRUE)
    }else{
      attr(X,"scaled:center") <- rep(0, ncol(X))
      attr(X,"scaled:scale") <- rep(1, ncol(X))
    }

    n <- length(y)
    if(Yscale) y <- scale(y, center=TRUE, scale=FALSE) # If use mean, don't scale

    nlsep <- function(par, X, Y)
    {
      theta <- par # lengthscale
      K <- cor.sep(X, theta=theta, nu=nu)
      Ki <- solve(K+diag(g,n))
      ldetK <- determinant(K, logarithm=TRUE)$modulus
      tau2hat <- drop(t(Y)%*%Ki%*%(Y)/n)
      ll <- - (n/2)*log(tau2hat) - (1/2)*ldetK
      return(drop(-ll))
    }

    outg <- optim(init, nlsep, method="L-BFGS-B", lower=lower, upper=upper, X=X, Y=y)

    K <- cor.sep(X, theta=outg$par, nu=nu)
    Ki <- solve(K+diag(g,n))
    tau2hat <- drop(t(y) %*% Ki %*% (y) / n)
    mu.hat <- 0

    return(list(theta = outg$par, nu=nu, g=g, Ki=Ki, mu.hat=mu.hat, X = X, y = y, tau2hat=tau2hat, Xscale=Xscale, Yscale=Yscale, constant=constant))
  }
}

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
