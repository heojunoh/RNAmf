#' predictive posterior mean and variance of the model with fidelity level 1 and 2.
#'
#' @description The function computes the posterior mean and variance of RNA models with two fidelity levels
#' by fitted model using \code{\link{RNAmf_two_level}}.
#'
#' @seealso \code{\link{RNAmf_two_level}} for prediction.
#'
#' @details From the model fitted by \code{\link{RNAmf_two_level}},
#' the posterior mean and variance are calculated based on the closed form expression derived by a recursive fashion.
#' The formulae depend on its kernel choices.
#' For Matern kernels, \code{\link{zetafun}} computes the part of the posterior variance.
#' For details, see Heo and Sung (2023+, <arXiv:2309.11772>).
#'
#' @param fit an object of class RNAmf.
#' @param x vector or matrix of new input locations to predict.
#'
#' @return A list predictive posterior mean and variance:
#' \itemize{
#'   \item \code{mu}: vector of predictive posterior mean.
#'   \item \code{sig2}: vector of predictive posterior variance.
#' }
#'
#' #' @importFrom plgp distance
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

predRNAmf_two_level <- function(fit, x){
  kernel <- fit$kernel
  constant <- fit$constant
  fit1 <- fit$fit1
  fit2 <- fit$fit2

  if(kernel=="sqex"){
    if(constant){
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      pred.fit <- pred.GP(fit1, x)
      x.mu <- pred.fit$mu # mean of f1(u)
      sig2 <- pred.fit$sig2 #*0

      ### calculate the closed form ###
      X2 <- matrix(fit2$X[,-(d+1)], ncol=d)
      w1.x2 <- fit2$X[,d+1]
      y2 <- fit2$y
      n <- length(y2)
      theta <- fit2$theta
      tau2hat <- fit2$tau2hat
      mu2 <- fit2$mu.hat

      Ci <- fit2$Ki
      a <- Ci %*% (y2 - mu2)

      ### scale new inputs ###
      x <- scale.inputs(x, attr(fit2$X,"scaled:center")[1:d], attr(fit2$X,"scaled:scale")[1:d])
      x.mu <- scale.inputs(x.mu, attr(fit2$X,"scaled:center")[d+1], attr(fit2$X,"scaled:scale")[d+1])
      sig2 <- sig2/attr(fit2$X,"scaled:scale")[d+1]^2

      # mean
      predy <- mu2 + (exp(-distance(t(t(x)/sqrt(theta[-(d+1)])), t(t(X2)/sqrt(theta[-(d+1)])))) *
                        1/sqrt(1+2*sig2/theta[d+1]) *
                        exp(-(drop(outer(x.mu, w1.x2, FUN="-")))^2/(theta[d+1]+2*sig2))) %*% a

      # var
      predsig2 <- c(rep(0, nrow(x)))
      for(i in 1: nrow(x)){ # each test point
        mat <- drop(exp(-distance(t(x[i,])/sqrt(theta[-(d+1)]), t(t(X2)/sqrt(theta[-(d+1)])))) %o%
                      exp(-distance(t(x[i,])/sqrt(theta[-(d+1)]), t(t(X2)/sqrt(theta[-(d+1)]))))) * # common components
          1/sqrt(1+4*sig2[i]/theta[d+1]) *
          exp(-(outer(w1.x2, w1.x2, FUN="+")/2 - x.mu[i])^2/(theta[d+1]/2+2*sig2[i])) *
          exp(-(outer(w1.x2, w1.x2, FUN="-"))^2/(2*theta[d+1]))

        predsig2[i] <- pmax(0, tau2hat - (predy[i]-mu2)^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }
    }else{
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      pred.fit <- pred.GP(fit1, x)
      x.mu <- pred.fit$mu # mean of f1(u)
      sig2 <- pred.fit$sig2

      ### calculate the closed form ###
      X2 <- matrix(fit2$X[,-(d+1)], ncol=d)
      w1.x2 <- fit2$X[,d+1]
      y2 <- fit2$y
      n <- length(y2)
      theta <- fit2$theta
      tau2hat <- fit2$tau2hat

      Ci <- fit2$Ki
      a <- Ci %*% (y2 + attr(y2, "scaled:center"))

      ### scale new inputs ###
      x <- scale.inputs(x, attr(fit2$X,"scaled:center")[1:d], attr(fit2$X,"scaled:scale")[1:d])
      x.mu <- scale.inputs(x.mu, attr(fit2$X,"scaled:center")[d+1], attr(fit2$X,"scaled:scale")[d+1])
      sig2 <- sig2/attr(fit2$X,"scaled:scale")[d+1]^2

      # mean
      predy <- (exp(-distance(t(t(x)/sqrt(theta[-(d+1)])), t(t(X2)/sqrt(theta[-(d+1)])))) *
                        1/sqrt(1+2*sig2/theta[d+1]) *
                        exp(-(drop(outer(x.mu, w1.x2, FUN="-")))^2/(theta[d+1]+2*sig2))) %*% a

      # var
      predsig2 <- c(rep(0, nrow(x)))
      for(i in 1: nrow(x)){ # each test point
        mat <- drop(exp(-distance(t(x[i,])/sqrt(theta[-(d+1)]), t(t(X2)/sqrt(theta[-(d+1)])))) %o%
                      exp(-distance(t(x[i,])/sqrt(theta[-(d+1)]), t(t(X2)/sqrt(theta[-(d+1)]))))) * # common components
          1/sqrt(1+4*sig2[i]/theta[d+1]) *
          exp(-(outer(w1.x2, w1.x2, FUN="+")/2 - x.mu[i])^2/(theta[d+1]/2+2*sig2[i])) *
          exp(-(outer(w1.x2, w1.x2, FUN="-"))^2/(2*theta[d+1]))

        predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }
    }
  }else if(kernel=="matern1.5"){
    if(constant){
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      pred.fit <- pred.matGP(fit1, x)
      x.mu <- pred.fit$mu # mean of f1(u)
      sig2 <- pred.fit$sig2 #*0

      ### calculate the closed form ###
      X2 <- matrix(fit2$X[,-(d+1)], ncol=d)
      w1.x2 <- fit2$X[,d+1]
      y2 <- fit2$y
      n <- length(y2)
      theta <- fit2$theta
      tau2hat <- fit2$tau2hat
      mu2 <- fit2$mu.hat

      Ci <- fit2$Ki
      a <- Ci %*% (y2 - mu2)

      ### scale new inputs ###
      x <- scale.inputs(x, attr(fit2$X,"scaled:center")[1:d], attr(fit2$X,"scaled:scale")[1:d])
      x.mu <- scale.inputs(x.mu, attr(fit2$X,"scaled:center")[d+1], attr(fit2$X,"scaled:scale")[d+1])
      sig2 <- sig2/attr(fit2$X,"scaled:scale")[d+1]^2

      # mean
      predy <- c(rep(0, nrow(x)))
      for(j in 1: nrow(x)){ # each test point
        mua <- x.mu[j] - sqrt(3)*sig2[j]/theta[d+1]
        mub <- x.mu[j] + sqrt(3)*sig2[j]/theta[d+1]

        lambda11 <- c(1, mua)
        lambda12 <- c(0, 1)
        lambda21 <- c(1, -mub)

        e1 <- cbind(matrix(1-sqrt(3)*w1.x2/theta[d+1]), sqrt(3)/theta[d+1])
        e2 <- cbind(matrix(1+sqrt(3)*w1.x2/theta[d+1]), sqrt(3)/theta[d+1])

        predy[j] <- mu2 + drop(t(t(cor.sep(t(x[j,]), X2, theta[-(d+1)], nu=1.5)) * # common but depends on kernel
                                   (exp((3*sig2[j] + 2*sqrt(3)*theta[d+1]*(w1.x2 - x.mu[j]))/(2*theta[d+1]^2)) *
                                      (e1 %*% lambda11 * pnorm((mua - w1.x2)/sqrt(sig2[j])) +
                                         e1 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2 - mua)^2/(2*sig2[j]))) +
                                      exp((3*sig2[j] - 2*sqrt(3)*theta[d+1]*(w1.x2 - x.mu[j]))/(2*theta[d+1]^2)) *
                                      (e2 %*% lambda21 * pnorm((-mub + w1.x2)/sqrt(sig2[j])) +
                                         e2 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2 - mub)^2/(2*sig2[j]))))) %*% a)
      }

      # var
      predsig2 <- c(rep(0, nrow(x)))
      for(i in 1: nrow(x)){
        zeta <- function(x,y){zetafun(w1=x, w2=y, m=x.mu[i], s=sig2[i], nu=1.5, theta=theta[d+1])}

        mat <- drop(t(cor.sep(t(x[i,]), X2, theta[-(d+1)], nu=1.5)) %o% t(cor.sep(t(x[i,]), X2, theta[-(d+1)], nu=1.5))) * # constant depends on kernel
          outer(w1.x2, w1.x2, FUN=Vectorize(zeta))

        predsig2[i] <- pmax(0, tau2hat - (predy[i]-mu2)^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }
    }else{
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      pred.fit <- pred.matGP(fit1, x)
      x.mu <- pred.fit$mu # mean of f1(u)
      sig2 <- pred.fit$sig2

      ### calculate the closed form ###
      X2 <- matrix(fit2$X[,-(d+1)], ncol=d)
      w1.x2 <- fit2$X[,d+1]
      y2 <- fit2$y
      n <- length(y2)
      theta <- fit2$theta
      tau2hat <- fit2$tau2hat

      Ci <- fit2$Ki
      a <- Ci %*% (y2 + attr(y2, "scaled:center"))

      ### scale new inputs ###
      x <- scale.inputs(x, attr(fit2$X,"scaled:center")[1:d], attr(fit2$X,"scaled:scale")[1:d])
      x.mu <- scale.inputs(x.mu, attr(fit2$X,"scaled:center")[d+1], attr(fit2$X,"scaled:scale")[d+1])
      sig2 <- sig2/attr(fit2$X,"scaled:scale")[d+1]^2

      # mean
      predy <- c(rep(0, nrow(x)))
      for(j in 1: nrow(x)){ # each test point
        mua <- x.mu[j] - sqrt(3)*sig2[j]/theta[d+1]
        mub <- x.mu[j] + sqrt(3)*sig2[j]/theta[d+1]

        lambda11 <- c(1, mua)
        lambda12 <- c(0, 1)
        lambda21 <- c(1, -mub)

        e1 <- cbind(matrix(1-sqrt(3)*w1.x2/theta[d+1]), sqrt(3)/theta[d+1])
        e2 <- cbind(matrix(1+sqrt(3)*w1.x2/theta[d+1]), sqrt(3)/theta[d+1])

        predy[j] <- drop(t(t(cor.sep(t(x[j,]), X2, theta[-(d+1)], nu=1.5)) * # common but depends on kernel
                                   (exp((3*sig2[j] + 2*sqrt(3)*theta[d+1]*(w1.x2 - x.mu[j]))/(2*theta[d+1]^2)) *
                                      (e1 %*% lambda11 * pnorm((mua - w1.x2)/sqrt(sig2[j])) +
                                         e1 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2 - mua)^2/(2*sig2[j]))) +
                                      exp((3*sig2[j] - 2*sqrt(3)*theta[d+1]*(w1.x2 - x.mu[j]))/(2*theta[d+1]^2)) *
                                      (e2 %*% lambda21 * pnorm((-mub + w1.x2)/sqrt(sig2[j])) +
                                         e2 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2 - mub)^2/(2*sig2[j]))))) %*% a)
      }

      # var
      predsig2 <- c(rep(0, nrow(x)))
      for(i in 1: nrow(x)){
        zeta <- function(x,y){zetafun(w1=x, w2=y, m=x.mu[i], s=sig2[i], nu=1.5, theta=theta[d+1])}

        mat <- drop(t(cor.sep(t(x[i,]), X2, theta[-(d+1)], nu=1.5)) %o% t(cor.sep(t(x[i,]), X2, theta[-(d+1)], nu=1.5))) * # constant depends on kernel
          outer(w1.x2, w1.x2, FUN=Vectorize(zeta))

        predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }
    }
  }else if(kernel=="matern2.5"){
    if(constant){
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      pred.fit <- pred.matGP(fit1, x)
      x.mu <- pred.fit$mu # mean of f1(u)
      sig2 <- pred.fit$sig2 #*0

      ### calculate the closed form ###
      X2 <- matrix(fit2$X[,-(d+1)], ncol=d)
      w1.x2 <- fit2$X[,d+1]
      y2 <- fit2$y
      n <- length(y2)
      theta <- fit2$theta
      tau2hat <- fit2$tau2hat
      mu2 <- fit2$mu.hat

      Ci <- fit2$Ki
      a <- Ci %*% (y2 - mu2)

      ### scale new inputs ###
      x <- scale.inputs(x, attr(fit2$X,"scaled:center")[1:d], attr(fit2$X,"scaled:scale")[1:d])
      x.mu <- scale.inputs(x.mu, attr(fit2$X,"scaled:center")[d+1], attr(fit2$X,"scaled:scale")[d+1])
      sig2 <- sig2/attr(fit2$X,"scaled:scale")[d+1]^2

      # mean
      predy <- c(rep(0, nrow(x)))
      for(j in 1: nrow(x)){ # each test point
        mua <- x.mu[j] - sqrt(5)*sig2[j]/theta[d+1]
        mub <- x.mu[j] + sqrt(5)*sig2[j]/theta[d+1]

        lambda11 <- c(1, mua, mua^2+sig2[j])
        lambda12 <- cbind(0, 1, matrix(mua+w1.x2))
        lambda21 <- c(1, -mub, mub^2+sig2[j])
        lambda22 <- cbind(0, 1, matrix(-mub-w1.x2))

        e1 <- cbind(matrix(1-sqrt(5)*w1.x2/theta[d+1]+5*w1.x2^2/(3*theta[d+1]^2)),
                    matrix(sqrt(5)/theta[d+1]-10*w1.x2/(3*theta[d+1]^2)),
                    5/(3*theta[d+1]^2))
        e2 <- cbind(matrix(1+sqrt(5)*w1.x2/theta[d+1]+5*w1.x2^2/(3*theta[d+1]^2)),
                    matrix(sqrt(5)/theta[d+1]+10*w1.x2/(3*theta[d+1]^2)),
                    5/(3*theta[d+1]^2))

        predy[j] <- mu2 + drop(t(t(cor.sep(t(x[j,]), X2, theta[-(d+1)], nu=2.5)) *
                                   (exp((5*sig2[j] + 2*sqrt(5)*theta[d+1]*(w1.x2 - x.mu[j]))/(2*theta[d+1]^2)) *
                                      (e1 %*% lambda11 * pnorm((mua - w1.x2)/sqrt(sig2[j])) +
                                         rowSums(e1 * lambda12) * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2 - mua)^2/(2*sig2[j]))) +
                                      exp((5*sig2[j] - 2*sqrt(5)*theta[d+1]*(w1.x2 - x.mu[j]))/(2*theta[d+1]^2)) *
                                      (e2 %*% lambda21 * pnorm((-mub + w1.x2)/sqrt(sig2[j])) +
                                         rowSums(e2 * lambda22) * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2 - mub)^2/(2*sig2[j]))))) %*% a)
      }

      # var
      predsig2 <- c(rep(0, nrow(x)))
      for(i in 1: nrow(x)){
        zeta <- function(x,y){zetafun(w1=x, w2=y, m=x.mu[i], s=sig2[i], nu=2.5, theta=theta[d+1])}

        mat <- drop(t(cor.sep(t(x[i,]), X2, theta[-(d+1)], nu=2.5)) %o% t(cor.sep(t(x[i,]), X2, theta[-(d+1)], nu=2.5))) * # constant depends on kernel
          outer(w1.x2, w1.x2, FUN=Vectorize(zeta))

        predsig2[i] <- pmax(0, tau2hat - (predy[i]-mu2)^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }
    }else{
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      pred.fit <- pred.matGP(fit1, x)
      x.mu <- pred.fit$mu # mean of f1(u)
      sig2 <- pred.fit$sig2

      ### calculate the closed form ###
      X2 <- matrix(fit2$X[,-(d+1)], ncol=d)
      w1.x2 <- fit2$X[,d+1]
      y2 <- fit2$y
      n <- length(y2)
      theta <- fit2$theta
      tau2hat <- fit2$tau2hat

      Ci <- fit2$Ki
      a <- Ci %*% (y2 + attr(y2, "scaled:center"))

      ### scale new inputs ###
      x <- scale.inputs(x, attr(fit2$X,"scaled:center")[1:d], attr(fit2$X,"scaled:scale")[1:d])
      x.mu <- scale.inputs(x.mu, attr(fit2$X,"scaled:center")[d+1], attr(fit2$X,"scaled:scale")[d+1])
      sig2 <- sig2/attr(fit2$X,"scaled:scale")[d+1]^2

      # mean
      predy <- c(rep(0, nrow(x)))

      for(j in 1: nrow(x)){ # each test point
        mua <- x.mu[j] - sqrt(5)*sig2[j]/theta[d+1]
        mub <- x.mu[j] + sqrt(5)*sig2[j]/theta[d+1]

        lambda11 <- c(1, mua, mua^2+sig2[j])
        lambda12 <- cbind(0, 1, matrix(mua+w1.x2))
        lambda21 <- c(1, -mub, mub^2+sig2[j])
        lambda22 <- cbind(0, 1, matrix(-mub-w1.x2))

        e1 <- cbind(matrix(1-sqrt(5)*w1.x2/theta[d+1]+5*w1.x2^2/(3*theta[d+1]^2)),
                    matrix(sqrt(5)/theta[d+1]-10*w1.x2/(3*theta[d+1]^2)),
                    5/(3*theta[d+1]^2))
        e2 <- cbind(matrix(1+sqrt(5)*w1.x2/theta[d+1]+5*w1.x2^2/(3*theta[d+1]^2)),
                    matrix(sqrt(5)/theta[d+1]+10*w1.x2/(3*theta[d+1]^2)),
                    5/(3*theta[d+1]^2))

        predy[j] <- drop(t(t(cor.sep(t(x[j,]), X2, theta[-(d+1)], nu=2.5)) *
                                   (exp((5*sig2[j] + 2*sqrt(5)*theta[d+1]*(w1.x2 - x.mu[j]))/(2*theta[d+1]^2)) *
                                      (e1 %*% lambda11 * pnorm((mua - w1.x2)/sqrt(sig2[j])) +
                                         rowSums(e1 * lambda12) * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2 - mua)^2/(2*sig2[j]))) +
                                      exp((5*sig2[j] - 2*sqrt(5)*theta[d+1]*(w1.x2 - x.mu[j]))/(2*theta[d+1]^2)) *
                                      (e2 %*% lambda21 * pnorm((-mub + w1.x2)/sqrt(sig2[j])) +
                                         rowSums(e2 * lambda22) * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2 - mub)^2/(2*sig2[j]))))) %*% a)
      }

      # var
      predsig2 <- c(rep(0, nrow(x)))
      for(i in 1: nrow(x)){
        zeta <- function(x,y){zetafun(w1=x, w2=y, m=x.mu[i], s=sig2[i], nu=2.5, theta=theta[d+1])}

        mat <- drop(t(cor.sep(t(x[i,]), X2, theta[-(d+1)], nu=2.5)) %o% t(cor.sep(t(x[i,]), X2, theta[-(d+1)], nu=2.5))) * # constant depends on kernel
          outer(w1.x2, w1.x2, FUN=Vectorize(zeta))

        predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }
    }
  }

  return(list(mu=predy, sig2=predsig2))
}

