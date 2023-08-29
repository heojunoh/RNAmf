#' predRNAmf2
#'
#' predictive posterior mean and variance of the model with fidelity level 1, 2 and 3.
#'
#' @param fit an object of class RNAmf2.
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
#'

predRNAmf2 <- function(fit, x){
  kernel <- fit$kernel
  constant <- fit$constant
  fit.RNAmf1 <- fit$fit.RNAmf1
  fit1 <- fit.RNAmf1$fit1
  fit2 <- fit.RNAmf1$fit2
  fit3 <- fit$fit3

  if(kernel=="sqex"){
    if(constant){
      pred.RNAmf1 <- predRNAmf(fit.RNAmf1, x)
      x.mu <- pred.RNAmf1$mu
      sig2 <- pred.RNAmf1$sig2

      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)

      ### calculate the closed form ###
      X3 <- matrix(fit3$X[,-(d+1)], ncol=d)
      w2.x3 <- fit3$X[,d+1]
      y3 <- fit3$y
      n <- length(y3)
      theta <- fit3$theta
      tau2hat <- fit3$tau2hat
      mu3 <- fit3$mu.hat

      Ci <- fit3$Ki
      a <- Ci %*% (y3 - mu3)

      ### scale new inputs ###
      x <- t((t(x)-attr(fit3$X,"scaled:center")[1:d])/attr(fit3$X,"scaled:scale")[1:d])
      x.mu <- t((t(x.mu)-attr(fit3$X,"scaled:center")[d+1])/attr(fit3$X,"scaled:scale")[d+1])
      sig2 <- sig2/attr(fit3$X,"scaled:scale")[d+1]^2

      # mean
      predy <- mu3 + (exp(-distance(t(t(x)/sqrt(theta[-(d+1)])), t(t(X3)/sqrt(theta[-(d+1)])))) *
                        1/sqrt(1+2*sig2/theta[d+1]) *
                        exp(-(drop(outer(x.mu, w2.x3, FUN="-")))^2/(theta[d+1]+2*sig2))) %*% a

      # var
      predsig2 <- c(rep(0, nrow(x)))
      for(i in 1: nrow(x)){ # each test point
        mat <- drop(exp(-distance(t(x[i,])/sqrt(theta[-(d+1)]), t(t(X3)/sqrt(theta[-(d+1)])))) %o%
                      exp(-distance(t(x[i,])/sqrt(theta[-(d+1)]), t(t(X3)/sqrt(theta[-(d+1)]))))) * # common components
          1/sqrt(1+4*sig2[i]/theta[d+1]) *
          exp(-(outer(w2.x3, w2.x3, FUN="+")/2 - x.mu[i])^2/(theta[d+1]/2+2*sig2[i])) *
          exp(-(outer(w2.x3, w2.x3, FUN="-"))^2/(2*theta[d+1]))

        predsig2[i] <- pmax(0, tau2hat - (predy[i]-mu3)^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }
    }else{
      pred.RNAmf1 <- predRNAmf(fit.RNAmf1, x)
      x.mu <- pred.RNAmf1$mu
      sig2 <- pred.RNAmf1$sig2

      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)

      ### calculate the closed form ###
      X3 <- matrix(fit3$X[,-(d+1)], ncol=d)
      w2.x3 <- fit3$X[,d+1]
      y3 <- fit3$y
      n <- length(y3)
      theta <- fit3$theta
      tau2hat <- fit3$tau2hat

      Ci <- fit3$Ki
      a <- Ci %*% (y3 + attr(y3, "scaled:center"))

      ### scale new inputs ###
      x <- t((t(x)-attr(fit3$X,"scaled:center")[1:d])/attr(fit3$X,"scaled:scale")[1:d])
      x.mu <- t((t(x.mu)-attr(fit3$X,"scaled:center")[d+1])/attr(fit3$X,"scaled:scale")[d+1])
      sig2 <- sig2/attr(fit3$X,"scaled:scale")[d+1]^2

      # mean
      predy <- (exp(-distance(t(t(x)/sqrt(theta[-(d+1)])), t(t(X3)/sqrt(theta[-(d+1)])))) *
                        1/sqrt(1+2*sig2/theta[d+1]) *
                        exp(-(drop(outer(x.mu, w2.x3, FUN="-")))^2/(theta[d+1]+2*sig2))) %*% a

      # var
      predsig2 <- c(rep(0, nrow(x)))
      for(i in 1: nrow(x)){ # each test point
        mat <- drop(exp(-distance(t(x[i,])/sqrt(theta[-(d+1)]), t(t(X3)/sqrt(theta[-(d+1)])))) %o%
                      exp(-distance(t(x[i,])/sqrt(theta[-(d+1)]), t(t(X3)/sqrt(theta[-(d+1)]))))) * # common components
          1/sqrt(1+4*sig2[i]/theta[d+1]) *
          exp(-(outer(w2.x3, w2.x3, FUN="+")/2 - x.mu[i])^2/(theta[d+1]/2+2*sig2[i])) *
          exp(-(outer(w2.x3, w2.x3, FUN="-"))^2/(2*theta[d+1]))

        predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }
    }
  }else if(kernel=="matern1.5"){
    if(constant){
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      pred.RNAmf1 <- predRNAmf(fit.RNAmf1, x)
      x.mu <- pred.RNAmf1$mu
      sig2 <- pred.RNAmf1$sig2

      ### calculate the closed form ###
      X3 <- matrix(fit3$X[,-(d+1)], ncol=d)
      w2.x3 <- fit3$X[,d+1]
      y3 <- fit3$y
      n <- length(y3)
      theta <- fit3$theta
      tau2hat <- fit3$tau2hat
      mu3 <- fit3$mu.hat

      Ci <- fit3$Ki
      a <- Ci %*% (y3 - mu3)

      ### scale new inputs ###
      x <- t((t(x)-attr(fit3$X,"scaled:center")[1:d])/attr(fit3$X,"scaled:scale")[1:d])
      x.mu <- t((t(x.mu)-attr(fit3$X,"scaled:center")[d+1])/attr(fit3$X,"scaled:scale")[d+1])
      sig2 <- sig2/attr(fit3$X,"scaled:scale")[d+1]^2

      # mean
      predy <- c(rep(0, nrow(x)))
      for(j in 1: nrow(x)){ # each test point
        mua <- x.mu[j] - sqrt(3)*sig2[j]/theta[d+1]
        mub <- x.mu[j] + sqrt(3)*sig2[j]/theta[d+1]

        lambda11 <- c(1, mua)
        lambda12 <- c(0, 1)
        lambda21 <- c(1, -mub)

        e1 <- cbind(matrix(1-sqrt(3)*w2.x3/theta[d+1]), sqrt(3)/theta[d+1])
        e2 <- cbind(matrix(1+sqrt(3)*w2.x3/theta[d+1]), sqrt(3)/theta[d+1])

        predy[j] <- mu3 + drop(t(t(cor.sep(t(x[j,]), X3, theta[-(d+1)], nu=1.5)) * # common but depends on kernel
                                   (exp((3*sig2[j] + 2*sqrt(3)*theta[d+1]*(w2.x3 - x.mu[j]))/(2*theta[d+1]^2)) *
                                      (e1 %*% lambda11 * pnorm((mua - w2.x3)/sqrt(sig2[j])) +
                                         e1 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3 - mua)^2/(2*sig2[j]))) +
                                      exp((3*sig2[j] - 2*sqrt(3)*theta[d+1]*(w2.x3 - x.mu[j]))/(2*theta[d+1]^2)) *
                                      (e2 %*% lambda21 * pnorm((-mub + w2.x3)/sqrt(sig2[j])) +
                                         e2 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3 - mub)^2/(2*sig2[j]))))) %*% a)
      }


      # var
      predsig2 <- c(rep(0, nrow(x)))
      for(i in 1: nrow(x)){
        zeta <- function(x,y){zetafun(w1=x, w2=y, m=x.mu[i], s=sig2[i], nu=1.5, theta=theta[d+1])}

        mat <- drop(t(cor.sep(t(x[i,]), X3, theta[-(d+1)], nu=1.5)) %o% t(cor.sep(t(x[i,]), X3, theta[-(d+1)], nu=1.5))) * # constant depends on kernel
          outer(w2.x3, w2.x3, FUN=Vectorize(zeta))

        predsig2[i] <- pmax(0, tau2hat - (predy[i]-mu3)^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }

    }else{
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      pred.RNAmf1 <- predRNAmf(fit.RNAmf1, x)
      x.mu <- pred.RNAmf1$mu
      sig2 <- pred.RNAmf1$sig2

      ### calculate the closed form ###
      X3 <- matrix(fit3$X[,-(d+1)], ncol=d)
      w2.x3 <- fit3$X[,d+1]
      y3 <- fit3$y
      n <- length(y3)
      theta <- fit3$theta
      tau2hat <- fit3$tau2hat

      Ci <- fit3$Ki
      a <- Ci %*% (y3 + attr(y3, "scaled:center"))

      ### scale new inputs ###
      x <- t((t(x)-attr(fit3$X,"scaled:center")[1:d])/attr(fit3$X,"scaled:scale")[1:d])
      x.mu <- t((t(x.mu)-attr(fit3$X,"scaled:center")[d+1])/attr(fit3$X,"scaled:scale")[d+1])
      sig2 <- sig2/attr(fit3$X,"scaled:scale")[d+1]^2

      # mean
      predy <- c(rep(0, nrow(x)))
      for(j in 1: nrow(x)){ # each test point
        mua <- x.mu[j] - sqrt(3)*sig2[j]/theta[d+1]
        mub <- x.mu[j] + sqrt(3)*sig2[j]/theta[d+1]

        lambda11 <- c(1, mua)
        lambda12 <- c(0, 1)
        lambda21 <- c(1, -mub)

        e1 <- cbind(matrix(1-sqrt(3)*w2.x3/theta[d+1]), sqrt(3)/theta[d+1])
        e2 <- cbind(matrix(1+sqrt(3)*w2.x3/theta[d+1]), sqrt(3)/theta[d+1])

        predy[j] <- drop(t(t(cor.sep(t(x[j,]), X3, theta[-(d+1)], nu=1.5)) * # common but depends on kernel
                             (exp((3*sig2[j] + 2*sqrt(3)*theta[d+1]*(w2.x3 - x.mu[j]))/(2*theta[d+1]^2)) *
                                (e1 %*% lambda11 * pnorm((mua - w2.x3)/sqrt(sig2[j])) +
                                   e1 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3 - mua)^2/(2*sig2[j]))) +
                                exp((3*sig2[j] - 2*sqrt(3)*theta[d+1]*(w2.x3 - x.mu[j]))/(2*theta[d+1]^2)) *
                                (e2 %*% lambda21 * pnorm((-mub + w2.x3)/sqrt(sig2[j])) +
                                   e2 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3 - mub)^2/(2*sig2[j]))))) %*% a)
      }

      # var
      predsig2 <- c(rep(0, nrow(x)))
      for(i in 1: nrow(x)){
        zeta <- function(x,y){zetafun(w1=x, w2=y, m=x.mu[i], s=sig2[i], nu=1.5, theta=theta[d+1])}

        mat <- drop(t(cor.sep(t(x[i,]), X3, theta[-(d+1)], nu=1.5)) %o% t(cor.sep(t(x[i,]), X3, theta[-(d+1)], nu=1.5))) * # constant depends on kernel
          outer(w2.x3, w2.x3, FUN=Vectorize(zeta))

        predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }

    }
  }else if(kernel=="matern2.5"){
    if(constant){
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      pred.RNAmf1 <- predRNAmf(fit.RNAmf1, x)
      x.mu <- pred.RNAmf1$mu
      sig2 <- pred.RNAmf1$sig2

      ### calculate the closed form ###
      X3 <- matrix(fit3$X[,-(d+1)], ncol=d)
      w2.x3 <- fit3$X[,d+1]
      y3 <- fit3$y
      n <- length(y3)
      theta <- fit3$theta
      tau2hat <- fit3$tau2hat
      mu3 <- fit3$mu.hat

      Ci <- fit3$Ki
      a <- Ci %*% (y3 - mu3)

      ### scale new inputs ###
      x <- t((t(x)-attr(fit3$X,"scaled:center")[1:d])/attr(fit3$X,"scaled:scale")[1:d])
      x.mu <- t((t(x.mu)-attr(fit3$X,"scaled:center")[d+1])/attr(fit3$X,"scaled:scale")[d+1])
      sig2 <- sig2/attr(fit3$X,"scaled:scale")[d+1]^2

      # mean
      predy <- c(rep(0, nrow(x)))
      for(j in 1: nrow(x)){ # each test point
        mua <- x.mu[j] - sqrt(5)*sig2[j]/theta[d+1]
        mub <- x.mu[j] + sqrt(5)*sig2[j]/theta[d+1]

        lambda11 <- c(1, mua, mua^2+sig2[j])
        lambda12 <- cbind(0, 1, matrix(mua+w2.x3))
        lambda21 <- c(1, -mub, mub^2+sig2[j])
        lambda22 <- cbind(0, 1, matrix(-mub-w2.x3))

        e1 <- cbind(matrix(1-sqrt(5)*w2.x3/theta[d+1]+5*w2.x3^2/(3*theta[d+1]^2)),
                    matrix(sqrt(5)/theta[d+1]-10*w2.x3/(3*theta[d+1]^2)),
                    5/(3*theta[d+1]^2))
        e2 <- cbind(matrix(1+sqrt(5)*w2.x3/theta[d+1]+5*w2.x3^2/(3*theta[d+1]^2)),
                    matrix(sqrt(5)/theta[d+1]+10*w2.x3/(3*theta[d+1]^2)),
                    5/(3*theta[d+1]^2))

        predy[j] <- mu3 + drop(t(t(cor.sep(t(x[j,]), X3, theta[-(d+1)], nu=2.5)) *
                                   (exp((5*sig2[j] + 2*sqrt(5)*theta[d+1]*(w2.x3 - x.mu[j]))/(2*theta[d+1]^2)) *
                                      (e1 %*% lambda11 * pnorm((mua - w2.x3)/sqrt(sig2[j])) +
                                         rowSums(e1 * lambda12) * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3 - mua)^2/(2*sig2[j]))) +
                                      exp((5*sig2[j] - 2*sqrt(5)*theta[d+1]*(w2.x3 - x.mu[j]))/(2*theta[d+1]^2)) *
                                      (e2 %*% lambda21 * pnorm((-mub + w2.x3)/sqrt(sig2[j])) +
                                         rowSums(e2 * lambda22) * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3 - mub)^2/(2*sig2[j]))))) %*% a)
      }

      # var
      predsig2 <- c(rep(0, nrow(x)))
      for(i in 1: nrow(x)){
        zeta <- function(x,y){zetafun(w1=x, w2=y, m=x.mu[i], s=sig2[i], nu=2.5, theta=theta[d+1])}

        mat <- drop(t(cor.sep(t(x[i,]), X3, theta[-(d+1)], nu=2.5)) %o% t(cor.sep(t(x[i,]), X3, theta[-(d+1)], nu=2.5))) * # constant depends on kernel
          outer(w2.x3, w2.x3, FUN=Vectorize(zeta))

        predsig2[i] <- pmax(0, tau2hat - (predy[i]-mu3)^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }
    }else{
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      pred.RNAmf1 <- predRNAmf(fit.RNAmf1, x)
      x.mu <- pred.RNAmf1$mu
      sig2 <- pred.RNAmf1$sig2

      ### calculate the closed form ###
      X3 <- matrix(fit3$X[,-(d+1)], ncol=d)
      w2.x3 <- fit3$X[,d+1]
      y3 <- fit3$y
      n <- length(y3)
      theta <- fit3$theta
      tau2hat <- fit3$tau2hat

      Ci <- fit3$Ki
      a <- Ci %*% (y3 + attr(y3, "scaled:center"))

      ### scale new inputs ###
      x <- t((t(x)-attr(fit3$X,"scaled:center")[1:d])/attr(fit3$X,"scaled:scale")[1:d])
      x.mu <- t((t(x.mu)-attr(fit3$X,"scaled:center")[d+1])/attr(fit3$X,"scaled:scale")[d+1])
      sig2 <- sig2/attr(fit3$X,"scaled:scale")[d+1]^2

      # mean
      predy <- c(rep(0, nrow(x)))
      for(j in 1: nrow(x)){ # each test point
        mua <- x.mu[j] - sqrt(5)*sig2[j]/theta[d+1]
        mub <- x.mu[j] + sqrt(5)*sig2[j]/theta[d+1]

        lambda11 <- c(1, mua, mua^2+sig2[j])
        lambda12 <- cbind(0, 1, matrix(mua+w2.x3))
        lambda21 <- c(1, -mub, mub^2+sig2[j])
        lambda22 <- cbind(0, 1, matrix(-mub-w2.x3))

        e1 <- cbind(matrix(1-sqrt(5)*w2.x3/theta[d+1]+5*w2.x3^2/(3*theta[d+1]^2)),
                    matrix(sqrt(5)/theta[d+1]-10*w2.x3/(3*theta[d+1]^2)),
                    5/(3*theta[d+1]^2))
        e2 <- cbind(matrix(1+sqrt(5)*w2.x3/theta[d+1]+5*w2.x3^2/(3*theta[d+1]^2)),
                    matrix(sqrt(5)/theta[d+1]+10*w2.x3/(3*theta[d+1]^2)),
                    5/(3*theta[d+1]^2))

        predy[j] <- drop(t(t(cor.sep(t(x[j,]), X3, theta[-(d+1)], nu=2.5)) *
                             (exp((5*sig2[j] + 2*sqrt(5)*theta[d+1]*(w2.x3 - x.mu[j]))/(2*theta[d+1]^2)) *
                                (e1 %*% lambda11 * pnorm((mua - w2.x3)/sqrt(sig2[j])) +
                                   rowSums(e1 * lambda12) * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3 - mua)^2/(2*sig2[j]))) +
                                exp((5*sig2[j] - 2*sqrt(5)*theta[d+1]*(w2.x3 - x.mu[j]))/(2*theta[d+1]^2)) *
                                (e2 %*% lambda21 * pnorm((-mub + w2.x3)/sqrt(sig2[j])) +
                                   rowSums(e2 * lambda22) * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3 - mub)^2/(2*sig2[j]))))) %*% a)
      }

      # var
      predsig2 <- c(rep(0, nrow(x)))
      for(i in 1: nrow(x)){
        zeta <- function(x,y){zetafun(w1=x, w2=y, m=x.mu[i], s=sig2[i], nu=2.5, theta=theta[d+1])}

        mat <- drop(t(cor.sep(t(x[i,]), X3, theta[-(d+1)], nu=2.5)) %o% t(cor.sep(t(x[i,]), X3, theta[-(d+1)], nu=2.5))) * # constant depends on kernel
          outer(w2.x3, w2.x3, FUN=Vectorize(zeta))

        predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }
    }
  }

  return(list(mu=predy, sig2=predsig2))
}

