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
      fit.RNAmf1 <- predRNAmf(fit.RNAmf1, x)
      x.mu <- fit.RNAmf1$mu
      sig2 <- fit.RNAmf1$sig2

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
      predy <- c(rep(0, nrow(x)))

      for(j in 1: nrow(x)){
        predv <- c(rep(1,n))
        for(i in 1:n){
          for(m in 1:d){
            predv[i] <- predv[i] * exp(-(x[j,m]-X3[i,m])^2/theta[m])
          }
          predv[i] <- predv[i] * 1/sqrt(1+2*sig2[j]/theta[d+1]) *
            exp(-(w2.x3[i]-x.mu[j])^2/(theta[d+1]+2*sig2[j]))
        }
        predy[j] <- mu3 + drop(predv%*%a)
      }

      # var
      predsig2 <- c(rep(0, nrow(x)))

      for(i in 1: nrow(x)){
        mat <- matrix(1, n, n)
        for(k in 1:n){
          for(l in 1:n){
            for(m in 1:d){
              mat[k,l] <- mat[k,l] * exp(-((x[i,m]-X3[k,m])^2+(x[i,m]-X3[l,m])^2)/theta[m])
            }
            mat[k,l] <- mat[k,l] * (a[k]*a[l] - tau2hat*Ci[k,l]) *
              1/sqrt(1+4*sig2[i]/theta[d+1]) *
              exp(-((w2.x3[k]+w2.x3[l])/2-x.mu[i])^2/(theta[d+1]/2+2*sig2[i])) *
              exp(-(w2.x3[k]-w2.x3[l])^2/(2*theta[d+1]))
          }}

        predsig2[i] <- pmax(0, tau2hat - (predy[i]-mu3)^2 + sum(mat))
      }
    }else{
      fit.RNAmf1 <- predRNAmf(fit.RNAmf1, x)
      x.mu <- fit.RNAmf1$mu
      sig2 <- fit.RNAmf1$sig2

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
      predy <- c(rep(0, nrow(x)))

      for(j in 1: nrow(x)){
        predv <- c(rep(1,n))
        for(i in 1:n){
          for(m in 1:d){
            predv[i] <- predv[i] * exp(-(x[j,m]-X3[i,m])^2/theta[m])
          }
          predv[i] <- predv[i] * 1/sqrt(1+2*sig2[j]/theta[d+1]) *
            exp(-(w2.x3[i]-x.mu[j])^2/(theta[d+1]+2*sig2[j]))
        }
        predy[j] <- drop(predv%*%a)
      }

      # var
      predsig2 <- c(rep(0, nrow(x)))

      for(i in 1: nrow(x)){
        mat <- matrix(1, n, n)
        for(k in 1:n){
          for(l in 1:n){
            for(m in 1:d){
              mat[k,l] <- mat[k,l] * exp(-((x[i,m]-X3[k,m])^2+(x[i,m]-X3[l,m])^2)/theta[m])
            }
            mat[k,l] <- mat[k,l] * (a[k]*a[l] - tau2hat*Ci[k,l]) *
              1/sqrt(1+4*sig2[i]/theta[d+1]) *
              exp(-((w2.x3[k]+w2.x3[l])/2-x.mu[i])^2/(theta[d+1]/2+2*sig2[i])) *
              exp(-(w2.x3[k]-w2.x3[l])^2/(2*theta[d+1]))
          }}

        predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + sum(mat))
      }
    }
  }else if(kernel=="matern1.5"){
    if(constant){
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      fit.RNAmf1 <- predRNAmf(fit.RNAmf1, x)
      x.mu <- fit.RNAmf1$mu
      sig2 <- fit.RNAmf1$sig2

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
        predv <- c(rep(1,n))
        for(i in 1:n){ # each row of train set
          for(m in 1:d){ # dim of train set
            predv[i] <- predv[i] * matern.kernel(sqrt(distance(t(t(x[j,m])/theta[m]), t(t(X3[i,m])/theta[m]))), nu=1.5) # common but depends on kernel
          } # depends on kernel structure

          mua <- x.mu[j] - sqrt(3)*sig2[j]/theta[d+1]
          mub <- x.mu[j] + sqrt(3)*sig2[j]/theta[d+1]

          lambda11 <- c(1, mua)
          lambda12 <- c(0, 1)
          lambda21 <- c(1, -mub)

          e1 <- c(1-sqrt(3)*w2.x3[i]/theta[d+1], sqrt(3)/theta[d+1])
          e2 <- c(1+sqrt(3)*w2.x3[i]/theta[d+1], sqrt(3)/theta[d+1])

          predv[i] <- predv[i] * (exp((3*sig2[j] + 2*sqrt(3)*theta[d+1]*(w2.x3[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                    (e1 %*% lambda11 * pnorm((mua - w2.x3[i])/sqrt(sig2[j])) +
                                       e1 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3[i] - mua)^2/(2*sig2[j]))) +
                                    exp((3*sig2[j] - 2*sqrt(3)*theta[d+1]*(w2.x3[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                    (e2 %*% lambda21 * pnorm((-mub + w2.x3[i])/sqrt(sig2[j])) +
                                       e2 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3[i] - mub)^2/(2*sig2[j]))))
        }
        predy[j] <- mu3 + drop(predv%*%a)
      }

      # var
      predsig2 <- c(rep(0, nrow(x)))

      for(i in 1: nrow(x)){
        mat <- matrix(1, n, n) # matrix J in Ming's paper
        for(k in 1:n){
          for(l in 1:n){
            for(m in 1:d){ # constant depends on kernel
              mat[k,l] <- mat[k,l] * matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X3[k,m])/theta[m]))), nu=1.5) *
                matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X3[l,m])/theta[m]))), nu=1.5)
            } # expected depends on kernel structure
            mat[k,l] <- mat[k,l] * zetafun(w1=w2.x3[k], w2=w2.x3[l], m=x.mu[i], s=sig2[i], nu=1.5, theta=theta[d+1])
          }}

        predsig2[i] <- pmax(0, tau2hat - (predy[i]-mu3)^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }
    }else{
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      fit.RNAmf1 <- predRNAmf(fit.RNAmf1, x)
      x.mu <- fit.RNAmf1$mu
      sig2 <- fit.RNAmf1$sig2

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
        predv <- c(rep(1,n))
        for(i in 1:n){ # each row of train set
          for(m in 1:d){ # dim of train set
            predv[i] <- predv[i] * matern.kernel(sqrt(distance(t(t(x[j,m])/theta[m]), t(t(X3[i,m])/theta[m]))), nu=1.5) # common but depends on kernel
          } # depends on kernel structure

          mua <- x.mu[j] - sqrt(3)*sig2[j]/theta[d+1]
          mub <- x.mu[j] + sqrt(3)*sig2[j]/theta[d+1]

          lambda11 <- c(1, mua)
          lambda12 <- c(0, 1)
          lambda21 <- c(1, -mub)

          e1 <- c(1-sqrt(3)*w2.x3[i]/theta[d+1], sqrt(3)/theta[d+1])
          e2 <- c(1+sqrt(3)*w2.x3[i]/theta[d+1], sqrt(3)/theta[d+1])

          predv[i] <- predv[i] * (exp((3*sig2[j] + 2*sqrt(3)*theta[d+1]*(w2.x3[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                    (e1 %*% lambda11 * pnorm((mua - w2.x3[i])/sqrt(sig2[j])) +
                                       e1 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3[i] - mua)^2/(2*sig2[j]))) +
                                    exp((3*sig2[j] - 2*sqrt(3)*theta[d+1]*(w2.x3[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                    (e2 %*% lambda21 * pnorm((-mub + w2.x3[i])/sqrt(sig2[j])) +
                                       e2 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3[i] - mub)^2/(2*sig2[j]))))
        }
        predy[j] <- drop(predv%*%a)
      }

      # var
      predsig2 <- c(rep(0, nrow(x)))

      for(i in 1: nrow(x)){
        mat <- matrix(1, n, n) # matrix J in Ming's paper
        for(k in 1:n){
          for(l in 1:n){
            for(m in 1:d){ # constant depends on kernel
              mat[k,l] <- mat[k,l] * matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X3[k,m])/theta[m]))), nu=1.5) *
                matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X3[l,m])/theta[m]))), nu=1.5)
            } # expected depends on kernel structure
            mat[k,l] <- mat[k,l] * zetafun(w1=w2.x3[k], w2=w2.x3[l], m=x.mu[i], s=sig2[i], nu=1.5, theta=theta[d+1])
          }}

        predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }
    }
  }else if(kernel=="matern2.5"){
    if(constant){
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      fit.RNAmf1 <- predRNAmf(fit.RNAmf1, x)
      x.mu <- fit.RNAmf1$mu
      sig2 <- fit.RNAmf1$sig2

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

      for(j in 1: nrow(x)){
        predv <- c(rep(1,n))
        for(i in 1:n){
          for(m in 1:d){
            predv[i] <- predv[i] * matern.kernel(sqrt(distance(t(t(x[j,m])/theta[m]), t(t(X3[i,m])/theta[m]))), nu=2.5) # common but depends on kernel
          } # depends on kernel structure

          mua <- x.mu[j] - sqrt(5)*sig2[j]/theta[d+1]
          mub <- x.mu[j] + sqrt(5)*sig2[j]/theta[d+1]

          lambda11 <- c(1, mua, mua^2+sig2[j])
          lambda12 <- c(0, 1, mua+w2.x3[i])
          lambda21 <- c(1, -mub, mub^2+sig2[j])
          lambda22 <- c(0, 1, -mub-w2.x3[i])

          e1 <- c(1-sqrt(5)*w2.x3[i]/theta[d+1]+5*w2.x3[i]^2/(3*theta[d+1]^2),
                  sqrt(5)/theta[d+1]-10*w2.x3[i]/(3*theta[d+1]^2),
                  5/(3*theta[d+1]^2))
          e2 <- c(1+sqrt(5)*w2.x3[i]/theta[d+1]+5*w2.x3[i]^2/(3*theta[d+1]^2),
                  sqrt(5)/theta[d+1]+10*w2.x3[i]/(3*theta[d+1]^2),
                  5/(3*theta[d+1]^2))

          predv[i] <- predv[i] * (exp((5*sig2[j] + 2*sqrt(5)*theta[d+1]*(w2.x3[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                    (e1 %*% lambda11 * pnorm((mua - w2.x3[i])/sqrt(sig2[j])) +
                                       e1 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3[i] - mua)^2/(2*sig2[j]))) +
                                    exp((5*sig2[j] - 2*sqrt(5)*theta[d+1]*(w2.x3[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                    (e2 %*% lambda21 * pnorm((-mub + w2.x3[i])/sqrt(sig2[j])) +
                                       e2 %*% lambda22 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3[i] - mub)^2/(2*sig2[j]))))
        }
        predy[j] <- mu3 + drop(predv%*%a)
      }

      # var
      predsig2 <- c(rep(0, nrow(x)))

      for(i in 1: nrow(x)){
        mat <- matrix(1, n, n)
        for(k in 1:n){
          for(l in 1:n){
            for(m in 1:d){ # common but depends on kernel
              mat[k,l] <- mat[k,l] * matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X3[k,m])/theta[m]))), nu=2.5) *
                matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X3[l,m])/theta[m]))), nu=2.5)
            } # expected depends on kernel structure
            mat[k,l] <- mat[k,l] * zetafun(w1=w2.x3[k], w2=w2.x3[l], m=x.mu[i], s=sig2[i], nu=2.5, theta=theta[d+1])
          }}

        predsig2[i] <- pmax(0, tau2hat - (predy[i]-mu3)^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }
    }else{
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      fit.RNAmf1 <- predRNAmf(fit.RNAmf1, x)
      x.mu <- fit.RNAmf1$mu
      sig2 <- fit.RNAmf1$sig2

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

      for(j in 1: nrow(x)){
        predv <- c(rep(1,n))
        for(i in 1:n){
          for(m in 1:d){
            predv[i] <- predv[i] * matern.kernel(sqrt(distance(t(t(x[j,m])/theta[m]), t(t(X3[i,m])/theta[m]))), nu=2.5) # common but depends on kernel
          } # depends on kernel structure

          mua <- x.mu[j] - sqrt(5)*sig2[j]/theta[d+1]
          mub <- x.mu[j] + sqrt(5)*sig2[j]/theta[d+1]

          lambda11 <- c(1, mua, mua^2+sig2[j])
          lambda12 <- c(0, 1, mua+w2.x3[i])
          lambda21 <- c(1, -mub, mub^2+sig2[j])
          lambda22 <- c(0, 1, -mub-w2.x3[i])

          e1 <- c(1-sqrt(5)*w2.x3[i]/theta[d+1]+5*w2.x3[i]^2/(3*theta[d+1]^2),
                  sqrt(5)/theta[d+1]-10*w2.x3[i]/(3*theta[d+1]^2),
                  5/(3*theta[d+1]^2))
          e2 <- c(1+sqrt(5)*w2.x3[i]/theta[d+1]+5*w2.x3[i]^2/(3*theta[d+1]^2),
                  sqrt(5)/theta[d+1]+10*w2.x3[i]/(3*theta[d+1]^2),
                  5/(3*theta[d+1]^2))

          predv[i] <- predv[i] * (exp((5*sig2[j] + 2*sqrt(5)*theta[d+1]*(w2.x3[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                    (e1 %*% lambda11 * pnorm((mua - w2.x3[i])/sqrt(sig2[j])) +
                                       e1 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3[i] - mua)^2/(2*sig2[j]))) +
                                    exp((5*sig2[j] - 2*sqrt(5)*theta[d+1]*(w2.x3[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                    (e2 %*% lambda21 * pnorm((-mub + w2.x3[i])/sqrt(sig2[j])) +
                                       e2 %*% lambda22 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w2.x3[i] - mub)^2/(2*sig2[j]))))
        }
        predy[j] <- drop(predv%*%a)
      }

      # var
      predsig2 <- c(rep(0, nrow(x)))

      for(i in 1: nrow(x)){
        mat <- matrix(1, n, n)
        for(k in 1:n){
          for(l in 1:n){
            for(m in 1:d){ # common but depends on kernel
              mat[k,l] <- mat[k,l] * matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X3[k,m])/theta[m]))), nu=2.5) *
                matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X3[l,m])/theta[m]))), nu=2.5)
            } # expected depends on kernel structure
            mat[k,l] <- mat[k,l] * zetafun(w1=w2.x3[k], w2=w2.x3[l], m=x.mu[i], s=sig2[i], nu=2.5, theta=theta[d+1])
          }}

        predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }
    }
  }

  return(list(mu=predy, sig2=predsig2))
}

