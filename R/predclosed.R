#' predclosed
#'
#' predictive posterior mean and variance of the model with fidelity level 1 and 2.
#'
#' @param fit an object of class closed.
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

predclosed <- function(fit, x){
  kernel <- fit$kernel
  constant <- fit$constant
  fit1 <- fit$fit1
  fit2 <- fit$fit2

  if(kernel=="sqex"){
    if(constant){
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      x.mu <- pred.GP(fit1, x)$mu # mean of f1(u)
      sig2 <- pred.GP(fit1, x)$sig2 #*0

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
      x <- t((t(x)-attr(fit2$X,"scaled:center")[1:d])/attr(fit2$X,"scaled:scale")[1:d])
      x.mu <- t((t(x.mu)-attr(fit2$X,"scaled:center")[d+1])/attr(fit2$X,"scaled:scale")[d+1])
      sig2 <- sig2/attr(fit2$X,"scaled:scale")[d+1]^2

      # mean
      predy <- c(rep(0, nrow(x)))

      for(j in 1: nrow(x)){ # each test point
        predv <- c(rep(1,n))
        for(i in 1:n){ # each row of train set
          for(m in 1:d){ # dim of train set
            predv[i] <- predv[i] * exp(-(x[j,m]-X2[i,m])^2/theta[m]) # common components
          } # depends on kernel structure
          predv[i] <- predv[i] * 1/sqrt(1+2*sig2[j]/theta[d+1]) *
            exp(-(w1.x2[i]-x.mu[j])^2/(theta[d+1]+2*sig2[j]))
        }
        predy[j] <- mu2 + drop(predv%*%a)
      }

      # var
      predsig2 <- c(rep(0, nrow(x)))

      for(i in 1: nrow(x)){ # each test point
        mat <- matrix(1, n, n)
        for(k in 1:n){ # each row of train set
          for(l in 1:n){ # dim of train set
            for(m in 1:d){
              mat[k,l] <- mat[k,l] * exp(-((x[i,m]-X2[k,m])^2+(x[i,m]-X2[l,m])^2)/theta[m]) # common components
            }
            mat[k,l] <- mat[k,l] * #(a[k]*a[l] - tau2hat*Ci[k,l]) *
              1/sqrt(1+4*sig2[i]/theta[d+1]) *
              exp(-((w1.x2[k]+w1.x2[l])/2-x.mu[i])^2/(theta[d+1]/2+2*sig2[i])) *
              exp(-(w1.x2[k]-w1.x2[l])^2/(2*theta[d+1]))
          }}

        predsig2[i] <- pmax(0, tau2hat - (predy[i]-mu2)^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }
    }else{
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      x.mu <- pred.GP(fit1, x)$mu # mean of f1(u)
      sig2 <- pred.GP(fit1, x)$sig2

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
      x <- t((t(x)-attr(fit2$X,"scaled:center")[1:d])/attr(fit2$X,"scaled:scale")[1:d])
      x.mu <- t((t(x.mu)-attr(fit2$X,"scaled:center")[d+1])/attr(fit2$X,"scaled:scale")[d+1])
      sig2 <- sig2/attr(fit2$X,"scaled:scale")[d+1]^2

      # mean
      predy <- c(rep(0, nrow(x)))

      for(j in 1: nrow(x)){
        predv <- c(rep(1,n))
        for(i in 1:n){
          for(m in 1:d){
            predv[i] <- predv[i] * exp(-(x[j,m]-X2[i,m])^2/theta[m])
          }
          predv[i] <- predv[i] * 1/sqrt(1+2*sig2[j]/theta[d+1]) *
            exp(-(w1.x2[i]-x.mu[j])^2/(theta[d+1]+2*sig2[j]))
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
              mat[k,l] <- mat[k,l] * exp(-((x[i,m]-X2[k,m])^2+(x[i,m]-X2[l,m])^2)/theta[m])
            }
            mat[k,l] <- mat[k,l] * #(a[k]*a[l] - tau2hat*Ci[k,l]) *
              1/sqrt(1+4*sig2[i]/theta[d+1]) *
              exp(-((w1.x2[k]+w1.x2[l])/2-x.mu[i])^2/(theta[d+1]/2+2*sig2[i])) *
              exp(-(w1.x2[k]-w1.x2[l])^2/(2*theta[d+1]))
          }}

        predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }
    }
  }else if(kernel=="matern1.5"){
    if(constant){
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      x.mu <- pred.matGP(fit1, x)$mu # mean of f1(u)
      sig2 <- pred.matGP(fit1, x)$sig2 #*0

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
      x <- t((t(x)-attr(fit2$X,"scaled:center")[1:d])/attr(fit2$X,"scaled:scale")[1:d])
      x.mu <- t((t(x.mu)-attr(fit2$X,"scaled:center")[d+1])/attr(fit2$X,"scaled:scale")[d+1])
      sig2 <- sig2/attr(fit2$X,"scaled:scale")[d+1]^2

      # mean
      predy <- c(rep(0, nrow(x)))

      for(j in 1: nrow(x)){ # each test point
        predv <- c(rep(1,n))
        for(i in 1:n){ # each row of train set
          for(m in 1:d){ # dim of train set
            predv[i] <- predv[i] * matern.kernel(sqrt(distance(t(t(x[j,m])/theta[m]), t(t(X2[i,m])/theta[m]))), nu=1.5) # common but depends on kernel
          } # depends on kernel structure

          mua <- x.mu[j] - sqrt(3)*sig2[j]/theta[d+1]
          mub <- x.mu[j] + sqrt(3)*sig2[j]/theta[d+1]

          lambda11 <- c(1, mua)
          lambda12 <- c(0, 1)
          lambda21 <- c(1, -mub)

          e1 <- c(1-sqrt(3)*w1.x2[i]/theta[d+1], sqrt(3)/theta[d+1])
          e2 <- c(1+sqrt(3)*w1.x2[i]/theta[d+1], sqrt(3)/theta[d+1])

          predv[i] <- predv[i] * (exp((3*sig2[j] + 2*sqrt(3)*theta[d+1]*(w1.x2[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                    (e1 %*% lambda11 * pnorm((mua - w1.x2[i])/sqrt(sig2[j])) +
                                       e1 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2[i] - mua)^2/(2*sig2[j]))) +
                                    exp((3*sig2[j] - 2*sqrt(3)*theta[d+1]*(w1.x2[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                    (e2 %*% lambda21 * pnorm((-mub + w1.x2[i])/sqrt(sig2[j])) +
                                       e2 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2[i] - mub)^2/(2*sig2[j]))))
        }
        predy[j] <- mu2 + drop(predv%*%a)
      }

      # var
      predsig2 <- c(rep(0, nrow(x)))

      for(i in 1: nrow(x)){
        mat <- matrix(1, n, n) # matrix J in Ming's paper
        for(k in 1:n){
          for(l in 1:n){
            for(m in 1:d){ # constant depends on kernel
              mat[k,l] <- mat[k,l] * matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X2[k,m])/theta[m]))), nu=1.5) *
                matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X2[l,m])/theta[m]))), nu=1.5)
            } # expected depends on kernel structure
            mat[k,l] <- mat[k,l] * zetafun(w1=w1.x2[k], w2=w1.x2[l], m=x.mu[i], s=sig2[i], nu=1.5, theta=theta[d+1])
          }}

        predsig2[i] <- pmax(0, tau2hat - (predy[i]-mu2)^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }

    }else{
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      x.mu <- pred.matGP(fit1, x)$mu # mean of f1(u)
      sig2 <- pred.matGP(fit1, x)$sig2

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
      x <- t((t(x)-attr(fit2$X,"scaled:center")[1:d])/attr(fit2$X,"scaled:scale")[1:d])
      x.mu <- t((t(x.mu)-attr(fit2$X,"scaled:center")[d+1])/attr(fit2$X,"scaled:scale")[d+1])
      sig2 <- sig2/attr(fit2$X,"scaled:scale")[d+1]^2

      # mean
      predy <- c(rep(0, nrow(x)))

      for(j in 1: nrow(x)){ # each test point
        predv <- c(rep(1,n))
        for(i in 1:n){ # each row of train set
          for(m in 1:d){ # dim of train set
            predv[i] <- predv[i] * matern.kernel(sqrt(distance(t(t(x[j,m])/theta[m]), t(t(X2[i,m])/theta[m]))), nu=1.5) # common but depends on kernel
          } # depends on kernel structure

          mua <- x.mu[j] - sqrt(3)*sig2[j]/theta[d+1]
          mub <- x.mu[j] + sqrt(3)*sig2[j]/theta[d+1]

          lambda11 <- c(1, mua)
          lambda12 <- c(0, 1)
          lambda21 <- c(1, -mub)

          e1 <- c(1-sqrt(3)*w1.x2[i]/theta[d+1], sqrt(3)/theta[d+1])
          e2 <- c(1+sqrt(3)*w1.x2[i]/theta[d+1], sqrt(3)/theta[d+1])

          predv[i] <- predv[i] * (exp((3*sig2[j] + 2*sqrt(3)*theta[d+1]*(w1.x2[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                    (e1 %*% lambda11 * pnorm((mua - w1.x2[i])/sqrt(sig2[j])) +
                                       e1 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2[i] - mua)^2/(2*sig2[j]))) +
                                    exp((3*sig2[j] - 2*sqrt(3)*theta[d+1]*(w1.x2[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                    (e2 %*% lambda21 * pnorm((-mub + w1.x2[i])/sqrt(sig2[j])) +
                                       e2 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2[i] - mub)^2/(2*sig2[j]))))
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
              mat[k,l] <- mat[k,l] * matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X2[k,m])/theta[m]))), nu=1.5) *
                matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X2[l,m])/theta[m]))), nu=1.5)
            } # expected depends on kernel structure
            mat[k,l] <- mat[k,l] * zetafun(w1=w1.x2[k], w2=w1.x2[l], m=x.mu[i], s=sig2[i], nu=1.5, theta=theta[d+1])
          }}

        predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }
    }
  }else if(kernel=="matern2.5"){
    if(constant){
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      x.mu <- pred.matGP(fit1, x)$mu # mean of f1(u)
      sig2 <- pred.matGP(fit1, x)$sig2 #*0

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
      x <- t((t(x)-attr(fit2$X,"scaled:center")[1:d])/attr(fit2$X,"scaled:scale")[1:d])
      x.mu <- t((t(x.mu)-attr(fit2$X,"scaled:center")[d+1])/attr(fit2$X,"scaled:scale")[d+1])
      sig2 <- sig2/attr(fit2$X,"scaled:scale")[d+1]^2

      # mean
      predy <- c(rep(0, nrow(x)))

      for(j in 1: nrow(x)){ # each test point
        predv <- c(rep(1,n))
        for(i in 1:n){ # each row of train set
          for(m in 1:d){ # dim of train set
            predv[i] <- predv[i] * matern.kernel(sqrt(distance(t(t(x[j,m])/theta[m]), t(t(X2[i,m])/theta[m]))), nu=2.5) # common but depends on kernel
          } # depends on kernel structure

          mua <- x.mu[j] - sqrt(5)*sig2[j]/theta[d+1]
          mub <- x.mu[j] + sqrt(5)*sig2[j]/theta[d+1]

          lambda11 <- c(1, mua, mua^2+sig2[j])
          lambda12 <- c(0, 1, mua+w1.x2[i])
          lambda21 <- c(1, -mub, mub^2+sig2[j])
          lambda22 <- c(0, 1, -mub-w1.x2[i])

          e1 <- c(1-sqrt(5)*w1.x2[i]/theta[d+1]+5*w1.x2[i]^2/(3*theta[d+1]^2),
                  sqrt(5)/theta[d+1]-10*w1.x2[i]/(3*theta[d+1]^2),
                  5/(3*theta[d+1]^2))
          e2 <- c(1+sqrt(5)*w1.x2[i]/theta[d+1]+5*w1.x2[i]^2/(3*theta[d+1]^2),
                  sqrt(5)/theta[d+1]+10*w1.x2[i]/(3*theta[d+1]^2),
                  5/(3*theta[d+1]^2))

          predv[i] <- predv[i] * (exp((5*sig2[j] + 2*sqrt(5)*theta[d+1]*(w1.x2[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                    (e1 %*% lambda11 * pnorm((mua - w1.x2[i])/sqrt(sig2[j])) +
                                       e1 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2[i] - mua)^2/(2*sig2[j]))) +
                                    exp((5*sig2[j] - 2*sqrt(5)*theta[d+1]*(w1.x2[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                    (e2 %*% lambda21 * pnorm((-mub + w1.x2[i])/sqrt(sig2[j])) +
                                       e2 %*% lambda22 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2[i] - mub)^2/(2*sig2[j]))))
        }
        predy[j] <- mu2 + drop(predv%*%a)
      }

      # var
      predsig2 <- c(rep(0, nrow(x)))

      for(i in 1: nrow(x)){
        mat <- matrix(1, n, n)
        for(k in 1:n){
          for(l in 1:n){
            for(m in 1:d){ # common but depends on kernel
              mat[k,l] <- mat[k,l] * matern.kernel(sqrt(distance(x[i,m]/theta[m], X2[k,m]/theta[m])), nu=2.5) *
                matern.kernel(sqrt(distance(x[i,m]/theta[m], X2[l,m]/theta[m])), nu=2.5)
            } # expected depends on kernel structure
            mat[k,l] <- mat[k,l] * zetafun(w1=w1.x2[k], w2=w1.x2[l], m=x.mu[i], s=sig2[i], nu=2.5, theta=theta[d+1])
          }}

        predsig2[i] <- pmax(0, tau2hat - (predy[i]-mu2)^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }
    }else{
      d <- ncol(fit1$X)
      x <- matrix(x, ncol=d)
      x.mu <- pred.matGP(fit1, x)$mu # mean of f1(u)
      sig2 <- pred.matGP(fit1, x)$sig2

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
      x <- t((t(x)-attr(fit2$X,"scaled:center")[1:d])/attr(fit2$X,"scaled:scale")[1:d])
      x.mu <- t((t(x.mu)-attr(fit2$X,"scaled:center")[d+1])/attr(fit2$X,"scaled:scale")[d+1])
      sig2 <- sig2/attr(fit2$X,"scaled:scale")[d+1]^2

      # mean
      predy <- c(rep(0, nrow(x)))

      for(j in 1: nrow(x)){ # each test point
        predv <- c(rep(1,n))
        for(i in 1:n){ # each row of train set
          for(m in 1:d){ # dim of train set
            predv[i] <- predv[i] * matern.kernel(sqrt(distance(x[j,m]/theta[m], X2[i,m]/theta[m])), nu=2.5) # common but depends on kernel
          } # depends on kernel structure

          mua <- x.mu[j] - sqrt(5)*sig2[j]/theta[d+1]
          mub <- x.mu[j] + sqrt(5)*sig2[j]/theta[d+1]

          lambda11 <- c(1, mua, mua^2+sig2[j])
          lambda12 <- c(0, 1, mua+w1.x2[i])
          lambda21 <- c(1, -mub, mub^2+sig2[j])
          lambda22 <- c(0, 1, -mub-w1.x2[i])

          e1 <- c(1-sqrt(5)*w1.x2[i]/theta[d+1]+5*w1.x2[i]^2/(3*theta[d+1]^2),
                  sqrt(5)/theta[d+1]-10*w1.x2[i]/(3*theta[d+1]^2),
                  5/(3*theta[d+1]^2))
          e2 <- c(1+sqrt(5)*w1.x2[i]/theta[d+1]+5*w1.x2[i]^2/(3*theta[d+1]^2),
                  sqrt(5)/theta[d+1]+10*w1.x2[i]/(3*theta[d+1]^2),
                  5/(3*theta[d+1]^2))

          predv[i] <- predv[i] * (exp((5*sig2[j] + 2*sqrt(5)*theta[d+1]*(w1.x2[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                    (e1 %*% lambda11 * pnorm((mua - w1.x2[i])/sqrt(sig2[j])) +
                                       e1 %*% lambda12 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2[i] - mua)^2/(2*sig2[j]))) +
                                    exp((5*sig2[j] - 2*sqrt(5)*theta[d+1]*(w1.x2[i] - x.mu[j]))/(2*theta[d+1]^2)) *
                                    (e2 %*% lambda21 * pnorm((-mub + w1.x2[i])/sqrt(sig2[j])) +
                                       e2 %*% lambda22 * sqrt(sig2[j])/sqrt(2*pi) * exp(-(w1.x2[i] - mub)^2/(2*sig2[j]))))
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
              mat[k,l] <- mat[k,l] * matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X2[k,m])/theta[m]))), nu=2.5) *
                matern.kernel(sqrt(distance(t(t(x[i,m])/theta[m]), t(t(X2[l,m])/theta[m]))), nu=2.5)
            } # expected depends on kernel structure
            mat[k,l] <- mat[k,l] * zetafun(w1=w1.x2[k], w2=w1.x2[l], m=x.mu[i], s=sig2[i], nu=2.5, theta=theta[d+1])
          }}

        predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a)%*%mat%*%a) - tau2hat*sum(diag(Ci%*%mat)))
      }
    }
  }

  return(list(mu=predy, sig2=predsig2))
}

