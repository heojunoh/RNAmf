#' predictive posterior mean and variance of the model with fidelity level 1, 2 and 3.
#'
#' @description The function computes the posterior mean and variance of RNA models with three fidelity levels
#' by fitted model using \code{\link{RNAmf_three_level}}.
#'
#' @seealso \code{\link{RNAmf_three_level}} for prediction.
#'
#' @details From the model fitted by \code{\link{RNAmf_three_level}},
#' the posterior mean and variance are calculated based on the closed form expression derived by a recursive fashion.
#' The formulae depend on its kernel choices.
#' For Matern kernels, \code{\link{zetafun}} computes the part of the posterior variance.
#' For details, see Heo and Sung (2023+, <arXiv:2309.11772>).
#'
#' @param fit an object of class RNAmf_three_level.
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

predRNAmf_three_level <- function(fit, x){
  kernel <- fit$kernel
  constant <- fit$constant
  fit.RNAmf_two_level <- fit$fit.RNAmf_two_level
  fit1 <- fit.RNAmf_two_level$fit1
  fit2 <- fit.RNAmf_two_level$fit2
  fit3 <- fit$fit3

  if(kernel=="sqex"){
    if(constant){
      pred.RNAmf_two_level <- predRNAmf_two_level(fit.RNAmf_two_level, x)
      x.mu <- pred.RNAmf_two_level$mu
      sig2 <- pred.RNAmf_two_level$sig2

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
      x <- scale.inputs(x, attr(fit3$X,"scaled:center")[1:d], attr(fit3$X,"scaled:scale")[1:d])
      x.mu <- scale.inputs(x.mu, attr(fit3$X,"scaled:center")[d+1], attr(fit3$X,"scaled:scale")[d+1])
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
      pred.RNAmf_two_level <- predRNAmf_two_level(fit.RNAmf_two_level, x)
      x.mu <- pred.RNAmf_two_level$mu
      sig2 <- pred.RNAmf_two_level$sig2

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
      x <- scale.inputs(x, attr(fit3$X,"scaled:center")[1:d], attr(fit3$X,"scaled:scale")[1:d])
      x.mu <- scale.inputs(x.mu, attr(fit3$X,"scaled:center")[d+1], attr(fit3$X,"scaled:scale")[d+1])
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
      pred.RNAmf_two_level <- predRNAmf_two_level(fit.RNAmf_two_level, x)
      x.mu <- pred.RNAmf_two_level$mu
      sig2 <- pred.RNAmf_two_level$sig2

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
      x <- scale.inputs(x, attr(fit3$X,"scaled:center")[1:d], attr(fit3$X,"scaled:scale")[1:d])
      x.mu <- scale.inputs(x.mu, attr(fit3$X,"scaled:center")[d+1], attr(fit3$X,"scaled:scale")[d+1])
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
      pred.RNAmf_two_level <- predRNAmf_two_level(fit.RNAmf_two_level, x)
      x.mu <- pred.RNAmf_two_level$mu
      sig2 <- pred.RNAmf_two_level$sig2

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
      x <- scale.inputs(x, attr(fit3$X,"scaled:center")[1:d], attr(fit3$X,"scaled:scale")[1:d])
      x.mu <- scale.inputs(x.mu, attr(fit3$X,"scaled:center")[d+1], attr(fit3$X,"scaled:scale")[d+1])
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
      pred.RNAmf_two_level <- predRNAmf_two_level(fit.RNAmf_two_level, x)
      x.mu <- pred.RNAmf_two_level$mu
      sig2 <- pred.RNAmf_two_level$sig2

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
      x <- scale.inputs(x, attr(fit3$X,"scaled:center")[1:d], attr(fit3$X,"scaled:scale")[1:d])
      x.mu <- scale.inputs(x.mu, attr(fit3$X,"scaled:center")[d+1], attr(fit3$X,"scaled:scale")[d+1])
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
      pred.RNAmf_two_level <- predRNAmf_two_level(fit.RNAmf_two_level, x)
      x.mu <- pred.RNAmf_two_level$mu
      sig2 <- pred.RNAmf_two_level$sig2

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
      x <- scale.inputs(x, attr(fit3$X,"scaled:center")[1:d], attr(fit3$X,"scaled:scale")[1:d])
      x.mu <- scale.inputs(x.mu, attr(fit3$X,"scaled:center")[d+1], attr(fit3$X,"scaled:scale")[d+1])
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

