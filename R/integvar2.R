#' integvar2
#'
#' computing integrated one-step ahead variance with three fidelity levels
#'
#' @param x vector or matrix of test data.
#' @param fit an object of class RNAmf2.
#' @param mc.sample a number of mc samples generated for this approach. Default is 10.
#' @return A list containing the integrated one-step ahead variance:
#' \itemize{
#'   \item \code{intvar1}: vector of integrated variance when each data point is added at level 1.
#'   \item \code{intvar2}: vector of integrated variance when each data point is added at level 2.
#'   \item \code{intvar3}: vector of integrated variance when each data point is added at level 3.
#' }
#' @importFrom plgp covar.sep
#' @importFrom stats rnorm
#' @export
#'

integvar2 <- function(x, fit, mc.sample=10){

  intvar1 <- c(rep(0, dim(t(x))[2])) # IMSPE candidates
  intvar2 <- c(rep(0, dim(t(x))[2])) # IMSPE candidates
  intvar3 <- c(rep(0, dim(t(x))[2])) # IMSPE candidates
  pseudointvar1 <- c(rep(0, mc.sample))
  pseudointvar2 <- c(rep(0, mc.sample))
  pseudointvar3 <- c(rep(0, mc.sample))

  fit1 <- f1 <- fit$fit.RNAmf1$fit1
  fit2 <- f2 <- fit$fit.RNAmf1$fit2
  fit3 <- f3 <- fit$fit3
  constant <- fit$constant
  kernel <- fit$kernel
  g <- fit1$g
  x <- matrix(x, nrow=dim(t(x))[2])

  x.center1 <- attr(fit1$X, "scaled:center")
  x.scale1 <- attr(fit1$X, "scaled:scale")
  y.center1 <- attr(fit1$y, "scaled:center")

  x.center2 <- attr(fit2$X, "scaled:center")
  x.scale2 <- attr(fit2$X, "scaled:scale")
  y.center2 <- attr(fit2$y, "scaled:center")

  x.center3 <- attr(fit3$X, "scaled:center")
  x.scale3 <- attr(fit3$X, "scaled:scale")
  y.center3 <- attr(fit3$y, "scaled:center")


  for(i in 1:length(intvar1)){

    newx <- matrix(x[i,], nrow=1)

    if(kernel=="sqex"){
      x1.sample <- rnorm(mc.sample, mean=pred.GP(f1, newx)$mu, sd=sqrt(pred.GP(f1, newx)$sig2))
    }else if(kernel=="matern1.5"){
      x1.sample <- rnorm(mc.sample, mean=pred.matGP(f1, newx)$mu, sd=sqrt(pred.matGP(f1, newx)$sig2))
    }else if(kernel=="matern2.5"){
      x1.sample <- rnorm(mc.sample, mean=pred.matGP(f1, newx)$mu, sd=sqrt(pred.matGP(f1, newx)$sig2))
    }

    newx1 <- matrix((newx-x.center1)/x.scale1, nrow=1)


    ### Choose level 1 ###
    ### update Ki1
    if(kernel=="sqex"){
      cov.newx1 <- covar.sep(X1=newx1, d=f1$theta, g=g)
      cov.Xnewx1 <- covar.sep(X1=f1$X, X2=newx1, d=f1$theta, g=0)
    }else if(kernel=="matern1.5"){
      cov.newx1 <- cor.sep(X=newx1, theta=f1$theta, nu=1.5)
      cov.Xnewx1 <- cor.sep(X=f1$X, x=newx1, theta=f1$theta, nu=1.5)
    }else if(kernel=="matern2.5"){
      cov.newx1 <- cor.sep(X=newx1, theta=f1$theta, nu=2.5)
      cov.Xnewx1 <- cor.sep(X=f1$X, x=newx1, theta=f1$theta, nu=2.5)
    }

    v.next1 <- drop(cov.newx1 - t(cov.Xnewx1) %*% f1$Ki %*% cov.Xnewx1)
    g.next1 <- - 1/drop(v.next1) * f1$Ki %*% cov.Xnewx1

    for(j in 1:mc.sample){

      fit1$Ki <- rbind(cbind(f1$Ki+g.next1%*%t(g.next1)*v.next1, g.next1),
                       cbind(t(g.next1), 1/drop(v.next1)))

      fit1$X <- rbind(f1$X, newx1)
      attr(fit1$X, "scaled:center") <- x.center1
      attr(fit1$X, "scaled:scale") <- x.scale1

      if(constant){
        fit1$y <- c(f1$y, x1.sample[j])
      }else{
        fit1$y <- c(f1$y, x1.sample[j]-y.center1)
        attr(fit1$y, "scaled:center") <- y.center1
      }

      fit1$tau2hat <- drop(t(fit1$y - fit1$mu.hat) %*% fit1$Ki %*% (fit1$y - fit1$mu.hat) / length(fit1$y))

      fit$fit.RNAmf1$fit1 <- fit1

      pseudointvar1[j] <- mean(predRNAmf2(fit, x)$sig2)


      ### Choose level 2 ###
      if(kernel=="sqex"){
        x2.sample <- pred.GP(fit2, cbind(newx, x1.sample[j]))$mu
      }else if(kernel=="matern1.5"){
        x2.sample <- pred.matGP(fit2, cbind(newx, x1.sample[j]))$mu
      }else if(kernel=="matern2.5"){
        x2.sample <- pred.matGP(fit2, cbind(newx, x1.sample[j]))$mu
      }


      ### update Ki2
      newx2 <- t((t(cbind(newx, x1.sample[j]))-x.center2)/x.scale2)

      if(kernel=="sqex"){
        cov.newx2 <- covar.sep(X1=newx2, d=f2$theta, g=g)
        cov.Xnewx2 <- covar.sep(X1=f2$X, X2=newx2, d=f2$theta, g=0)
      }else if(kernel=="matern1.5"){
        cov.newx2 <- cor.sep(X=newx2, theta=f2$theta, nu=1.5)
        cov.Xnewx2 <- cor.sep(X=f2$X, x=newx2, theta=f2$theta, nu=1.5)
      }else if(kernel=="matern2.5"){
        cov.newx2 <- cor.sep(X=newx2, theta=f2$theta, nu=2.5)
        cov.Xnewx2 <- cor.sep(X=f2$X, x=newx2, theta=f2$theta, nu=2.5)
      }

      v.next2 <- drop(cov.newx2 - t(cov.Xnewx2) %*% f2$Ki %*% cov.Xnewx2)
      g.next2 <- - 1/drop(v.next2) * f2$Ki %*% cov.Xnewx2

      fit2$Ki <- rbind(cbind(f2$Ki+g.next2%*%t(g.next2)*v.next2, g.next2),
                       cbind(t(g.next2), 1/drop(v.next2)))

      fit2$X <- rbind(f2$X, newx2)
      attr(fit2$X, "scaled:center") <- x.center2
      attr(fit2$X, "scaled:scale") <- x.scale2

      if(constant){
        fit2$y <- c(f2$y, x2.sample)
      }else{
        fit2$y <- c(f2$y, x2.sample-y.center2)
        attr(fit2$y, "scaled:center") <- y.center2
      }

      fit2$tau2hat <- drop(t(fit2$y - fit2$mu.hat) %*% fit2$Ki %*% (fit2$y - fit2$mu.hat) / length(fit2$y))

      fit$fit.RNAmf1$fit2 <- fit2

      pseudointvar2[j] <- mean(predRNAmf2(fit, x)$sig2)


      ### Choose level 3 ###
      if(kernel=="sqex"){
        x3.sample <- pred.GP(fit3, cbind(newx, x2.sample))$mu
      }else if(kernel=="matern1.5"){
        x3.sample <- pred.matGP(fit3, cbind(newx, x2.sample))$mu
      }else if(kernel=="matern2.5"){
        x3.sample <- pred.matGP(fit3, cbind(newx, x2.sample))$mu
      }


      ### update Ki3
      newx3 <- t((t(cbind(newx, x2.sample))-x.center3)/x.scale3)

      if(kernel=="sqex"){
        cov.newx3 <- covar.sep(X1=newx3, d=f3$theta, g=g)
        cov.Xnewx3 <- covar.sep(X1=f3$X, X2=newx3, d=f3$theta, g=0)
      }else if(kernel=="matern1.5"){
        cov.newx3 <- cor.sep(X=newx3, theta=f3$theta, nu=1.5)
        cov.Xnewx3 <- cor.sep(X=f3$X, x=newx3, theta=f3$theta, nu=1.5)
      }else if(kernel=="matern2.5"){
        cov.newx3 <- cor.sep(X=newx3, theta=f3$theta, nu=2.5)
        cov.Xnewx3 <- cor.sep(X=f3$X, x=newx3, theta=f3$theta, nu=2.5)
      }

      v.next3 <- drop(cov.newx3 - t(cov.Xnewx3) %*% f3$Ki %*% cov.Xnewx3)
      g.next3 <- - 1/drop(v.next3) * f3$Ki %*% cov.Xnewx3

      fit3$Ki <- rbind(cbind(f3$Ki+g.next3%*%t(g.next3)*v.next3, g.next3),
                       cbind(t(g.next3), 1/drop(v.next3)))

      fit3$X <- rbind(f3$X, newx3)
      attr(fit3$X, "scaled:center") <- x.center3
      attr(fit3$X, "scaled:scale") <- x.scale3

      if(constant){
        fit3$y <- c(f3$y, x3.sample)
      }else{
        fit3$y <- c(f3$y, x3.sample-y.center3)
        attr(fit3$y, "scaled:center") <- y.center3
      }

      fit3$tau2hat <- drop(t(fit3$y - fit3$mu.hat) %*% fit3$Ki %*% (fit3$y - fit3$mu.hat) / length(fit3$y))

      fit$fit3 <- fit3

      pseudointvar3[j] <- mean(predRNAmf2(fit, x)$sig2)

      fit$fit.RNAmf1$fit1 <- fit1 <- f1
      fit$fit.RNAmf1$fit2 <- fit2 <- f2
      fit$fit3 <- fit3 <- f3
    }

    intvar1[i] <- mean(pseudointvar1)
    intvar2[i] <- mean(pseudointvar2)
    intvar3[i] <- mean(pseudointvar3)
  }

  return(list(intvar1=intvar1, intvar2=intvar2, intvar3=intvar3))
}
