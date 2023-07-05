#' integvar2
#'
#' computing integrated one-step ahead variance with three fidelity levels
#'
#' @param x vector or matrix of test data.
#' @param fit an object of class closed2.
#' @param mc.sample a number of mc samples generated for this approach. Default is 10000.
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

integvar2 <- function(x, fit, mc.sample=10000){

  intvar1 <- c(rep(0, dim(t(x))[2])) # IMSPE candidates
  intvar2 <- c(rep(0, dim(t(x))[2])) # IMSPE candidates
  intvar3 <- c(rep(0, dim(t(x))[2])) # IMSPE candidates

  fit1 <- fit$closed1$fit1
  fit2 <- fit$closed1$fit2
  fit3 <- fit$fit3
  constant <- fit$constant
  kernel <- fit$kernel
  x <- matrix(x, nrow=dim(t(x))[2])


  for(i in 1:length(intvar1)){
    ### Choose level 1 ###
    newx <- matrix(x[i,], nrow=1)

    if(kernel=="sqex"){
      x1.sample <- rnorm(mc.sample, mean=pred.GP(fit1, newx)$mu, sd=sqrt(pred.GP(fit1, newx)$sig2))
    }else if(kernel=="matern1.5"){
      x1.sample <- rnorm(mc.sample, mean=pred.matGP(fit1, newx)$mu, sd=sqrt(pred.matGP(fit1, newx)$sig2))
    }else if(kernel=="matern2.5"){
      x1.sample <- rnorm(mc.sample, mean=pred.matGP(fit1, newx)$mu, sd=sqrt(pred.matGP(fit1, newx)$sig2))
    }

    mu.cand1 <- mean(x1.sample) # f1(newx)

    x.center1 <- attr(fit1$X, "scaled:center")
    x.scale1 <- attr(fit1$X, "scaled:scale")
    y.center1 <- attr(fit1$y, "scaled:center")

    newx <- matrix((newx-attr(fit1$X,"scaled:center"))/attr(fit1$X,"scaled:scale"), nrow=1)

    ### update Ki
    if(kernel=="sqex"){
      v.next <- drop(covar.sep(X1=newx, d=fit1$theta, g=0) -
                       t(covar.sep(X1=fit1$X, X2=newx, d=fit1$theta, g=0)) %*%
                       fit1$Ki %*%
                       covar.sep(X1=fit1$X, X2=newx, d=fit1$theta, g=0))
      g.next <- - drop(solve(v.next)) * fit1$Ki %*% covar.sep(X1=fit1$X, X2=newx, d=fit1$theta, g=0)
    }else if(kernel=="matern1.5"){
      v.next <- drop(cor.sep(X=newx, theta=fit1$theta, nu=1.5) -
                       t(cor.sep(X=fit1$X, x=newx, theta=fit1$theta, nu=1.5)) %*%
                       fit1$Ki %*%
                       cor.sep(X=fit1$X, x=newx, theta=fit1$theta, nu=1.5))
      g.next <- - drop(solve(v.next)) * fit1$Ki %*% cor.sep(X=fit1$X, x=newx, theta=fit1$theta, nu=1.5)
    }else if(kernel=="matern2.5"){
      v.next <- drop(cor.sep(X=newx, theta=fit1$theta, nu=2.5) -
                       t(cor.sep(X=fit1$X, x=newx, theta=fit1$theta, nu=2.5)) %*%
                       fit1$Ki %*%
                       cor.sep(X=fit1$X, x=newx, theta=fit1$theta, nu=2.5))
      g.next <- - drop(solve(v.next)) * fit1$Ki %*% cor.sep(X=fit1$X, x=newx, theta=fit1$theta, nu=2.5)
    }

    fit1$Ki <- rbind(cbind(fit1$Ki+g.next%*%t(g.next)*v.next, g.next),
                     cbind(t(g.next), solve(v.next)))

    fit1$X <- rbind(fit1$X, newx)
    attr(fit1$X, "scaled:center") <- x.center1
    attr(fit1$X, "scaled:scale") <- x.scale1

    if(constant){
      fit1$y <- c(fit1$y, mu.cand1)
    }else{
      fit1$y <- c(fit1$y, mu.cand1-attr(fit1$y,"scaled:center"))
      attr(fit1$y, "scaled:center") <- y.center1
    }

    fit1$tau2hat <- drop(t(fit1$y) %*% fit1$Ki %*% fit1$y / length(fit1$y))

    fit$closed1$fit1 <- fit1

    intvar1[i] <- mean(predclosed2(fit, x)$sig2)
  }


  for(i in 1:length(intvar2)){
    ### Choose level 2 ###
    newx <- matrix(x[i,], nrow=1)

    if(kernel=="sqex"){
      x1.sample <- rnorm(mc.sample, mean=pred.GP(fit1, newx)$mu, sd=sqrt(pred.GP(fit1, newx)$sig2))
      mu.cand1 <- mean(x1.sample) # f1(newx)

      x2.sample <- rnorm(mc.sample,
                         mean=pred.GP(fit2, cbind(newx, mu.cand1))$mu,
                         sd=sqrt(pred.GP(fit2, cbind(newx, mu.cand1))$sig2))
    }else if(kernel=="matern1.5"){
      x1.sample <- rnorm(mc.sample, mean=pred.matGP(fit1, newx)$mu, sd=sqrt(pred.matGP(fit1, newx)$sig2))
      mu.cand1 <- mean(x1.sample) # f1(newx)

      x2.sample <- rnorm(mc.sample,
                         mean=pred.matGP(fit2, cbind(newx, mu.cand1))$mu,
                         sd=sqrt(pred.matGP(fit2, cbind(newx, mu.cand1))$sig2))
    }else if(kernel=="matern2.5"){
      x1.sample <- rnorm(mc.sample, mean=pred.matGP(fit1, newx)$mu, sd=sqrt(pred.matGP(fit1, newx)$sig2))
      mu.cand1 <- mean(x1.sample) # f1(newx)

      x2.sample <- rnorm(mc.sample,
                         mean=pred.matGP(fit2, cbind(newx, mu.cand1))$mu,
                         sd=sqrt(pred.matGP(fit2, cbind(newx, mu.cand1))$sig2))
    }

    mu.cand2 <- mean(x2.sample) # f2(newx)

    x.center1 <- attr(fit1$X, "scaled:center")
    x.scale1 <- attr(fit1$X, "scaled:scale")
    y.center1 <- attr(fit1$y, "scaled:center")

    x.center2 <- attr(fit2$X, "scaled:center")
    x.scale2 <- attr(fit2$X, "scaled:scale")
    y.center2 <- attr(fit2$y, "scaled:center")

    newx1 <- matrix((newx-attr(fit1$X,"scaled:center"))/attr(fit1$X,"scaled:scale"), nrow=1)
    newx2 <- t((t(cbind(newx, mean(x1.sample)))-attr(fit2$X,"scaled:center"))/attr(fit2$X,"scaled:scale"))

    ### update Ki1
    if(kernel=="sqex"){
      v.next1 <- drop(covar.sep(X1=newx1, d=fit1$theta, g=0) -
                        t(covar.sep(X1=fit1$X, X2=newx1, d=fit1$theta, g=0)) %*%
                        fit1$Ki %*%
                        covar.sep(X1=fit1$X, X2=newx1, d=fit1$theta, g=0))
      g.next1 <- - drop(solve(v.next1)) * fit1$Ki %*% covar.sep(X1=fit1$X, X2=newx1, d=fit1$theta, g=0)
    }else if(kernel=="matern1.5"){
      v.next1 <- drop(cor.sep(X=newx1, theta=fit1$theta, nu=1.5) -
                        t(cor.sep(X=fit1$X, x=newx1, theta=fit1$theta, nu=1.5)) %*%
                        fit1$Ki %*%
                        cor.sep(X=fit1$X, x=newx1, theta=fit1$theta, nu=1.5))
      g.next1 <- - drop(solve(v.next1)) * fit1$Ki %*% cor.sep(X=fit1$X, x=newx1, theta=fit1$theta, nu=1.5)
    }else if(kernel=="matern2.5"){
      v.next1 <- drop(cor.sep(X=newx1, theta=fit1$theta, nu=2.5) -
                        t(cor.sep(X=fit1$X, x=newx1, theta=fit1$theta, nu=2.5)) %*%
                        fit1$Ki %*%
                        cor.sep(X=fit1$X, x=newx1, theta=fit1$theta, nu=2.5))
      g.next1 <- - drop(solve(v.next1)) * fit1$Ki %*% cor.sep(X=fit1$X, x=newx1, theta=fit1$theta, nu=2.5)
    }

    fit1$Ki <- rbind(cbind(fit1$Ki+g.next1%*%t(g.next1)*v.next1, g.next1),
                     cbind(t(g.next1), solve(v.next1)))

    fit1$X <- rbind(fit1$X, newx1)
    attr(fit1$X, "scaled:center") <- x.center1
    attr(fit1$X, "scaled:scale") <- x.scale1

    if(constant){
      fit1$y <- c(fit1$y, mu.cand1)
    }else{
      fit1$y <- c(fit1$y, mu.cand1-attr(fit1$y,"scaled:center"))
      attr(fit1$y, "scaled:center") <- y.center1
    }

    fit1$tau2hat <- drop(t(fit1$y) %*% fit1$Ki %*% fit1$y / length(fit1$y))

    ### update Ki2
    if(kernel=="sqex"){
      v.next2 <- drop(covar.sep(X1=newx2, d=fit2$theta, g=0) -
                        t(covar.sep(X1=fit2$X, X2=newx2, d=fit2$theta, g=0)) %*%
                        fit2$Ki %*%
                        covar.sep(X1=fit2$X, X2=newx2, d=fit2$theta, g=0))
      g.next2 <- - drop(solve(v.next2)) * fit2$Ki %*% covar.sep(X1=fit2$X, X2=newx2, d=fit2$theta, g=0)
    }else if(kernel=="matern1.5"){
      v.next2 <- drop(cor.sep(X=newx2, theta=fit2$theta, nu=1.5) -
                        t(cor.sep(X=fit2$X, x=newx2, theta=fit2$theta, nu=1.5)) %*%
                        fit2$Ki %*%
                        cor.sep(X=fit2$X, x=newx2, theta=fit2$theta, nu=1.5))
      g.next2 <- - drop(solve(v.next2)) * fit2$Ki %*% cor.sep(X=fit2$X, x=newx2, theta=fit2$theta, nu=1.5)
    }else if(kernel=="matern2.5"){
      v.next2 <- drop(cor.sep(X=newx2, theta=fit2$theta, nu=2.5) -
                        t(cor.sep(X=fit2$X, x=newx2, theta=fit2$theta, nu=2.5)) %*%
                        fit2$Ki %*%
                        cor.sep(X=fit2$X, x=newx2, theta=fit2$theta, nu=2.5))
      g.next2 <- - drop(solve(v.next2)) * fit2$Ki %*% cor.sep(X=fit2$X, x=newx2, theta=fit2$theta, nu=2.5)
    }

    fit2$Ki <- rbind(cbind(fit2$Ki+g.next2%*%t(g.next2)*v.next2, g.next2),
                     cbind(t(g.next2), solve(v.next2)))

    fit2$X <- rbind(fit2$X, newx2)
    attr(fit2$X, "scaled:center") <- x.center2
    attr(fit2$X, "scaled:scale") <- x.scale2

    if(constant){
      fit2$y <- c(fit2$y, mu.cand2)
    }else{
      fit2$y <- c(fit2$y, mu.cand2-attr(fit2$y,"scaled:center"))
      attr(fit2$y, "scaled:center") <- y.center2
    }

    fit2$tau2hat <- drop(t(fit2$y) %*% fit2$Ki %*% fit2$y / length(fit2$y))

    fit$closed1$fit1 <- fit1
    fit$closed1$fit2 <- fit2

    intvar2[i] <- mean(predclosed2(fit, x)$sig2)
  }


  for(i in 1:length(intvar3)){
    ### Choose level 3 ###
    newx <- matrix(x[i,], nrow=1)

    if(kernel=="sqex"){
      x1.sample <- rnorm(mc.sample, mean=pred.GP(fit1, newx)$mu, sd=sqrt(pred.GP(fit1, newx)$sig2))
      mu.cand1 <- mean(x1.sample) # f1(newx)

      x2.sample <- rnorm(mc.sample,
                         mean=pred.GP(fit2, cbind(newx, mu.cand1))$mu,
                         sd=sqrt(pred.GP(fit2, cbind(newx, mu.cand1))$sig2))
      mu.cand2 <- mean(x2.sample) # f2(newx)

      x3.sample <- rnorm(mc.sample,
                         mean=pred.GP(fit3, cbind(newx, mu.cand2))$mu,
                         sd=sqrt(pred.GP(fit3, cbind(newx, mu.cand2))$sig2))
    }else if(kernel=="matern1.5"){
      x1.sample <- rnorm(mc.sample, mean=pred.matGP(fit1, newx)$mu, sd=sqrt(pred.matGP(fit1, newx)$sig2))
      mu.cand1 <- mean(x1.sample) # f1(newx)

      x2.sample <- rnorm(mc.sample,
                         mean=pred.matGP(fit2, cbind(newx, mu.cand1))$mu,
                         sd=sqrt(pred.matGP(fit2, cbind(newx, mu.cand1))$sig2))
      mu.cand2 <- mean(x2.sample) # f2(newx)

      x3.sample <- rnorm(mc.sample,
                         mean=pred.matGP(fit3, cbind(newx, mu.cand2))$mu,
                         sd=sqrt(pred.matGP(fit3, cbind(newx, mu.cand2))$sig2))
    }else if(kernel=="matern2.5"){
      x1.sample <- rnorm(mc.sample, mean=pred.matGP(fit1, newx)$mu, sd=sqrt(pred.matGP(fit1, newx)$sig2))
      mu.cand1 <- mean(x1.sample) # f1(newx)

      x2.sample <- rnorm(mc.sample,
                         mean=pred.matGP(fit2, cbind(newx, mu.cand1))$mu,
                         sd=sqrt(pred.matGP(fit2, cbind(newx, mu.cand1))$sig2))
      mu.cand2 <- mean(x2.sample) # f2(newx)

      x3.sample <- rnorm(mc.sample,
                         mean=pred.matGP(fit3, cbind(newx, mu.cand2))$mu,
                         sd=sqrt(pred.matGP(fit3, cbind(newx, mu.cand2))$sig2))
    }

    mu.cand3 <- mean(x3.sample) # f2(newx)

    x.center1 <- attr(fit1$X, "scaled:center")
    x.scale1 <- attr(fit1$X, "scaled:scale")
    y.center1 <- attr(fit1$y, "scaled:center")

    x.center2 <- attr(fit2$X, "scaled:center")
    x.scale2 <- attr(fit2$X, "scaled:scale")
    y.center2 <- attr(fit2$y, "scaled:center")

    x.center3 <- attr(fit3$X, "scaled:center")
    x.scale3 <- attr(fit3$X, "scaled:scale")
    y.center3 <- attr(fit3$y, "scaled:center")

    newx1 <- matrix((newx-attr(fit1$X,"scaled:center"))/attr(fit1$X,"scaled:scale"), nrow=1)
    newx2 <- t((t(cbind(newx, mean(x1.sample)))-attr(fit2$X,"scaled:center"))/attr(fit2$X,"scaled:scale"))
    newx3 <- t((t(cbind(newx, mean(x2.sample)))-attr(fit3$X,"scaled:center"))/attr(fit3$X,"scaled:scale"))

    ### update Ki1
    if(kernel=="sqex"){
      v.next1 <- drop(covar.sep(X1=newx1, d=fit1$theta, g=0) -
                        t(covar.sep(X1=fit1$X, X2=newx1, d=fit1$theta, g=0)) %*%
                        fit1$Ki %*%
                        covar.sep(X1=fit1$X, X2=newx1, d=fit1$theta, g=0))
      g.next1 <- - drop(solve(v.next1)) * fit1$Ki %*% covar.sep(X1=fit1$X, X2=newx1, d=fit1$theta, g=0)
    }else if(kernel=="matern1.5"){
      v.next1 <- drop(cor.sep(X=newx1, theta=fit1$theta, nu=1.5) -
                        t(cor.sep(X=fit1$X, x=newx1, theta=fit1$theta, nu=1.5)) %*%
                        fit1$Ki %*%
                        cor.sep(X=fit1$X, x=newx1, theta=fit1$theta, nu=1.5))
      g.next1 <- - drop(solve(v.next1)) * fit1$Ki %*% cor.sep(X=fit1$X, x=newx1, theta=fit1$theta, nu=1.5)
    }else if(kernel=="matern2.5"){
      v.next1 <- drop(cor.sep(X=newx1, theta=fit1$theta, nu=2.5) -
                        t(cor.sep(X=fit1$X, x=newx1, theta=fit1$theta, nu=2.5)) %*%
                        fit1$Ki %*%
                        cor.sep(X=fit1$X, x=newx1, theta=fit1$theta, nu=2.5))
      g.next1 <- - drop(solve(v.next1)) * fit1$Ki %*% cor.sep(X=fit1$X, x=newx1, theta=fit1$theta, nu=2.5)
    }

    fit1$Ki <- rbind(cbind(fit1$Ki+g.next1%*%t(g.next1)*v.next1, g.next1),
                     cbind(t(g.next1), solve(v.next1)))

    fit1$X <- rbind(fit1$X, newx1)
    attr(fit1$X, "scaled:center") <- x.center1
    attr(fit1$X, "scaled:scale") <- x.scale1

    if(constant){
      fit1$y <- c(fit1$y, mu.cand1)
    }else{
      fit1$y <- c(fit1$y, mu.cand1-attr(fit1$y,"scaled:center"))
      attr(fit1$y, "scaled:center") <- y.center1
    }

    fit1$tau2hat <- drop(t(fit1$y) %*% fit1$Ki %*% fit1$y / length(fit1$y))

    ### update Ki2
    if(kernel=="sqex"){
      v.next2 <- drop(covar.sep(X1=newx2, d=fit2$theta, g=0) -
                        t(covar.sep(X1=fit2$X, X2=newx2, d=fit2$theta, g=0)) %*%
                        fit2$Ki %*%
                        covar.sep(X1=fit2$X, X2=newx2, d=fit2$theta, g=0))
      g.next2 <- - drop(solve(v.next2)) * fit2$Ki %*% covar.sep(X1=fit2$X, X2=newx2, d=fit2$theta, g=0)
    }else if(kernel=="matern1.5"){
      v.next2 <- drop(cor.sep(X=newx2, theta=fit2$theta, nu=1.5) -
                        t(cor.sep(X=fit2$X, x=newx2, theta=fit2$theta, nu=1.5)) %*%
                        fit2$Ki %*%
                        cor.sep(X=fit2$X, x=newx2, theta=fit2$theta, nu=1.5))
      g.next2 <- - drop(solve(v.next2)) * fit2$Ki %*% cor.sep(X=fit2$X, x=newx2, theta=fit2$theta, nu=1.5)
    }else if(kernel=="matern2.5"){
      v.next2 <- drop(cor.sep(X=newx2, theta=fit2$theta, nu=2.5) -
                        t(cor.sep(X=fit2$X, x=newx2, theta=fit2$theta, nu=2.5)) %*%
                        fit2$Ki %*%
                        cor.sep(X=fit2$X, x=newx2, theta=fit2$theta, nu=2.5))
      g.next2 <- - drop(solve(v.next2)) * fit2$Ki %*% cor.sep(X=fit2$X, x=newx2, theta=fit2$theta, nu=2.5)
    }

    fit2$Ki <- rbind(cbind(fit2$Ki+g.next2%*%t(g.next2)*v.next2, g.next2),
                     cbind(t(g.next2), solve(v.next2)))

    fit2$X <- rbind(fit2$X, newx2)
    attr(fit2$X, "scaled:center") <- x.center2
    attr(fit2$X, "scaled:scale") <- x.scale2

    if(constant){
      fit2$y <- c(fit2$y, mu.cand2)
    }else{
      fit2$y <- c(fit2$y, mu.cand2-attr(fit2$y,"scaled:center"))
      attr(fit2$y, "scaled:center") <- y.center2
    }

    fit2$tau2hat <- drop(t(fit2$y) %*% fit2$Ki %*% fit2$y / length(fit2$y))

    ### update Ki3
    if(kernel=="sqex"){
      v.next3 <- drop(covar.sep(X1=newx3, d=fit3$theta, g=0) -
                        t(covar.sep(X1=fit3$X, X2=newx3, d=fit3$theta, g=0)) %*%
                        fit3$Ki %*%
                        covar.sep(X1=fit3$X, X2=newx3, d=fit3$theta, g=0))
      g.next3 <- - drop(solve(v.next3)) * fit3$Ki %*% covar.sep(X1=fit3$X, X2=newx3, d=fit3$theta, g=0)
    }else if(kernel=="matern1.5"){
      v.next3 <- drop(cor.sep(X=newx3, theta=fit3$theta, nu=1.5) -
                        t(cor.sep(X=fit3$X, x=newx3, theta=fit3$theta, nu=1.5)) %*%
                        fit3$Ki %*%
                        cor.sep(X=fit3$X, x=newx3, theta=fit3$theta, nu=1.5))
      g.next3 <- - drop(solve(v.next3)) * fit3$Ki %*% cor.sep(X=fit3$X, x=newx3, theta=fit3$theta, nu=1.5)
    }else if(kernel=="matern2.5"){
      v.next3 <- drop(cor.sep(X=newx3, theta=fit3$theta, nu=2.5) -
                        t(cor.sep(X=fit3$X, x=newx3, theta=fit3$theta, nu=2.5)) %*%
                        fit3$Ki %*%
                        cor.sep(X=fit3$X, x=newx3, theta=fit3$theta, nu=2.5))
      g.next3 <- - drop(solve(v.next3)) * fit3$Ki %*% cor.sep(X=fit3$X, x=newx3, theta=fit3$theta, nu=2.5)
    }

    fit3$Ki <- rbind(cbind(fit3$Ki+g.next3%*%t(g.next3)*v.next3, g.next3),
                     cbind(t(g.next3), solve(v.next3)))

    fit3$X <- rbind(fit3$X, newx3)
    attr(fit3$X, "scaled:center") <- x.center3
    attr(fit3$X, "scaled:scale") <- x.scale3

    if(constant){
      fit3$y <- rbind(fit3$y, mu.cand3)
    }else{
      fit3$y <- rbind(fit3$y, mu.cand3-attr(fit3$y,"scaled:center"))
      attr(fit3$y, "scaled:center") <- y.center3
    }

    fit3$tau2hat <- drop(t(fit3$y) %*% fit3$Ki %*% fit3$y / length(fit3$y))


    fit$closed1$fit1 <- fit1
    fit$closed1$fit2 <- fit2
    fit$fit3 <- fit3

    intvar3[i] <- mean(predclosed2(fit, x)$sig2)
  }

  return(list(intvar1=intvar1, intvar2=intvar2, intvar3=intvar3))
}
