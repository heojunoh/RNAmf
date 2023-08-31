#' obj.ALC_two_level_1
#'
#' object to optimize the next point by ALC criterion updating at level 1 with two levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param Xref vector or matrix of reference data.
#' @param fit an object of class RNAmf.
#' @return A mean of the predictive posterior variance at Xref.
#' @importFrom plgp covar.sep
#' @importFrom stats rnorm
#' @importFrom maximin maximin
#' @export
#'

obj.ALC_two_level_1 <- function(Xcand, Xref, fit){

  fit1 <- f1 <- fit$fit1
  fit2 <- f2 <- fit$fit2
  constant <- fit$constant
  kernel <- fit$kernel
  g <- fit1$g

  x.center1 <- attr(fit1$X, "scaled:center")
  x.scale1 <- attr(fit1$X, "scaled:scale")
  y.center1 <- attr(fit1$y, "scaled:center")

  x.center2 <- attr(fit2$X, "scaled:center")
  x.scale2 <- attr(fit2$X, "scaled:scale")
  y.center2 <- attr(fit2$y, "scaled:center")


  ### update with new point ###
  newx <- matrix(Xcand, nrow=1)

  if(kernel=="sqex"){
    x1.sample <- rnorm(1, mean=pred.GP(f1, newx)$mu, sd=sqrt(pred.GP(f1, newx)$sig2))
  }else if(kernel=="matern1.5"){
    x1.sample <- rnorm(1, mean=pred.matGP(f1, newx)$mu, sd=sqrt(pred.matGP(f1, newx)$sig2))
  }else if(kernel=="matern2.5"){
    x1.sample <- rnorm(1, mean=pred.matGP(f1, newx)$mu, sd=sqrt(pred.matGP(f1, newx)$sig2))
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


  fit1$Ki <- rbind(cbind(f1$Ki+g.next1%*%t(g.next1)*v.next1, g.next1),
                   cbind(t(g.next1), 1/drop(v.next1)))

  fit1$X <- rbind(f1$X, newx1)
  attr(fit1$X, "scaled:center") <- x.center1
  attr(fit1$X, "scaled:scale") <- x.scale1

  if(constant){
    fit1$y <- c(f1$y, x1.sample)
  }else{
    fit1$y <- c(f1$y, x1.sample-y.center1)
    attr(fit1$y, "scaled:center") <- y.center1
  }

  fit1$tau2hat <- drop(t(fit1$y - fit1$mu.hat) %*% fit1$Ki %*% (fit1$y - fit1$mu.hat) / length(fit1$y))

  fit$fit1 <- fit1

  mean(predRNAmf(fit, Xref)$sig2) # to minimize the deduced variance. To maximize, -mean
}


#' obj.ALC_two_level_2
#'
#' object to optimize the next point by ALC criterion updating at level 2 with two levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param Xref vector or matrix of reference data.
#' @param fit an object of class RNAmf.
#' @return A mean of the predictive posterior variance at Xref.
#' @importFrom plgp covar.sep
#' @importFrom stats rnorm
#' @importFrom maximin maximin
#' @export
#'

obj.ALC_two_level_2 <- function(Xcand, Xref, fit){

  fit1 <- f1 <- fit$fit1
  fit2 <- f2 <- fit$fit2
  constant <- fit$constant
  kernel <- fit$kernel
  g <- fit1$g

  x.center1 <- attr(fit1$X, "scaled:center")
  x.scale1 <- attr(fit1$X, "scaled:scale")
  y.center1 <- attr(fit1$y, "scaled:center")

  x.center2 <- attr(fit2$X, "scaled:center")
  x.scale2 <- attr(fit2$X, "scaled:scale")
  y.center2 <- attr(fit2$y, "scaled:center")


  ### update with new point ###
  newx <- matrix(Xcand, nrow=1)

  if(kernel=="sqex"){
    x1.sample <- rnorm(1, mean=pred.GP(f1, newx)$mu, sd=sqrt(pred.GP(f1, newx)$sig2))
  }else if(kernel=="matern1.5"){
    x1.sample <- rnorm(1, mean=pred.matGP(f1, newx)$mu, sd=sqrt(pred.matGP(f1, newx)$sig2))
  }else if(kernel=="matern2.5"){
    x1.sample <- rnorm(1, mean=pred.matGP(f1, newx)$mu, sd=sqrt(pred.matGP(f1, newx)$sig2))
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


  fit1$Ki <- rbind(cbind(f1$Ki+g.next1%*%t(g.next1)*v.next1, g.next1),
                   cbind(t(g.next1), 1/drop(v.next1)))

  fit1$X <- rbind(f1$X, newx1)
  attr(fit1$X, "scaled:center") <- x.center1
  attr(fit1$X, "scaled:scale") <- x.scale1

  if(constant){
    fit1$y <- c(f1$y, x1.sample)
  }else{
    fit1$y <- c(f1$y, x1.sample-y.center1)
    attr(fit1$y, "scaled:center") <- y.center1
  }

  fit1$tau2hat <- drop(t(fit1$y - fit1$mu.hat) %*% fit1$Ki %*% (fit1$y - fit1$mu.hat) / length(fit1$y))

  fit$fit1 <- fit1


  ### Choose level 2 ###
  if(kernel=="sqex"){
    pred2 <- pred.GP(fit2, cbind(newx, x1.sample))
    x2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
  }else if(kernel=="matern1.5"){
    pred2 <- pred.matGP(fit2, cbind(newx, x1.sample))
    x2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
  }else if(kernel=="matern2.5"){
    pred2 <- pred.matGP(fit2, cbind(newx, x1.sample))
    x2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
  }

  ### update Ki2
  newx2 <- t((t(cbind(newx, x1.sample))-x.center2)/x.scale2)

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

  fit$fit2 <- fit2

  mean(predRNAmf(fit, Xref)$sig2) # to minimize the deduced variance. To maximize, -mean
}


#' ALC_two_level_optm
#'
#' find the next point by ALC criterion with two fidelity levels and update the model
#'
#' @param Xref vector or matrix of reference data.
#' @param fit an object of class RNAmf.
#' @param cost a vector of the costs for each level of fidelity.
#' @param funcs list of functions for each level of fidelity.
#' @return A list containing the integrated one-step ahead variance:
#' \itemize{
#'   \item \code{fit}: fitted model after acquire chosen point.
#'   \item \code{intvar1}: vector of integrated variance when each data point is added at level 1.
#'   \item \code{intvar2}: vector of integrated variance when each data point is added at level 2.
#'   \item \code{cost}: vector of the costs for each level of fidelity.
#'   \item \code{Xcand}: candidate data points considered to be acquired.
#'   \item \code{chosen}: list of chosen level, location, and point.
#' }
#' @importFrom plgp covar.sep
#' @importFrom stats rnorm
#' @importFrom maximin maximin
#' @export
#'

ALC_two_level_optm <- function(Xref, fit, cost, funcs){

  if (length(cost)!=2) stop("The length of cost should be 2")
  if (cost[1] >= cost[2]) stop("If the cost for high-fidelity is cheaper, just acquire the high-fidelity")


  Icurrent <- mean(predRNAmf(fit, Xref)$sig2)

  fit1 <- f1 <- fit$fit1
  fit2 <- f2 <- fit$fit2
  constant <- fit$constant
  kernel <- fit$kernel
  g <- fit1$g

  x.center1 <- attr(fit1$X, "scaled:center")
  x.scale1 <- attr(fit1$X, "scaled:scale")
  y.center1 <- attr(fit1$y, "scaled:center")

  x.center2 <- attr(fit2$X, "scaled:center")
  x.scale2 <- attr(fit2$X, "scaled:scale")
  y.center2 <- attr(fit2$y, "scaled:center")

  if(ncol(fit1$X)==1){ # Xcand
    Xcand <- maximin.1d(Xorig=t(t(fit1$X) * x.scale1 + x.center1))
  }else{
    Xcand <- maximin(nrow(fit1$X), ncol(fit1$X), T=10*nrow(fit1$X),
                     Xorig=t(t(fit1$X) * x.scale1 + x.center1))$Xf
  }

  xoptm1 <- xoptm2 <- matrix(NA, nrow=nrow(Xcand), ncol=ncol(Xcand))
  value1 <- value2 <- rep(NA, nrow(Xcand))

  ### Optimize at level 1 ###
  for(i in 1:nrow(Xcand)){
    out <- optim(Xcand[i,], obj.ALC_two_level_1, method="L-BFGS-B", lower=0, upper=1, fit=fit, Xref=Xref)
    xoptm1[i,] <- out$par
    value1[i] <- out$value
  }

  ### Find the optimal point at level 1 ###
  xnext1 <- xoptm1[which.min(value1),]


  ### Optimize at level 2 ###
  for(i in 1:nrow(Xcand)){
    out <- optim(Xcand[i,], obj.ALC_two_level_2, method="L-BFGS-B", lower=0, upper=1, fit=fit, Xref=Xref)
    xoptm2[i,] <- out$par
    value2[i] <- out$value
  }

  ### Find the optimal point at level 2 ###
  xnext2 <- xoptm2[which.min(value2),]


  ### Select the next level and point ###
  ALCvalue <- c(Icurrent - min(value1), Icurrent - min(value2))/c(cost[1], cost[1]+cost[2])

  chosen <- list("level"=which.max(ALCvalue), # next level
                 "location"=c(which.min(value1), which.min(value2))[which.max(ALCvalue)], # next location
                 "Xnext"=Xcand[c(which.min(value1), which.min(value2))[which.max(ALCvalue)],]) # next point


  ### Update the model ###
  newx <- matrix(chosen$Xnext, nrow=1)
  level <- chosen$level

  X1 <- t(t(fit1$X)*attr(fit1$X,"scaled:scale")+attr(fit1$X,"scaled:center"))
  X2 <- matrix(t(t(fit2$X)*attr(fit2$X,"scaled:scale")+attr(fit2$X,"scaled:center"))[,-ncol(fit2$X)], ncol=ncol(fit2$X)-1)

  if(constant){
    y1 <- fit1$y
    y2 <- fit2$y
  }else{
    y1 <- fit1$y+attr(fit1$y,"scaled:center")
    y2 <- fit2$y+attr(fit2$y,"scaled:center")
  }


  ### Choose level 1 ###
  if(level == 1){
    y1.select <- funcs[[1]](newx)

    X1 <- rbind(X1, newx)
    y1 <- c(y1, y1.select)
  }

  ### Choose level 2 ###
  if(level == 2){
    y1.select <- funcs[[1]](newx)
    y2.select <- funcs[[2]](newx)

    X1 <- rbind(X1, newx)
    y1 <- c(y1, y1.select)
    X2 <- rbind(X2, newx)
    y2 <- c(y2, y2.select)
  }

  fit <- RNAmf(X1, y1, X2, y2, kernel=kernel, constant=constant)


  return(list(fit=fit, intvar1=value1, intvar2=value2, cost=cost, Xcand=Xcand, chosen=chosen))
}
