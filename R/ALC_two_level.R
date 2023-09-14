#' obj.ALC_two_level_1
#'
#' object to optimize the next point by ALC criterion updating at level 1 with two levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param Xref vector or matrix of reference data.
#' @param fit an object of class RNAmf.
#' @param mc.sample a number of mc samples generated for this approach. Default is 100.
#' @param parallel logical indicating whether to run parallel or not. Default is FALSE.
#' @param ncore the number of core for parallel. Default is 1.
#' @return A mean of the deduced variance at Xref.
#' @importFrom plgp covar.sep
#' @importFrom stats rnorm
#' @importFrom maximin maximin
#' @export
#'

obj.ALC_two_level_1 <- function(Xcand, Xref, fit, mc.sample, parallel=FALSE, ncore=1){

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

  Xcand <- matrix(Xcand, nrow=1)
  if(kernel=="sqex"){
    y1.sample <- rnorm(mc.sample, mean=pred.GP(f1, Xcand)$mu, sd=sqrt(pred.GP(f1, Xcand)$sig2))
  }else if(kernel=="matern1.5"){
    y1.sample <- rnorm(mc.sample, mean=pred.matGP(f1, Xcand)$mu, sd=sqrt(pred.matGP(f1, Xcand)$sig2))
  }else if(kernel=="matern2.5"){
    y1.sample <- rnorm(mc.sample, mean=pred.matGP(f1, Xcand)$mu, sd=sqrt(pred.matGP(f1, Xcand)$sig2))
  }

  Xcand <- matrix((Xcand-x.center1)/x.scale1, nrow=1)

  ### Choose level 1 ###
  ### update Ki1
  if(kernel=="sqex"){
    cov.newx1 <- covar.sep(X1=Xcand, d=f1$theta, g=g)
    cov.Xnewx1 <- covar.sep(X1=f1$X, X2=Xcand, d=f1$theta, g=0)
  }else if(kernel=="matern1.5"){
    cov.newx1 <- cor.sep(X=Xcand, theta=f1$theta, nu=1.5)
    cov.Xnewx1 <- cor.sep(X=f1$X, x=Xcand, theta=f1$theta, nu=1.5)
  }else if(kernel=="matern2.5"){
    cov.newx1 <- cor.sep(X=Xcand, theta=f1$theta, nu=2.5)
    cov.Xnewx1 <- cor.sep(X=f1$X, x=Xcand, theta=f1$theta, nu=2.5)
  }
  v.next1 <- drop(cov.newx1 - t(cov.Xnewx1) %*% f1$Ki %*% cov.Xnewx1)
  g.next1 <- - 1/drop(v.next1) * f1$Ki %*% cov.Xnewx1


  fit1$Ki <- rbind(cbind(f1$Ki+g.next1%*%t(g.next1)*v.next1, g.next1),
                   cbind(t(g.next1), 1/drop(v.next1)))

  fit1$X <- rbind(f1$X, Xcand)
  attr(fit1$X, "scaled:center") <- x.center1
  attr(fit1$X, "scaled:scale") <- x.scale1


  ALC.out <- rep(0, mc.sample)
  if(parallel){
    ALC.out <- foreach(i = 1:mc.sample, .combine=c) %dopar% {
      fit.tmp <- fit
      if(constant){
        fit1$y <- c(f1$y, y1.sample[i])
      }else{
        fit1$y <- c(f1$y, y1.sample[i]-y.center1)
        attr(fit1$y, "scaled:center") <- y.center1
      }
      fit1$tau2hat <- drop(t(fit1$y - fit1$mu.hat) %*% fit1$Ki %*% (fit1$y - fit1$mu.hat) / length(fit1$y))

      fit.tmp$fit1 <- fit1

      return(mean(predRNAmf(fit.tmp, Xref)$sig2)) # to minimize the deduced variance. To maximize, -mean
    }
  }else{
    for(i in 1:mc.sample){
      fit.tmp <- fit
      if(constant){
        fit1$y <- c(f1$y, y1.sample[i])
      }else{
        fit1$y <- c(f1$y, y1.sample[i]-y.center1)
        attr(fit1$y, "scaled:center") <- y.center1
      }
      fit1$tau2hat <- drop(t(fit1$y - fit1$mu.hat) %*% fit1$Ki %*% (fit1$y - fit1$mu.hat) / length(fit1$y))

      fit.tmp$fit1 <- fit1

      ALC.out[i] <- mean(predRNAmf(fit.tmp, Xref)$sig2) # to minimize the deduced variance. To maximize, -mean
    }
  }
  return(mean(ALC.out))
}


#' obj.ALC_two_level_2
#'
#' object to optimize the next point by ALC criterion updating at level 2 with two levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param Xref vector or matrix of reference data.
#' @param fit an object of class RNAmf.
#' @param mc.sample a number of mc samples generated for this approach. Default is 100.
#' @param parallel logical indicating whether to run parallel or not. Default is FALSE.
#' @param ncore the number of core for parallel. Default is 1.
#' @return A mean of the deduced variance at Xref.
#' @importFrom plgp covar.sep
#' @importFrom stats rnorm
#' @importFrom maximin maximin
#' @export
#'

obj.ALC_two_level_2 <- function(Xcand, Xref, fit, mc.sample, parallel=FALSE, ncore=1){

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

  Xcand <- matrix(Xcand, nrow=1)
  if(kernel=="sqex"){
    y1.sample <- rnorm(mc.sample, mean=pred.GP(f1, Xcand)$mu, sd=sqrt(pred.GP(f1, Xcand)$sig2))
  }else if(kernel=="matern1.5"){
    y1.sample <- rnorm(mc.sample, mean=pred.matGP(f1, Xcand)$mu, sd=sqrt(pred.matGP(f1, Xcand)$sig2))
  }else if(kernel=="matern2.5"){
    y1.sample <- rnorm(mc.sample, mean=pred.matGP(f1, Xcand)$mu, sd=sqrt(pred.matGP(f1, Xcand)$sig2))
  }
  Xcand <- matrix((Xcand-x.center1)/x.scale1, nrow=1)

  ### Choose level 1 ###
  ### update Ki1
  if(kernel=="sqex"){
    cov.newx1 <- covar.sep(X1=Xcand, d=f1$theta, g=g)
    cov.Xnewx1 <- covar.sep(X1=f1$X, X2=Xcand, d=f1$theta, g=0)
  }else if(kernel=="matern1.5"){
    cov.newx1 <- cor.sep(X=Xcand, theta=f1$theta, nu=1.5)
    cov.Xnewx1 <- cor.sep(X=f1$X, x=Xcand, theta=f1$theta, nu=1.5)
  }else if(kernel=="matern2.5"){
    cov.newx1 <- cor.sep(X=Xcand, theta=f1$theta, nu=2.5)
    cov.Xnewx1 <- cor.sep(X=f1$X, x=Xcand, theta=f1$theta, nu=2.5)
  }
  v.next1 <- drop(cov.newx1 - t(cov.Xnewx1) %*% f1$Ki %*% cov.Xnewx1)
  g.next1 <- - 1/drop(v.next1) * f1$Ki %*% cov.Xnewx1


  fit1$Ki <- rbind(cbind(f1$Ki+g.next1%*%t(g.next1)*v.next1, g.next1),
                   cbind(t(g.next1), 1/drop(v.next1)))

  fit1$X <- rbind(f1$X, Xcand)
  attr(fit1$X, "scaled:center") <- x.center1
  attr(fit1$X, "scaled:scale") <- x.scale1

  if(parallel){
    ALC.out <- foreach(i = 1:mc.sample, .combine=c) %dopar% {
      fit.tmp <- fit
      if(constant){
        fit1$y <- c(f1$y, y1.sample[i])
      }else{
        fit1$y <- c(f1$y, y1.sample[i]-y.center1)
        attr(fit1$y, "scaled:center") <- y.center1
      }

      fit1$tau2hat <- drop(t(fit1$y - fit1$mu.hat) %*% fit1$Ki %*% (fit1$y - fit1$mu.hat) / length(fit1$y))

      fit.tmp$fit1 <- fit1

      ### Choose level 2 ###
      if(kernel=="sqex"){
        pred2 <- pred.GP(fit2, cbind(Xcand, y1.sample[i]))
        y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      }else if(kernel=="matern1.5"){
        pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
        y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      }else if(kernel=="matern2.5"){
        pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
        y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      }

      ### update Ki2
      newx2 <- t((t(cbind(Xcand, y1.sample[i]))-x.center2)/x.scale2)

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
        fit2$y <- c(f2$y, y2.sample)
      }else{
        fit2$y <- c(f2$y, y2.sample-y.center2)
        attr(fit2$y, "scaled:center") <- y.center2
      }

      fit2$tau2hat <- drop(t(fit2$y - fit2$mu.hat) %*% fit2$Ki %*% (fit2$y - fit2$mu.hat) / length(fit2$y))

      fit.tmp$fit2 <- fit2

      return(mean(predRNAmf(fit.tmp, Xref)$sig2)) # to minimize the deduced variance. To maximize, -mean
    }
  }else{
    ALC.out <- rep(0, mc.sample)
    for(i in 1:mc.sample){
      fit.tmp <- fit
      if(constant){
        fit1$y <- c(f1$y, y1.sample[i])
      }else{
        fit1$y <- c(f1$y, y1.sample[i]-y.center1)
        attr(fit1$y, "scaled:center") <- y.center1
      }

      fit1$tau2hat <- drop(t(fit1$y - fit1$mu.hat) %*% fit1$Ki %*% (fit1$y - fit1$mu.hat) / length(fit1$y))

      fit.tmp$fit1 <- fit1

      ### Choose level 2 ###
      if(kernel=="sqex"){
        pred2 <- pred.GP(fit2, cbind(Xcand, y1.sample[i]))
        y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      }else if(kernel=="matern1.5"){
        pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
        y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      }else if(kernel=="matern2.5"){
        pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
        y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      }

      ### update Ki2
      newx2 <- t((t(cbind(Xcand, y1.sample[i]))-x.center2)/x.scale2)

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
        fit2$y <- c(f2$y, y2.sample)
      }else{
        fit2$y <- c(f2$y, y2.sample-y.center2)
        attr(fit2$y, "scaled:center") <- y.center2
      }

      fit2$tau2hat <- drop(t(fit2$y - fit2$mu.hat) %*% fit2$Ki %*% (fit2$y - fit2$mu.hat) / length(fit2$y))

      fit.tmp$fit2 <- fit2

      ALC.out[i] <- mean(predRNAmf(fit.tmp, Xref)$sig2) # to minimize the deduced variance. To maximize, -mean
    }
  }

  return(mean(ALC.out))
}

#' ALC_two_level
#'
#' find the next point by ALC criterion with two fidelity levels and update the model
#'
#' @param Xref vector or matrix of reference data.
#' @param fit an object of class RNAmf.
#' @param mc.sample a number of mc samples generated for this approach. Default is 100.
#' @param cost a vector of the costs for each level of fidelity.
#' @param funcs list of functions for each level of fidelity.
#' @param n.start a number of candidate point. Default is 10*d.
#' @param parallel logical indicating whether to run parallel or not. Default is FALSE.
#' @param ncore the number of core for parallel. Default is 1.
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
#' @importFrom lhs randomLHS
#' @importFrom lhs maximinLHS
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom doParallel stopImplicitCluster
#' @export
#'

ALC_two_level <- function(Xref=NULL, fit, mc.sample=100, cost, funcs, n.start, parallel=FALSE, ncore=1){

  if(length(cost)!=2) stop("The length of cost should be 2")
  if(cost[1] >= cost[2]) stop("If the cost for high-fidelity is cheaper, just acquire the high-fidelity")
  if(is.null(Xref)) Xref <- randomLHS(dim(fit$fit1$X)[1], dim(fit$fit1$X)[2])
  if(missing(n.start)) n.start <- 10 * dim(fit$fit1$X)[2]
  if(parallel) registerDoParallel(ncore)

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

  ### Generate the candidate set ###
  Xcand <- maximinLHS(n.start, ncol(fit1$X))

  ### Calculate the deduced variance ###
  cat("running starting points: \n")
  time.start <- proc.time()[3]
  if(parallel){
    pseudointvar <- foreach(i = 1:nrow(Xcand), .combine=cbind) %dopar% {
      newx <- matrix(Xcand[i,], nrow=1)

      return(c(obj.ALC_two_level_1(newx, Xref, fit, mc.sample),
               obj.ALC_two_level_2(newx, Xref, fit, mc.sample)))
    }
    intvar1 <- pseudointvar[,1]
    intvar2 <- pseudointvar[,2]
  }else{
    intvar1 <- c(rep(0, nrow(Xcand))) # IMSPE candidates
    intvar2 <- c(rep(0, nrow(Xcand))) # IMSPE candidates
    for(i in 1:nrow(Xcand)){
      print(paste(i, nrow(Xcand), sep="/"))
      newx <- matrix(Xcand[i,], nrow=1)

      intvar1[i] <- obj.ALC_two_level_1(newx, Xref, fit, mc.sample)
      intvar2[i] <- obj.ALC_two_level_2(newx, Xref, fit, mc.sample)
    }
  }
  print(proc.time()[3]- time.start)

  ### Find the next point ###
  cat("running optim for level 1: \n")
  time.start <- proc.time()[3]
  X.start <- matrix(Xcand[which.min(intvar1),], nrow=1)
  optim.out <- optim(X.start, obj.ALC_two_level_1, method="L-BFGS-B", lower=0, upper=1, fit=fit, Xref=Xref, mc.sample=mc.sample, parallel=parallel, ncore=ncore)
  Xnext.1 <- optim.out$par
  ALC.1 <- optim.out$value
  print(proc.time()[3]- time.start)

  cat("running optim for level 2: \n")
  time.start <- proc.time()[3]
  X.start <- matrix(Xcand[which.min(intvar2),], nrow=1)
  optim.out <- optim(X.start, obj.ALC_two_level_2, method="L-BFGS-B", lower=0, upper=1, fit=fit, Xref=Xref, mc.sample=mc.sample, parallel=parallel, ncore=ncore)
  Xnext.2 <- optim.out$par
  ALC.2 <- optim.out$value
  print(proc.time()[3]- time.start)

  ALCvalue <- c(Icurrent - ALC.1, Icurrent - ALC.2)/c(cost[1], cost[1]+cost[2])
  if(ALCvalue[2] > ALCvalue[1]){
    level <- 2
    Xnext <- Xnext.2
  }else{
    level <- 1
    Xnext <- Xnext.1
  }

  chosen <- list("level"=level, # next level
                 "Xnext"=Xnext) # next point


  ### Update the model ###
  newx <- matrix(chosen$Xnext, nrow=1)
  # if(newx[1,1]==0) newx[1,1] <- sqrt(.Machine$double.eps) # To prevent yielding infinity for Park function
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

  if(parallel)  stopImplicitCluster()

  return(list(fit=fit, intvar1=intvar1, intvar2=intvar2, cost=cost, Xcand=Xcand, chosen=chosen))
}
