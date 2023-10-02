#' object to optimize the point by ALM criterion updating at level 1 with two levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param fit an object of class RNAmf.
#' @return A mean of the negative predictive posterior variance at Xcand.
#' @noRd
#'

obj.ALM_level_1 <- function(Xcand, fit){

  newx <- matrix(Xcand, nrow=1)
  fit1 <- fit$fit1
  kernel <- fit$kernel

  ### calculate the posterior predictive variance ###
  if(kernel=="sqex"){
    predsig2 <- pred.GP(fit1, newx)$sig2
  }else if(kernel=="matern1.5"){
    predsig2 <- pred.matGP(fit1, newx)$sig2
  }else if(kernel=="matern2.5"){
    predsig2 <- pred.matGP(fit1, newx)$sig2
  }

  -predsig2 # to maximize the current variance.
}


#' object to optimize the point by ALM criterion updating at level 2 with two levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param fit an object of class RNAmf.
#' @return A mean of the negative predictive posterior variance at Xcand.
#' @noRd
#'

obj.ALM_level_2 <- function(Xcand, fit){
  newx <- matrix(Xcand, nrow=1)
  -predRNAmf_two_level(fit, newx)$sig2 # to maximize the current variance.
}


#' object to optimize the point by ALM criterion updating at level 3 with three levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param fit an object of class RNAmf.
#' @return A mean of the negative predictive posterior variance at Xcand.
#' @noRd
#'

obj.ALM_level_3 <- function(Xcand, fit){
  newx <- matrix(Xcand, nrow=1)
  -predRNAmf_three_level(fit, newx)$sig2 # to maximize the current variance.
}


#' find the next point by ALM criterion with two fidelity levels and update the model
#'
#' @description The function acquires the new point in two fidelity levels setting
#' by the Active learning MacKay(ALM) criterion.
#' The function generates the candidate set by maximin design,
#' and compute the posterior predictive variance at each level.
#' Then, it optimizes the ALM criterion starting at the candidate point with the maximum posterior predictive variance.
#' After optimization, the function calculates the ratio between the posterior predictive variance and the simulation cost.
#' Finally, it acquires the point and the level maximizing the ALM criterion.
#'
#' @param fit an object of class RNAmf.
#' @param cost a vector of the costs for each level of fidelity.
#' @param funcs list of functions for each level of fidelity.
#' @param n.start a number of candidate point. Default is 10*d.
#' @param parallel logical indicating whether to run parallel or not. Default is FALSE.
#' @param ncore the number of core for parallel.
#' @return A list containing the integrated one-step ahead variance:
#' \itemize{
#'   \item \code{fit}: fitted model after acquire chosen point.
#'   \item \code{predsig2.low}: vector of predictive variance of optimized data point at level 1.
#'   \item \code{predsig2}: vector of predictive variance of optimized data point at level 2.
#'   \item \code{cost}: vector of the costs for each level of fidelity.
#'   \item \code{Xcand}: candidate data points considered to be acquired.
#'   \item \code{chosen}: list of chosen level, location, and point.
#' }
#' @importFrom plgp covar.sep
#' @importFrom lhs maximinLHS
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom doParallel stopImplicitCluster
#' @export
#'

ALM_two_level <- function(fit, cost, funcs, n.start, parallel=FALSE, ncore=1){

  if(length(cost)!=2) stop("The length of cost should be 2")
  if(cost[1] >= cost[2]) stop("If the cost for high-fidelity is cheaper, just acquire the high-fidelity")
  if(missing(n.start)) n.start <- 10 * dim(fit$fit1$X)[2]
  if(parallel) registerDoParallel(ncore)

  fit1 <- fit$fit1
  fit2 <- fit$fit2
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

  ### Calculate the current variance at each level ###
  cat("running starting points: \n")
  time.start <- proc.time()[3]
  if(parallel){
    optm.mat <- foreach(i = 1:nrow(Xcand), .combine=cbind) %dopar% {
      newx <- matrix(Xcand[i,], nrow=1)
      return(c(-obj.ALM_level_1(newx, fit=fit),
               -obj.ALM_level_2(newx, fit=fit)))
    }
  }else{
    optm.mat <- cbind(c(rep(0, nrow(Xcand))), c(rep(0, nrow(Xcand))))
    for(i in 1:nrow(Xcand)){
      print(paste(i, nrow(Xcand), sep="/"))
      newx <- matrix(Xcand[i,], nrow=1)

      optm.mat[1,i] <- -obj.ALM_level_1(newx, fit=fit)
      optm.mat[2,i] <- -obj.ALM_level_2(newx, fit=fit)
    }
  }
  print(proc.time()[3]- time.start)

  ### Find the next point ###
  cat("running optim for level 1: \n")
  time.start <- proc.time()[3]
  X.start <- matrix(Xcand[which.max(optm.mat[1,]),], nrow=1)
  optim.out <- optim(X.start, obj.ALM_level_1, method="L-BFGS-B", lower=0, upper=1, fit=fit)
  Xnext.1 <- optim.out$par
  ALM.1 <- -optim.out$value
  print(proc.time()[3]- time.start)

  cat("running optim for level 2: \n")
  time.start <- proc.time()[3]
  X.start <- matrix(Xcand[which.max(optm.mat[2,]),], nrow=1)
  optim.out <- optim(X.start, obj.ALM_level_2, method="L-BFGS-B", lower=0, upper=1, fit=fit)
  Xnext.2 <- optim.out$par
  ALM.2 <- -optim.out$value
  print(proc.time()[3]- time.start)

  ALMvalue <- c(ALM.1, ALM.2)/c(cost[1], cost[1]+cost[2])
  if(ALMvalue[2] > ALMvalue[1]){
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
  level <- chosen$level

  X1 <- scale.inputs(fit1$X, x.center1, x.scale1, back=TRUE)
  X2 <- matrix(scale.inputs(fit2$X, x.center2, x.scale2, back=TRUE)[,-ncol(fit2$X)], ncol=ncol(fit2$X)-1)

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
    if(checknested(scale.inputs(fit1$X, x.center1, x.scale1, back=TRUE), chosen$Xnext)==TRUE &
       checknested(scale.inputs(fit2$X, x.center2, x.scale2, back=TRUE), chosen$Xnext)==FALSE){
      y2.select <- funcs[[2]](newx)

      X2 <- rbind(X2, newx)
      y2 <- c(y2, y2.select)
    }else{
      y1.select <- funcs[[1]](newx)
      y2.select <- funcs[[2]](newx)

      X1 <- rbind(X1, newx)
      y1 <- c(y1, y1.select)
      X2 <- rbind(X2, newx)
      y2 <- c(y2, y2.select)
    }
  }

  fit <- RNAmf_two_level(X1, y1, X2, y2, kernel=kernel, constant=constant)

  if(parallel)  stopImplicitCluster()

  return(list(fit=fit, predsig2.low=optm.mat[1,], predsig2=optm.mat[2,], cost=cost, Xcand=Xcand, chosen=chosen))
}


#' find the next point by ALM criterion with three fidelity levels and update the model
#'
#' @description The function acquires the new point in three fidelity levels setting
#' by the Active learning MacKay(ALM) criterion.
#' The function generates the candidate set by maximin design,
#' and compute the posterior predictive variance at each level.
#' Then, it optimizes the ALM criterion starting at the candidate point with the maximum posterior predictive variance.
#' After optimization, the function calculates the ratio between the posterior predictive variance and the simulation cost.
#' Finally, it acquires the point and the level maximizing the ALM criterion.
#'
#' @param fit an object of class RNAmf.
#' @param cost a vector of the costs for each level of fidelity.
#' @param funcs list of functions for each level of fidelity.
#' @param n.start a number of candidate point. Default is 10*d.
#' @param parallel logical indicating whether to run parallel or not. Default is FALSE.
#' @param ncore the number of core for parallel.
#' @return A list containing the integrated one-step ahead variance:
#' \itemize{
#'   \item \code{fit}: fitted model after acquire chosen point.
#'   \item \code{predsig2.low}: vector of predictive variance of optimized data point at level 1.
#'   \item \code{predsig2.med}: vector of predictive variance of optimized data point at level 2.
#'   \item \code{predsig2.high}: vector of predictive variance of optimized data point at level 3.
#'   \item \code{cost}: vector of the costs for each level of fidelity.
#'   \item \code{Xcand}: candidate data points considered to be acquired.
#'   \item \code{chosen}: list of chosen level, location, and point.
#' }
#' @importFrom plgp covar.sep
#' @importFrom lhs maximinLHS
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom doParallel stopImplicitCluster
#' @export
#'

ALM_three_level <- function(fit, cost, funcs, n.start, parallel=FALSE, ncore=1){

  if(length(cost)!=3) stop("The length of cost should be 3")
  if(cost[1] >= cost[2] | cost[2] >= cost[3]) stop("If the cost for high-fidelity is cheaper, just acquire the high-fidelity")
  if(missing(n.start)) n.start <- 10 * dim(fit$fit.RNAmf_two_level$fit1$X)[2]
  if(parallel) registerDoParallel(ncore)

  fit_two_level <- fit$fit.RNAmf_two_level
  fit1 <- fit_two_level$fit1
  fit2 <- fit_two_level$fit2
  fit3 <- fit$fit3
  constant <- fit$constant
  kernel <- fit$kernel
  g <- fit1$g

  x.center1 <- attr(fit1$X, "scaled:center")
  x.scale1 <- attr(fit1$X, "scaled:scale")
  y.center1 <- attr(fit1$y, "scaled:center")

  x.center2 <- attr(fit2$X, "scaled:center")
  x.scale2 <- attr(fit2$X, "scaled:scale")
  y.center2 <- attr(fit2$y, "scaled:center")

  x.center3 <- attr(fit3$X, "scaled:center")
  x.scale3 <- attr(fit3$X, "scaled:scale")
  y.center3 <- attr(fit3$y, "scaled:center")


  ### Generate the candidate set ###
  Xcand <- maximinLHS(n.start, ncol(fit1$X))

  ### Calculate the current variance at each level ###
  cat("running starting points: \n")
  time.start <- proc.time()[3]
  if(parallel){
    optm.mat <- foreach(i = 1:nrow(Xcand), .combine=cbind) %dopar% {
      newx <- matrix(Xcand[i,], nrow=1)

      return(c(-obj.ALM_level_1(newx, fit=fit_two_level),
               -obj.ALM_level_2(newx, fit=fit_two_level),
               -obj.ALM_level_3(newx, fit=fit)))
    }
  }else{
    optm.mat <- cbind(c(rep(0, nrow(Xcand))), c(rep(0, nrow(Xcand))), c(rep(0, nrow(Xcand))))
    for(i in 1:nrow(Xcand)){
      print(paste(i, nrow(Xcand), sep="/"))
      newx <- matrix(Xcand[i,], nrow=1)

      optm.mat[1,i] <- -obj.ALM_level_1(newx, fit=fit_two_level)
      optm.mat[2,i] <- -obj.ALM_level_2(newx, fit=fit_two_level)
      optm.mat[3,i] <- -obj.ALM_level_3(newx, fit=fit)
    }
  }
  print(proc.time()[3]- time.start)

  ### Find the next point ###
  cat("running optim for level 1: \n")
  time.start <- proc.time()[3]
  X.start <- matrix(Xcand[which.max(optm.mat[1,]),], nrow=1)
  optim.out <- optim(X.start, obj.ALM_level_1, method="L-BFGS-B", lower=0, upper=1, fit=fit_two_level)
  Xnext.1 <- optim.out$par
  ALM.1 <- -optim.out$value
  print(proc.time()[3]- time.start)

  cat("running optim for level 2: \n")
  time.start <- proc.time()[3]
  X.start <- matrix(Xcand[which.max(optm.mat[2,]),], nrow=1)
  optim.out <- optim(X.start, obj.ALM_level_2, method="L-BFGS-B", lower=0, upper=1, fit=fit_two_level)
  Xnext.2 <- optim.out$par
  ALM.2 <- -optim.out$value
  print(proc.time()[3]- time.start)

  cat("running optim for level 3: \n")
  time.start <- proc.time()[3]
  X.start <- matrix(Xcand[which.max(optm.mat[3,]),], nrow=1)
  optim.out <- optim(X.start, obj.ALM_level_3, method="L-BFGS-B", lower=0, upper=1, fit=fit)
  Xnext.3 <- optim.out$par
  ALM.3 <- -optim.out$value
  print(proc.time()[3]- time.start)

  ALMvalue <- c(ALM.1, ALM.2, ALM.3)/c(cost[1], cost[1]+cost[2], cost[1]+cost[2]+cost[3])
  if(ALMvalue[3] > ALMvalue[2]){
    level <- 3
    Xnext <- Xnext.3
  }else if(ALMvalue[2] > ALMvalue[1]){
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
  level <- chosen$level

  X1 <- scale.inputs(fit1$X, x.center1, x.scale1, back=TRUE)
  X2 <- matrix(scale.inputs(fit2$X, x.center2, x.scale2, back=TRUE)[,-ncol(fit2$X)], ncol=ncol(fit2$X)-1)
  X3 <- matrix(scale.inputs(fit3$X, x.center3, x.scale3, back=TRUE)[,-ncol(fit3$X)], ncol=ncol(fit3$X)-1)

  if(constant){
    y1 <- fit1$y
    y2 <- fit2$y
    y3 <- fit3$y
  }else{
    y1 <- fit1$y+attr(fit1$y,"scaled:center")
    y2 <- fit2$y+attr(fit2$y,"scaled:center")
    y3 <- fit3$y+attr(fit3$y,"scaled:center")
  }


  ### Choose level 1 ###
  if(level == 1){
    y1.select <- funcs[[1]](newx)

    X1 <- rbind(X1, newx)
    y1 <- c(y1, y1.select)
  }

  ### Choose level 2 ###
  if(level == 2){
    if(checknested(scale.inputs(fit1$X, x.center1, x.scale1, back=TRUE), chosen$Xnext)==TRUE &
       checknested(scale.inputs(fit2$X, x.center2, x.scale2, back=TRUE), chosen$Xnext)==FALSE){
      y2.select <- funcs[[2]](newx)

      X2 <- rbind(X2, newx)
      y2 <- c(y2, y2.select)
    }else{
      y1.select <- funcs[[1]](newx)
      y2.select <- funcs[[2]](newx)

      X1 <- rbind(X1, newx)
      y1 <- c(y1, y1.select)
      X2 <- rbind(X2, newx)
      y2 <- c(y2, y2.select)
    }
  }

  ### Choose level 3 ###
  if(level == 3){
    if(checknested(scale.inputs(fit1$X, x.center1, x.scale1, back=TRUE), chosen$Xnext)==TRUE &
       checknested(scale.inputs(fit2$X, x.center2, x.scale2, back=TRUE), chosen$Xnext)==TRUE &
       checknested(scale.inputs(fit3$X, x.center3, x.scale3, back=TRUE), chosen$Xnext)==FALSE){
      y3.select <- funcs[[3]](newx)

      X3 <- rbind(X3, newx)
      y3 <- c(y3, y3.select)
    }else if(checknested(scale.inputs(fit1$X, x.center1, x.scale1, back=TRUE), chosen$Xnext)==TRUE &
             checknested(scale.inputs(fit2$X, x.center2, x.scale2, back=TRUE), chosen$Xnext)==FALSE &
             checknested(scale.inputs(fit3$X, x.center3, x.scale3, back=TRUE), chosen$Xnext)==FALSE){
      y2.select <- funcs[[2]](newx)
      y3.select <- funcs[[3]](newx)

      X2 <- rbind(X2, newx)
      y2 <- c(y2, y2.select)
      X3 <- rbind(X3, newx)
      y3 <- c(y3, y3.select)
    }else{
      y1.select <- funcs[[1]](newx)
      y2.select <- funcs[[2]](newx)
      y3.select <- funcs[[3]](newx)

      X1 <- rbind(X1, newx)
      y1 <- c(y1, y1.select)
      X2 <- rbind(X2, newx)
      y2 <- c(y2, y2.select)
      X3 <- rbind(X3, newx)
      y3 <- c(y3, y3.select)
    }
  }

  fit <- RNAmf_three_level(X1, y1, X2, y2, X3, y3, kernel=kernel, constant=constant)

  if(parallel)  stopImplicitCluster()

  return(list(fit=fit, predsig2.low=optm.mat[1,], predsig2.med=optm.mat[2,], predsig2.high=optm.mat[3,], cost=cost, Xcand=Xcand, chosen=chosen))
}
