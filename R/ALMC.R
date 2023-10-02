#' find the next point by ALMC criterion with two fidelity levels and update the model
#'
#' @description The function acquires the new point in two fidelity levels setting
#' by the hybrid approach, Active learning MacKay-Cohn(ALMC) criterion.
#' The function generates the candidate set by maximin design,
#' and compute the posterior predictive variance at the highest level.
#' Then, it optimizes the ALC criterion starting at the candidate point with the maximum posterior predictive variance.
#' After optimization, the function calculates the ratio between the reduction in posterior variances and the simulation cost.
#' Finally, it acquires the point and the level maximizing the ALC criterion.
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

ALMC_two_level <- function(Xref=NULL, fit, mc.sample=100, cost, funcs, n.start, parallel=FALSE, ncore=1){

  if(length(cost)!=2) stop("The length of cost should be 2")
  if(cost[1] >= cost[2]) stop("If the cost for high-fidelity is cheaper, just acquire the high-fidelity")
  if(is.null(Xref)) Xref <- randomLHS(dim(fit$fit1$X)[1], dim(fit$fit1$X)[2])
  if(missing(n.start)) n.start <- 10 * dim(fit$fit1$X)[2]
  if(parallel) registerDoParallel(ncore)

  Icurrent <- mean(predRNAmf_two_level(fit, Xref)$sig2)

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

  ### Find the next point ###
  cat("running optim: \n")
  time.start <- proc.time()[3]
  if(parallel){
    optm.mat <- foreach(i = 1:nrow(Xcand), .combine=cbind) %dopar% {
      newx <- matrix(Xcand[i,], nrow=1)
      optim.ALM <- optim(newx, obj.ALM_level_2, method="L-BFGS-B", lower=0, upper=1, fit=fit)

      return(c(-optim.ALM$value, optim.ALM$par))
    }
  }else{
    optm.mat <- cbind(c(rep(0, nrow(Xcand))), c(rep(0, nrow(Xcand))))
    for(i in 1:nrow(Xcand)){
      print(paste(i, nrow(Xcand), sep="/"))
      newx <- matrix(Xcand[i,], nrow=1)
      optim.ALM <- optim(newx, obj.ALM_level_2, method="L-BFGS-B", lower=0, upper=1, fit=fit)

      optm.mat[,i] <- c(-optim.ALM$value, optim.ALM$par)
    }
  }
  print(proc.time()[3]- time.start)

  Xnext <- matrix(optm.mat[-1, which.max(optm.mat[1,])], nrow=1)

  ### Calculate the deduced variance ###
  cat("calculating deduced variance for level 1: \n")
  time.start <- proc.time()[3]
  ALMC.1 <- obj.ALC_level_1(Xnext, Xref=Xref, fit=fit, mc.sample=mc.sample, parallel=parallel, ncore=ncore)
  print(proc.time()[3]- time.start)

  cat("calculating deduced variance for level 2: \n")
  time.start <- proc.time()[3]
  ALMC.2 <- obj.ALC_level_2(Xnext, Xref=Xref, fit=fit, mc.sample=mc.sample, parallel=parallel, ncore=ncore)
  print(proc.time()[3]- time.start)

  ALMCvalue <- c(Icurrent - ALMC.1, Icurrent - ALMC.2)/c(cost[1], cost[1]+cost[2])
  if(ALMCvalue[2] > ALMCvalue[1]){
    level <- 2
  }else{
    level <- 1
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

  return(list(fit=fit, cost=cost, Xcand=Xcand, chosen=chosen))
}


#' find the next point by ALMC criterion with three fidelity levels and update the model
#'
#' @description The function acquires the new point in three fidelity levels setting
#' by the hybrid approach, Active learning MacKay-Cohn(ALMC) criterion.
#' The function generates the candidate set by maximin design,
#' and compute the posterior predictive variance at the highest level.
#' Then, it optimizes the ALC criterion starting at the candidate point with the maximum posterior predictive variance.
#' After optimization, the function calculates the ratio between the reduction in posterior variances and the simulation cost.
#' Finally, it acquires the point and the level maximizing the ALC criterion.
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

ALMC_three_level <- function(Xref=NULL, fit, mc.sample=100, cost, funcs, n.start, parallel=FALSE, ncore=1){

  if(length(cost)!=3) stop("The length of cost should be 3")
  if(cost[1] >= cost[2] | cost[2] >= cost[3]) stop("If the cost for high-fidelity is cheaper, just acquire the high-fidelity")
  if(is.null(Xref)) Xref <- randomLHS(dim(fit$fit.RNAmf_two_level$fit1$X)[1], dim(fit$fit.RNAmf_two_level$fit1$X)[2])
  if(missing(n.start)) n.start <- 10 * dim(fit$fit.RNAmf_two_level$fit1$X)[2]
  if(parallel) registerDoParallel(ncore)

  Icurrent <- mean(predRNAmf_three_level(fit, Xref)$sig2)

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

  ### Find the next point ###
  cat("running optim: \n")
  time.start <- proc.time()[3]
  if(parallel){
    optm.mat <- foreach(i = 1:nrow(Xcand), .combine=cbind) %dopar% {
      newx <- matrix(Xcand[i,], nrow=1)
      optim.ALM <- optim(newx, obj.ALM_level_3, method="L-BFGS-B", lower=0, upper=1, fit=fit)

      return(c(-optim.ALM$value, optim.ALM$par))
    }
  }else{
    optm.mat <- cbind(c(rep(0, nrow(Xcand))), c(rep(0, nrow(Xcand))))
    for(i in 1:nrow(Xcand)){
      print(paste(i, nrow(Xcand), sep="/"))
      newx <- matrix(Xcand[i,], nrow=1)
      optim.ALM <- optim(newx, obj.ALM_level_3, method="L-BFGS-B", lower=0, upper=1, fit=fit)

      optm.mat[,i] <- c(-optim.ALM$value, optim.ALM$par)
    }
  }
  print(proc.time()[3]- time.start)

  Xnext <- matrix(optm.mat[-1, which.max(optm.mat[1,])], nrow=1)

  ### Calculate the deduced variance ###
  cat("calculating deduced variance for level 1: \n")
  time.start <- proc.time()[3]
  ALMC.1 <- obj.ALC_level_1(Xnext, Xref=Xref, fit=fit_two_level, mc.sample=mc.sample, parallel=parallel, ncore=ncore)
  print(proc.time()[3]- time.start)

  cat("calculating deduced variance for level 2: \n")
  time.start <- proc.time()[3]
  ALMC.2 <- obj.ALC_level_2(Xnext, Xref=Xref, fit=fit_two_level, mc.sample=mc.sample, parallel=parallel, ncore=ncore)
  print(proc.time()[3]- time.start)

  cat("calculating deduced variance for level 3: \n")
  time.start <- proc.time()[3]
  ALMC.3 <- obj.ALC_level_3(Xnext, Xref=Xref, fit=fit, mc.sample=mc.sample, parallel=parallel, ncore=ncore)
  print(proc.time()[3]- time.start)

  ALMCvalue <- c(Icurrent - ALMC.1, Icurrent - ALMC.2, Icurrent - ALMC.3)/c(cost[1], cost[1]+cost[2], cost[1]+cost[2]+cost[3])
  if(ALMCvalue[3] > ALMCvalue[2]){
    level <- 3
  }else if(ALMCvalue[2] > ALMCvalue[1]){
    level <- 2
  }else{
    level <- 1
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

  return(list(fit=fit, cost=cost, Xcand=Xcand, chosen=chosen))
}
