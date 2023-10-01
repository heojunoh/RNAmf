#' find the next point by ALMC criterion with two fidelity levels and update the model
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

  predsig2 <- predRNAmf_two_level(fit, Xcand)$sig2

  ### Find the next point ###
  cat("running optim: \n")
  time.start <- proc.time()[3]
  if(parallel){
    optm.mat <- foreach(i = 1:nrow(Xcand), .combine=cbind) %dopar% {
      newx <- matrix(Xcand[i,], nrow=1)
      optim.ALM <- optim(newx, obj.ALM_two_level_2, method="L-BFGS-B", lower=0, upper=1, fit=fit)

      return(c(-optim.ALM$value, optim.ALM$par))
    }
  }else{
    optm.mat <- cbind(c(rep(0, nrow(Xcand))), c(rep(0, nrow(Xcand))))
    for(i in 1:nrow(Xcand)){
      print(paste(i, nrow(Xcand), sep="/"))
      newx <- matrix(Xcand[i,], nrow=1)
      optim.ALM <- optim(newx, obj.ALM_two_level_2, method="L-BFGS-B", lower=0, upper=1, fit=fit)

      optm.mat[,i] <- c(-optim.ALM$value, optim.ALM$par)
    }
  }
  print(proc.time()[3]- time.start)

  Xnext <- matrix(optm.mat[-1, which.max(optm.mat[1,])], nrow=1)

  ### Calculate the deduced variance ###
  cat("calculating deduced variance for level 1: \n")
  time.start <- proc.time()[3]
  ALMC.1 <- obj.ALC_two_level_1(Xnext, Xref=Xref, fit=fit, mc.sample=mc.sample, parallel=parallel, ncore=ncore)
  print(proc.time()[3]- time.start)

  cat("calculating deduced variance for level 2: \n")
  time.start <- proc.time()[3]
  ALMC.2 <- obj.ALC_two_level_2(Xnext, Xref=Xref, fit=fit, mc.sample=mc.sample, parallel=parallel, ncore=ncore)
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
    y1.select <- funcs[[1]](newx)
    y2.select <- funcs[[2]](newx)

    X1 <- rbind(X1, newx)
    y1 <- c(y1, y1.select)
    X2 <- rbind(X2, newx)
    y2 <- c(y2, y2.select)
  }

  fit <- RNAmf_two_level(X1, y1, X2, y2, kernel=kernel, constant=constant)

  if(parallel)  stopImplicitCluster()

  return(list(fit=fit, cost=cost, Xcand=Xcand, chosen=chosen))
}
