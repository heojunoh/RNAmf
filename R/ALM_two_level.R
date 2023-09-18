#' obj.ALM_two_level_1
#'
#' object to optimize the point by ALM criterion updating at level 1 with two levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param fit an object of class RNAmf.
#' @return A mean of the negative predictive posterior variance at Xcand.
#' @export
#'

obj.ALM_two_level_1 <- function(Xcand, fit){

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


#' obj.ALM_two_level_2
#'
#' object to optimize the point by ALM criterion updating at level 2 with two levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param fit an object of class RNAmf.
#' @return A mean of the negative predictive posterior variance at Xcand.
#' @export
#'

obj.ALM_two_level_2 <- function(Xcand, fit){
  newx <- matrix(Xcand, nrow=1)
  -predRNAmf(fit, newx)$sig2 # to maximize the current variance.
}


#' ALM_two_level_optm
#'
#' find the next point by ALM criterion with two fidelity levels and update the model
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

  if (length(cost)!=2) stop("The length of cost should be 2")
  if (cost[1] >= cost[2]) stop("If the cost for high-fidelity is cheaper, just acquire the high-fidelity")
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
      return(c(-obj.ALM_two_level_1(newx, fit=fit),
               -obj.ALM_two_level_2(newx, fit=fit)))
    }
  }else{
    optm.mat <- cbind(c(rep(0, nrow(Xcand))), c(rep(0, nrow(Xcand))))
    for(i in 1:nrow(Xcand)){
      print(paste(i, nrow(Xcand), sep="/"))
      newx <- matrix(Xcand[i,], nrow=1)

      optm.mat[1,i] <- -obj.ALC_two_level_1(newx, fit=fit)
      optm.mat[2,i] <- -obj.ALC_two_level_2(newx, fit=fit)
    }
  }
  print(proc.time()[3]- time.start)

  ### Find the next point ###
  cat("running optim for level 1: \n")
  time.start <- proc.time()[3]
  X.start <- matrix(Xcand[which.max(optm.mat[1,]),], nrow=1)
  optim.out <- optim(X.start, obj.ALM_two_level_1, method="L-BFGS-B", lower=0, upper=1, fit=fit)
  Xnext.1 <- optim.out$par
  ALM.1 <- -optim.out$value
  print(proc.time()[3]- time.start)

  cat("running optim for level 2: \n")
  time.start <- proc.time()[3]
  X.start <- matrix(Xcand[which.max(optm.mat[2,]),], nrow=1)
  optim.out <- optim(X.start, obj.ALM_two_level_2, method="L-BFGS-B", lower=0, upper=1, fit=fit)
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

    # ### Blade ###
    # d1 <- data.frame(newx*0.5+0.25, rep(0.05, 1)) # scale X to [-1,1]
    # write.csv(d1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
    # run_matlab_script("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE,
    #                   splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
    #                   intern = TRUE)
    # d2 <- read.table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
    # y1.select <- d2$V4


    X1 <- rbind(X1, newx)
    y1 <- c(y1, y1.select)
  }

  ### Choose level 2 ###
  if(level == 2){
    y1.select <- funcs[[1]](newx)
    y2.select <- funcs[[2]](newx)

    # ### Blade ###
    # d1 <- data.frame(newx*0.5+0.25, rep(0.05, 1)) # scale X to [-1,1]
    # write.csv(d1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
    # run_matlab_script("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE,
    #                   splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
    #                   intern = TRUE)
    # d2 <- read.table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
    # y1.select <- d2$V4
    #
    # d1 <- data.frame(newx*0.5+0.25, rep(0.0125, 1)) # scale X to [-1,1]
    # write.csv(d1, "/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
    # run_matlab_script("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE,
    #                   splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
    #                   intern = TRUE)
    # d2 <- read.table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
    # y2.select <- d2$V4

    X1 <- rbind(X1, newx)
    y1 <- c(y1, y1.select)
    X2 <- rbind(X2, newx)
    y2 <- c(y2, y2.select)
  }

  fit <- RNAmf(X1, y1, X2, y2, kernel=kernel, constant=constant)

  if(parallel)  stopImplicitCluster()

  return(list(fit=fit, predsig2.low=optm.mat[1,], predsig2=optm.mat[2,], cost=cost, Xcand=Xcand, chosen=chosen))
}
