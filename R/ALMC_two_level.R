#' ALMC_two_level
#'
#' find the next point by ALMC criterion with two fidelity levels and update the model
#'
#' @param Xref vector or matrix of reference data.
#' @param fit an object of class RNAmf.
#' @param mc.sample a number of mc samples generated for this approach. Default is 100.
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
#' @importFrom lhs randomLHS
#' @importFrom maximin maximin
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @export
#'

ALMC_two_level <- function(Xref, fit, mc.sample=100, cost, funcs){

  if (length(cost)!=2) stop("The length of cost should be 2")
  if (cost[1] >= cost[2]) stop("If the cost for high-fidelity is cheaper, just acquire the high-fidelity")
  if(is.null(Xref)) Xref <- randomLHS(dim(fit$fit1$X)[1], dim(fit$fit1$X)[2])

  registerDoParallel(5)

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
  predsig2 <- predRNAmf(fit, Xcand)$sig2
  intvar1 <- c(0) # IMSPE candidates
  intvar2 <- c(0) # IMSPE candidates

  newx <- matrix(Xcand[which.max(predsig2),], nrow=1)

  if(kernel=="sqex"){
    x1.sample <- rnorm(mc.sample, mean=pred.GP(f1, newx)$mu, sd=sqrt(pred.GP(f1, newx)$sig2))
  }else if(kernel=="matern1.5"){
    x1.sample <- rnorm(mc.sample, mean=pred.matGP(f1, newx)$mu, sd=sqrt(pred.matGP(f1, newx)$sig2))
  }else if(kernel=="matern2.5"){
    x1.sample <- rnorm(mc.sample, mean=pred.matGP(f1, newx)$mu, sd=sqrt(pred.matGP(f1, newx)$sig2))
  }

  ### MC ###
  pseudointvar <- foreach(j = 1:mc.sample, .combine=rbind) %dopar% {
    # print(j)
    ### Optimize at level 1 ###
    out1 <- optim(newx, obj.ALC_two_level_1, method="L-BFGS-B", lower=0, upper=1, fit=fit, Xref=Xref, y1.sample=x1.sample[j])
    ### Optimize at level 2 ###
    out2 <- optim(newx, obj.ALC_two_level_2, method="L-BFGS-B", lower=0, upper=1, fit=fit, Xref=Xref, y1.sample=x1.sample[j])

    return(c(out1$value, out2$value))
  }

  intvar1 <- mean(pseudointvar[,1])
  intvar2 <- mean(pseudointvar[,2])


  ### Select the next level and point ###
  ALMCvalue <- c(Icurrent - intvar1, Icurrent - intvar2)/c(cost[1], cost[1]+cost[2])

  if(which.max(ALMCvalue)==1){
    location <- optim(newx, obj.ALC_two_level_1, method="L-BFGS-B", lower=0, upper=1, fit=fit, Xref=Xref, y1.sample=mean(x1.sample))$par
  }else if(which.max(ALMCvalue)==2){
    location <- optim(newx, obj.ALC_two_level_2, method="L-BFGS-B", lower=0, upper=1, fit=fit, Xref=Xref, y1.sample=mean(x1.sample))$par
  }

  chosen <- list("level"=which.max(ALMCvalue), # next level
                 "location"=which.max(predsig2), # next location
                 "Xnext"=location) # next point


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


  return(list(fit=fit, intvar1=intvar1, intvar2=intvar2, cost=cost, Xcand=Xcand, chosen=chosen))
}
