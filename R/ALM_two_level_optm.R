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
#' @importFrom maximin maximin
#' @export
#'

ALM_two_level_optm <- function(fit, cost, funcs){

  if (length(cost)!=2) stop("The length of cost should be 2")
  if (cost[1] >= cost[2]) stop("If the cost for high-fidelity is cheaper, just acquire the high-fidelity")

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
    out <- optim(Xcand[i,], obj.ALM_two_level_1, method="L-BFGS-B", lower=0, upper=1, fit=fit)
    xoptm1[i,] <- out$par
    value1[i] <- -out$value
  }

  ### Find the optimal point at level 1 ###
  xnext1 <- xoptm1[which.max(value1),]


  ### Optimize at level 2 ###
  for(i in 1:nrow(Xcand)){
    out <- optim(Xcand[i,], obj.ALM_two_level_2, method="L-BFGS-B", lower=0, upper=1, fit=fit)
    xoptm2[i,] <- out$par
    value2[i] <- -out$value
  }

  ### Find the optimal point at level 2 ###
  xnext2 <- xoptm2[which.max(value2),]


  ### Find the next point ###
  ALMvalue <- c(value1[which.max(value1)], value2[which.max(value2)])/c(cost[1], cost[1]+cost[2])

  chosen <- list("level"=which.max(ALMvalue), # next level
                 "location"=c(which.max(value1), which.max(value2))[which.max(ALMvalue)], # next location
                 "Xnext"=Xcand[c(which.max(value1), which.max(value2))[which.max(ALMvalue)],]) # next point


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


  return(list(fit=fit, predsig2.low=value1, predsig2=value2, cost=cost, Xcand=Xcand, chosen=chosen))
}
