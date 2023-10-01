#' scale inputs before fit the model
#'
#' @param X vector or matrix of input locations to be scaled.
#' @param back logical indicating for scale back to the original vector or matrix (back=TRUE) or scale the original vector or matrix (back=FALSE). Default is FALSE.
#'
#' @return A scaled X or original X.
#'
#' #' @importFrom plgp distance
#' @noRd
#'

scale.inputs <- function(X, center=NULL, scale=NULL, back=FALSE){

  if(back){
    if(is.null(center) | is.null(scale)) stop("center and scale are required to scale back")
    X <- t(t(X) * scale + center)
  }else{
    if(is.null(center)) center <- attr(scale(X),"scaled:center")
    if(is.null(scale)) scale <- attr(scale(X),"scaled:scale")

    X <- t((t(X)-center)/scale)
    attr(X, "scaled:center") <- center
    attr(X, "scaled:scale") <- scale
  }
  return(X)
}

#' Check the design is nested
#'
#' @param X1 vector or matrix of input locations at lower fidelity
#' @param X2 vector or matrix of input locations at higher fidelity
#'
#' @return A logical indicating if X2 is nested or not.
#'
#' @noRd
#'

checknested <- function(XX1, XX2){
  checknest <- c()
  for(i in 1:nrow(XX2)){
    checknest <- c(checknest, suppressWarnings(any(apply(XX1, 1, function(xx){all.equal(XX2[i,], xx, tolerance=sqrt(.Machine$double.eps))}))))
  }
  checknest[is.na(checknest)] <- FALSE
  all(checknest)
}




