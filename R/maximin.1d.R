#' maximin.1d
#'
#' maximin design for 1-dimensional dataset.
#'
#' @param Xorig The existing design.
#' @param T a number of iterations. Default is 10000*length(Xorig)
#' @return A 1d matrix of data with maximin design.
#' @importFrom plgp distance
#' @importFrom stats runif
#' @noRd
#'

maximin.1d <- function(Xorig, T=10*length(Xorig)) {

  if(is.vector(Xorig)) Xorig <- matrix(Xorig, ncol=1)
  if(ncol(Xorig) != 1) stop("This function is for 1d")

  n <- nrow(Xorig)
  p <- ncol(Xorig)

  X <- matrix(runif(n*p), ncol=p) ## initial design
  d <- distance(X)
  d <- d[upper.tri(d)]
  md <- min(d)
  if(!is.null(Xorig)) { ## new code
    md2 <- min(distance(X, Xorig))
    if(md2 < md) md <- md2
  }
  for(t in 1:T) {
    row <- sample(1:n, 1)
    xold <- X[row,] ## random row selection
    X[row,] <- runif(p) ## random new row
    d <- distance(X)
    d <- d[upper.tri(d)]
    mdprime <- min(d)
    if(!is.null(Xorig)) { ## new code
      mdprime2 <- min(distance(X, Xorig))
      if(mdprime2 < mdprime)
        mdprime <- mdprime2
    }
    if(mdprime > md) {
      md <- mdprime ## accept
    } else { X[row,] <- xold } ## reject
  }
  return(X)
}

#' maximin.multid
#'
#' maximin design for multi-dimensional dataset.
#'
#' @param Xorig The existing design.
#' @param T a number of iterations. Default is 10000*length(Xorig)
#' @return A multi-dimensional matrix of data with maximin design.
#' @importFrom plgp distance
#' @importFrom stats runif
#' @noRd
#'

maximin.multid <- function(Xorig, T=10*nrow(Xorig)) {

  n <- nrow(Xorig)
  p <- ncol(Xorig)

  X <- matrix(runif(n*p), ncol=p) ## initial design
  d <- distance(X)
  d <- d[upper.tri(d)]
  md <- min(d)
  if(!is.null(Xorig)) { ## new code
    md2 <- min(distance(X, Xorig))
    if(md2 < md) md <- md2
  }
  for(t in 1:T) {
    row <- sample(1:n, 1)
    xold <- X[row,] ## random row selection
    X[row,] <- runif(p) ## random new row
    d <- distance(X)
    d <- d[upper.tri(d)]
    mdprime <- min(d)
    if(!is.null(Xorig)) { ## new code
      mdprime2 <- min(distance(X, Xorig))
      if(mdprime2 < mdprime)
        mdprime <- mdprime2
    }
    if(mdprime > md) {
      md <- mdprime ## accept
    } else { X[row,] <- xold } ## reject
  }
  return(X)
}
