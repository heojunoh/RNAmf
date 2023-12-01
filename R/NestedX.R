#' subsetting matrix to be nested
#'
#' @param X1 matrix of the design at the lower level
#' @param X2 matrix of the design at the higher level
#' @return A list containing the design and the length of its subset
#' \itemize{
#'   \item \code{X}: matrix of the design at the lower level.
#'   \item \code{le}: length of the indices of subset.
#' }
#' @noRd
#'

subsetX <- function(X1 = NULL, X2 = NULL) {
  d <- dim(as.matrix(X2))[2] # d2
  n2 <- dim(as.matrix(X2))[1] # n2
  n1 <- dim(as.matrix(X1))[1] # n1

  dist <- 0
  for (i in 1:d) {
    grid <- expand.grid(X2[, i], X1[, i])
    dist <- dist + (grid[, 1] - grid[, 2])^2
  }
  dist.mat <- matrix(dist, n2, n1)
  indice <- max.col(-(dist.mat)) # find the minimum distance of column at each row

  X1 <- matrix(X1[-indice, ], ncol=d)
  X1 <- rbind(X1, matrix(X2, ncol=d))

  return(list(X = X1, le = length(indice)))
}

#' Constructing the nested design sets for RNA model.
#'
#' @description The function constructs the nested design sets with two fidelity levels
#' \eqn{\mathcal{X}_2 \subseteq \mathcal{X}_{1}} for \code{\link{RNAmf_two_level}} or
#' three fidelity levels \eqn{\mathcal{X}_3 \subseteq \mathcal{X}_2 \subseteq \mathcal{X}_{1}}
#' for \code{\link{RNAmf_three_level}}.
#'
#' @details The procedure replace the points of lower level design \eqn{\mathcal{X}_{l-1}}
#' to the closest points of higher level design \eqn{\mathcal{X}_{l}}.
#' The length of the \eqn{\mathcal{X}_{l-1}} could be larger than the user specified.
#' For details, see "\href{http://cran.nexr.com/web/packages/MuFiCokriging/MuFiCokriging.pdf}{\code{NestedDesign}}".
#'
#' @references
#' L. Le Gratiet and J. Garnier (2014). Recursive co-kriging model for design of computer experiments
#' with multiple levels of fidelity. \emph{International Journal for Uncertainty Quantification}, 4(5), 365-386;
#' doi:10.1615/Int.J.UncertaintyQuantification.2014006914
#'
#' @param n vector of the number of design points at each fidelity level \eqn{l}. Thus, the vector must have a positive value \eqn{n_1, n_2} or \eqn{n_1, n_2, n_3} where \eqn{n_1 > n_2 > n_3}.
#' @param d constant of the dimension of the design.
#' @return A list containing the design at each level, i.e., \eqn{\mathcal{X}_{1}, \mathcal{X}_{2}} or \eqn{\mathcal{X}_{1}, \mathcal{X}_{2}, \mathcal{X}_{3}}.
#' @export
#' @examples
#' ### number of design points ###
#' n1 <- 30
#' n2 <- 15
#'
#' ### dimension of the design ###
#' d <- 2
#'
#' ### fix seed to reproduce the result ###
#' set.seed(1)
#'
#' ### generate the nested design ###
#' NX <- NestedX(c(n1, n2), d)
#'
#' ### visualize nested design ###
#' plot(NX[[1]], col="red", pch=1, xlab="x1", ylab="x2")
#' points(NX[[2]], col="blue", pch=4)
#'


NestedX <- function(n, d) { # n; vector, d; dim

  if (is.unsorted(rev(n), strictly = TRUE)) {
    stop("The number of design at each level must be descending order \n")
  }
  if (!all(n > 0)) {
    stop("The number of design at each level must be positive \n")
  }
  if (length(d) != 1) {
    stop("The dimension of design at each level must be same \n")
  }

  ### generating initial designs ###
  level <- length(n)
  if (level < 2) {
    stop("The level of design should be larger than 1 \n")
  } else if (level == 2) {
    X1 <- maximinLHS(n[1], d)
    X2 <- maximinLHS(n[2], d)
    X <- list(X1, X2)
  } else if (level == 3) {
    X1 <- maximinLHS(n[1], d)
    X2 <- maximinLHS(n[2], d)
    X3 <- maximinLHS(n[3], d)
    X <- list(X1, X2, X3)
  }

  ### subsetting designs ###
  indices <- list()
  for (i in (level - 1):1) {
    SB <- subsetX(matrix(X[[i]], ncol=d), matrix(X[[i + 1]], ncol=d))
    X[[i]] <- SB$X
    n <- dim(SB$X)[1]
    indices[[i]] <- seq(n - SB$le + 1, n, by = 1)
  }

  ### assigning the design, indices, and the number of data at each level ###
  n <- c()
  for (i in 1:(level - 1)) {
    n[i] <- length(indices[[i]])
  }
  X.nested <- list()
  X.nested$X <- X[[1]]
  X.nested$ind <- indices
  X.nested$n <- n


  if (level == 2) {
    X <- list(X.nested$X, matrix(X.nested$X[X.nested$ind[[1]], ], ncol = d))
  } else if (level == 3) {
    X <- list(
      X.nested$X, matrix(X.nested$X[X.nested$ind[[1]], ], ncol = d),
      matrix(X.nested$X[X.nested$ind[[1]][X.nested$ind[[2]]], ], ncol = d)
    )
  }

  return(X = X)
}
