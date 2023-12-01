#' prediction of the RNAmf emulator with 2 or 3 fidelity levels.
#'
#' @description The function computes the posterior mean and variance of RNA models with two or three fidelity levels
#' by fitted model using \code{\link{RNAmf_two_level}} or \code{\link{RNAmf_three_level}}.
#'
#' @seealso \code{\link{RNAmf_two_level}} or \code{\link{RNAmf_three_level}} for the model.
#'
#' @details From the model fitted by \code{\link{RNAmf_two_level}} or \code{\link{RNAmf_three_level}},
#' the posterior mean and variance are calculated based on the closed form expression derived by a recursive fashion.
#' The formulas depend on its kernel choices.
#' For details, see Heo and Sung (2023+, <arXiv:2309.11772>).
#'
#' @param object a class \code{RNAmf} object fitted by \code{\link{RNAmf_two_level}} or \code{\link{RNAmf_three_level}}.
#' @param x vector or matrix of new input locations to predict.
#' @param ... for compatibility with generic method \code{predict}.
#'
#' @return
#' \itemize{
#'   \item \code{mu}: vector of predictive posterior mean.
#'   \item \code{sig2}: vector of predictive posterior variance.
#'   \item \code{time}: a scalar of the time for the computation.
#' }
#'
#' @importFrom plgp distance
#' @importFrom stats predict
#' @rdname predict
#' @title predict
#' @method predict RNAmf
#' @export
#' @examples
#' ### two levels example ###
#' library(lhs)
#'
#' ### Perdikaris function ###
#' f1 <- function(x) {
#'   sin(8 * pi * x)
#' }
#'
#' f2 <- function(x) {
#'   (x - sqrt(2)) * (sin(8 * pi * x))^2
#' }
#'
#' ### training data ###
#' n1 <- 13
#' n2 <- 8
#'
#' ### fix seed to reproduce the result ###
#' set.seed(1)
#'
#' ### generate initial nested design ###
#' X <- NestedX(c(n1, n2), 1)
#' X1 <- X[[1]]
#' X2 <- X[[2]]
#'
#' ### n1 and n2 might be changed from NestedX ###
#' ### assign n1 and n2 again ###
#' n1 <- nrow(X1)
#' n2 <- nrow(X2)
#'
#' y1 <- f1(X1)
#' y2 <- f2(X2)
#'
#' ### n=100 uniform test data ###
#' x <- seq(0, 1, length.out = 100)
#'
#' ### fit an RNAmf ###
#' fit.RNAmf <- RNAmf_two_level(X1, y1, X2, y2, kernel = "sqex")
#'
#' ### predict ###
#' predy <- predict(fit.RNAmf, x)$mu
#' predsig2 <- predict(fit.RNAmf, x)$sig2
#'
#' ### RMSE ###
#' print(sqrt(mean((predy - f2(x))^2)))
#'
#' ### visualize the emulation performance ###
#' plot(x, predy,
#'   type = "l", lwd = 2, col = 3, # emulator and confidence interval
#'   ylim = c(-2, 1)
#' )
#' lines(x, predy + 1.96 * sqrt(predsig2 * length(y2) / (length(y2) - 2)), col = 3, lty = 2)
#' lines(x, predy - 1.96 * sqrt(predsig2 * length(y2) / (length(y2) - 2)), col = 3, lty = 2)
#'
#' curve(f2(x), add = TRUE, col = 1, lwd = 2, lty = 2) # high fidelity function
#'
#' points(X1, y1, pch = 1, col = "red") # low-fidelity design
#' points(X2, y2, pch = 4, col = "blue") # high-fidelity design
#'
#' ### three levels example ###
#' ### Branin function ###
#' branin <- function(xx, l){
#'   x1 <- xx[1]
#'   x2 <- xx[2]
#'   if(l == 1){
#'     10*sqrt((-1.275*(1.2*x1+0.4)^2/pi^2+5*(1.2*x1+0.4)/pi+(1.2*x2+0.4)-6)^2 +
#'     (10-5/(4*pi))*cos((1.2*x1+0.4))+ 10) + 2*(1.2*x1+1.9) - 3*(3*(1.2*x2+2.4)-1) - 1 - 3*x2 + 1
#'   }else if(l == 2){
#'     10*sqrt((-1.275*(x1+2)^2/pi^2+5*(x1+2)/pi+(x2+2)-6)^2 +
#'     (10-5/(4*pi))*cos((x1+2))+ 10) + 2*(x1-0.5) - 3*(3*x2-1) - 1
#'   }else if(l == 3){
#'     (-1.275*x1^2/pi^2+5*x1/pi+x2-6)^2 + (10-5/(4*pi))*cos(x1)+ 10
#'   }
#' }
#'
#' output.branin <- function(x, l){
#'   factor_range <- list("x1" = c(-5, 10), "x2" = c(0, 15))
#'
#'   for(i in 1:length(factor_range)) x[i] <- factor_range[[i]][1] + x[i] * diff(factor_range[[i]])
#'   branin(x[1:2], l)
#' }
#'
#' ### training data ###
#' n1 <- 20; n2 <- 15; n3 <- 10
#'
#' ### fix seed to reproduce the result ###
#' set.seed(1)
#'
#' ### generate initial nested design ###
#' X <- NestedX(c(n1, n2, n3), 2)
#' X1 <- X[[1]]
#' X2 <- X[[2]]
#' X3 <- X[[3]]
#'
#' ### n1, n2 and n3 might be changed from NestedX ###
#' ### assign n1, n2 and n3 again ###
#' n1 <- nrow(X1)
#' n2 <- nrow(X2)
#' n3 <- nrow(X3)
#'
#' y1 <- apply(X1,1,output.branin, l=1)
#' y2 <- apply(X2,1,output.branin, l=2)
#' y3 <- apply(X3,1,output.branin, l=3)
#'
#' ### n=10000 grid test data ###
#' x <- as.matrix(expand.grid(seq(0, 1, length.out = 100),seq(0, 1, length.out = 100)))
#'
#' ### fit an RNAmf ###
#' fit.RNAmf <- RNAmf_three_level(X1, y1, X2, y2, X3, y3, kernel = "sqex")
#'
#' ### predict ###
#' pred.RNAmf <- predict(fit.RNAmf, x)
#' predy <- pred.RNAmf$mu
#' predsig2 <- pred.RNAmf$sig2
#'
#' ### RMSE ###
#' print(sqrt(mean((predy - apply(x,1,output.branin, l=3))^2)))
#'
#' ### visualize the emulation performance ###
#' x1 <- x2 <- seq(0, 1, length.out = 100)
#' par(mfrow=c(1,2))
#' image(x1, x2, matrix(apply(x,1,output.branin, l=3), ncol=100),
#' zlim=c(0,310), main="Branin function")
#' image(x1, x2, matrix(predy, ncol=100),
#' zlim=c(0,310), main="RNAmf prediction")
#'
#' ### predictive variance ###
#' print(predsig2)
#'
predict.RNAmf <- function(object, x, ...) {
  t1 <- proc.time()[3]
  ### check the object ###
  if (!inherits(object, "RNAmf")) {
    stop("The object is not of class \"RNAmf\" \n")
  }

  ### prediction ###
  if (object$level == 2) { # level 2
    kernel <- object$kernel
    constant <- object$constant
    fit1 <- object$fit1
    fit2 <- object$fit2

    if (kernel == "sqex") {
      if (constant) {
        d <- ncol(fit1$X)
        x <- matrix(x, ncol = d)
        pred.fit <- pred.GP(fit1, x)
        x.mu <- pred.fit$mu # mean of f1(u)
        sig2 <- pred.fit$sig2 #* 0

        ### calculate the closed form ###
        X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
        w1.x2 <- fit2$X[, d + 1]
        y2 <- fit2$y
        n <- length(y2)
        theta <- fit2$theta
        tau2hat <- fit2$tau2hat
        mu2 <- fit2$mu.hat

        Ci <- fit2$Ki
        a <- Ci %*% (y2 - mu2)

        ### scale new inputs ###
        x <- scale.inputs(x, attr(fit2$X, "scaled:center")[1:d], attr(fit2$X, "scaled:scale")[1:d])
        x.mu <- scale.inputs(x.mu, attr(fit2$X, "scaled:center")[d + 1], attr(fit2$X, "scaled:scale")[d + 1])
        sig2 <- sig2 / attr(fit2$X, "scaled:scale")[d + 1]^2

        # mean
        predy <- mu2 + (exp(-distance(t(t(x) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))) *
          1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
          exp(-(drop(outer(x.mu, w1.x2, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))) %*% a

        # var
        predsig2 <- c(rep(0, nrow(x)))
        for (i in 1:nrow(x)) { # each test point
          mat <- drop(exp(-distance(t(x[i, ]) / sqrt(theta[-(d + 1)]), t(t(X2) / sqrt(theta[-(d + 1)])))) %o%
            exp(-distance(t(x[i, ]) / sqrt(theta[-(d + 1)]), t(t(X2) / sqrt(theta[-(d + 1)]))))) * # common components
            1 / sqrt(1 + 4 * sig2[i] / theta[d + 1]) *
            exp(-(outer(w1.x2, w1.x2, FUN = "+") / 2 - x.mu[i])^2 / (theta[d + 1] / 2 + 2 * sig2[i])) *
            exp(-(outer(w1.x2, w1.x2, FUN = "-"))^2 / (2 * theta[d + 1]))

          predsig2[i] <- pmax(0, tau2hat - (predy[i] - mu2)^2 + drop(t(a) %*% mat %*% a) - tau2hat * sum(diag(Ci %*% mat)))
        }
      } else {
        d <- ncol(fit1$X)
        x <- matrix(x, ncol = d)
        pred.fit <- pred.GP(fit1, x)
        x.mu <- pred.fit$mu # mean of f1(u)
        sig2 <- pred.fit$sig2

        ### calculate the closed form ###
        X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
        w1.x2 <- fit2$X[, d + 1]
        y2 <- fit2$y
        n <- length(y2)
        theta <- fit2$theta
        tau2hat <- fit2$tau2hat

        Ci <- fit2$Ki
        a <- Ci %*% (y2 + attr(y2, "scaled:center"))

        ### scale new inputs ###
        x <- scale.inputs(x, attr(fit2$X, "scaled:center")[1:d], attr(fit2$X, "scaled:scale")[1:d])
        x.mu <- scale.inputs(x.mu, attr(fit2$X, "scaled:center")[d + 1], attr(fit2$X, "scaled:scale")[d + 1])
        sig2 <- sig2 / attr(fit2$X, "scaled:scale")[d + 1]^2

        # mean
        predy <- (exp(-distance(t(t(x) / sqrt(theta[-(d + 1)])), t(t(X2) / sqrt(theta[-(d + 1)])))) *
          1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
          exp(-(drop(outer(x.mu, w1.x2, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))) %*% a

        # var
        predsig2 <- c(rep(0, nrow(x)))
        for (i in 1:nrow(x)) { # each test point
          mat <- drop(exp(-distance(t(x[i, ]) / sqrt(theta[-(d + 1)]), t(t(X2) / sqrt(theta[-(d + 1)])))) %o%
            exp(-distance(t(x[i, ]) / sqrt(theta[-(d + 1)]), t(t(X2) / sqrt(theta[-(d + 1)]))))) * # common components
            1 / sqrt(1 + 4 * sig2[i] / theta[d + 1]) *
            exp(-(outer(w1.x2, w1.x2, FUN = "+") / 2 - x.mu[i])^2 / (theta[d + 1] / 2 + 2 * sig2[i])) *
            exp(-(outer(w1.x2, w1.x2, FUN = "-"))^2 / (2 * theta[d + 1]))

          predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a) %*% mat %*% a) - tau2hat * sum(diag(Ci %*% mat)))
        }
      }
    } else if (kernel == "matern1.5") {
      if (constant) {
        d <- ncol(fit1$X)
        x <- matrix(x, ncol = d)
        pred.fit <- pred.matGP(fit1, x)
        x.mu <- pred.fit$mu # mean of f1(u)
        sig2 <- pred.fit$sig2 #* 0

        ### calculate the closed form ###
        X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
        w1.x2 <- fit2$X[, d + 1]
        y2 <- fit2$y
        n <- length(y2)
        theta <- fit2$theta
        tau2hat <- fit2$tau2hat
        mu2 <- fit2$mu.hat

        Ci <- fit2$Ki
        a <- Ci %*% (y2 - mu2)

        ### scale new inputs ###
        x <- scale.inputs(x, attr(fit2$X, "scaled:center")[1:d], attr(fit2$X, "scaled:scale")[1:d])
        x.mu <- scale.inputs(x.mu, attr(fit2$X, "scaled:center")[d + 1], attr(fit2$X, "scaled:scale")[d + 1])
        sig2 <- sig2 / attr(fit2$X, "scaled:scale")[d + 1]^2

        # mean
        predy <- c(rep(0, nrow(x)))
        for (j in 1:nrow(x)) { # each test point
          mua <- x.mu[j] - sqrt(3) * sig2[j] / theta[d + 1]
          mub <- x.mu[j] + sqrt(3) * sig2[j] / theta[d + 1]

          lambda11 <- c(1, mua)
          lambda12 <- c(0, 1)
          lambda21 <- c(1, -mub)

          e1 <- cbind(matrix(1 - sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])
          e2 <- cbind(matrix(1 + sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])

          predy[j] <- mu2 + drop(t(t(cor.sep(t(x[j, ]), X2, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
            (exp((3 * sig2[j] + 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu[j])) / (2 * theta[d + 1]^2)) *
              (e1 %*% lambda11 * pnorm((mua - w1.x2) / sqrt(sig2[j])) +
                e1 %*% lambda12 * sqrt(sig2[j]) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2[j]))) +
              exp((3 * sig2[j] - 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu[j])) / (2 * theta[d + 1]^2)) *
                (e2 %*% lambda21 * pnorm((-mub + w1.x2) / sqrt(sig2[j])) +
                  e2 %*% lambda12 * sqrt(sig2[j]) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2[j]))))) %*% a)
        }

        # var
        predsig2 <- c(rep(0, nrow(x)))
        for (i in 1:nrow(x)) {
          zeta <- function(x, y) {
            zetafun(w1 = x, w2 = y, m = x.mu[i], s = sig2[i], nu = 1.5, theta = theta[d + 1])
          }

          mat <- drop(t(cor.sep(t(x[i, ]), X2, theta[-(d + 1)], nu = 1.5)) %o% t(cor.sep(t(x[i, ]), X2, theta[-(d + 1)], nu = 1.5))) * # constant depends on kernel
            outer(w1.x2, w1.x2, FUN = Vectorize(zeta))

          predsig2[i] <- pmax(0, tau2hat - (predy[i] - mu2)^2 + drop(t(a) %*% mat %*% a) - tau2hat * sum(diag(Ci %*% mat)))
        }
      } else {
        d <- ncol(fit1$X)
        x <- matrix(x, ncol = d)
        pred.fit <- pred.matGP(fit1, x)
        x.mu <- pred.fit$mu # mean of f1(u)
        sig2 <- pred.fit$sig2

        ### calculate the closed form ###
        X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
        w1.x2 <- fit2$X[, d + 1]
        y2 <- fit2$y
        n <- length(y2)
        theta <- fit2$theta
        tau2hat <- fit2$tau2hat

        Ci <- fit2$Ki
        a <- Ci %*% (y2 + attr(y2, "scaled:center"))

        ### scale new inputs ###
        x <- scale.inputs(x, attr(fit2$X, "scaled:center")[1:d], attr(fit2$X, "scaled:scale")[1:d])
        x.mu <- scale.inputs(x.mu, attr(fit2$X, "scaled:center")[d + 1], attr(fit2$X, "scaled:scale")[d + 1])
        sig2 <- sig2 / attr(fit2$X, "scaled:scale")[d + 1]^2

        # mean
        predy <- c(rep(0, nrow(x)))
        for (j in 1:nrow(x)) { # each test point
          mua <- x.mu[j] - sqrt(3) * sig2[j] / theta[d + 1]
          mub <- x.mu[j] + sqrt(3) * sig2[j] / theta[d + 1]

          lambda11 <- c(1, mua)
          lambda12 <- c(0, 1)
          lambda21 <- c(1, -mub)

          e1 <- cbind(matrix(1 - sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])
          e2 <- cbind(matrix(1 + sqrt(3) * w1.x2 / theta[d + 1]), sqrt(3) / theta[d + 1])

          predy[j] <- drop(t(t(cor.sep(t(x[j, ]), X2, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
            (exp((3 * sig2[j] + 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu[j])) / (2 * theta[d + 1]^2)) *
              (e1 %*% lambda11 * pnorm((mua - w1.x2) / sqrt(sig2[j])) +
                e1 %*% lambda12 * sqrt(sig2[j]) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2[j]))) +
              exp((3 * sig2[j] - 2 * sqrt(3) * theta[d + 1] * (w1.x2 - x.mu[j])) / (2 * theta[d + 1]^2)) *
                (e2 %*% lambda21 * pnorm((-mub + w1.x2) / sqrt(sig2[j])) +
                  e2 %*% lambda12 * sqrt(sig2[j]) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2[j]))))) %*% a)
        }

        # var
        predsig2 <- c(rep(0, nrow(x)))
        for (i in 1:nrow(x)) {
          zeta <- function(x, y) {
            zetafun(w1 = x, w2 = y, m = x.mu[i], s = sig2[i], nu = 1.5, theta = theta[d + 1])
          }

          mat <- drop(t(cor.sep(t(x[i, ]), X2, theta[-(d + 1)], nu = 1.5)) %o% t(cor.sep(t(x[i, ]), X2, theta[-(d + 1)], nu = 1.5))) * # constant depends on kernel
            outer(w1.x2, w1.x2, FUN = Vectorize(zeta))

          predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a) %*% mat %*% a) - tau2hat * sum(diag(Ci %*% mat)))
        }
      }
    } else if (kernel == "matern2.5") {
      if (constant) {
        d <- ncol(fit1$X)
        x <- matrix(x, ncol = d)
        pred.fit <- pred.matGP(fit1, x)
        x.mu <- pred.fit$mu # mean of f1(u)
        sig2 <- pred.fit$sig2 #* 0

        ### calculate the closed form ###
        X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
        w1.x2 <- fit2$X[, d + 1]
        y2 <- fit2$y
        n <- length(y2)
        theta <- fit2$theta
        tau2hat <- fit2$tau2hat
        mu2 <- fit2$mu.hat

        Ci <- fit2$Ki
        a <- Ci %*% (y2 - mu2)

        ### scale new inputs ###
        x <- scale.inputs(x, attr(fit2$X, "scaled:center")[1:d], attr(fit2$X, "scaled:scale")[1:d])
        x.mu <- scale.inputs(x.mu, attr(fit2$X, "scaled:center")[d + 1], attr(fit2$X, "scaled:scale")[d + 1])
        sig2 <- sig2 / attr(fit2$X, "scaled:scale")[d + 1]^2

        # mean
        predy <- c(rep(0, nrow(x)))
        for (j in 1:nrow(x)) { # each test point
          mua <- x.mu[j] - sqrt(5) * sig2[j] / theta[d + 1]
          mub <- x.mu[j] + sqrt(5) * sig2[j] / theta[d + 1]

          lambda11 <- c(1, mua, mua^2 + sig2[j])
          lambda12 <- cbind(0, 1, matrix(mua + w1.x2))
          lambda21 <- c(1, -mub, mub^2 + sig2[j])
          lambda22 <- cbind(0, 1, matrix(-mub - w1.x2))

          e1 <- cbind(
            matrix(1 - sqrt(5) * w1.x2 / theta[d + 1] + 5 * w1.x2^2 / (3 * theta[d + 1]^2)),
            matrix(sqrt(5) / theta[d + 1] - 10 * w1.x2 / (3 * theta[d + 1]^2)),
            5 / (3 * theta[d + 1]^2)
          )
          e2 <- cbind(
            matrix(1 + sqrt(5) * w1.x2 / theta[d + 1] + 5 * w1.x2^2 / (3 * theta[d + 1]^2)),
            matrix(sqrt(5) / theta[d + 1] + 10 * w1.x2 / (3 * theta[d + 1]^2)),
            5 / (3 * theta[d + 1]^2)
          )

          predy[j] <- mu2 + drop(t(t(cor.sep(t(x[j, ]), X2, theta[-(d + 1)], nu = 2.5)) *
            (exp((5 * sig2[j] + 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu[j])) / (2 * theta[d + 1]^2)) *
              (e1 %*% lambda11 * pnorm((mua - w1.x2) / sqrt(sig2[j])) +
                rowSums(e1 * lambda12) * sqrt(sig2[j]) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2[j]))) +
              exp((5 * sig2[j] - 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu[j])) / (2 * theta[d + 1]^2)) *
                (e2 %*% lambda21 * pnorm((-mub + w1.x2) / sqrt(sig2[j])) +
                  rowSums(e2 * lambda22) * sqrt(sig2[j]) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2[j]))))) %*% a)
        }

        # var
        predsig2 <- c(rep(0, nrow(x)))
        for (i in 1:nrow(x)) {
          zeta <- function(x, y) {
            zetafun(w1 = x, w2 = y, m = x.mu[i], s = sig2[i], nu = 2.5, theta = theta[d + 1])
          }

          mat <- drop(t(cor.sep(t(x[i, ]), X2, theta[-(d + 1)], nu = 2.5)) %o% t(cor.sep(t(x[i, ]), X2, theta[-(d + 1)], nu = 2.5))) * # constant depends on kernel
            outer(w1.x2, w1.x2, FUN = Vectorize(zeta))

          predsig2[i] <- pmax(0, tau2hat - (predy[i] - mu2)^2 + drop(t(a) %*% mat %*% a) - tau2hat * sum(diag(Ci %*% mat)))
        }
      } else {
        d <- ncol(fit1$X)
        x <- matrix(x, ncol = d)
        pred.fit <- pred.matGP(fit1, x)
        x.mu <- pred.fit$mu # mean of f1(u)
        sig2 <- pred.fit$sig2

        ### calculate the closed form ###
        X2 <- matrix(fit2$X[, -(d + 1)], ncol = d)
        w1.x2 <- fit2$X[, d + 1]
        y2 <- fit2$y
        n <- length(y2)
        theta <- fit2$theta
        tau2hat <- fit2$tau2hat

        Ci <- fit2$Ki
        a <- Ci %*% (y2 + attr(y2, "scaled:center"))

        ### scale new inputs ###
        x <- scale.inputs(x, attr(fit2$X, "scaled:center")[1:d], attr(fit2$X, "scaled:scale")[1:d])
        x.mu <- scale.inputs(x.mu, attr(fit2$X, "scaled:center")[d + 1], attr(fit2$X, "scaled:scale")[d + 1])
        sig2 <- sig2 / attr(fit2$X, "scaled:scale")[d + 1]^2

        # mean
        predy <- c(rep(0, nrow(x)))

        for (j in 1:nrow(x)) { # each test point
          mua <- x.mu[j] - sqrt(5) * sig2[j] / theta[d + 1]
          mub <- x.mu[j] + sqrt(5) * sig2[j] / theta[d + 1]

          lambda11 <- c(1, mua, mua^2 + sig2[j])
          lambda12 <- cbind(0, 1, matrix(mua + w1.x2))
          lambda21 <- c(1, -mub, mub^2 + sig2[j])
          lambda22 <- cbind(0, 1, matrix(-mub - w1.x2))

          e1 <- cbind(
            matrix(1 - sqrt(5) * w1.x2 / theta[d + 1] + 5 * w1.x2^2 / (3 * theta[d + 1]^2)),
            matrix(sqrt(5) / theta[d + 1] - 10 * w1.x2 / (3 * theta[d + 1]^2)),
            5 / (3 * theta[d + 1]^2)
          )
          e2 <- cbind(
            matrix(1 + sqrt(5) * w1.x2 / theta[d + 1] + 5 * w1.x2^2 / (3 * theta[d + 1]^2)),
            matrix(sqrt(5) / theta[d + 1] + 10 * w1.x2 / (3 * theta[d + 1]^2)),
            5 / (3 * theta[d + 1]^2)
          )

          predy[j] <- drop(t(t(cor.sep(t(x[j, ]), X2, theta[-(d + 1)], nu = 2.5)) *
            (exp((5 * sig2[j] + 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu[j])) / (2 * theta[d + 1]^2)) *
              (e1 %*% lambda11 * pnorm((mua - w1.x2) / sqrt(sig2[j])) +
                rowSums(e1 * lambda12) * sqrt(sig2[j]) / sqrt(2 * pi) * exp(-(w1.x2 - mua)^2 / (2 * sig2[j]))) +
              exp((5 * sig2[j] - 2 * sqrt(5) * theta[d + 1] * (w1.x2 - x.mu[j])) / (2 * theta[d + 1]^2)) *
                (e2 %*% lambda21 * pnorm((-mub + w1.x2) / sqrt(sig2[j])) +
                  rowSums(e2 * lambda22) * sqrt(sig2[j]) / sqrt(2 * pi) * exp(-(w1.x2 - mub)^2 / (2 * sig2[j]))))) %*% a)
        }

        # var
        predsig2 <- c(rep(0, nrow(x)))
        for (i in 1:nrow(x)) {
          zeta <- function(x, y) {
            zetafun(w1 = x, w2 = y, m = x.mu[i], s = sig2[i], nu = 2.5, theta = theta[d + 1])
          }

          mat <- drop(t(cor.sep(t(x[i, ]), X2, theta[-(d + 1)], nu = 2.5)) %o% t(cor.sep(t(x[i, ]), X2, theta[-(d + 1)], nu = 2.5))) * # constant depends on kernel
            outer(w1.x2, w1.x2, FUN = Vectorize(zeta))

          predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a) %*% mat %*% a) - tau2hat * sum(diag(Ci %*% mat)))
        }
      }
    }
  } else if (object$level == 3) { # level 3
    kernel <- object$kernel
    constant <- object$constant
    fit.RNAmf_two_level <- object$fit.RNAmf_two_level
    fit1 <- fit.RNAmf_two_level$fit1
    fit2 <- fit.RNAmf_two_level$fit2
    fit3 <- object$fit3

    if (kernel == "sqex") {
      if (constant) {
        pred.RNAmf_two_level <- predict(fit.RNAmf_two_level, x)
        x.mu <- pred.RNAmf_two_level$mu
        sig2 <- pred.RNAmf_two_level$sig2

        d <- ncol(fit1$X)
        x <- matrix(x, ncol = d)

        ### calculate the closed form ###
        X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
        w2.x3 <- fit3$X[, d + 1]
        y3 <- fit3$y
        n <- length(y3)
        theta <- fit3$theta
        tau2hat <- fit3$tau2hat
        mu3 <- fit3$mu.hat

        Ci <- fit3$Ki
        a <- Ci %*% (y3 - mu3)

        ### scale new inputs ###
        x <- scale.inputs(x, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
        x.mu <- scale.inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
        sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

        # mean
        predy <- mu3 + (exp(-distance(t(t(x) / sqrt(theta[-(d + 1)])), t(t(X3) / sqrt(theta[-(d + 1)])))) *
          1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
          exp(-(drop(outer(x.mu, w2.x3, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))) %*% a

        # var
        predsig2 <- c(rep(0, nrow(x)))
        for (i in 1:nrow(x)) { # each test point
          mat <- drop(exp(-distance(t(x[i, ]) / sqrt(theta[-(d + 1)]), t(t(X3) / sqrt(theta[-(d + 1)])))) %o%
            exp(-distance(t(x[i, ]) / sqrt(theta[-(d + 1)]), t(t(X3) / sqrt(theta[-(d + 1)]))))) * # common components
            1 / sqrt(1 + 4 * sig2[i] / theta[d + 1]) *
            exp(-(outer(w2.x3, w2.x3, FUN = "+") / 2 - x.mu[i])^2 / (theta[d + 1] / 2 + 2 * sig2[i])) *
            exp(-(outer(w2.x3, w2.x3, FUN = "-"))^2 / (2 * theta[d + 1]))

          predsig2[i] <- pmax(0, tau2hat - (predy[i] - mu3)^2 + drop(t(a) %*% mat %*% a) - tau2hat * sum(diag(Ci %*% mat)))
        }
      } else {
        pred.RNAmf_two_level <- predict(fit.RNAmf_two_level, x)
        x.mu <- pred.RNAmf_two_level$mu
        sig2 <- pred.RNAmf_two_level$sig2

        d <- ncol(fit1$X)
        x <- matrix(x, ncol = d)

        ### calculate the closed form ###
        X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
        w2.x3 <- fit3$X[, d + 1]
        y3 <- fit3$y
        n <- length(y3)
        theta <- fit3$theta
        tau2hat <- fit3$tau2hat

        Ci <- fit3$Ki
        a <- Ci %*% (y3 + attr(y3, "scaled:center"))

        ### scale new inputs ###
        x <- scale.inputs(x, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
        x.mu <- scale.inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
        sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

        # mean
        predy <- (exp(-distance(t(t(x) / sqrt(theta[-(d + 1)])), t(t(X3) / sqrt(theta[-(d + 1)])))) *
          1 / sqrt(1 + 2 * sig2 / theta[d + 1]) *
          exp(-(drop(outer(x.mu, w2.x3, FUN = "-")))^2 / (theta[d + 1] + 2 * sig2))) %*% a

        # var
        predsig2 <- c(rep(0, nrow(x)))
        for (i in 1:nrow(x)) { # each test point
          mat <- drop(exp(-distance(t(x[i, ]) / sqrt(theta[-(d + 1)]), t(t(X3) / sqrt(theta[-(d + 1)])))) %o%
            exp(-distance(t(x[i, ]) / sqrt(theta[-(d + 1)]), t(t(X3) / sqrt(theta[-(d + 1)]))))) * # common components
            1 / sqrt(1 + 4 * sig2[i] / theta[d + 1]) *
            exp(-(outer(w2.x3, w2.x3, FUN = "+") / 2 - x.mu[i])^2 / (theta[d + 1] / 2 + 2 * sig2[i])) *
            exp(-(outer(w2.x3, w2.x3, FUN = "-"))^2 / (2 * theta[d + 1]))

          predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a) %*% mat %*% a) - tau2hat * sum(diag(Ci %*% mat)))
        }
      }
    } else if (kernel == "matern1.5") {
      if (constant) {
        d <- ncol(fit1$X)
        x <- matrix(x, ncol = d)
        pred.RNAmf_two_level <- predict(fit.RNAmf_two_level, x)
        x.mu <- pred.RNAmf_two_level$mu
        sig2 <- pred.RNAmf_two_level$sig2

        ### calculate the closed form ###
        X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
        w2.x3 <- fit3$X[, d + 1]
        y3 <- fit3$y
        n <- length(y3)
        theta <- fit3$theta
        tau2hat <- fit3$tau2hat
        mu3 <- fit3$mu.hat

        Ci <- fit3$Ki
        a <- Ci %*% (y3 - mu3)

        ### scale new inputs ###
        x <- scale.inputs(x, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
        x.mu <- scale.inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
        sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

        # mean
        predy <- c(rep(0, nrow(x)))
        for (j in 1:nrow(x)) { # each test point
          mua <- x.mu[j] - sqrt(3) * sig2[j] / theta[d + 1]
          mub <- x.mu[j] + sqrt(3) * sig2[j] / theta[d + 1]

          lambda11 <- c(1, mua)
          lambda12 <- c(0, 1)
          lambda21 <- c(1, -mub)

          e1 <- cbind(matrix(1 - sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])
          e2 <- cbind(matrix(1 + sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])

          predy[j] <- mu3 + drop(t(t(cor.sep(t(x[j, ]), X3, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
            (exp((3 * sig2[j] + 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu[j])) / (2 * theta[d + 1]^2)) *
              (e1 %*% lambda11 * pnorm((mua - w2.x3) / sqrt(sig2[j])) +
                e1 %*% lambda12 * sqrt(sig2[j]) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2[j]))) +
              exp((3 * sig2[j] - 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu[j])) / (2 * theta[d + 1]^2)) *
                (e2 %*% lambda21 * pnorm((-mub + w2.x3) / sqrt(sig2[j])) +
                  e2 %*% lambda12 * sqrt(sig2[j]) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2[j]))))) %*% a)
        }


        # var
        predsig2 <- c(rep(0, nrow(x)))
        for (i in 1:nrow(x)) {
          zeta <- function(x, y) {
            zetafun(w1 = x, w2 = y, m = x.mu[i], s = sig2[i], nu = 1.5, theta = theta[d + 1])
          }

          mat <- drop(t(cor.sep(t(x[i, ]), X3, theta[-(d + 1)], nu = 1.5)) %o% t(cor.sep(t(x[i, ]), X3, theta[-(d + 1)], nu = 1.5))) * # constant depends on kernel
            outer(w2.x3, w2.x3, FUN = Vectorize(zeta))

          predsig2[i] <- pmax(0, tau2hat - (predy[i] - mu3)^2 + drop(t(a) %*% mat %*% a) - tau2hat * sum(diag(Ci %*% mat)))
        }
      } else {
        d <- ncol(fit1$X)
        x <- matrix(x, ncol = d)
        pred.RNAmf_two_level <- predict(fit.RNAmf_two_level, x)
        x.mu <- pred.RNAmf_two_level$mu
        sig2 <- pred.RNAmf_two_level$sig2

        ### calculate the closed form ###
        X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
        w2.x3 <- fit3$X[, d + 1]
        y3 <- fit3$y
        n <- length(y3)
        theta <- fit3$theta
        tau2hat <- fit3$tau2hat

        Ci <- fit3$Ki
        a <- Ci %*% (y3 + attr(y3, "scaled:center"))

        ### scale new inputs ###
        x <- scale.inputs(x, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
        x.mu <- scale.inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
        sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

        # mean
        predy <- c(rep(0, nrow(x)))
        for (j in 1:nrow(x)) { # each test point
          mua <- x.mu[j] - sqrt(3) * sig2[j] / theta[d + 1]
          mub <- x.mu[j] + sqrt(3) * sig2[j] / theta[d + 1]

          lambda11 <- c(1, mua)
          lambda12 <- c(0, 1)
          lambda21 <- c(1, -mub)

          e1 <- cbind(matrix(1 - sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])
          e2 <- cbind(matrix(1 + sqrt(3) * w2.x3 / theta[d + 1]), sqrt(3) / theta[d + 1])

          predy[j] <- drop(t(t(cor.sep(t(x[j, ]), X3, theta[-(d + 1)], nu = 1.5)) * # common but depends on kernel
            (exp((3 * sig2[j] + 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu[j])) / (2 * theta[d + 1]^2)) *
              (e1 %*% lambda11 * pnorm((mua - w2.x3) / sqrt(sig2[j])) +
                e1 %*% lambda12 * sqrt(sig2[j]) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2[j]))) +
              exp((3 * sig2[j] - 2 * sqrt(3) * theta[d + 1] * (w2.x3 - x.mu[j])) / (2 * theta[d + 1]^2)) *
                (e2 %*% lambda21 * pnorm((-mub + w2.x3) / sqrt(sig2[j])) +
                  e2 %*% lambda12 * sqrt(sig2[j]) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2[j]))))) %*% a)
        }

        # var
        predsig2 <- c(rep(0, nrow(x)))
        for (i in 1:nrow(x)) {
          zeta <- function(x, y) {
            zetafun(w1 = x, w2 = y, m = x.mu[i], s = sig2[i], nu = 1.5, theta = theta[d + 1])
          }

          mat <- drop(t(cor.sep(t(x[i, ]), X3, theta[-(d + 1)], nu = 1.5)) %o% t(cor.sep(t(x[i, ]), X3, theta[-(d + 1)], nu = 1.5))) * # constant depends on kernel
            outer(w2.x3, w2.x3, FUN = Vectorize(zeta))

          predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a) %*% mat %*% a) - tau2hat * sum(diag(Ci %*% mat)))
        }
      }
    } else if (kernel == "matern2.5") {
      if (constant) {
        d <- ncol(fit1$X)
        x <- matrix(x, ncol = d)
        pred.RNAmf_two_level <- predict(fit.RNAmf_two_level, x)
        x.mu <- pred.RNAmf_two_level$mu
        sig2 <- pred.RNAmf_two_level$sig2

        ### calculate the closed form ###
        X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
        w2.x3 <- fit3$X[, d + 1]
        y3 <- fit3$y
        n <- length(y3)
        theta <- fit3$theta
        tau2hat <- fit3$tau2hat
        mu3 <- fit3$mu.hat

        Ci <- fit3$Ki
        a <- Ci %*% (y3 - mu3)

        ### scale new inputs ###
        x <- scale.inputs(x, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
        x.mu <- scale.inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
        sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

        # mean
        predy <- c(rep(0, nrow(x)))
        for (j in 1:nrow(x)) { # each test point
          mua <- x.mu[j] - sqrt(5) * sig2[j] / theta[d + 1]
          mub <- x.mu[j] + sqrt(5) * sig2[j] / theta[d + 1]

          lambda11 <- c(1, mua, mua^2 + sig2[j])
          lambda12 <- cbind(0, 1, matrix(mua + w2.x3))
          lambda21 <- c(1, -mub, mub^2 + sig2[j])
          lambda22 <- cbind(0, 1, matrix(-mub - w2.x3))

          e1 <- cbind(
            matrix(1 - sqrt(5) * w2.x3 / theta[d + 1] + 5 * w2.x3^2 / (3 * theta[d + 1]^2)),
            matrix(sqrt(5) / theta[d + 1] - 10 * w2.x3 / (3 * theta[d + 1]^2)),
            5 / (3 * theta[d + 1]^2)
          )
          e2 <- cbind(
            matrix(1 + sqrt(5) * w2.x3 / theta[d + 1] + 5 * w2.x3^2 / (3 * theta[d + 1]^2)),
            matrix(sqrt(5) / theta[d + 1] + 10 * w2.x3 / (3 * theta[d + 1]^2)),
            5 / (3 * theta[d + 1]^2)
          )

          predy[j] <- mu3 + drop(t(t(cor.sep(t(x[j, ]), X3, theta[-(d + 1)], nu = 2.5)) *
            (exp((5 * sig2[j] + 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu[j])) / (2 * theta[d + 1]^2)) *
              (e1 %*% lambda11 * pnorm((mua - w2.x3) / sqrt(sig2[j])) +
                rowSums(e1 * lambda12) * sqrt(sig2[j]) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2[j]))) +
              exp((5 * sig2[j] - 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu[j])) / (2 * theta[d + 1]^2)) *
                (e2 %*% lambda21 * pnorm((-mub + w2.x3) / sqrt(sig2[j])) +
                  rowSums(e2 * lambda22) * sqrt(sig2[j]) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2[j]))))) %*% a)
        }

        # var
        predsig2 <- c(rep(0, nrow(x)))
        for (i in 1:nrow(x)) {
          zeta <- function(x, y) {
            zetafun(w1 = x, w2 = y, m = x.mu[i], s = sig2[i], nu = 2.5, theta = theta[d + 1])
          }

          mat <- drop(t(cor.sep(t(x[i, ]), X3, theta[-(d + 1)], nu = 2.5)) %o% t(cor.sep(t(x[i, ]), X3, theta[-(d + 1)], nu = 2.5))) * # constant depends on kernel
            outer(w2.x3, w2.x3, FUN = Vectorize(zeta))

          predsig2[i] <- pmax(0, tau2hat - (predy[i] - mu3)^2 + drop(t(a) %*% mat %*% a) - tau2hat * sum(diag(Ci %*% mat)))
        }
      } else {
        d <- ncol(fit1$X)
        x <- matrix(x, ncol = d)
        pred.RNAmf_two_level <- predict(fit.RNAmf_two_level, x)
        x.mu <- pred.RNAmf_two_level$mu
        sig2 <- pred.RNAmf_two_level$sig2

        ### calculate the closed form ###
        X3 <- matrix(fit3$X[, -(d + 1)], ncol = d)
        w2.x3 <- fit3$X[, d + 1]
        y3 <- fit3$y
        n <- length(y3)
        theta <- fit3$theta
        tau2hat <- fit3$tau2hat

        Ci <- fit3$Ki
        a <- Ci %*% (y3 + attr(y3, "scaled:center"))

        ### scale new inputs ###
        x <- scale.inputs(x, attr(fit3$X, "scaled:center")[1:d], attr(fit3$X, "scaled:scale")[1:d])
        x.mu <- scale.inputs(x.mu, attr(fit3$X, "scaled:center")[d + 1], attr(fit3$X, "scaled:scale")[d + 1])
        sig2 <- sig2 / attr(fit3$X, "scaled:scale")[d + 1]^2

        # mean
        predy <- c(rep(0, nrow(x)))
        for (j in 1:nrow(x)) { # each test point
          mua <- x.mu[j] - sqrt(5) * sig2[j] / theta[d + 1]
          mub <- x.mu[j] + sqrt(5) * sig2[j] / theta[d + 1]

          lambda11 <- c(1, mua, mua^2 + sig2[j])
          lambda12 <- cbind(0, 1, matrix(mua + w2.x3))
          lambda21 <- c(1, -mub, mub^2 + sig2[j])
          lambda22 <- cbind(0, 1, matrix(-mub - w2.x3))

          e1 <- cbind(
            matrix(1 - sqrt(5) * w2.x3 / theta[d + 1] + 5 * w2.x3^2 / (3 * theta[d + 1]^2)),
            matrix(sqrt(5) / theta[d + 1] - 10 * w2.x3 / (3 * theta[d + 1]^2)),
            5 / (3 * theta[d + 1]^2)
          )
          e2 <- cbind(
            matrix(1 + sqrt(5) * w2.x3 / theta[d + 1] + 5 * w2.x3^2 / (3 * theta[d + 1]^2)),
            matrix(sqrt(5) / theta[d + 1] + 10 * w2.x3 / (3 * theta[d + 1]^2)),
            5 / (3 * theta[d + 1]^2)
          )

          predy[j] <- drop(t(t(cor.sep(t(x[j, ]), X3, theta[-(d + 1)], nu = 2.5)) *
            (exp((5 * sig2[j] + 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu[j])) / (2 * theta[d + 1]^2)) *
              (e1 %*% lambda11 * pnorm((mua - w2.x3) / sqrt(sig2[j])) +
                rowSums(e1 * lambda12) * sqrt(sig2[j]) / sqrt(2 * pi) * exp(-(w2.x3 - mua)^2 / (2 * sig2[j]))) +
              exp((5 * sig2[j] - 2 * sqrt(5) * theta[d + 1] * (w2.x3 - x.mu[j])) / (2 * theta[d + 1]^2)) *
                (e2 %*% lambda21 * pnorm((-mub + w2.x3) / sqrt(sig2[j])) +
                  rowSums(e2 * lambda22) * sqrt(sig2[j]) / sqrt(2 * pi) * exp(-(w2.x3 - mub)^2 / (2 * sig2[j]))))) %*% a)
        }

        # var
        predsig2 <- c(rep(0, nrow(x)))
        for (i in 1:nrow(x)) {
          zeta <- function(x, y) {
            zetafun(w1 = x, w2 = y, m = x.mu[i], s = sig2[i], nu = 2.5, theta = theta[d + 1])
          }

          mat <- drop(t(cor.sep(t(x[i, ]), X3, theta[-(d + 1)], nu = 2.5)) %o% t(cor.sep(t(x[i, ]), X3, theta[-(d + 1)], nu = 2.5))) * # constant depends on kernel
            outer(w2.x3, w2.x3, FUN = Vectorize(zeta))

          predsig2[i] <- pmax(0, tau2hat - predy[i]^2 + drop(t(a) %*% mat %*% a) - tau2hat * sum(diag(Ci %*% mat)))
        }
      }
    }
  } else {
    stop("level is missing")
  }

  return(list(mu = predy, sig2 = predsig2, time = proc.time()[3] - t1))
}
