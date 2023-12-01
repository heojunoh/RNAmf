#' object to optimize the next point by ALC criterion updating at level 1 with two levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param Xref vector or matrix of reference data.
#' @param fit an object of class RNAmf.
#' @param mc.sample a number of mc samples generated for this approach. Default is 100.
#' @param parallel logical indicating whether to run parallel or not. Default is FALSE.
#' @param ncore the number of core for parallel. Default is 1.
#' @return A mean of the deduced variance at Xref.
#' @importFrom plgp covar.sep
#' @importFrom stats rnorm
#' @noRd
#'

obj.ALC_level_1 <- function(Xcand, Xref, fit, mc.sample, parallel = FALSE, ncore = 1) {
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

  Xcand <- matrix(Xcand, nrow = 1)
  if (kernel == "sqex") {
    y1.sample <- rnorm(mc.sample, mean = pred.GP(f1, Xcand)$mu, sd = sqrt(pred.GP(f1, Xcand)$sig2))
  } else if (kernel == "matern1.5") {
    y1.sample <- rnorm(mc.sample, mean = pred.matGP(f1, Xcand)$mu, sd = sqrt(pred.matGP(f1, Xcand)$sig2))
  } else if (kernel == "matern2.5") {
    y1.sample <- rnorm(mc.sample, mean = pred.matGP(f1, Xcand)$mu, sd = sqrt(pred.matGP(f1, Xcand)$sig2))
  }

  Xcand <- scale.inputs(Xcand, x.center1, x.scale1)

  ### Choose level 1 ###
  ### update Ki1
  if (kernel == "sqex") {
    cov.newx1 <- covar.sep(X1 = Xcand, d = f1$theta, g = g)
    cov.Xnewx1 <- covar.sep(X1 = f1$X, X2 = Xcand, d = f1$theta, g = 0)
  } else if (kernel == "matern1.5") {
    cov.newx1 <- cor.sep(X = Xcand, theta = f1$theta, nu = 1.5)
    cov.Xnewx1 <- cor.sep(X = f1$X, x = Xcand, theta = f1$theta, nu = 1.5)
  } else if (kernel == "matern2.5") {
    cov.newx1 <- cor.sep(X = Xcand, theta = f1$theta, nu = 2.5)
    cov.Xnewx1 <- cor.sep(X = f1$X, x = Xcand, theta = f1$theta, nu = 2.5)
  }
  v.next1 <- drop(cov.newx1 - t(cov.Xnewx1) %*% f1$Ki %*% cov.Xnewx1)
  g.next1 <- -1 / drop(v.next1) * f1$Ki %*% cov.Xnewx1


  fit1$Ki <- rbind(
    cbind(f1$Ki + g.next1 %*% t(g.next1) * v.next1, g.next1),
    cbind(t(g.next1), 1 / drop(v.next1))
  )

  fit1$X <- rbind(f1$X, Xcand)
  attr(fit1$X, "scaled:center") <- x.center1
  attr(fit1$X, "scaled:scale") <- x.scale1


  ALC.out <- rep(0, mc.sample)
  if (parallel) {
    ALC.out <- foreach(i = 1:mc.sample, .combine = c) %dopar% {
      fit.tmp <- fit
      if (constant) {
        fit1$y <- c(f1$y, y1.sample[i])
      } else {
        fit1$y <- c(f1$y, y1.sample[i] - y.center1)
        attr(fit1$y, "scaled:center") <- y.center1
      }
      fit1$tau2hat <- drop(t(fit1$y - fit1$mu.hat) %*% fit1$Ki %*% (fit1$y - fit1$mu.hat) / length(fit1$y))

      fit.tmp$fit1 <- fit1

      return(mean(predict(fit.tmp, Xref)$sig2)) # to minimize the deduced variance. To maximize, -mean
    }
  } else {
    for (i in 1:mc.sample) {
      fit.tmp <- fit
      if (constant) {
        fit1$y <- c(f1$y, y1.sample[i])
      } else {
        fit1$y <- c(f1$y, y1.sample[i] - y.center1)
        attr(fit1$y, "scaled:center") <- y.center1
      }
      fit1$tau2hat <- drop(t(fit1$y - fit1$mu.hat) %*% fit1$Ki %*% (fit1$y - fit1$mu.hat) / length(fit1$y))

      fit.tmp$fit1 <- fit1

      ALC.out[i] <- mean(predict(fit.tmp, Xref)$sig2) # to minimize the deduced variance. To maximize, -mean
    }
  }
  return(mean(ALC.out))
}


#' object to optimize the next point by ALC criterion updating at level 2 with two levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param Xref vector or matrix of reference data.
#' @param fit an object of class RNAmf.
#' @param mc.sample a number of mc samples generated for this approach. Default is 100.
#' @param parallel logical indicating whether to run parallel or not. Default is FALSE.
#' @param ncore the number of core for parallel. Default is 1.
#' @return A mean of the deduced variance at Xref.
#' @importFrom plgp covar.sep
#' @importFrom stats rnorm
#' @noRd
#'

obj.ALC_level_2 <- function(Xcand, Xref, fit, mc.sample, parallel = FALSE, ncore = 1) {
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

  Xcand <- matrix(Xcand, nrow = 1)
  if (kernel == "sqex") {
    y1.sample <- rnorm(mc.sample, mean = pred.GP(f1, Xcand)$mu, sd = sqrt(pred.GP(f1, Xcand)$sig2))
  } else if (kernel == "matern1.5") {
    y1.sample <- rnorm(mc.sample, mean = pred.matGP(f1, Xcand)$mu, sd = sqrt(pred.matGP(f1, Xcand)$sig2))
  } else if (kernel == "matern2.5") {
    y1.sample <- rnorm(mc.sample, mean = pred.matGP(f1, Xcand)$mu, sd = sqrt(pred.matGP(f1, Xcand)$sig2))
  }
  Xcand <- scale.inputs(Xcand, x.center1, x.scale1)

  ### Choose level 1 ###
  ### update Ki1
  if (kernel == "sqex") {
    cov.newx1 <- covar.sep(X1 = Xcand, d = f1$theta, g = g)
    cov.Xnewx1 <- covar.sep(X1 = f1$X, X2 = Xcand, d = f1$theta, g = 0)
  } else if (kernel == "matern1.5") {
    cov.newx1 <- cor.sep(X = Xcand, theta = f1$theta, nu = 1.5)
    cov.Xnewx1 <- cor.sep(X = f1$X, x = Xcand, theta = f1$theta, nu = 1.5)
  } else if (kernel == "matern2.5") {
    cov.newx1 <- cor.sep(X = Xcand, theta = f1$theta, nu = 2.5)
    cov.Xnewx1 <- cor.sep(X = f1$X, x = Xcand, theta = f1$theta, nu = 2.5)
  }
  v.next1 <- drop(cov.newx1 - t(cov.Xnewx1) %*% f1$Ki %*% cov.Xnewx1)
  g.next1 <- -1 / drop(v.next1) * f1$Ki %*% cov.Xnewx1


  fit1$Ki <- rbind(
    cbind(f1$Ki + g.next1 %*% t(g.next1) * v.next1, g.next1),
    cbind(t(g.next1), 1 / drop(v.next1))
  )

  fit1$X <- rbind(f1$X, Xcand)
  attr(fit1$X, "scaled:center") <- x.center1
  attr(fit1$X, "scaled:scale") <- x.scale1

  if (parallel) {
    ALC.out <- foreach(i = 1:mc.sample, .combine = c) %dopar% {
      fit.tmp <- fit
      if (constant) {
        fit1$y <- c(f1$y, y1.sample[i])
      } else {
        fit1$y <- c(f1$y, y1.sample[i] - y.center1)
        attr(fit1$y, "scaled:center") <- y.center1
      }

      fit1$tau2hat <- drop(t(fit1$y - fit1$mu.hat) %*% fit1$Ki %*% (fit1$y - fit1$mu.hat) / length(fit1$y))

      fit.tmp$fit1 <- fit1

      ### Choose level 2 ###
      if (kernel == "sqex") {
        pred2 <- pred.GP(fit2, cbind(Xcand, y1.sample[i]))
        y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      } else if (kernel == "matern1.5") {
        pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
        y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      } else if (kernel == "matern2.5") {
        pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
        y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      }

      ### update Ki2
      newx2 <- scale.inputs(cbind(Xcand, y1.sample[i]), x.center2, x.scale2)

      if (kernel == "sqex") {
        cov.newx2 <- covar.sep(X1 = newx2, d = f2$theta, g = g)
        cov.Xnewx2 <- covar.sep(X1 = f2$X, X2 = newx2, d = f2$theta, g = 0)
      } else if (kernel == "matern1.5") {
        cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 1.5)
        cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 1.5)
      } else if (kernel == "matern2.5") {
        cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 2.5)
        cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 2.5)
      }
      v.next2 <- drop(cov.newx2 - t(cov.Xnewx2) %*% f2$Ki %*% cov.Xnewx2)
      g.next2 <- -1 / drop(v.next2) * f2$Ki %*% cov.Xnewx2

      fit2$Ki <- rbind(
        cbind(f2$Ki + g.next2 %*% t(g.next2) * v.next2, g.next2),
        cbind(t(g.next2), 1 / drop(v.next2))
      )

      fit2$X <- rbind(f2$X, newx2)
      attr(fit2$X, "scaled:center") <- x.center2
      attr(fit2$X, "scaled:scale") <- x.scale2

      if (constant) {
        fit2$y <- c(f2$y, y2.sample)
      } else {
        fit2$y <- c(f2$y, y2.sample - y.center2)
        attr(fit2$y, "scaled:center") <- y.center2
      }

      fit2$tau2hat <- drop(t(fit2$y - fit2$mu.hat) %*% fit2$Ki %*% (fit2$y - fit2$mu.hat) / length(fit2$y))

      fit.tmp$fit2 <- fit2

      return(mean(predict(fit.tmp, Xref)$sig2)) # to minimize the deduced variance. To maximize, -mean
    }
  } else {
    ALC.out <- rep(0, mc.sample)
    for (i in 1:mc.sample) {
      fit.tmp <- fit
      if (constant) {
        fit1$y <- c(f1$y, y1.sample[i])
      } else {
        fit1$y <- c(f1$y, y1.sample[i] - y.center1)
        attr(fit1$y, "scaled:center") <- y.center1
      }

      fit1$tau2hat <- drop(t(fit1$y - fit1$mu.hat) %*% fit1$Ki %*% (fit1$y - fit1$mu.hat) / length(fit1$y))

      fit.tmp$fit1 <- fit1

      ### Choose level 2 ###
      if (kernel == "sqex") {
        pred2 <- pred.GP(fit2, cbind(Xcand, y1.sample[i]))
        y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      } else if (kernel == "matern1.5") {
        pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
        y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      } else if (kernel == "matern2.5") {
        pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
        y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      }

      ### update Ki2
      newx2 <- scale.inputs(cbind(Xcand, y1.sample[i]), x.center2, x.scale2)

      if (kernel == "sqex") {
        cov.newx2 <- covar.sep(X1 = newx2, d = f2$theta, g = g)
        cov.Xnewx2 <- covar.sep(X1 = f2$X, X2 = newx2, d = f2$theta, g = 0)
      } else if (kernel == "matern1.5") {
        cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 1.5)
        cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 1.5)
      } else if (kernel == "matern2.5") {
        cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 2.5)
        cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 2.5)
      }
      v.next2 <- drop(cov.newx2 - t(cov.Xnewx2) %*% f2$Ki %*% cov.Xnewx2)
      g.next2 <- -1 / drop(v.next2) * f2$Ki %*% cov.Xnewx2

      fit2$Ki <- rbind(
        cbind(f2$Ki + g.next2 %*% t(g.next2) * v.next2, g.next2),
        cbind(t(g.next2), 1 / drop(v.next2))
      )

      fit2$X <- rbind(f2$X, newx2)
      attr(fit2$X, "scaled:center") <- x.center2
      attr(fit2$X, "scaled:scale") <- x.scale2

      if (constant) {
        fit2$y <- c(f2$y, y2.sample)
      } else {
        fit2$y <- c(f2$y, y2.sample - y.center2)
        attr(fit2$y, "scaled:center") <- y.center2
      }

      fit2$tau2hat <- drop(t(fit2$y - fit2$mu.hat) %*% fit2$Ki %*% (fit2$y - fit2$mu.hat) / length(fit2$y))

      fit.tmp$fit2 <- fit2

      ALC.out[i] <- mean(predict(fit.tmp, Xref)$sig2) # to minimize the deduced variance. To maximize, -mean
    }
  }

  return(mean(ALC.out))
}


#' object to optimize the next point by ALC criterion updating at level 3 with three levels of fidelity
#'
#' @param Xcand candidate data point to be optimized.
#' @param Xref vector or matrix of reference data.
#' @param fit an object of class RNAmf.
#' @param mc.sample a number of mc samples generated for this approach. Default is 100.
#' @param parallel logical indicating whether to run parallel or not. Default is FALSE.
#' @param ncore the number of core for parallel. Default is 1.
#' @return A mean of the deduced variance at Xref.
#' @importFrom plgp covar.sep
#' @importFrom stats rnorm
#' @noRd
#'

obj.ALC_level_3 <- function(Xcand, Xref, fit, mc.sample, parallel = FALSE, ncore = 1) {
  fit_two_level <- fit$fit.RNAmf_two_level
  fit1 <- f1 <- fit_two_level$fit1
  fit2 <- f2 <- fit_two_level$fit2
  fit3 <- f3 <- fit$fit3
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

  Xcand <- matrix(Xcand, nrow = 1)
  if (kernel == "sqex") {
    y1.sample <- rnorm(mc.sample, mean = pred.GP(f1, Xcand)$mu, sd = sqrt(pred.GP(f1, Xcand)$sig2))
  } else if (kernel == "matern1.5") {
    y1.sample <- rnorm(mc.sample, mean = pred.matGP(f1, Xcand)$mu, sd = sqrt(pred.matGP(f1, Xcand)$sig2))
  } else if (kernel == "matern2.5") {
    y1.sample <- rnorm(mc.sample, mean = pred.matGP(f1, Xcand)$mu, sd = sqrt(pred.matGP(f1, Xcand)$sig2))
  }
  Xcand <- scale.inputs(Xcand, x.center1, x.scale1)

  ### Choose level 1 ###
  ### update Ki1
  if (kernel == "sqex") {
    cov.newx1 <- covar.sep(X1 = Xcand, d = f1$theta, g = g)
    cov.Xnewx1 <- covar.sep(X1 = f1$X, X2 = Xcand, d = f1$theta, g = 0)
  } else if (kernel == "matern1.5") {
    cov.newx1 <- cor.sep(X = Xcand, theta = f1$theta, nu = 1.5)
    cov.Xnewx1 <- cor.sep(X = f1$X, x = Xcand, theta = f1$theta, nu = 1.5)
  } else if (kernel == "matern2.5") {
    cov.newx1 <- cor.sep(X = Xcand, theta = f1$theta, nu = 2.5)
    cov.Xnewx1 <- cor.sep(X = f1$X, x = Xcand, theta = f1$theta, nu = 2.5)
  }
  v.next1 <- drop(cov.newx1 - t(cov.Xnewx1) %*% f1$Ki %*% cov.Xnewx1)
  g.next1 <- -1 / drop(v.next1) * f1$Ki %*% cov.Xnewx1


  fit1$Ki <- rbind(
    cbind(f1$Ki + g.next1 %*% t(g.next1) * v.next1, g.next1),
    cbind(t(g.next1), 1 / drop(v.next1))
  )

  fit1$X <- rbind(f1$X, Xcand)
  attr(fit1$X, "scaled:center") <- x.center1
  attr(fit1$X, "scaled:scale") <- x.scale1

  if (parallel) {
    ALC.out <- foreach(i = 1:mc.sample, .combine = c) %dopar% {
      fit.tmp <- fit
      if (constant) {
        fit1$y <- c(f1$y, y1.sample[i])
      } else {
        fit1$y <- c(f1$y, y1.sample[i] - y.center1)
        attr(fit1$y, "scaled:center") <- y.center1
      }

      fit1$tau2hat <- drop(t(fit1$y - fit1$mu.hat) %*% fit1$Ki %*% (fit1$y - fit1$mu.hat) / length(fit1$y))

      fit.tmp$fit1 <- fit1

      ### Choose level 2 ###
      if (kernel == "sqex") {
        pred2 <- pred.GP(fit2, cbind(Xcand, y1.sample[i]))
        y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      } else if (kernel == "matern1.5") {
        pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
        y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      } else if (kernel == "matern2.5") {
        pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
        y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      }

      ### update Ki2
      newx2 <- scale.inputs(cbind(Xcand, y1.sample[i]), x.center2, x.scale2)

      if (kernel == "sqex") {
        cov.newx2 <- covar.sep(X1 = newx2, d = f2$theta, g = g)
        cov.Xnewx2 <- covar.sep(X1 = f2$X, X2 = newx2, d = f2$theta, g = 0)
      } else if (kernel == "matern1.5") {
        cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 1.5)
        cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 1.5)
      } else if (kernel == "matern2.5") {
        cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 2.5)
        cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 2.5)
      }
      v.next2 <- drop(cov.newx2 - t(cov.Xnewx2) %*% f2$Ki %*% cov.Xnewx2)
      g.next2 <- -1 / drop(v.next2) * f2$Ki %*% cov.Xnewx2

      fit2$Ki <- rbind(
        cbind(f2$Ki + g.next2 %*% t(g.next2) * v.next2, g.next2),
        cbind(t(g.next2), 1 / drop(v.next2))
      )

      fit2$X <- rbind(f2$X, newx2)
      attr(fit2$X, "scaled:center") <- x.center2
      attr(fit2$X, "scaled:scale") <- x.scale2

      if (constant) {
        fit2$y <- c(f2$y, y2.sample)
      } else {
        fit2$y <- c(f2$y, y2.sample - y.center2)
        attr(fit2$y, "scaled:center") <- y.center2
      }

      fit2$tau2hat <- drop(t(fit2$y - fit2$mu.hat) %*% fit2$Ki %*% (fit2$y - fit2$mu.hat) / length(fit2$y))

      fit.tmp$fit2 <- fit2

      ### Choose level 3 ###
      if (kernel == "sqex") {
        pred3 <- pred.GP(fit2, cbind(Xcand, y2.sample[i]))
        y3.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      } else if (kernel == "matern1.5") {
        pred3 <- pred.matGP(fit2, cbind(Xcand, y2.sample[i]))
        y3.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      } else if (kernel == "matern2.5") {
        pred3 <- pred.matGP(fit2, cbind(Xcand, y2.sample[i]))
        y3.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      }

      ### update Ki2
      newx3 <- scale.inputs(cbind(Xcand, y2.sample[i]), x.center3, x.scale3)

      if (kernel == "sqex") {
        cov.newx3 <- covar.sep(X1 = newx3, d = f3$theta, g = g)
        cov.Xnewx3 <- covar.sep(X1 = f3$X, X2 = newx3, d = f3$theta, g = 0)
      } else if (kernel == "matern1.5") {
        cov.newx3 <- cor.sep(X = newx3, theta = f3$theta, nu = 1.5)
        cov.Xnewx3 <- cor.sep(X = f3$X, x = newx3, theta = f3$theta, nu = 1.5)
      } else if (kernel == "matern2.5") {
        cov.newx3 <- cor.sep(X = newx3, theta = f3$theta, nu = 2.5)
        cov.Xnewx3 <- cor.sep(X = f3$X, x = newx3, theta = f3$theta, nu = 2.5)
      }
      v.next3 <- drop(cov.newx3 - t(cov.Xnewx3) %*% f3$Ki %*% cov.Xnewx3)
      g.next3 <- -1 / drop(v.next3) * f3$Ki %*% cov.Xnewx3

      fit3$Ki <- rbind(
        cbind(f3$Ki + g.next3 %*% t(g.next3) * v.next3, g.next3),
        cbind(t(g.next3), 1 / drop(v.next3))
      )

      fit3$X <- rbind(f3$X, newx3)
      attr(fit3$X, "scaled:center") <- x.center3
      attr(fit3$X, "scaled:scale") <- x.scale3

      if (constant) {
        fit3$y <- c(f3$y, y3.sample)
      } else {
        fit3$y <- c(f3$y, y3.sample - y.center3)
        attr(fit3$y, "scaled:center") <- y.center3
      }

      fit3$tau2hat <- drop(t(fit3$y - fit3$mu.hat) %*% fit3$Ki %*% (fit3$y - fit3$mu.hat) / length(fit3$y))

      fit.tmp$fit3 <- fit3

      return(mean(predict(fit.tmp, Xref)$sig2)) # to minimize the deduced variance. To maximize, -mean
    }
  } else {
    ALC.out <- rep(0, mc.sample)
    for (i in 1:mc.sample) {
      fit.tmp <- fit
      if (constant) {
        fit1$y <- c(f1$y, y1.sample[i])
      } else {
        fit1$y <- c(f1$y, y1.sample[i] - y.center1)
        attr(fit1$y, "scaled:center") <- y.center1
      }

      fit1$tau2hat <- drop(t(fit1$y - fit1$mu.hat) %*% fit1$Ki %*% (fit1$y - fit1$mu.hat) / length(fit1$y))

      fit.tmp$fit1 <- fit1

      ### Choose level 2 ###
      if (kernel == "sqex") {
        pred2 <- pred.GP(fit2, cbind(Xcand, y1.sample[i]))
        y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      } else if (kernel == "matern1.5") {
        pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
        y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      } else if (kernel == "matern2.5") {
        pred2 <- pred.matGP(fit2, cbind(Xcand, y1.sample[i]))
        y2.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      }

      ### update Ki2
      newx2 <- scale.inputs(cbind(Xcand, y1.sample[i]), x.center2, x.scale2)

      if (kernel == "sqex") {
        cov.newx2 <- covar.sep(X1 = newx2, d = f2$theta, g = g)
        cov.Xnewx2 <- covar.sep(X1 = f2$X, X2 = newx2, d = f2$theta, g = 0)
      } else if (kernel == "matern1.5") {
        cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 1.5)
        cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 1.5)
      } else if (kernel == "matern2.5") {
        cov.newx2 <- cor.sep(X = newx2, theta = f2$theta, nu = 2.5)
        cov.Xnewx2 <- cor.sep(X = f2$X, x = newx2, theta = f2$theta, nu = 2.5)
      }
      v.next2 <- drop(cov.newx2 - t(cov.Xnewx2) %*% f2$Ki %*% cov.Xnewx2)
      g.next2 <- -1 / drop(v.next2) * f2$Ki %*% cov.Xnewx2

      fit2$Ki <- rbind(
        cbind(f2$Ki + g.next2 %*% t(g.next2) * v.next2, g.next2),
        cbind(t(g.next2), 1 / drop(v.next2))
      )

      fit2$X <- rbind(f2$X, newx2)
      attr(fit2$X, "scaled:center") <- x.center2
      attr(fit2$X, "scaled:scale") <- x.scale2

      if (constant) {
        fit2$y <- c(f2$y, y2.sample)
      } else {
        fit2$y <- c(f2$y, y2.sample - y.center2)
        attr(fit2$y, "scaled:center") <- y.center2
      }

      fit2$tau2hat <- drop(t(fit2$y - fit2$mu.hat) %*% fit2$Ki %*% (fit2$y - fit2$mu.hat) / length(fit2$y))

      fit.tmp$fit2 <- fit2

      ### Choose level 3 ###
      if (kernel == "sqex") {
        pred3 <- pred.GP(fit2, cbind(Xcand, y2.sample[i]))
        y3.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      } else if (kernel == "matern1.5") {
        pred3 <- pred.matGP(fit2, cbind(Xcand, y2.sample[i]))
        y3.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      } else if (kernel == "matern2.5") {
        pred3 <- pred.matGP(fit2, cbind(Xcand, y2.sample[i]))
        y3.sample <- rnorm(1, pred2$mu, sqrt(pred2$sig2))
      }

      ### update Ki2
      newx3 <- scale.inputs(cbind(Xcand, y2.sample[i]), x.center3, x.scale3)

      if (kernel == "sqex") {
        cov.newx3 <- covar.sep(X1 = newx3, d = f3$theta, g = g)
        cov.Xnewx3 <- covar.sep(X1 = f3$X, X2 = newx3, d = f3$theta, g = 0)
      } else if (kernel == "matern1.5") {
        cov.newx3 <- cor.sep(X = newx3, theta = f3$theta, nu = 1.5)
        cov.Xnewx3 <- cor.sep(X = f3$X, x = newx3, theta = f3$theta, nu = 1.5)
      } else if (kernel == "matern2.5") {
        cov.newx3 <- cor.sep(X = newx3, theta = f3$theta, nu = 2.5)
        cov.Xnewx3 <- cor.sep(X = f3$X, x = newx3, theta = f3$theta, nu = 2.5)
      }
      v.next3 <- drop(cov.newx3 - t(cov.Xnewx3) %*% f3$Ki %*% cov.Xnewx3)
      g.next3 <- -1 / drop(v.next3) * f3$Ki %*% cov.Xnewx3

      fit3$Ki <- rbind(
        cbind(f3$Ki + g.next3 %*% t(g.next3) * v.next3, g.next3),
        cbind(t(g.next3), 1 / drop(v.next3))
      )

      fit3$X <- rbind(f3$X, newx3)
      attr(fit3$X, "scaled:center") <- x.center3
      attr(fit3$X, "scaled:scale") <- x.scale3

      if (constant) {
        fit3$y <- c(f3$y, y3.sample)
      } else {
        fit3$y <- c(f3$y, y3.sample - y.center3)
        attr(fit3$y, "scaled:center") <- y.center3
      }

      fit3$tau2hat <- drop(t(fit3$y - fit3$mu.hat) %*% fit3$Ki %*% (fit3$y - fit3$mu.hat) / length(fit3$y))

      fit.tmp$fit3 <- fit3

      ALC.out[i] <- mean(predict(fit.tmp, Xref)$sig2) # to minimize the deduced variance. To maximize, -mean
    }
  }

  return(mean(ALC.out))
}


#' @title find the next point by ALC criterion
#'
#' @description The function acquires the new point by the Active learning Cohn (ALC) criterion.
#' It calculates the ALC criterion
#' \eqn{\frac{\Delta \sigma_L^{2}(l,\bm{x})}{\sum^l_{j=1}C_j} =
#' \frac{\int_{\Omega} \sigma_L^{*2}(\bm{\xi})-\tilde{\sigma}_L^{*2}(\bm{\xi};l,\bm{x}){\rm{d}}\bm{\xi}}{\sum^l_{j=1}C_j}},
#' where \eqn{f_L} is the highest-fidelity simulation code,
#' \eqn{\sigma_L^{*2}(\bm{\xi})} is the posterior variance of \eqn{f_L(\bm{\xi})},
#' \eqn{C_j} is the simulation cost at fidelity level \eqn{j},
#' and \eqn{\tilde{\sigma}_L^{*2}(\bm{\xi};l,\bm{x})} is the posterior variance
#' based on the augmented design combining the current design and a new input location \eqn{\bm{x}}
#' at each fidelity level lower than or equal to \eqn{l}.
#' The integration is approximated by MC integration using uniform reference samples.
#'
#' @details \code{Xref} plays a role of \eqn{\bm{\xi}} to approximate the integration.
#' To impute the posterior variance based on the augmented design \eqn{\tilde{\sigma}_L^{*2}(\bm{\xi};l,\bm{x})},
#' MC approximation is used.
#' Due to the nested assumption, imputing \eqn{y^{[s]}_{n_s+1}} for each \eqn{1\leq s\leq l} by drawing samples
#' from the posterior distribution of \eqn{f_s(\bm{x}^{[s]}_{n_s+1})}
#' based on the current design allows to compute \eqn{\tilde{\sigma}_L^{*2}(\bm{\xi};l,\bm{x})}.
#' Inverse of covariance matrix is computed by the Sherman-Morrison formula.
#' For details, see Heo and Sung (2023+, <arXiv:2309.11772>).
#'
#' To search for the next acquisition \eqn{\bm{x^*}} by maximizing AL criterion,
#' the gradient-based optimization can be used by \code{optim=TRUE}.
#' Firstly, \eqn{\tilde{\sigma}_L^{*2}(\bm{\xi};l,\bm{x})} is computed on the
#' \eqn{5 \times d} number of points.
#' After that, the point minimizing \eqn{\tilde{\sigma}_L^{*2}(\bm{\xi};l,\bm{x})}
#' serves as a starting point of optimization by \code{L-BFGS-B} method.
#' Otherwise, when \code{optim=FALSE}, AL criterion is optimized only on \code{Xcand}.
#'
#' The point is selected by maximizing the ALC criterion:
#' \eqn{\text{argmax}_{l\in\{1,\ldots,L\}; \bm{x} \in \Omega}
#' \frac{\Delta \sigma_L^{2}(l,\bm{x})}{\sum^l_{j=1}C_j}}.
#'
#'
#' @param Xref vector or matrix of reference location to approximate the integral of ALC. If \code{Xref=NULL}, \eqn{100 \times d} points are generated by Latin hypercube design. Default is \code{NULL}.
#' @param Xcand vector or matrix of candidate set which could be added into the current design only when \code{optim=FALSE}. \code{Xcand} is the set of the points where ALC criterion is evaluated. If \code{Xcand=NULL}, \code{Xref} is used. Default is \code{NULL}. See details.
#' @param fit object of class \code{RNAmf}.
#' @param mc.sample a number of mc samples generated for the imputation through MC approximation. Default is \code{100}.
#' @param cost vector of the costs for each level of fidelity. If \code{cost=NULL}, total costs at all levels would be 1. \code{cost} is encouraged to have a ascending order of positive value. Default is \code{NULL}.
#' @param optim logical indicating whether to optimize AL criterion by \code{optim}'s gradient-based \code{L-BFGS-B} method. If \code{optim=TRUE}, \eqn{5 \times d} starting points are generated by Latin hypercube design for optimization. If \code{optim=FALSE}, AL criterion is optimized on the \code{Xcand}. Default is \code{TRUE}.
#' @param parallel logical indicating whether to compute the AL criterion in parallel or not. If \code{parallel=TRUE}, parallel computation is utilized. Default is \code{FALSE}.
#' @param ncore a number of core for parallel. It is only used if \code{parallel=TRUE}. Default is 1.
#' @return
#' \itemize{
#'   \item \code{ALC}: list of ALC criterion integrated on \code{Xref} when each data point on \code{Xcand} is added at each level \eqn{l} if \code{optim=FALSE}. If \code{optim=TRUE}, \code{ALC} returns \code{NULL}.
#'   \item \code{cost}: a copy of \code{cost}.
#'   \item \code{Xcand}: a copy of \code{Xcand}.
#'   \item \code{chosen}: list of chosen level and point.
#'   \item \code{time}: a scalar of the time for the computation.
#' }
#' @importFrom plgp covar.sep
#' @importFrom stats rnorm
#' @importFrom lhs randomLHS
#' @importFrom lhs maximinLHS
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom doParallel stopImplicitCluster
#' @usage ALC_RNAmf(Xref = NULL, Xcand = NULL, fit, mc.sample = 100,
#' cost = NULL, optim = TRUE, parallel = FALSE, ncore = 1)
#' @export
#' @examples
#' \dontrun{
#' library(lhs)
#' library(doParallel)
#' library(foreach)
#'
#' ### simulation costs ###
#' cost <- c(1, 3)
#'
#' ### 1-d Perdikaris function in Perdikaris, et al. (2017) ###
#' # low-fidelity function
#' f1 <- function(x) {
#'   sin(8 * pi * x)
#' }
#'
#' # high-fidelity function
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
#' ### active learning with optim=TRUE ###
#' alc.RNAmf.optim <- ALC_RNAmf(
#'   Xref = x, Xcand = x, fit.RNAmf, cost = cost,
#'   optim = TRUE, parallel = TRUE, ncore = 10
#' )
#' alc.RNAmf.optim$time # computation time of optim=TRUE
#'
#' ### active learning with optim=FALSE ###
#' alc.RNAmf <- ALC_RNAmf(
#'   Xref = x, Xcand = x, fit.RNAmf, cost = cost,
#'   optim = FALSE, parallel = TRUE, ncore = 10
#' )
#' alc.RNAmf$time # computation time of optim=FALSE
#'
#' ### visualize ALC ###
#' par(mfrow = c(1, 2))
#' plot(x, alc.RNAmf$ALC$ALC1,
#'   type = "l", lty = 2,
#'   xlab = "x", ylab = "ALC criterion augmented at the low-fidelity level",
#'   ylim = c(min(c(alc.RNAmf$ALC$ALC1, alc.RNAmf$ALC$ALC2)),
#'            max(c(alc.RNAmf$ALC$ALC1, alc.RNAmf$ALC$ALC2)))
#' )
#' plot(x, alc.RNAmf$ALC$ALC2,
#'   type = "l", lty = 2,
#'   xlab = "x", ylab = "ALC criterion augmented at the high-fidelity level",
#'   ylim = c(min(c(alc.RNAmf$ALC$1, alc.RNAmf$ALC$2)),
#'            max(c(alc.RNAmf$ALC$1, alc.RNAmf$ALC$2)))
#' )
#' points(alc.RNAmf$chosen$Xnext,
#'   alc.RNAmf$ALC$2[which(x == drop(alc.RNAmf$chosen$Xnext))],
#'   pch = 16, cex = 1, col = "red"
#' )}
#'
ALC_RNAmf <- function(Xref = NULL, Xcand = NULL, fit, mc.sample = 100, cost = NULL, optim = TRUE, parallel = FALSE, ncore = 1) {
  t1 <- proc.time()[3]
  ### check the object ###
  if (!inherits(fit, "RNAmf")) {
    stop("The object is not of class \"RNAmf\" \n")
  }
  if (length(cost) != fit$level) stop("The length of cost should be the level of object")

  ### ALC ###
  if (fit$level == 2) { # level 2
    if (!is.null(cost) & cost[1] >= cost[2]) {
      warning("If the cost for high-fidelity is cheaper, acquire the high-fidelity")
    } else if (is.null(cost)) {
      cost <- c(1, 0)
    }
    if (is.null(Xref)) Xref <- randomLHS(100 * dim(fit$fit1$X)[2], dim(fit$fit1$X)[2])
    if (parallel) registerDoParallel(ncore)

    Icurrent <- mean(predict(fit, Xref)$sig2)

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
    if (optim){ # optim = TRUE
      Xcand <- randomLHS(5*ncol(fit1$X), ncol(fit1$X))
    }else{ # optim = FALSE
      if (is.null(Xcand)){
        Xcand <- Xref
      }
    }
    # if (is.null(Xcand)) Xcand <- maximinLHS(n.start, ncol(fit1$X))
    # Xcand <- matrix(Xcand, ncol = ncol(fit1$X))
    # if (ncol(Xcand) != dim(fit$fit1$X)[2]) stop("The dimension of candidate set should be equal to the dimension of the design")

    ### Calculate the deduced variance ###
    cat("running starting points: \n")
    time.start <- proc.time()[3]
    if (parallel) {
      pseudointvar <- foreach(i = 1:nrow(Xcand), .combine = cbind) %dopar% {
        newx <- matrix(Xcand[i, ], nrow = 1)

        return(c(
          obj.ALC_level_1(newx, Xref, fit, mc.sample),
          obj.ALC_level_2(newx, Xref, fit, mc.sample)
        ))
      }
      intvar1 <- pseudointvar[1, ]
      intvar2 <- pseudointvar[2, ]
    } else {
      intvar1 <- c(rep(0, nrow(Xcand))) # IMSPE candidates
      intvar2 <- c(rep(0, nrow(Xcand))) # IMSPE candidates
      for (i in 1:nrow(Xcand)) {
        print(paste(i, nrow(Xcand), sep = "/"))
        newx <- matrix(Xcand[i, ], nrow = 1)

        intvar1[i] <- obj.ALC_level_1(newx, Xref, fit, mc.sample)
        intvar2[i] <- obj.ALC_level_2(newx, Xref, fit, mc.sample)
      }
    }
    print(proc.time()[3] - time.start)

    ### Find the next point ###
    if (optim) {
      cat("running optim for level 1: \n")
      time.start <- proc.time()[3]
      X.start <- matrix(Xcand[which.min(intvar1), ], nrow = 1)
      optim.out <- optim(X.start, obj.ALC_level_1, method = "L-BFGS-B", lower = 0, upper = 1, fit = fit, Xref = Xref, mc.sample = mc.sample, parallel = parallel, ncore = ncore)
      Xnext.1 <- optim.out$par
      ALC.1 <- optim.out$value
      print(proc.time()[3] - time.start)

      cat("running optim for level 2: \n")
      time.start <- proc.time()[3]
      X.start <- matrix(Xcand[which.min(intvar2), ], nrow = 1)
      optim.out <- optim(X.start, obj.ALC_level_2, method = "L-BFGS-B", lower = 0, upper = 1, fit = fit, Xref = Xref, mc.sample = mc.sample, parallel = parallel, ncore = ncore)
      Xnext.2 <- optim.out$par
      ALC.2 <- optim.out$value
      print(proc.time()[3] - time.start)

      ALCvalue <- c(Icurrent - ALC.1, Icurrent - ALC.2) / c(cost[1], cost[1] + cost[2])
      if (ALCvalue[2] > ALCvalue[1]) {
        level <- 2
        Xnext <- Xnext.2
      } else {
        level <- 1
        Xnext <- Xnext.1
      }
    } else {
      ALCvalue <- c(Icurrent - which.min(intvar1), Icurrent - which.min(intvar2)) / c(cost[1], cost[1] + cost[2])
      if (ALCvalue[2] > ALCvalue[1]) {
        level <- 2
        Xnext <- matrix(Xcand[which.min(intvar2), ], nrow = 1)
      } else {
        level <- 1
        Xnext <- matrix(Xcand[which.min(intvar1), ], nrow = 1)
      }
    }

    chosen <- list(
      "level" = level, # next level
      "Xnext" = Xnext
    ) # next point


    # ### Update the model ###
    # newx <- matrix(chosen$Xnext, nrow=1)
    # level <- chosen$level
    #
    # X1 <- scale.inputs(fit1$X, x.center1, x.scale1, back=TRUE)
    # X2 <- matrix(scale.inputs(fit2$X, x.center2, x.scale2, back=TRUE)[,-ncol(fit2$X)], ncol=ncol(fit2$X)-1)
    #
    # if(constant){
    #   y1 <- fit1$y
    #   y2 <- fit2$y
    # }else{
    #   y1 <- fit1$y+attr(fit1$y,"scaled:center")
    #   y2 <- fit2$y+attr(fit2$y,"scaled:center")
    # }
    #
    #
    # ### Choose level 1 ###
    # if(level == 1){
    #   y1.select <- funcs[[1]](newx)
    #
    #   X1 <- rbind(X1, newx)
    #   y1 <- c(y1, y1.select)
    # }
    #
    # ### Choose level 2 ###
    # if(level == 2){
    #   if(checknested(scale.inputs(fit1$X, x.center1, x.scale1, back=TRUE), chosen$Xnext)==TRUE &
    #      checknested(scale.inputs(fit2$X, x.center2, x.scale2, back=TRUE), chosen$Xnext)==FALSE){
    #     y2.select <- funcs[[2]](newx)
    #
    #     X2 <- rbind(X2, newx)
    #     y2 <- c(y2, y2.select)
    #   }else{
    #     y1.select <- funcs[[1]](newx)
    #     y2.select <- funcs[[2]](newx)
    #
    #     X1 <- rbind(X1, newx)
    #     y1 <- c(y1, y1.select)
    #     X2 <- rbind(X2, newx)
    #     y2 <- c(y2, y2.select)
    #   }
    # }
    #
    # fit <- RNAmf_two_level(X1, y1, X2, y2, kernel=kernel, constant=constant)
    ALC <- list(ALC1 = intvar1 / cost[1], ALC2 = intvar2 / (cost[1] + cost[2]))
  } else if (fit$level == 3) { # level 3

    if (!is.null(cost) & (cost[1] >= cost[2] | cost[2] >= cost[3])) {
      warning("If the cost for high-fidelity is cheaper, acquire the high-fidelity")
    } else if (is.null(cost)) {
      cost <- c(1, 0, 0)
    }
    if (is.null(Xref)) Xref <- randomLHS(100 * dim(fit$fit.RNAmf_two_level$fit1$X)[2], dim(fit$fit.RNAmf_two_level$fit1$X)[2])
    if (parallel) registerDoParallel(ncore)

    Icurrent <- mean(predict(fit, Xref)$sig2)

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
    if (optim){ # optim = TRUE
      Xcand <- randomLHS(5*ncol(fit1$X), ncol(fit1$X))
    }else{ # optim = FALSE
      if (is.null(Xcand)){
        Xcand <- Xref
      }
    }
    # if (is.null(Xcand)) Xcand <- maximinLHS(n.start, ncol(fit1$X))
    # Xcand <- matrix(Xcand, ncol = ncol(fit1$X))
    # if (ncol(Xcand) != dim(fit$fit1$X)[2]) stop("The dimension of candidate set should be equal to the dimension of the design")

    ### Calculate the deduced variance ###
    cat("running starting points: \n")
    time.start <- proc.time()[3]
    if (parallel) {
      pseudointvar <- foreach(i = 1:nrow(Xcand), .combine = cbind) %dopar% {
        newx <- matrix(Xcand[i, ], nrow = 1)

        return(c(
          obj.ALC_level_1(newx, Xref, fit_two_level, mc.sample),
          obj.ALC_level_2(newx, Xref, fit_two_level, mc.sample),
          obj.ALC_level_3(newx, Xref, fit, mc.sample)
        ))
      }
      intvar1 <- pseudointvar[1, ]
      intvar2 <- pseudointvar[2, ]
      intvar2 <- pseudointvar[3, ]
    } else {
      intvar1 <- c(rep(0, nrow(Xcand))) # IMSPE candidates
      intvar2 <- c(rep(0, nrow(Xcand))) # IMSPE candidates
      intvar3 <- c(rep(0, nrow(Xcand))) # IMSPE candidates
      for (i in 1:nrow(Xcand)) {
        print(paste(i, nrow(Xcand), sep = "/"))
        newx <- matrix(Xcand[i, ], nrow = 1)

        intvar1[i] <- obj.ALC_level_1(newx, Xref, fit_two_level, mc.sample)
        intvar2[i] <- obj.ALC_level_2(newx, Xref, fit_two_level, mc.sample)
        intvar3[i] <- obj.ALC_level_3(newx, Xref, fit, mc.sample)
      }
    }
    print(proc.time()[3] - time.start)

    ### Find the next point ###
    if (optim) {
      cat("running optim for level 1: \n")
      time.start <- proc.time()[3]
      X.start <- matrix(Xcand[which.min(intvar1), ], nrow = 1)
      optim.out <- optim(X.start, obj.ALC_level_1, method = "L-BFGS-B", lower = 0, upper = 1, fit = fit_two_level, Xref = Xref, mc.sample = mc.sample, parallel = parallel, ncore = ncore)
      Xnext.1 <- optim.out$par
      ALC.1 <- optim.out$value
      print(proc.time()[3] - time.start)

      cat("running optim for level 2: \n")
      time.start <- proc.time()[3]
      X.start <- matrix(Xcand[which.min(intvar2), ], nrow = 1)
      optim.out <- optim(X.start, obj.ALC_level_2, method = "L-BFGS-B", lower = 0, upper = 1, fit = fit_two_level, Xref = Xref, mc.sample = mc.sample, parallel = parallel, ncore = ncore)
      Xnext.2 <- optim.out$par
      ALC.2 <- optim.out$value
      print(proc.time()[3] - time.start)

      cat("running optim for level 3: \n")
      time.start <- proc.time()[3]
      X.start <- matrix(Xcand[which.min(intvar3), ], nrow = 1)
      optim.out <- optim(X.start, obj.ALC_level_3, method = "L-BFGS-B", lower = 0, upper = 1, fit = fit, Xref = Xref, mc.sample = mc.sample, parallel = parallel, ncore = ncore)
      Xnext.3 <- optim.out$par
      ALC.3 <- optim.out$value
      print(proc.time()[3] - time.start)

      ALCvalue <- c(Icurrent - ALC.1, Icurrent - ALC.2, Icurrent - ALC.3) / c(cost[1], cost[1] + cost[2], cost[1] + cost[2] + cost[3])
      if (ALCvalue[3] > ALCvalue[2]) {
        level <- 3
        Xnext <- Xnext.3
      } else if (ALCvalue[2] > ALCvalue[1]) {
        level <- 2
        Xnext <- Xnext.2
      } else {
        level <- 1
        Xnext <- Xnext.1
      }
    } else {
      ALCvalue <- c(Icurrent - which.min(intvar1), Icurrent - which.min(intvar2), Icurrent - which.min(intvar3)) / c(cost[1], cost[1] + cost[2], cost[1] + cost[2] + cost[3])
      if (ALCvalue[3] > ALCvalue[2]) {
        level <- 3
        Xnext <- matrix(Xcand[which.min(intvar3), ], nrow = 1)
      } else if (ALCvalue[2] > ALCvalue[1]) {
        level <- 2
        Xnext <- matrix(Xcand[which.min(intvar2), ], nrow = 1)
      } else {
        level <- 1
        Xnext <- matrix(Xcand[which.min(intvar1), ], nrow = 1)
      }
    }

    chosen <- list(
      "level" = level, # next level
      "Xnext" = Xnext
    ) # next point


    # ### Update the model ###
    # newx <- matrix(chosen$Xnext, nrow=1)
    # level <- chosen$level
    #
    # X1 <- scale.inputs(fit1$X, x.center1, x.scale1, back=TRUE)
    # X2 <- matrix(scale.inputs(fit2$X, x.center2, x.scale2, back=TRUE)[,-ncol(fit2$X)], ncol=ncol(fit2$X)-1)
    # X3 <- matrix(scale.inputs(fit3$X, x.center3, x.scale3, back=TRUE)[,-ncol(fit3$X)], ncol=ncol(fit3$X)-1)
    #
    # if(constant){
    #   y1 <- fit1$y
    #   y2 <- fit2$y
    #   y3 <- fit3$y
    # }else{
    #   y1 <- fit1$y+attr(fit1$y,"scaled:center")
    #   y2 <- fit2$y+attr(fit2$y,"scaled:center")
    #   y3 <- fit3$y+attr(fit3$y,"scaled:center")
    # }
    #
    #
    # ### Choose level 1 ###
    # if(level == 1){
    #   y1.select <- funcs[[1]](newx)
    #
    #   X1 <- rbind(X1, newx)
    #   y1 <- c(y1, y1.select)
    # }
    #
    # ### Choose level 2 ###
    # if(level == 2){
    #   if(checknested(scale.inputs(fit1$X, x.center1, x.scale1, back=TRUE), chosen$Xnext)==TRUE &
    #      checknested(scale.inputs(fit2$X, x.center2, x.scale2, back=TRUE), chosen$Xnext)==FALSE){
    #     y2.select <- funcs[[2]](newx)
    #
    #     X2 <- rbind(X2, newx)
    #     y2 <- c(y2, y2.select)
    #   }else{
    #     y1.select <- funcs[[1]](newx)
    #     y2.select <- funcs[[2]](newx)
    #
    #     X1 <- rbind(X1, newx)
    #     y1 <- c(y1, y1.select)
    #     X2 <- rbind(X2, newx)
    #     y2 <- c(y2, y2.select)
    #   }
    # }
    #
    # ### Choose level 3 ###
    # if(level == 3){
    #   if(checknested(scale.inputs(fit1$X, x.center1, x.scale1, back=TRUE), chosen$Xnext)==TRUE &
    #      checknested(scale.inputs(fit2$X, x.center2, x.scale2, back=TRUE), chosen$Xnext)==TRUE &
    #      checknested(scale.inputs(fit3$X, x.center3, x.scale3, back=TRUE), chosen$Xnext)==FALSE){
    #     y3.select <- funcs[[3]](newx)
    #
    #     X3 <- rbind(X3, newx)
    #     y3 <- c(y3, y3.select)
    #   }else if(checknested(scale.inputs(fit1$X, x.center1, x.scale1, back=TRUE), chosen$Xnext)==TRUE &
    #            checknested(scale.inputs(fit2$X, x.center2, x.scale2, back=TRUE), chosen$Xnext)==FALSE &
    #            checknested(scale.inputs(fit3$X, x.center3, x.scale3, back=TRUE), chosen$Xnext)==FALSE){
    #     y2.select <- funcs[[2]](newx)
    #     y3.select <- funcs[[3]](newx)
    #
    #     X2 <- rbind(X2, newx)
    #     y2 <- c(y2, y2.select)
    #     X3 <- rbind(X3, newx)
    #     y3 <- c(y3, y3.select)
    #   }else{
    #     y1.select <- funcs[[1]](newx)
    #     y2.select <- funcs[[2]](newx)
    #     y3.select <- funcs[[3]](newx)
    #
    #     X1 <- rbind(X1, newx)
    #     y1 <- c(y1, y1.select)
    #     X2 <- rbind(X2, newx)
    #     y2 <- c(y2, y2.select)
    #     X3 <- rbind(X3, newx)
    #     y3 <- c(y3, y3.select)
    #   }
    # }
    #
    # fit <- RNAmf_three_level(X1, y1, X2, y2, X3, y3, kernel=kernel, constant=constant)
    ALC <- list(ALC1 = intvar1 / cost[1], ALC2 = intvar2 / (cost[1] + cost[2]), ALC3 = intvar3 / (cost[1] + cost[2] + cost[3]))
  } else {
    stop("level is missing")
  }

  if (parallel) stopImplicitCluster()
  if (optim) ALC <- NULL

  return(list(ALC = ALC, cost = cost, Xcand = Xcand, chosen = chosen, time = proc.time()[3] - t1))
}
