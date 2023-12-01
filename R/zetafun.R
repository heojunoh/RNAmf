#' computing zeta component in the closed form posterior variance of matern kernel.
#'
#' @param w1 numerical value of \eqn{f_1(x_k)}.
#' @param w2 numerical value of \eqn{f_1(x_l)}.
#' @param m numerical value of \eqn{\mu^*_1(x)}.
#' @param s numerical value of \eqn{\sigma^{*2}_1(x)}.
#' @param nu numerical value of smoothness hyperparameter. It should be 1.5 or 2.5.
#' @param theta vector of lengthscale hyperparameter.
#'
#' @return calculated value of zeta component.
#'
#' @importFrom plgp distance
#' @importFrom stats pnorm
#' @noRd
#'

zetafun <- function(w1, w2, m, s, nu, theta) {
  if (nu == 1.5) {
    if (w1 < w2) {
      muc <- m - 2 * sqrt(3) * s / theta
      mud <- m + 2 * sqrt(3) * s / theta

      lambda31 <- c(1, muc, muc^2 + s)
      lambda32 <- c(0, 1, muc + w2)
      lambda41 <- c(1, m, m^2 + s)
      lambda42 <- c(0, 1, m + w1)
      lambda43 <- c(0, 1, m + w2)
      lambda51 <- c(1, -mud, mud^2 + s)
      lambda52 <- c(0, 1, -mud - w1)

      e3 <- c(
        1 + (3 * w1 * w2 - sqrt(3) * theta * (w1 + w2)) / theta^2,
        (2 * sqrt(3) * theta - 3 * (w1 + w2)) / theta^2,
        3 / theta^2
      )
      e4 <- c(
        1 + (-3 * w1 * w2 + sqrt(3) * theta * (w2 - w1)) / theta^2,
        3 * (w1 + w2) / theta^2,
        -3 / theta^2
      )
      e5 <- c(
        1 + (3 * w1 * w2 + sqrt(3) * theta * (w1 + w2)) / theta^2,
        (2 * sqrt(3) * theta + 3 * (w1 + w2)) / theta^2,
        3 / theta^2
      )


      return(exp((6 * s + sqrt(3) * theta * (w1 + w2 - 2 * m)) / theta^2) *
        (e3 %*% lambda31 * pnorm((muc - w2) / sqrt(s)) +
          e3 %*% lambda32 * sqrt(s) / sqrt(2 * pi) * exp(-(w2 - muc)^2 / (2 * s))) +
        exp(-(sqrt(3) * (w2 - w1) / theta)) *
          (e4 %*% lambda41 * (pnorm((w2 - m) / sqrt(s)) - pnorm((w1 - m) / sqrt(s))) +
            e4 %*% lambda42 * sqrt(s) / sqrt(2 * pi) * exp(-(w1 - m)^2 / (2 * s)) -
            e4 %*% lambda43 * sqrt(s) / sqrt(2 * pi) * exp(-(w2 - m)^2 / (2 * s))) +
        exp((6 * s - sqrt(3) * theta * (w1 + w2 - 2 * m)) / theta^2) *
          (e5 %*% lambda51 * pnorm((w1 - mud) / sqrt(s)) +
            e5 %*% lambda52 * sqrt(s) / sqrt(2 * pi) * exp(-(w1 - mud)^2 / (2 * s))))
    } else { # b <= a
      muc <- m - 2 * sqrt(3) * s / theta
      mud <- m + 2 * sqrt(3) * s / theta

      lambda31 <- c(1, muc, muc^2 + s)
      lambda32 <- c(0, 1, muc + w1)
      lambda41 <- c(1, m, m^2 + s)
      lambda42 <- c(0, 1, m + w2)
      lambda43 <- c(0, 1, m + w1)
      lambda51 <- c(1, -mud, mud^2 + s)
      lambda52 <- c(0, 1, -mud - w2)

      e3 <- c(
        1 + (3 * w1 * w2 - sqrt(3) * theta * (w1 + w2)) / theta^2,
        (2 * sqrt(3) * theta - 3 * (w1 + w2)) / theta^2,
        3 / theta^2
      )
      e4 <- c(
        1 + (-3 * w1 * w2 + sqrt(3) * theta * (w1 - w2)) / theta^2,
        3 * (w1 + w2) / theta^2,
        -3 / theta^2
      )
      e5 <- c(
        1 + (3 * w1 * w2 + sqrt(3) * theta * (w1 + w2)) / theta^2,
        (2 * sqrt(3) * theta + 3 * (w1 + w2)) / theta^2,
        3 / theta^2
      )


      return(exp((6 * s + sqrt(3) * theta * (w1 + w2 - 2 * m)) / theta^2) *
        (e3 %*% lambda31 * pnorm((muc - w1) / sqrt(s)) +
          e3 %*% lambda32 * sqrt(s) / sqrt(2 * pi) * exp(-(w1 - muc)^2 / (2 * s))) +
        exp(-(sqrt(3) * (w1 - w2) / theta)) *
          (e4 %*% lambda41 * (pnorm((w1 - m) / sqrt(s)) - pnorm((w2 - m) / sqrt(s))) +
            e4 %*% lambda42 * sqrt(s) / sqrt(2 * pi) * exp(-(w2 - m)^2 / (2 * s)) -
            e4 %*% lambda43 * sqrt(s) / sqrt(2 * pi) * exp(-(w1 - m)^2 / (2 * s))) +
        exp((6 * s - sqrt(3) * theta * (w1 + w2 - 2 * m)) / theta^2) *
          (e5 %*% lambda51 * pnorm((w2 - mud) / sqrt(s)) +
            e5 %*% lambda52 * sqrt(s) / sqrt(2 * pi) * exp(-(w2 - mud)^2 / (2 * s))))
    }
  } else if (nu == 2.5) {
    if (w1 < w2) {
      muc <- m - 2 * sqrt(5) * s / theta
      mud <- m + 2 * sqrt(5) * s / theta

      lambda31 <- c(1, muc, muc^2 + s, muc^3 + 3 * muc * s, muc^4 + 6 * muc^2 * s + 3 * s^2)
      lambda32 <- c(0, 1, muc + w2, muc^2 + 2 * s + w2^2 + muc * w2, muc^3 + w2^3 + muc^2 * w2 + muc * w2^2 + 3 * s * w2 + 5 * muc * s)
      lambda41 <- c(1, m, m^2 + s, m^3 + 3 * m * s, m^4 + 6 * m^2 * s + 3 * s^2)
      lambda42 <- c(0, 1, m + w1, m^2 + 2 * s + w1^2 + m * w1, m^3 + w1^3 + w1 * m^2 + m * w1^2 + 3 * s * w1 + 5 * s * m)
      lambda43 <- c(0, 1, m + w2, m^2 + 2 * s + w2^2 + m * w2, m^3 + w2^3 + w2 * m^2 + m * w2^2 + 3 * s * w2 + 5 * s * m)
      lambda51 <- c(1, -mud, mud^2 + s, -mud^3 - 3 * mud * s, mud^4 + 6 * mud^2 * s + 3 * s^2)
      lambda52 <- c(0, 1, -mud - w1, mud^2 + 2 * s + w1^2 + mud * w1, -mud^3 - w1^3 - mud^2 * w1 - mud * w1^2 - 3 * s * w1 - 5 * mud * s)

      e3 <- c(
        1 + (25 * w1^2 * w2^2 - 3 * sqrt(5) * (3 * theta^3 + 5 * theta * w1 * w2) * (w1 + w2) + 15 * theta^2 * (w1^2 + w2^2 + 3 * w1 * w2)) / (9 * theta^4),
        (18 * sqrt(5) * theta^3 + 15 * sqrt(5) * theta * (w1^2 + w2^2) - 75 * theta^2 * (w1 + w2) - 50 * w1 * w2 * (w1 + w2) + 60 * sqrt(5) * theta * w1 * w2) / (9 * theta^4),
        5 * (5 * w1^2 + 5 * w2^2 + 15 * theta^2 - 9 * sqrt(5) * theta * (w1 + w2) + 20 * w1 * w2) / (9 * theta^4),
        10 * (3 * sqrt(5) * theta - 5 * (w1 + w2)) / (9 * theta^4),
        25 / (9 * theta^4)
      )
      e4 <- c(
        1 + (25 * w1^2 * w2^2 + 3 * sqrt(5) * (3 * theta^3 - 5 * theta * w1 * w2) * (w2 - w1) + 15 * theta^2 * (w1^2 + w2^2 - 3 * w1 * w2)) / (9 * theta^4),
        5 * (3 * sqrt(5) * theta * (w2^2 - w1^2) + 3 * theta^2 * (w1 + w2) - 10 * w1 * w2 * (w1 + w2)) / (9 * theta^4),
        5 * (5 * w1^2 + 5 * w2^2 - 3 * theta^2 - 3 * sqrt(5) * theta * (w2 - w1) + 20 * w1 * w2) / (9 * theta^4),
        -50 * (w1 + w2) / (9 * theta^4),
        25 / (9 * theta^4)
      )
      e5 <- c(
        1 + (25 * w1^2 * w2^2 + 3 * sqrt(5) * (3 * theta^3 + 5 * theta * w1 * w2) * (w1 + w2) + 15 * theta^2 * (w1^2 + w2^2 + 3 * w1 * w2)) / (9 * theta^4),
        (18 * sqrt(5) * theta^3 + 15 * sqrt(5) * theta * (w1^2 + w2^2) + 75 * theta^2 * (w1 + w2) + 50 * w1 * w2 * (w1 + w2) + 60 * sqrt(5) * theta * w1 * w2) / (9 * theta^4),
        5 * (5 * w1^2 + 5 * w2^2 + 15 * theta^2 + 9 * sqrt(5) * theta * (w1 + w2) + 20 * w1 * w2) / (9 * theta^4),
        10 * (3 * sqrt(5) * theta + 5 * (w1 + w2)) / (9 * theta^4),
        25 / (9 * theta^4)
      )


      return(exp((10 * s + sqrt(5) * theta * (w1 + w2 - 2 * m)) / theta^2) *
        (e3 %*% lambda31 * pnorm((muc - w2) / sqrt(s)) +
          e3 %*% lambda32 * sqrt(s) / sqrt(2 * pi) * exp(-(w2 - muc)^2 / (2 * s))) +
        exp(-(sqrt(5) * (w2 - w1) / theta)) *
          (e4 %*% lambda41 * (pnorm((w2 - m) / sqrt(s)) - pnorm((w1 - m) / sqrt(s))) +
            e4 %*% lambda42 * sqrt(s) / sqrt(2 * pi) * exp(-(w1 - m)^2 / (2 * s)) -
            e4 %*% lambda43 * sqrt(s) / sqrt(2 * pi) * exp(-(w2 - m)^2 / (2 * s))) +
        exp((10 * s - sqrt(5) * theta * (w1 + w2 - 2 * m)) / theta^2) *
          (e5 %*% lambda51 * pnorm((w1 - mud) / sqrt(s)) +
            e5 %*% lambda52 * sqrt(s) / sqrt(2 * pi) * exp(-(w1 - mud)^2 / (2 * s))))
    } else {
      muc <- m - 2 * sqrt(5) * s / theta
      mud <- m + 2 * sqrt(5) * s / theta

      lambda31 <- c(1, muc, muc^2 + s, muc^3 + 3 * muc * s, muc^4 + 6 * muc^2 * s + 3 * s^2)
      lambda32 <- c(0, 1, muc + w1, muc^2 + 2 * s + w1^2 + muc * w1, muc^3 + w1^3 + muc^2 * w1 + muc * w1^2 + 3 * s * w1 + 5 * muc * s)
      lambda41 <- c(1, m, m^2 + s, m^3 + 3 * m * s, m^4 + 6 * m^2 * s + 3 * s^2)
      lambda42 <- c(0, 1, m + w2, m^2 + 2 * s + w2^2 + m * w2, m^3 + w2^3 + w2 * m^2 + m * w2^2 + 3 * s * w2 + 5 * s * m)
      lambda43 <- c(0, 1, m + w1, m^2 + 2 * s + w1^2 + m * w1, m^3 + w1^3 + w1 * m^2 + m * w1^2 + 3 * s * w1 + 5 * s * m)
      lambda51 <- c(1, -mud, mud^2 + s, -mud^3 - 3 * mud * s, mud^4 + 6 * mud^2 * s + 3 * s^2)
      lambda52 <- c(0, 1, -mud - w2, mud^2 + 2 * s + w2^2 + mud * w2, -mud^3 - w2^3 - mud^2 * w2 - mud * w2^2 - 3 * s * w2 - 5 * mud * s)

      e3 <- c(
        1 + (25 * w1^2 * w2^2 - 3 * sqrt(5) * (3 * theta^3 + 5 * theta * w1 * w2) * (w1 + w2) + 15 * theta^2 * (w1^2 + w2^2 + 3 * w1 * w2)) / (9 * theta^4),
        (18 * sqrt(5) * theta^3 + 15 * sqrt(5) * theta * (w1^2 + w2^2) - 75 * theta^2 * (w1 + w2) - 50 * w1 * w2 * (w1 + w2) + 60 * sqrt(5) * theta * w1 * w2) / (9 * theta^4),
        5 * (5 * w1^2 + 5 * w2^2 + 15 * theta^2 - 9 * sqrt(5) * theta * (w1 + w2) + 20 * w1 * w2) / (9 * theta^4),
        10 * (3 * sqrt(5) * theta - 5 * (w1 + w2)) / (9 * theta^4),
        25 / (9 * theta^4)
      )
      e4 <- c(
        1 + (25 * w1^2 * w2^2 + 3 * sqrt(5) * (3 * theta^3 - 5 * theta * w1 * w2) * (w1 - w2) + 15 * theta^2 * (w1^2 + w2^2 - 3 * w1 * w2)) / (9 * theta^4),
        5 * (3 * sqrt(5) * theta * (w1^2 - w2^2) + 3 * theta^2 * (w1 + w2) - 10 * w1 * w2 * (w1 + w2)) / (9 * theta^4),
        5 * (5 * w1^2 + 5 * w2^2 - 3 * theta^2 - 3 * sqrt(5) * theta * (w1 - w2) + 20 * w1 * w2) / (9 * theta^4),
        -50 * (w1 + w2) / (9 * theta^4),
        25 / (9 * theta^4)
      )
      e5 <- c(
        1 + (25 * w1^2 * w2^2 + 3 * sqrt(5) * (3 * theta^3 + 5 * theta * w1 * w2) * (w1 + w2) + 15 * theta^2 * (w1^2 + w2^2 + 3 * w1 * w2)) / (9 * theta^4),
        (18 * sqrt(5) * theta^3 + 15 * sqrt(5) * theta * (w1^2 + w2^2) + 75 * theta^2 * (w1 + w2) + 50 * w1 * w2 * (w1 + w2) + 60 * sqrt(5) * theta * w1 * w2) / (9 * theta^4),
        5 * (5 * w1^2 + 5 * w2^2 + 15 * theta^2 + 9 * sqrt(5) * theta * (w1 + w2) + 20 * w1 * w2) / (9 * theta^4),
        10 * (3 * sqrt(5) * theta + 5 * (w1 + w2)) / (9 * theta^4),
        25 / (9 * theta^4)
      )


      return(exp((10 * s + sqrt(5) * theta * (w1 + w2 - 2 * m)) / theta^2) *
        (e3 %*% lambda31 * pnorm((muc - w1) / sqrt(s)) +
          e3 %*% lambda32 * sqrt(s) / sqrt(2 * pi) * exp(-(w1 - muc)^2 / (2 * s))) +
        exp(-(sqrt(5) * (w1 - w2) / theta)) *
          (e4 %*% lambda41 * (pnorm((w1 - m) / sqrt(s)) - pnorm((w2 - m) / sqrt(s))) +
            e4 %*% lambda42 * sqrt(s) / sqrt(2 * pi) * exp(-(w2 - m)^2 / (2 * s)) -
            e4 %*% lambda43 * sqrt(s) / sqrt(2 * pi) * exp(-(w1 - m)^2 / (2 * s))) +
        exp((10 * s - sqrt(5) * theta * (w1 + w2 - 2 * m)) / theta^2) *
          (e5 %*% lambda51 * pnorm((w2 - mud) / sqrt(s)) +
            e5 %*% lambda52 * sqrt(s) / sqrt(2 * pi) * exp(-(w2 - mud)^2 / (2 * s))))
    }
  }
}
