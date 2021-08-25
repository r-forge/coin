
psmirnov <- function(obs, m, n = length(obs) - m, c) {

    obs <- sort(obs)
    TIES <- c(diff(obs) != 0, TRUE)
    cat("Ties =", TIES, "\n")
    k <- diag <- 1
    u <- 0
    cat("m n =", m, n, "\n")
    repeat {
      cat("k =", k, "\n")
      u <- c(u, 1 + u[length(u)])
      v <- k - u
      diag_bit <- (u <= m) & (v <= n) & (u >= 0) & (v >= 0)
      u <- u[diag_bit]
      v <- v[diag_bit]
      cat("u =", u, "\n")
      cat("v =", v, "\n")
      d <- abs(u / m - v / n)
      diag <- (c(diag, 0) + c(0, diag))[diag_bit]
      cat("diag =", diag, "\n")
      if (TIES[k])
          diag <- diag * (c > d)
      cat("diag =", diag, "\n")
      k <- k + 1
      if ((m + n) < k) break
    }
    cat("final diag =", diag, "\n")
    1 - diag / exp(lgamma(m + n + 1) - lgamma(m + 1) - lgamma(n + 1))
}

psmirnov(1:7, m = 3, c = 1/2)
psmirnov(1:12, m = 5, c = 3 / 7)
psmirnov(c(1, 2, 2, 3, 3, 1, 2, 3, 3, 4, 5, 6), m = 5, c = 3 / 7)
