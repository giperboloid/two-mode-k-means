source("two_mode_k_means.R")
library(plotly)

Simulation <- function(n, m, K, L, ESD) {
  # Function generates P, Q, V, E and X then it passes all the required
  # parameters to the funciton Two_mode_k_means(n, m, K, L, P, Q, X).
  # As a result it receives VAF, new P and Q for generating new X and 
  # building a plot.
  #
  # Args:
  #   n: Number of rows.
  #   m: Number of columns.
  #   K: Number of row clusters.
  #   L: Number of column clusters.
  #   ESD: Error standard deviation
  #
  # Returns:
  #   Nothing. Function prints a VAF criterion and builds a heat map.
  P = RandMembMatrix(n, K)
  Q = RandMembMatrix(m, L)
  quantiles <- qnorm(seq(1 / (K * L + 1), 1 - 1 / (K * L + 1), 1 / (K * L + 1)))
  V <- matrix(sample(quantiles), nrow = K, ncol = L)
  E <- matrix(rnorm(n * m), n, m) 
  X <- P %*% V %*% Conj(t.default(Q)) + E * ESD
 
  return.list <- Two_mode_k_means(n, m, K, L, P, Q, X)
  cat("VAF = ", return.list$vaf)
  
  new.X <- GenerateNewX(m, n, K, L, X, return.list$P, return.list$Q)
  
  BuildPlot(n, new.X)
}

BuildPlot <- function(n, X){
  p <- plot_ly(
    x = seq(1, n, 1), y = seq(1, 1, 1),
    z = X, type = "heatmap"
  )
  p
}

RandMembMatrix <- function(nrow, ncol) {
  m <- matrix(0L, nrow = nrow, ncol = ncol)
  m[cbind(sequence(nrow), sample(ncol, nrow, TRUE))] <- 1L
  return (m)
}

GenerateNewX <- function(m, n, K, L, X, P, Q) {
  c = 1
  new.X = X
  for (i in 1 : K) {
    for (j in 1 : n) {
      if (P[j, i] == 1) {
        new.X[c,] = X[j,]
        c = c + 1
      }
    }
  }
  c = 1
  X = new.X
  for (i in 1 : L) {
    for (j in 1 : m) {
      if (Q[j, i] == 1) {
        new.X[, c] = X[, j]
        c = c + 1
      }
    }
  }
  return (new.X)
}

