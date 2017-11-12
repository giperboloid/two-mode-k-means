Two_mode_k_means <- function(n, m, K, L, P, Q, X) {
  # Improves an initial partition by minimizing of f(P, Q, V) using 
  # formulas for changing:    
  #   P: Cluster membership matrix of the rows with K the number
  #      of row clusters, p[i, k] = 1 if row i belongs to row cluster k,
  #      and p[i, k] = 0 otherwise.
  #   Q: Cluster membership matrix of the columns with L the number
  #      of column clusters q[j, l] = 1 if column j belongs to column cluster l,
  #      and q[j, l] = 0 otherwise.
  #   V: Matrix with cluster centers for row cluster k and
  #      column cluster l.
  #
  # Args:
  #   n: Number of rows.
  #   m: Number of columns.
  #   K: Number of row clusters.
  #   L: Number of column clusters.
  #   P: Described above; size: [n, K].
  #   Q: Described above; size: [m, L].
  #   X: Two-mode data matrix; size: [n, m].
  #
  # Returns:
  #   List with the elements: 
  #     P: Updated P. 
  #     Q: Updated Q.
  #     vaf: Variance accounted for (VAF) criterion, which is comparable to
  #          the R^2 - measure used in regression analysis. The optimal value of VAF 
  #          ranges from 0 to 1 , and maximizing VAF corresponds to minimizing f(P, Q, V).
  V <- UpdateV(X, P, Q)
  old.V = 0
  while (any(abs(as.vector(old.V) - as.vector(V)) > 10 ^ -12)) {
    old.V = V
    
    P <- UpdateP(n, K, L, X, Q, V)
    P <- FillEmptyClustersInP(n, m, K, P, Q, V) 

    V <- UpdateV(X, P, Q)
    
    Q <- UpdateQ(m, K, L, X, P, V)
    Q <- FillEmptyClustersInQ(n, m, L, P, Q, V) 

    V <- UpdateV(X, P, Q)
  }
  
  vaf <- 1 - sum((X - P %*% V %*% Conj(t(Q)))^2) / sum(colSums((X - sum(colSums(X)) / (n * m)) ^ 2))
  return.list <- list("vaf" = vaf,"P" = P, "Q" = Q)
  return (return.list)
}

UpdateV <- function(X, P, Q) {
  V <- Conj(t(P)) %*% X %*% Q / (colSums(P) %*% t(colSums(Q)))
  return (V)
}

UpdateP <- function(n, K, L, X, Q, V) {
  SL <- colSums(Q) ^ 0.5
  Xsl <- X %*% Q %*% diag(1 / SL);
  Vsl <- diag(SL) %*% Conj(t(V))
  d <- matrix(0, n, K)
  Vsl <- Conj(t(as.vector(Vsl)))
  t <- array(Vsl[matrix(1, n),], c(n, L, K)) 
  for (k in 1 : K) {
    d[, k] = rowSums((Xsl - t[,, k]) ^ 2)
  }
  P <- matrix(0, n, K)
  min.index <- apply(d, 1, which.min)
  for (i in 1 : n) {
    P[i, min.index[i]] <- 1
  }
  return (P)
}

FillEmptyClustersInP <- function(n, m, K, P, Q, V) {
  if (any(colSums(P) == 0)) {
    temp.sum = colSums(P)
    QV <- Q %*% Conj(t(V))
    QV <- Conj(t(as.vector(QV)))
    ones = matrix(1, n)
    t <- array(QV[ones,], c(n, m, K))
    for (k in 1 : K) {
      d[,k] <- rowSums((X - t[,,k])^2)
    }
    distance <- rowSums(P * d)
    for (k in 1 : K) {
      if (temp.sum[k]) {
        max <- -1
        for (i in 1 : n) {
          if (distance[i] > max && sum(sum(ones %*% P[i,] * P)) > 1) {
            max.index <- i
            max <- distance[i]
          }
        }
        distance[max.index] <- 0
        P[max.index,] <- 0
        P[max.index, k] <- 1
      }
    }
  }
  return (P)
}

UpdateQ <- function(m, K, L, X, P, V) {
  SL = colSums(P) ^ 0.5
  Xsl = Conj(t(X)) %*% P %*% diag(1 / SL);
  Vsl = diag(SL) %*% V
  d = matrix(0, m, L)
  Vsl = Conj(t(as.vector(Vsl)))
  t = array(Vsl[matrix(1, m),], c(m, K, L))
  for (l in 1 : L) {
    d[, l] = rowSums((Xsl - t[,,l]) ^ 2)
  }
  Q = matrix(0, m, L)
  min.index = apply(d, 1, which.min)
  for (i in 1 : m) {
    Q[i, min.index[i]] = 1
  }
  return (Q)
}

FillEmptyClustersInQ <- function(n, m, L, P, Q, V) {
  if (any(colSums(Q) == 0)) {
    temp.sum = colSums(Q)
    PV <- Conj(t(as.vector(P %*% V)))
    ones = matrix(1, m)
    t = array(PV[ones,], c(m, n, L))
    for ( l in 1 : L) {
      d[, l] <- rowSums((Conj(t(X)) - t[,,l]) ^ 2)
    }
    distance <- rowSums(Q * d)
    for (l in 1 : L) {
      if (temp.sum[l]) {
        max <- -1
        for (j in 1 : m) {
          if (distance[j] > max && sum(sum(ones %*% Q[j,] * Q)) > 1) {
            max.index <- j
            max <- distance[j]
          }
        }
        distance[max.index] <- 0
        Q[max.index,] <- 0
        Q[max.index,l] <- 1
      }
    }
  }
  return (Q)
}