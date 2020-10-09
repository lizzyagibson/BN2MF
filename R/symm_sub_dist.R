#library(pracma)
#library(far)
#library(matlib)

#U <- matrix(rnorm(1:50), ncol = 5)
#V <- matrix(rnorm(1:60), ncol = 6)
#S <- U + 0.1

### Symmetric Subspace Distance
symm_subspace_dist <- function(U, V) {
  
  if (nrow(U) != max(nrow(U), ncol(U))) {U <- t(U)}
  if (nrow(V) != max(nrow(V), ncol(V))) {V <- t(V)}
  
  qrU <- qr.Q(qr(U))
  qrV <- qr.Q(qr(V))
  
  m <- ncol(U)
  n <- ncol(V)
  
  dUV <- sqrt( max(m,n) - sum((t(qrU) %*% qrV)^2) )
  ratio <- dUV/sqrt( max(m,n))
  ratio
}

symm_subspace_dist(U,V)
symm_subspace_dist(V,U)
symm_subspace_dist(S,U)

symm_subspace_dist(simo, simO)
simo = cbind(simO[,50:1])

