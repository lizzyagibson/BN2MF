#library(pracma)
#library(far)
#library(matlib)

#U <- matrix(rnorm(1:50), ncol = 5)
#V <- matrix(rnorm(1:60), ncol = 6)
#S <- U + 0.1

## base
qr(U) # upper triangle is R, lower triangle is info on Q...
## qy.qr(): return the results of the matrix multiplications: Q %*% y, 
## where Q is the order-nrow(x) orthogonal (or unitary) transformation represented by qr.
qr.Q(qr(U))
## Normal!!
apply(qr.Q(qr(U)), 2, function(x) sqrt(sum(x^2)))

svd(U)$u
apply(svd(U)$u, 2, function(x) sqrt(sum(x^2)))

## Orthogonal!! + Normal = ORTHONORMAL
t(qr.Q(qr(U))) %*% (qr.Q(qr(U)))

## Orthogonal!! + Normal = ORTHONORMAL
t(svd(U)$u) %*% (svd(U)$u)


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

