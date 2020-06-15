install.packages("matlib")

library(pracma)
library(far)
library(matlib)

U <- matrix(rnorm(1:50), ncol = 5)
V <- matrix(rnorm(1:60), ncol = 6)

## pracma
# gramSchmidt(U, tol = .Machine$double.eps^0.5)

## far
# orthonormalization(U, basis=TRUE, norm=TRUE)

## matlib
# GramSchmidt(U, normalize = TRUE, verbose = FALSE,
#            tol = sqrt(.Machine$double.eps))

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
  
  qrU <- qr.Q(qr(U))
  qrV <- qr.Q(qr(V))
  
  svdU <- svd(U)$u
  svdV <- svd(V)$u
  
  m <- ncol(U)
  n <- ncol(V)
  
  sqrt( max(m,n) - sum((t(svdU) %*% svdV)^2) )
  dUV <- sqrt( max(m,n) - sum((t(qrU) %*% qrV)^2) )
  
  ratio <- if (dUV/sqrt( max(m,n)) <= 1/2) {
    paste(round(dUV/sqrt( max(m,n)),2), "<= 1/2, subspaces are similar")} else {
      paste(round(dUV/sqrt( max(m,n)),2), "> 1/2, subspaces are NOT similar")}
  
  list(symm_subspace_dist = dUV, similarity_ratio = ratio)

  }


