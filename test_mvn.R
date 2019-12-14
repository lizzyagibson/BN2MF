dmvn(x, 0, 1, log = TRUE, ncores = 1, isChol = FALSE)

dmvnorm(x[i,j], mean = mu[u,j], log =TRUE)

N <- 100
d <- 5
mu <- 1:d
X <- t(t(matrix(rnorm(N*d), N, d)) + mu)
tmp <- matrix(rnorm(d^2), d, d)
mcov <- tcrossprod(tmp, tmp)  + diag(0.5, d)
myChol <- chol(mcov)

sum(dmvn(X, mu, diag(rep(1, times = ncol(X))), log = TRUE))

sum(dmvnorm(X, mu, diag(rep(1, times = ncol(X))), log = TRUE))

        