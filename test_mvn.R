N <- 100
d <- 5
mu <- 1:d
X <- t(t(matrix(rnorm(N*d), N, d)) + mu)
tmp <- matrix(rnorm(d^2), d, d)
mcov <- tcrossprod(tmp, tmp)  + diag(0.5, d)
myChol <- chol(mcov)

sum(dmvn(X, mu, diag(rep(1, times = ncol(X))), log = TRUE))
sum(dmvnorm(X, mu, diag(rep(1, times = ncol(X))), log = TRUE))
# Same!

diag(rep(1, times = ncol(X)))
chol(diag(rep(1, times = ncol(X))))
# Same!
        