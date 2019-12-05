
library(MASS)
library(LaplacesDemon)

# NOTE: sigma and lambda represent variance parameters. On
# paper these are designated using sigma^2 and lambda^2.
# x is expected to be a matrix, with ncol = number of dimensions
# and nrow = sample size.
# beta is an k * p matrix, where k = number of cluster and
# p = the number of dimension the data lives in.

assignment_conditional <- function (beta, theta, x, sigma) {
  k = nrow(beta)
  n = nrow(x)
  dim = ncol(x)
  z = rep(NA, n)

  for (j in 1:n) {
    p = rep(NA, k)
    for (i in 1:k) p[i] = dmvnorm(x[j, ], mean = beta[i, ],
                                  sigma = sigma * diag(dim)) * theta[i]
    p = p / sum(p)
    z[j] = sample(1:k, size = 1, prob = p)
  }

  z
}

proportions_conditional <- function (z, alpha) {
  alpha_posterior = as.vector(table(z)) + alpha
  theta = rdirichlet(1, alpha_posterior)

  as.vector(theta)
}

components_conditional <- function (z, x, sigma, mu, lambda) {
  n_cluster = as.vector(table(z))
  k = length(n_cluster)
  dim = ncol(x)
  # beta = rep(NA, k)
  beta = matrix(NA, nrow = k, ncol = dim)

  for (i in 1:k) {
    n_k = n_cluster[i]
    x_bar_k = colMeans(x[z == i, ])
    mu_hat = ((n_k * x_bar_k / sigma) + mu / lambda) / (n_k / sigma + 1 / lambda)
    lambda_hat = (1 / (n_k / sigma + 1 / lambda))^2
    beta[i, ] = rmvnorm(1, mean = mu_hat, sigma = lambda_hat * diag(dim))
  }
  beta
}

log_joint <- function(theta, beta, z, x, alpha, sigma, mu, lambda) {
  lp = log(ddirichlet(theta, alpha))
  k = nrow(beta)
  dim = ncol(x)

  for (i in 1:k) lp = lp + 
    dmvnorm(beta[i, ], mean = mu, sigma = lambda * diag(dim), log = TRUE)

  z_multi = as.vector(table(z))
  lp = lp + dmultinom(z_multi, prob = theta, log = TRUE)

  for (i in 1:n) lp = lp + dmvnorm(x[i, ], mean = beta[z[i], ],
                                   sigma = sigma * diag(dim), log = TRUE)

  lp
}

###############################################################################
## Gibb's sampler

gibbs_sampler <- function (init, x, sigma, mu, lambda,
                           n_sample, n_burn = 0, refresh = 100) {
  # The user manually specifies the number of iterations during the
  # warmup and the sample phases. We can use diagnostic tools to pick
  # a lag and select how much of the output we retain.
  # 
  # Arguments
  #  init: a list which contains initial values for theta, beta, and z.
  #  n_burn: number of burn_in samples.
  #  x: observed data (in matrix form).
  #  mu: mean hyper parameter for the Gaussian prior on beta
  #  lambda: variance hyper parameter for the Gaussian prior on beta
  #  refresh: tells the code when to print an iteration number.

  n_iter = n_burn + n_sample + 1
  n = nrow(x)
  dim = ncol(x)
  k = length(init$theta)

  # containers to save results from sampling
  z = matrix(NA, nrow = n_iter, ncol = n)
  theta = matrix(NA, nrow = n_iter, ncol = k)
  beta = array(NA, c(k, dim, n_iter))  # tensor
  lp = rep(NA, n_iter)
  
  # generated quantities (generate one observation per chain)
  z_pred = rep(NA, n_iter)
  x_pred = matrix(NA, nrow = n_iter, ncol = dim)

  # iteration "0"
  theta[1, ] = init$theta
  beta[, , 1] = init$beta
  z[1, ] = init$z
  lp[1] = log_joint(theta[1, ], beta[, , 1], z[1, ], 
                    x, alpha, sigma, mu, lambda)
  # generated quantities
  z_pred[1] = rcat(n = 1, p = theta[1, ])
  x_pred[1, ] = rmvnorm(n = 1, mean = beta[z[1], , 1],
                        sigma = sigma * diag(dim))

  # sample using the full conditional
  for (i in 2:n_iter) {
    if (i %% refresh == 0) cat("Iteration:", i, "/", n_iter, "\n")
    z[i, ] = assignment_conditional(beta[, , i - 1], theta[i - 1, ], x, sigma)
    theta[i, ] = proportions_conditional(z[i, ], alpha)
    beta[, , i] = components_conditional(z[i, ], x, sigma, mu, lambda)

    lp[i] = log_joint(theta[i, ], beta[, , i], z[i, ],
                      x, alpha, sigma, mu, lambda)

    # generated quantities
    z_pred[i] = rcat(n = 1, p = theta[i, ])
    x_pred[i, ] = rmvnorm(n = 1, mean = beta[z_pred[i], , i],
                          sigma = sigma * diag(dim))
  }

  list(theta = theta, beta = beta, z = z, lp = lp,
       z_pred = z_pred, x_pred = x_pred)
}
