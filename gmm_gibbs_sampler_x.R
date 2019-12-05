rm(list = ls())
gc()
set.seed(1954)

library(ggplot2)
library(mvtnorm)

## test density functions
# hyper-priors
alpha <- c(0.55, 0.20, 0.25)
lambda <- 1.2
mu <- c(0.1, 0.1)
sigma <- 1

# simulate x's from a mixture of Gaussians.
n <- 100
k <- 3
theta_exact <- c(0.5, 0.3, 0.2)
beta_exact <- matrix(NA, nrow = 3, ncol = 2)
beta_exact[, 1] <- c(-5, 0, 5)
beta_exact[, 2] <- c(0, -5, 5)
# beta_exact <- c(-5, 0, 5)
z_exact <- sample(1:3, n, replace = TRUE, prob = theta_exact)

x <- matrix(NA, nrow = n, ncol = 2)
for (i in 1:n) x[i, ] <- rmvnorm(1, mean = beta_exact[z_exact[i], ], 
                                 sigma = sigma * diag(2))

# Given our simulation data, z should retrieve z_exact.
# On the other hand, theta will not be exact, because our sample is small.
# Beta should be fairly exact, given how spread out the centroids are.
z <- assignment_conditional(beta_exact, theta_exact, x, sigma)
theta <- proportions_conditional(z, alpha)
beta <- components_conditional(z, x, sigma, mu, lambda)
lp <- log_joint(theta, beta, z, x, alpha, sigma, mu, lambda)

# draft variable
theta0 <- c(0.33, 0.33, 0.34)
beta0 <- matrix(0, nrow = 3, ncol = 2)
z0 <- sample(1:3, n, replace = TRUE, prob = theta0)

init <- list(theta = theta0,
             beta = beta0,
             z = z0)

n_burn <- 9
n_sample <- 100

samples <- gibbs_sampler(init, x, sigma, mu, lambda, n_burn, n_sample)

# Make a few plots to test for convergence. The simulated converges
# very well. Note that during the "sample using the full conditonal",
# it is crutial to use the most updated parameter, not simply the parameter
# from the previous iteration.

n_iter = n_burn + n_sample + 1
plot <- ggplot(data = data.frame(x = 1:n_iter, y = samples$lp),
               aes(x = 1:n_iter, y = samples$lp)) +
  geom_point() + xlab("iteration") + ylab("lp")

plot(x = 1:n_iter, samples$lp)
plot(x = 1:n_iter, samples$theta[, 1])
plot(x = 1:n_iter, samples$theta[, 2])
plot(x = 1:n_iter, samples$z[, 1])
plot(x = 1:n_iter, samples$beta[1, 1,])

# Diagnostics for autocorrelation
lp_acf <- acf(samples$lp[n_burn:n_iter])
theta_acf <- acf(samples$theta[n_burn:n_iter, 1])  # does not mix as well
beta_acf <- acf(samples$beta[1, 1, n_burn:n_iter])  # mixes very quickly
z_acf <- acf(samples$z[, 1])  # (return an error message)

# Difficult to read what may be the optimal lag.
# Eyeball the result and use lag <- 5. Keep it in mind for a potential
# source of errors.
lag <- 5
theta1_sampled <- samples$theta[seq(n_burn + 1, n_iter, lag), 1]
theta2_sampled <- samples$theta[seq(n_burn + 1, n_iter, lag), 2]
z1_sampled <- samples$z[seq(n_burn + 1, n_iter, lag), 1]

hist(theta2_sampled)
hist(z1_sampled)  # does not work when there is only one data type.
