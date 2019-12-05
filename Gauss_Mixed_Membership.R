##################################
## Gaussian Mixed Membership Model
## Gibb's Sampler
## 12/4/2019
##################################

library(mvtnorm)
library(tidyverse)
library(LaplacesDemon)

##
## Definitions:
##
## x = exposure matrix (n*p)
   ## sigma = hyperparameter
   ## x_ij ~ N(mu_z_i, sigma)
## mu = pattern distributions over chemicals (k*p)
   ## lambda = hyperparameter 
   ## mu_k ~ N(0, I*lambda)
## z = pattern assignments (each chemical concentration per person assigned to a topic) (n*p*k)
   ## z_ij ~ Categorical(theta_i)
## theta = individual distributions over patterns (k*n)
   ## alpha = hyperparameter
   ## theta_i ~ Dirichlet(alpha)
## n_k = # of chemicals assigned to pattern k
## x_bar_k = average concentration within a pattern
## k = # of clusters
##

# NOTE: sigma and lambda represent variance parameters. On
# paper these are designated using sigma^2 and lambda^2.

## Model inputs:
##
## x, k, alpha, m, lambda, sigma

x <- matrix(rnorm(70), nrow = 10)
n <- nrow(dat)
k <- 3
p <- ncol(dat)
alpha <- rep(1, times = k) # Uniform prior over cluster simplex
z <- sample(1:3, n, replace = TRUE)
m <- 0
lambda <- 1
sigma <- 1

################
## Pattern Means
################

mu_update <- function (z, x, sigma, m, lambda) {
  #  Arguments
  #  z: assignments from last step
  #  x: observed data
  #  sigma: variance hyperparameter for observed data
  #  m: mean hyperparameter for Gaussian prior on mu
  #  lambda: variance hyperparameter for Gaussian prior on mu
  
  n_kk <- as.vector(table(z)) # n_k = # of chemicals assigned to pattern k
  k <- length(n_kk) # number of patterns
  p <- ncol(x) # number of chemicals
  mu <- matrix(NA, nrow = k, ncol = p) # empty matrix to fill in for mu
  
  # update: Draw mu_k ~ N(0, I*lambda) for each pattern
  for (i in 1:k) {
    n_k <- n_kk[i] # Number of chemicals assigned to pattern k
    x_bar_k <- colMeans(x[z == i, ]) # Eq. 5.21
    m_hat <- (x_bar_k * ((n_k / sigma) / (n_k / sigma + 1 / lambda))) # Eq. 5.23
    
    lambda_hat <- (1 / (n_k / sigma + 1 / lambda)) # Eq. 5.24
    
    mu[i, ] <- rmvnorm(1, mean = m_hat, sigma = lambda_hat * diag(p)) 
    # RNG for the multivariate normal distribution with mean and covariance matrix sigma.
  }
  mu
}

##############
## Assignments
##############

z_update <- function (mu, theta, x, sigma) {
  k = nrow(mu)
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
#################
## Gibb's Sampler
#################

harness(
      InitGibbState <- function(){ # Initialize theta and z
         theta_init <- rdirichlet(n, alpha) # (n*k)
             z_init <- apply(theta_init, 1, function(w) rcat(p, p = w)) # (p*n)
                       # matrix(1/k, nrow = n, ncol = p)
      }
    
    , TransitionProposal=function(x){rnorm(1) }
    , ApplyTransition=function(state,proposal){state + proposal} # For Gibbs Sampler, transition == TRUE
    , ShouldWeTerminate=function(step,state,proposal){(step > 10)} # For now, 10 iterations
    )
