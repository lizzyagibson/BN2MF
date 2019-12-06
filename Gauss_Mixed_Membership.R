##################################
## Gaussian Mixed Membership Model
## Gibb's Sampler
## Lizzy Gibson
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
n <- nrow(x)
k <- 3
p <- ncol(x)
alpha <- rep(1, times = k) # Uniform prior over cluster simplex
z <- matrix(rep(sample(1:3, p, replace = TRUE), times = n), nrow = n, ncol = p) 
  # if this random sample doesn't include all k, results will be wrong dim !!!
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
  
  n_k <- apply(z, 2, table) # n_k = # of chemicals assigned to pattern k
  k   <- nrow(n_k) # number of patterns
  p   <- ncol(x) # number of chemicals
  mu  <- matrix(NA, nrow = k, ncol = p) # empty matrix to fill in for mu
  
  # update: Draw mu_k ~ N(0, I*lambda) for each pattern
  for (i in 1:k) {
    
    n_ki <- n_k[i,] # Number of chemicals assigned to pattern k
    
    x_new <- x
    x_new[z != i] <- 0 # Set chemical values not assigned to this pattern as zero, then take column ave
    x_bar_k <- colMeans(x_new) # Eq. 5.21

    m_hat <- (x_bar_k * ((n_ki / sigma) / (n_ki / sigma + 1 / lambda))) # Eq. 5.23
    lambda_hat <- (1 / (n_ki / sigma + 1 / lambda)) # Eq. 5.24
    mu[i, ] <- rmvnorm(1, mean = m_hat, sigma = lambda_hat * diag(p)) 
    # RNG for the multivariate normal distribution with mean and covariance matrix sigma.
  }
  mu
}

mu <- mu_update(z, x, sigma, m, lambda)

# should this sum to 1 over the vocabulary (all chemicals)?

###################################
## Individual Pattern Distributions
###################################

theta_update <- function (z, alpha) {
  #  Arguments
  #  z: assignments from last step
  #  alpha: hyperparameter on theta
  
  alpha_new <- t(apply(z, 1, table) + alpha) # alpha_new (k*n)
  theta     <- rdirichlet(n, alpha_new) # (k*n)
  
  theta
}

theta <- theta_update(z, alpha)

##############
## Assignments
##############

z_update <- function (x, mu, theta, sigma) {
  #  Arguments
  #  x: observed data
  #  mu: pattern means from last step
  #  theta: individual scores from last step
  #  sigma: variance hyperparameter for observed data
  
  k <- nrow(mu)
  n <- nrow(x)
  p <- ncol(x)
  z <- matrix(NA, nrow = n, ncol = p) # empty matrix to fill in for z

  for (i in 1:n) {
    for (j in 1:p) {
      
      prob_n <- matrix(NA, nrow = p, ncol = k) # Probability matrix for each individual [i is fixed]
    
      for (u in 1:k) {
        prob_n[j,u] <- theta[i,u] * dmvnorm(x[i,j], mean = mu[u,j]) # Eq. 5.17
        # density function for the multivariate normal distribution with mean and covariance matrix sigma.
      }
        prob_n        <- prob_n / rowSums(prob_n) # this normalizes chemical proportions over all patterns
        z[i,j]      <- sample(1:k, size = 1, prob = prob_n[j,]) # this chooses most likely pattern assignment
                       #rcat(1, p = prob_n[j,])
    }
  }
  z
}

z <- z_update(x, mu, theta, sigma)

#########################
## Log Joint aka Deviance
## Equal to the log posterior up to an additive constant
## To assess convergence
#########################

log_joint <- function(theta, mu, z, x, alpha, sigma, m, lambda) {
  
  k  <- nrow(mu)
  p  <- ncol(x)
  n  <- nrow(x)
  lp <- 0
  
  for (i in 1:n) {
    lp <- lp + log(ddirichlet(theta[i,], alpha)) # Shouldn't this be 1 number?
  }
  
  for (i in 1:k) {
    lp <- lp + dmvnorm(mu[i, ], mean = rep(m, times = p), sigma = lambda * diag(p), log = TRUE)
  }
    
  z_multi <- apply(z, 2, table) 
  
  table(z[,1])
  tabulate(z[1,])
  
  for (i in 1:n) {
    z_multi <- tabulate(z[1,])
    lp <- lp + dmultinom(z_multi, p = theta[i,], log = TRUE)
  }
  
  for (i in 1:n) {
    for (j in 1:p) {
    lp <- lp + dmvnorm(x[i,], mean = mu[z[i,j], ], sigma = sigma * diag(p), log = TRUE)
  }}
  lp
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


