##################################
## Gaussian Mixed Membership Model
## Gibb's Sampler
## Lizzy Gibson
## 12/4/2019
##################################

source("./core_loop.R")
library(mvtnorm)
library(tidyverse)
library(LaplacesDemon)
library(mvnfast)

##
## Definitions:
##
## x = exposure matrix (n*p)
   ## sigma_2 = hyperparameter
   ## x_ij ~ N(mu_z_i, sigma_2)
## mu = pattern distributions over chemicals (k*p)
   ## lambda_2 = hyperparameter 
   ## mu_k ~ N(0, I*lambda_2)
## z = pattern assignments (each chemical concentration per person assigned to a topic) (n*p*k)
   ## z_ij ~ Categorical(theta_i)
## theta = individual distributions over patterns (k*n)
   ## alpha = hyperparameter
   ## theta_i ~ Dirichlet(alpha)
## n_k = # of chemicals assigned to pattern k
## x_bar_k = average concentration within a pattern
## k = # of clusters
##

## Model inputs:
##
## x, k, alpha, m, lambda_2, sigma_2

### TEST
x <- matrix(rnorm(70), nrow = 10)
n <- nrow(x)
k <- 3
# p <- ncol(x)
# alpha <- rep(1, times = k) # Uniform prior over cluster simplex
z <- matrix(rep(sample(1:3, ncol(x), replace = TRUE), times = nrow(x)), nrow = nrow(x), ncol = ncol(x)) 
#   # if this random sample doesn't include all k, results will be wrong dim !!!
m <- 0
lambda_2 <- 1
sigma_2 <- 1
# theta <- rdirichlet(n, alpha) # (n*k)

# Cholesky is of an identity matrix is the identity matrix
# Matrix version of square root
=======
# x <- matrix(rnorm(70), nrow = 10)
# n <- nrow(x)
# k <- 3
# p <- ncol(x)
# alpha <- rep(1, times = k) # Uniform prior over cluster simplex
# z <- matrix(rep(sample(1:3, ncol(dat), replace = TRUE), times = nrow(dat)), nrow = nrow(dat), ncol = ncol(dat)) 
#   # if this random sample doesn't include all k, results will be wrong dim !!!
# m <- 0
# lambda_2 <- 1
# sigma_2 <- 1
# theta <- rdirichlet(n, alpha) # (n*k)

################
## Pattern Means
################

mu_update <- function (x, k, z, sigma_2, m, lambda_2, mu) {

mu_update <- function (x, k, z, sigma_2, m, lambda_2) {
  #  Arguments
  #  k: number of patterns
  #  z: assignments from last step
  #  x: observed data
  #  sigma_2: variance hyperparameter for observed data
  #  m: mean hyperparameter for Gaussian prior on mu
  #  lambda_2: variance hyperparameter for Gaussian prior on mu
  
  n_k <- apply(z, 2, tabulate, nbins = k) # n_k = # of chemicals assigned to pattern k
  n <- nrow(x)
  p <- ncol(x)
  alpha <- rep(1, times = k) # Uniform prior over cluster simplex

  theta <- rdirichlet(n, alpha) # (n*k)
  z <- t(apply(theta, 1, function(w) rcat(p, p = w))) # (p*n)
  
  n_k <- apply(z, 2, tabulate, nbins = k) # n_k = # of chemicals assigned to pattern k
  p   <- ncol(x) # number of chemicals
  mu  <- matrix(NA, nrow = k, ncol = p) # empty matrix to fill in for mu
  
  # update: Draw mu_k ~ N(0, I*lambda_2) for each pattern
  for (i in 1:k) {
    
    n_ki <- n_k[i,] # Number of chemicals assigned to pattern k
    
    x_new <- x
    x_new[z != i] <- 0 # Set chemical values not assigned to this pattern as zero, then take column ave
    x_bar_k <- colMeans(x_new) # Eq. 5.21

    m_hat <- x_bar_k * ((n_ki / sigma_2) / ((n_ki / sigma_2) + (1 / lambda_2))) # Eq. 5.23
    lambda_2_hat <- 1 / ((n_ki / sigma_2) + (1 / lambda_2)) # Eq. 5.24
    mu[i, ] <- abs(rmvnorm(1, mean = m_hat, sigma = lambda_2_hat * diag(p)))
    # RNG for the multivariate normal distribution with mean and covariance matrix sigma_2.
  }
  mu
}

# mu <- mu_update(x = x, k = 3, z= z, sigma_2 = 1, m = 0, lambda_2 = 1, mu)

###################################
## Individual Pattern Distributions
###################################

theta_update <- function (x, k, z) {
  #  Arguments
  #  z: assignments from last step
  #  alpha: hyperparameter on theta
  n         <- nrow(x) 
  alpha     <- rep(1, times = k)
  alpha_new <- t(apply(z, 1, tabulate, nbins = k) + alpha) # alpha_new (n*k)
  theta     <- rdirichlet(n, alpha_new) # (n*k)
  
  theta
}

# theta <- theta_update(x = x, k = 3, z = z)

##############
## Assignments
##############

z_update <- function (x, mu, theta, sigma_2, z) {

  #  Arguments
  #  x: observed data
  #  mu: pattern means from last step
  #  theta: individual scores from last step
  #  sigma_2: variance hyperparameter for observed data
  
  k <- nrow(mu)
  n <- nrow(x)
  p <- ncol(x)
  
  for (i in 1:n) {
    
    prob_n <- matrix(NA, nrow = p, ncol = k) # Probability matrix for each individual [i is fixed]
    
    for (j in 1:p) {
    
      for (u in 1:k) {

        prob_n[j,u] <- theta[i,u] * dmvn(x[i,j], mu = mu[u,j], sigma = sigma_2, log = TRUE) # Eq. 5.17
        # density function for the multivariate normal distribution with mean and covariance matrix sigma_2.
      }
        prob_n      <- prob_n / rowSums(prob_n) # this normalizes chemical proportions over all patterns
        z[i,j]      <- rcat(1, p = prob_n[j,]) # this chooses most likely pattern assignment
    }
  }
  z
}

# z <- z_update(x = x, mu = mu, theta = theta, sigma_2 = sigma_2, z = z)

#########################
## Log Joint aka Deviance
## Equal to the log posterior up to an additive constant
## To assess convergence
#########################

log_joint <- function(theta, mu, z, x, sigma_2, m, lambda_2) {
  
  k     <- nrow(mu)
  p     <- ncol(x)
  n     <- nrow(x)
  lp    <- 0
  alpha <- rep(1, times = nrow(mu))
  # joint p(theta, mu, z, x)
  
  # p(theta)
  lp <- lp + sum(apply(theta, 1, ddirichlet, alpha = alpha, log = TRUE))
  
  # p(mu)
  lp <-  lp + sum(apply(mu, 1, dmvn, mu = rep(m, times = p), sigma = diag(p), log = TRUE))
  
  for (i in 1:n) { # p(z)
    z_multi <- tabulate(z[i,], nbins = k) # Count number of chemicals assigned to each pattern in each person
    lp <- lp + dmultinom(z_multi, prob = theta[i,], log = TRUE)
  }

  for (i in 1:n) { # p(x)
    for (j in 1:p) {

    lp <- lp + dmvn(x[i,], mu = mu[z[i,j], ], sigma = diag(p), log = TRUE)
    }}

  lp
}

# lp <- log_joint(theta, mu, z, x, sigma_2 = 1, m = 0, lambda_2 = 1) 

#################
## Gibb's Sampler
#################

GMMM <- function(x, k, sigma_2, m, lambda_2){
            harness(
                  InitState = function(){ # Initialize theta and z
                      alpha      <- rep(1, times = k)
                      theta_init <- rdirichlet(nrow(x), alpha) # (n*k)
                      z_init     <- t(apply(theta_init, 1, function(w) rcat(ncol(x), p = w))) # (p*n)
                                   # matrix(1/k, nrow = n, ncol = p)
                      mu_init    <- rmvn(k, mu = rep(m, times = ncol(x)), sigma = diag(ncol(x)))
                      list(theta = theta_init, z = z_init, mu = mu_init)
                  }
          
                , TransitionProposal = function(previousState){ # Update Theta, Mu, Z
                      propose_mu    <- mu_update(x, k, previousState$z, sigma_2, m, lambda_2)
                      propose_z     <- z_update(x, propose_mu, previousState$theta, sigma_2)
                      propose_theta <- theta_update(x, k, propose_z)
                      deviance      <- log_joint(propose_theta, propose_mu, propose_z, x, sigma_2, m, lambda_2)
                      list(mu = propose_mu, z = propose_z, theta = propose_theta, deviance = deviance)
                  }
                
                , ApplyTransition = function(previousState, proposal){ # For Gibbs Sampler, transition == TRUE
                      proposal
                }
                
                , ShouldWeTerminate = function(step,state,proposal){# For now, 100 iterations
                  (step > 1000)
                  } 
                )
      }
                  } 

x <- matrix(rnorm(70), nrow = 10)
ex_run <- GMMM(x, k = 3, sigma_2 = 1, m = 0, lambda_2 = 1)
# ex_run$FinalState

# x <- matrix(rnorm(50000), nrow = 1000)
start_time <- Sys.time()
ex_run <- GMMM(x, k = 3, sigma_2 = 1, m = 0, lambda_2 = 1)
end_time <- Sys.time()
end_time - start_time

ex_run$TheStateRecord[grep("dev", names(ex_run$TheStateRecord))] %>% as_vector() %>% 
  cbind(., 1:ex_run$StepCount) %>% 
  as_tibble() %>% rename(deviance = 1, iteration = 2) %>% 
  ggplot(aes(x = iteration, y = deviance)) +
  geom_line() +
  theme_minimal() +
  labs(x = "Iterations", y = "Deviance")
