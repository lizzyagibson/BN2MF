library(tidyverse)
library(rstan)
library(rstanarm)
library(bayesplot)
library(brms)
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

sim <- read.csv("./Sims/sim_data_500.csv")

stan_dat <- list(
  K = 5,
  N = nrow(sim),
  P = ncol(sim),
  X = sim,
  mu1 = rep(2, times = ncol(sim)),
  sigma1 = rep(5, times = ncol(sim)),
  mu2 = rep(2, times = 5),
  sigma2 = rep(5, times = 5)
)
stan_dat

fit_mmm <- stan(
  file = "./Stan/mmm.stan",  # Stan program
  data = stan_dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 2,              # number of cores (could use one per chain)
  verbose = TRUE
)

#save(fit_mmm, file = "mmm_out.RDA")
load(file = "mmm_out.RDA")

vi_model <- rstan::stan_model(
  file = "./Stan/mmm.stan")

fit_vi <- rstan::vb(
  object = vi_model,  # Stan program
  data = stan_dat     # named list of data
)

print(fit_vi)

la <- rstan::extract(fit_vi, permuted = FALSE) # return a list of arrays 
### return an array of three dimensions: iterations, chains, parameters 

### use S3 functions on stanfit objects
a2 <- as.array(fit_vi)
m <- as.matrix(fit_vi)
dim(m)
d <- as.data.frame(fit_vi)

rstanarm::launch_shinystan(fit)

format(bayesplot::available_ppc()) # omit ppc_
format(bayesplot::available_mcmc()) # omit mcmc_

mcmc_rhat(rhat(fit_vi))

library(tidyverse)
names(read_csv("./Sims/sim_data_500.csv"))
