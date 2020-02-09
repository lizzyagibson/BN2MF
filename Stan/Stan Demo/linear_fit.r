rm(list = ls())
gc()

set.seed(1954)

# Make sure to adjust your working and lib directory.

library(rstan)  # (here used version 2.19.2)
library(tidyverse)
library(plyr)

rstan_options(auto_write = TRUE)

## read data
data <- read_rdump("linear.data.r")

## create initial estimates
init <- function() {
  list(sigma = rgamma(1, 1),
       beta = rnorm(1, mean = 1, sd = 1))
}

## run Stan
fit <- stan(file = "linear_demo.stan",
            data = data,
            init = init)

# Remark: the result is stored in the fit object.
summary(fit)[1]

# Do graphical checks.
pars = c("beta", "sigma", "lp__")
traceplot(fit, inc_warmup = TRUE, pars = pars)
traceplot(fit, pars = pars)

stan_dens(fit, separate_chains = TRUE, pars = pars)
pairs(fit, pars = pars)


# Posterior predictive checks
data_pred <- data.frame(data$x, data$y)
names(data_pred) <- c("x", "y")


pred <- as.data.frame(fit, pars = "y_pred") %>%
  gather(factor_key = TRUE) %>% 
  group_by(key) %>%
  dplyr::summarize(lb = quantile(value, probs = 0.05),
            median = quantile(value, probs = 0.5),
            ub = quantile(value, probs = 0.95)) %>%
  bind_cols(data_pred)

ggplot(pred, aes(x = x, y = y)) + geom_point() + 
  geom_line(aes(x = x, y = median)) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.25)
