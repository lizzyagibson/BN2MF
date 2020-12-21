# Testing annealed VCI
# on scores

library(R.matlab)
library(tidyverse)

#####
# Low Noise example
# DGP main simulations
#####

# Load simulations aka TRUTH
load("./Results/Main/sim_dgp_rep1.RDA")
true_patterns = sim_dgp_rep1[1,]$true_patterns[[1]]
true_scores = sim_dgp_rep1[1,]$true_scores[[1]]

# Load some descriptive stats from matlab
t0_grid <- readMat("./MATLAB/test_t0/t0_grid_out.mat")[[1]] %>% as_tibble()
colnames(t0_grid) = c("t0", "rank", "score", "iter", "mean_diff", "mean_T", "mean_reg")

# Load solution matrices with no annealing
ewa_reg <- readMat("./MATLAB/test_t0/ewa_reg.mat")[[1]]
varwa_reg <- readMat("./MATLAB/test_t0/varwa_reg.mat")[[1]]
eh_reg <- readMat("./MATLAB/test_t0/eh_reg.mat")[[1]]
varh_reg <- readMat("./MATLAB/test_t0/varh_reg.mat")[[1]]

tnot = c(seq(0, 0.01, by = .0005), seq(0.01,0.05, by = .001))

# load solution matrices with annealing across t0 values
t0 <- tibble()
for (i in 1:length(tnot)) {

  ewa <- readMat(paste0("./MATLAB/test_t0/ewa_t0", tnot[i], ".mat"))[[1]]
  ewa <- as_tibble(ewa) %>% nest(ewa = everything())
  
  varwa <- readMat(paste0("./MATLAB/test_t0/varwa_t0", tnot[i], ".mat"))[[1]]
  varwa <- as_tibble(varwa) %>% nest(varwa = everything())
  
  eh <- readMat(paste0("./MATLAB/test_t0/eh_t0", tnot[i], ".mat"))[[1]]
  eh <- as_tibble(eh) %>% nest(eh = everything())

  varh <- readMat(paste0("./MATLAB/test_t0/varh_t0", tnot[i], ".mat"))[[1]]
  varh <- as_tibble(varh) %>% nest(varh = everything())
  
  all_t0 = bind_cols(ewa, varwa, eh, varh)
  t0 = bind_rows(t0, all_t0)
  
}

# Clean up solution matrix
t0_all <- bind_cols(t0_grid, t0)
t0_4 <- t0_all %>% 
  filter(rank == 4) %>% 
  mutate(ewa = map(ewa, as.matrix),
         varwa = map(varwa, as.matrix),
         eh = map(eh, as.matrix),
         varh = map(varh, as.matrix),
         upperwa = map2(ewa, varwa, function(x,y) x + (1.96*sqrt(y))),
         lowerwa= map2(ewa, varwa, function(x,y) x - (1.96*sqrt(y))))

#####
# Quick Viz
#####

t0_4 %>% 
  arrange(desc(mean_diff))
# biggest var at 0.017

ggplot(t0_4, aes(x = t0, y = mean_diff)) + geom_point()
ggplot(t0_all, aes(x = t0, y = rank)) + geom_line()

#####
# Single Example
#####

# take the biggest var and see how often the truth is within
widest <- t0_4 %>% 
  arrange(desc(mean_diff)) %>% 
  slice(1)

upperwa <- widest$upperwa[[1]]
ewa <- widest$ewa[[1]]
lowerwa <- widest$lowerwa[[1]]

# need factor correspondence to rearrage solution
source("./Results/factor_correspondence.R")

# rearrange EWA
fc_out <- factor_correspondence(true_scores, ewa)

upperwa <- upperwa %*% fc_out$permutation_matrix
lowerwa <- lowerwa %*% fc_out$permutation_matrix
ewa <- fc_out$rearranged

sum(upperwa >= true_scores & lowerwa <= true_scores)/(1000*4)
# This is a bad proportion

# Try it L2 normalized
ewa_denom <- apply(ewa, 1, function(x) sqrt(sum(x^2)))
upperwa_norm <- upperwa/ewa_denom
# apply(upperwa_norm, 1, function(x) sqrt(sum(x^2)))
lowerwa_norm <- lowerwa/ewa_denom
# apply(lowerwa_norm, 1, function(x) sqrt(sum(x^2)))

norm_scores_denom <- apply(true_scores, 1, function(x) sqrt(sum(x^2)))
norm_scores <- true_scores/norm_scores_denom
# apply(norm_scores, 1, function(x) sqrt(sum(x^2)))

sum(upperwa_norm >= norm_scores & lowerwa_norm <= norm_scores)/(1000*4)

head(upperwa_norm)
head(norm_scores)
head(lowerwa_norm)

#####
# Repeat across dataframe
#####

t0_4 = t0_4 %>% 
  mutate(fc_out = map(ewa, function(x) factor_correspondence(true_scores, x)),
         upperwa = map2(upperwa, fc_out, function(x,y) x %*% y$permutation_matrix),
         lowerwa = map2(lowerwa, fc_out, function(x,y) x %*% y$permutation_matrix),
         ewa = map(fc_out, function(x) x$rearranged))

# Un-normalized
get_prop <- function(upperwa, lowerwa){
  sum(upperwa >= true_scores & lowerwa <= true_scores)/(1000*4)
}

t0_4 = t0_4 %>% 
  mutate(prop = map2(upperwa, lowerwa, get_prop))

t0_4 %>% 
  dplyr::select(t0, score, iter, prop) %>% 
  unnest(prop) %>% 
  arrange(desc(prop))

# Try it L2 normalized
get_norm <- function(ewa, dat) {
  denom <- apply(ewa, 1, function(x) sqrt(sum(x^2)))
  dat_norm <- dat/denom
  dat_norm
}

t0_4 = t0_4 %>% 
  mutate(upperwa_norm = map2(ewa, upperwa, get_norm),
         lowerwa_norm = map2(ewa, lowerwa, get_norm))

get_prop_norm <- function(upperwa, lowerwa){
  sum(upperwa >= norm_scores & lowerwa <= norm_scores)/(1000*4)
}

t0_4 = t0_4 %>% 
  mutate(prop_norm = map2(upperwa_norm, lowerwa_norm, get_prop_norm))

t0_4 %>% 
  dplyr::select(t0, score, iter, mean_diff, prop, prop_norm) %>% 
  unnest(c(prop, prop_norm)) %>% 
  arrange(desc(prop_norm))

# compare with no annealing
upper_reg <- ewa_reg + 1.96*sqrt(varwa_reg)
lower_reg <- ewa_reg - 1.96*sqrt(varwa_reg)

sum(upper_reg >= true_scores & lower_reg <= true_scores)/(1000*4)

upper_reg_norm = get_norm(ewa_reg, upper_reg)
lower_reg_norm = get_norm(ewa_reg, lower_reg)

sum(upper_reg_norm >= norm_scores & lower_reg_norm <= norm_scores)/(1000*4)

