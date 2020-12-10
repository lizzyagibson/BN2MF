library(R.matlab)
library(tidyverse)

# Load simulations aka TRUTH
load("/Users/lizzy/BN2MF/Sims/sim_dgp_rep1.RDA")
truth = sim_dgp_rep1[1,]$chem[[1]]
true_patterns = sim_dgp_rep1[1,]$true_patterns[[1]]
true_scores = sim_dgp_rep1[1,]$true_scores[[1]]

# Load some descriptive stats from matlab
t0_grid <- readMat("/Users/lizzy/BN2MF/MATLAB/test_t0/t0_grid_out.mat")[[1]] %>% as_tibble()
colnames(t0_grid) = c("t0", "rank", "score", "iter", "mean_diff", "mean_T", "mean_reg")
t0_grid

# Load solution matrices with no annealing
ewa_reg <- readMat("/Users/lizzy/BN2MF/MATLAB/test_t0/ewa_reg.mat")[[1]]
varwa_reg <- readMat("/Users/lizzy/BN2MF/MATLAB/test_t0/varwa_reg.mat")[[1]]
eh_reg <- readMat("/Users/lizzy/BN2MF/MATLAB/test_t0/eh_reg.mat")[[1]]
varh_reg <- readMat("/Users/lizzy/BN2MF/MATLAB/test_t0/varh_reg.mat")[[1]]
pred_reg <- ewa_reg %*% eh_reg

# load solution matrices with annealing across t0 values
t0 <- tibble()
for (i in 1:62) {

  ewa <- readMat(paste0("/Users/lizzy/BN2MF/MATLAB/test_t0/ewa_t0", tnot[i], ".mat"))[[1]]
  ewa <- as_tibble(ewa) %>% nest(ewa = everything())
  
  varwa <- readMat(paste0("/Users/lizzy/BN2MF/MATLAB/test_t0/varwa_t0", tnot[i], ".mat"))[[1]]
  varwa <- as_tibble(varwa) %>% nest(varwa = everything())
  
  eh <- readMat(paste0("/Users/lizzy/BN2MF/MATLAB/test_t0/eh_t0", tnot[i], ".mat"))[[1]]
  eh <- as_tibble(eh) %>% nest(eh = everything())

  varh <- readMat(paste0("/Users/lizzy/BN2MF/MATLAB/test_t0/varh_t0", tnot[i], ".mat"))[[1]]
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
         pred = map2(ewa, eh, function(x,y) x %*% y),
         norm = map(pred, function(x) norm(truth - x, "F")/norm(truth, "F"))) %>% 
  unnest(norm) %>% 
  mutate(upperwa = map2(ewa, varwa, function(x,y) x + (1.96*sqrt(y))),
         lowerwa= map2(ewa, varwa, function(x,y) x - (1.96*sqrt(y))))

# Check values
t0_4 %>% 
  arrange(norm)
#smallest norm at zero

t0_4 %>% 
  arrange(desc(mean_diff))
# biggest var at 0.017

ggplot(t0_4, aes(x = t0, y = norm)) + geom_point()
ggplot(t0_4, aes(x = t0, y = mean_diff)) + geom_point()
ggplot(t0_all, aes(x = t0, y = rank)) + geom_line()

# take the biggest var and see how often the truth is within
widest <- t0_4 %>% 
  arrange(desc(mean_diff)) %>% 
  slice(1)

upperwa <- widest$upperwa[[1]]
ewa <- widest$ewa[[1]]
lowerwa <- widest$lowerwa[[1]]

# need factor correspondence to rearrage solution
source("./R/factor_correspondence.R")

# rearrange EWA
fc_out <- factor_correspondence(true_scores, ewa)

upperwa <- upperwa %*% fc_out$permutation_matrix
lowerwa <- lowerwa %*% fc_out$permutation_matrix
ewa <- fc_out$rearranged

head(upperwa)
head(ewa)
head(true_scores)
head(lowerwa)

sum(upperwa >= true_scores & lowerwa <= true_scores)/(1000*4)
# This is a bad proportion

# Try it L2 normalized
upperwa_denom <- apply(upperwa, 1, function(x) sqrt(sum(x^2)))
upperwa_norm <- upperwa/upperwa_denom
# apply(upperwa_norm, 1, function(x) sqrt(sum(x^2)))

lowerwa_denom <- apply(lowerwa, 1, function(x) sqrt(sum(x^2)))
lowerwa_norm <- lowerwa/lowerwa_denom
# apply(lowerwa_norm, 1, function(x) sqrt(sum(x^2)))

norm_scores_denom <- apply(true_scores, 1, function(x) sqrt(sum(x^2)))
norm_scores <- true_scores/norm_scores_denom
# apply(norm_scores, 1, function(x) sqrt(sum(x^2)))

head(upperwa_norm)
head(norm_scores)
head(lowerwa_norm)

sum(upperwa_norm >= norm_scores & lowerwa_norm <= norm_scores)/(1000*4)

# Try it L1 normalized 
# The variance of X/n is equal to the variance of X divided by n^2, 
# or (np(1-p))/n² = (p(1-p))/n
upperwa_denom1 <- apply(upperwa, 1, function(x) (sum(abs(x))))
upperwa_norm1 <- upperwa/upperwa_denom1
# apply(upperwa_norm1, 1, function(x) sum(x))

lowerwa_denom1 <- apply(lowerwa, 1, function(x) sum(abs(x)))
lowerwa_norm1 <- lowerwa/lowerwa_denom1
# apply(lowerwa_norm1, 1, function(x) sum(abs(x)))

norm_scores_denom1 <- apply(true_scores, 1, function(x) sum(abs(x)))
norm_scores1 <- true_scores/norm_scores_denom1
# apply(norm_scores1, 1, function(x) sum(abs(x)))

head(upperwa_norm1)
head(norm_scores1)
head(lowerwa_norm1)

sum(upperwa_norm1 >= norm_scores1 & lowerwa_norm1 <= norm_scores1)/(1000*4)
