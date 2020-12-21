# Testing annealed VCI
# on scores

library(R.matlab)
library(tidyverse)

# need factor correspondence to rearrage solution
source("./Results/factor_correspondence.R")

#####
# Low Noise example
# DGP main simulations
#####

# Load simulations aka TRUTH
load("./Results/Main/sim_dgp_rep1.RDA")
sim_dgp_rep1
patterns <- as_tibble(sim_dgp_rep1$true_patterns[[1]])
scores <- as_tibble(sim_dgp_rep1$true_scores[[1]])

# load annealed resuts
t0 <- tibble()
for (i in 1:nrow(sim_dgp_rep1)) {
  
  ewa <- readMat(paste0("./Results/Main/bnmf_dgp_out/dgp_ewa_", i, ".mat"))[[1]]
  ewa <- as_tibble(ewa) %>% nest(ewa = everything())
  
  varwa <- readMat(paste0("./Results/Main/bnmf_dgp_out/dgp_varwa_", i, ".mat"))[[1]]
  varwa <- as_tibble(varwa) %>% nest(varwa = everything())
  
  eh <- readMat(paste0("./Results/Main/bnmf_dgp_out/dgp_eh_", i, ".mat"))[[1]]
  eh <- as_tibble(eh) %>% nest(eh = everything())
  
  alphah <- readMat(paste0("./Results/Main/bnmf_dgp_out/dgp_alphah_", i, ".mat"))[[1]]
  alphah <- as_tibble(alphah) %>% nest(alphah = everything())
  
  betah <- readMat(paste0("./Results/Main/bnmf_dgp_out/dgp_betah_", i, ".mat"))[[1]]
  betah <- as_tibble(betah) %>% nest(betah = everything())
  
  all_t0 = bind_cols(ewa, varwa, eh, alphah, betah)
  t0 = bind_rows(t0, all_t0)
}

# Clean up solution matrix
t0_all <- bind_cols(sim_dgp_rep1, t0)
t0_4 <- t0_all %>% 
  mutate(rank = map(ewa, ncol)) %>% 
  unnest(rank) %>% 
  filter(rank == 4) %>% 
  mutate(ewa = map(ewa, as.matrix),
         varwa = map(varwa, as.matrix),
         eh = map(eh, as.matrix),
         upperwa = map2(ewa, varwa, function(x,y) x + (1.96*sqrt(y))),
         lowerwa= map2(ewa, varwa, function(x,y) x - (1.96*sqrt(y))),
         pred = map2(ewa, eh, function(x,y) x %*% y))

#####
# Repeat across dataframe
#####

t0_4 = t0_4 %>% 
  mutate(fc_out = map2(true_scores, ewa, factor_correspondence),
         upperwa = map2(upperwa, fc_out, function(x,y) x %*% y$permutation_matrix),
         lowerwa = map2(lowerwa, fc_out, function(x,y) x %*% y$permutation_matrix),
         ewa = map(fc_out, function(x) x$rearranged)) # ewa rearranged eh is not!

# Un-normalized
get_prop <- function(upperwa, lowerwa, true_scores){
  sum(upperwa >= true_scores & lowerwa <= true_scores)/(1000*4)
}

t0_4 = t0_4 %>% 
  mutate(prop = pmap(list(upperwa, lowerwa, true_scores), get_prop))

t0_4 %>% 
  dplyr::select(seed, data, prop) %>% 
  unnest(prop) %>% 
  arrange(desc(prop))

# Try it L2 normalized
get_norm <- function(ewa, dat) {
  #denom <- apply(ewa, 1, function(x) sqrt(sum(x^2)))
  denom <- apply(ewa, 1, function(x) (sum(abs(x))))
  dat_norm <- dat/denom
  dat_norm
}

t0_4 = t0_4 %>% 
  mutate(upperwa_norm = map2(ewa, upperwa, get_norm),
         lowerwa_norm = map2(ewa, lowerwa, get_norm),
         norm_scores = map(true_scores, function(x) get_norm(x,x)))

get_prop_norm <- function(upperwa, lowerwa, norm_scores){
  sum(upperwa >= norm_scores & lowerwa <= norm_scores)/(1000*4)
}

t0_4 = t0_4 %>% 
  mutate(prop_norm = pmap(list(upperwa_norm, lowerwa_norm, norm_scores), get_prop_norm))

# save(t0_4, file = "./Results/Main/dgp_vci.RDA")
load("./Results/Main/dgp_vci.RDA")

t0_4 %>% 
  dplyr::select(seed, data, prop, prop_norm) %>% 
  unnest(c(prop, prop_norm)) %>% 
  group_by(data) %>% 
  summarise(qs = quantile(prop_norm, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(prop_norm),
            max = max(prop_norm),
            min = min(prop_norm)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, min, `0.25`, `0.5`, mean, `0.75`, max)

prop <- tibble()

for (i in 1:300) {
  add_prop = readMat(paste0("./Results/Main/prop_out/save_prop", i, ".mat"))[[1]] %>% as_tibble()  
  add_prop
  prop = bind_rows(prop, add_prop)
}
prop = prop %>% rename(row = V1, prop = V2) %>% 
  mutate(seed = rep(1:100, 3),
         type = c(rep("Distinct", 100), rep("Overlapping", 100), rep("Correlated", 100)))

prop %>% 
  group_by(type) %>% 
  summarise(qs = quantile(prop, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(prop),
            max = max(prop),
            min = min(prop)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(type, min, `0.25`, `0.5`, mean, `0.75`, max)

## Relative error
t0_relerror = t0_4 %>% 
  dplyr::select(seed, chem, sim, data, pred) %>% 
  mutate(l2_true = map2(chem, pred, function(x,y) norm(x-y, "F")/norm(x, "F")),
         l2_sim  = map2(sim,  pred, function(x,y) norm(x-y, "F")/norm(x, "F"))) %>% 
  unnest(c(l2_sim, l2_true))

t0_relerror %>%
  mutate(data = fct_relevel(data, "Distinct", "Overlapping", "Correlated")) %>% 
  ggplot(aes(x = data, y = l2_true)) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  scale_y_log10() +
  theme(legend.position = "none") + 
  labs(y = "BN2MF Relative Predictive Error")

