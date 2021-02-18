
## Packages
library(tidyverse)
library(R.matlab)

## Get functions
source("./Results/compare_functions.R")

sep_m1 = tibble() # to combine with sep_r
sep_m2 = tibble() # to combine with bootsraps

how_many = 50
for (j in 1:how_many) {
  
  bn2mf_loadings = readMat(paste0("./trouble_out/q_eh_NOTscaled", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data = everything()) %>% rename(bn2mf_loadings = data)

  bn2mf_scores = readMat(paste0("./trouble_out/q_ewa_NOTscaled", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(bn2mf_scores = data)

  eh_scaled = readMat(paste0("./trouble_out/q_eh_scaled", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data = everything()) %>% rename(eh_scaled = data)

  ewa_scaled = readMat(paste0("./trouble_out/q_ewa_scaled", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(ewa_scaled = data)

  upperWA = readMat(paste0("./trouble_out/q_upperWA_", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(upper_WA = data)

  lowerWA = readMat(paste0("./trouble_out/q_lowerWA_", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(lower_WA = data)

  upperH = readMat(paste0("./trouble_out/q_upperH_", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(upper_H = data)

  lowerH = readMat(paste0("./trouble_out/q_lowerH_", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(lower_H = data)

  m_out1   <- bind_cols(bn2mf_loadings, bn2mf_scores,eh_scaled, ewa_scaled, lowerWA, upperWA, lowerH, upperH)
  
  sep_m1   <- bind_rows(sep_m1, m_out1)
  print(j)
}
sep_m1 = sep_m1 %>% mutate(noise = "1")

for (j in 1:how_many) {
  
  bn2mf_loadings = readMat(paste0("./trouble_out/q_eh_NOTscaled", j, "_0.5.mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data = everything()) %>% rename(bn2mf_loadings = data)
  
  bn2mf_scores = readMat(paste0("./trouble_out/q_ewa_NOTscaled", j, "_0.5.mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(bn2mf_scores = data)
  
  eh_scaled = readMat(paste0("./trouble_out/q_eh_scaled", j, "_0.5.mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data = everything()) %>% rename(eh_scaled = data)
  
  ewa_scaled = readMat(paste0("./trouble_out/q_ewa_scaled", j, "_0.5.mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(ewa_scaled = data)
  
  upperWA = readMat(paste0("./trouble_out/q_upperWA_", j, "_0.5.mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(upper_WA = data)
  
  lowerWA = readMat(paste0("./trouble_out/q_lowerWA_", j, "_0.5.mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(lower_WA = data)
  
  upperH = readMat(paste0("./trouble_out/q_upperH_", j, "_0.5.mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(upper_H = data)
  
  lowerH = readMat(paste0("./trouble_out/q_lowerH_", j, "_0.5.mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(lower_H = data)
  
  m_out2   <- bind_cols(bn2mf_loadings, bn2mf_scores,eh_scaled, ewa_scaled, lowerWA, upperWA, lowerH, upperH)
  
  sep_m2   <- bind_rows(sep_m2, m_out2)
  print(j)
}
sep_m2 = sep_m2 %>% mutate(noise = "0.5")

sep_m <- bind_rows(sep_m1, sep_m2) %>%
         mutate(bn2mf_loadings = map(bn2mf_loadings, as.matrix),
                bn2mf_scores = map(bn2mf_scores, as.matrix),
                bn2mf_pred = map2(bn2mf_scores, bn2mf_loadings, function(x,y) x %*% y),
                bn2mf_rank = map(bn2mf_scores, ncol))
sep_m

load("./Sims/sim_trouble.RDA")
sim_dist1 = sim_dist

load("./Sims/sim_trouble_0.5.RDA")
sim_dist2 = sim_dist

#### Normalize truth ####
dist = bind_rows(sim_dist1, sim_dist2) %>% 
  mutate(denom = map(true_patterns, function(x) apply(x, 1, sum)),
         patterns_scaled = map2(true_patterns, denom, function(x,y) x/y),
         scores_scaled = map2(true_scores, denom, function(x,y) as.matrix(x) %*% diag(y)),
         id = 1:nrow(.))

sep = bind_cols(dist, sep_m)
sep

#### Tables ####
sep %>%
  dplyr::select(type, noise, bn2mf_rank) %>% unnest(bn2mf_rank) %>% table()
  
metrics = sep %>% 
  unnest(bn2mf_rank) %>% 
  filter(bn2mf_rank == 4) %>% 
  mutate_at(vars(11:18), ~map(., as.matrix)) %>% 
  mutate(iqr  = map2(upper_WA, lower_WA, function(x, y) mean(x-y)),
         prop = pmap(list(scores_scaled, lower_WA, upper_WA),
                     function(x,y,z) sum(x >= y & x <= z)/(nrow(x)*ncol(x)) )) %>%
  unnest(c(prop, iqr)) 

#### Results ####
summary(metrics$iqr)
summary(metrics$prop)

metrics %>%
  group_by(type, noise) %>% 
  summarise(qs = quantile(prop, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  unnest(type) %>% 
  mutate_if(is.numeric, round, 2)
