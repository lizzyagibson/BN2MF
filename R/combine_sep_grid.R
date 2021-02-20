#### Packages ####
library(tidyverse)
source("./Results/compare_functions.R")
library(furrr)
plan(strategy = multiprocess)

#### Load Data ####

# BN2MF output
load("./Sep/sep_out_grid.RDA")
sep_m

# Simulations
load("./Sims/sim_sep_grid.RDA")
sim_sep

# R output of comparison models
load("./Sep/combined_models.RDA")
sep_r

#### Normalize truth ####
sim_sep = sim_sep %>% 
      mutate(denom = map(true_patterns, function(x) apply(x, 1, sum)),
         true_patternsscaled = map2(true_patterns, denom, function(x,y) x/y),
         true_scoresscaled = map2(true_scores, denom, function(x,y) as.matrix(x) %*% diag(y)),
         id = 1:nrow(.)) %>% 
  mutate_at(vars(1:3), as.factor)

sep_bn2mf = bind_cols(sim_sep, sep_m)

#### VCI ####
vci = sep_bn2mf %>% 
  unnest(bn2mf_rank) %>% 
  filter(bn2mf_rank == 4) %>% 
  mutate_at(vars(12:18), ~map(., as.matrix)) %>% 
  mutate(iqr  = map2(upperWA, lowerWA, function(x, y) mean(x-y)),
         prop = pmap(list(true_scoresscaled, lowerWA, upperWA),
                     function(x,y,z) sum(x >= y & x <= z)/(nrow(x)*ncol(x)) )) %>%
  unnest(c(prop, iqr)) 

prop_table = vci %>%
              group_by(sep_num, noise_level) %>% 
              summarise(qs = quantile(prop, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75)) %>% 
              pivot_wider(names_from = "prob",
                          values_from = "qs") %>% 
              mutate_if(is.numeric, round, 2) %>% 
              rename(median = `0.5`)

prop_table %>% 
  ggplot(aes(x = sep_num, y = noise_level, fill = median)) +
  geom_tile() +
  geom_text(aes(label = round(median, 2)), size = 3.5, col = "coral") + 
  scale_x_discrete(limits = rev) +
  labs(x = "Number of distinct chemicals per pattern",
       y = "Noise level (as proportion of true SD)",
       fill = "Median coverage") +
  theme_minimal() +
  theme(legend.position = "bottom") + 
  scale_fill_distiller(palette = "YlGnBu", direction = 1)

#### ERROR ####

sep_all = sep_r %>% 
          mutate_at(vars(1:3), as.factor) %>% 
          full_join(., sep_bn2mf) %>% 
   dplyr::select(!grep("(scaled|upper|lower|denom)", colnames(.))) %>% 
   dplyr::select(true_scores, id, everything())
sep_all

#### Rank ####
sep_rank = sep_all %>% 
  dplyr::select(seed, sep_num, noise_level, grep("rank", colnames(.))) %>% 
  pivot_longer(pca_rank:bn2mf_rank,
               names_to = c("model", "drop"),
               names_sep = "_",
               values_to = "rank") %>% 
  unnest(rank) %>% 
  mutate(rank = ifelse(is.na(rank), 0, rank),
         rank_bin = ifelse(rank ==4, "right", "wrong"))

sep_rank %>% 
  group_by(model, rank_bin) %>% 
  summarise(n = n())

##### Metrics ####
metrics = sep_all %>%
  dplyr::select(!grep("rank", colnames(.))) %>% 
    pivot_longer(pca_pred:bn2mf_pred,
                 names_to = c("model", "matrix"),
                 names_sep = "_") %>% 
  mutate(truth = case_when(matrix == "loadings" ~ true_patterns,
                           matrix == "scores" ~ true_scores,
                           matrix == "pred" ~ chem),
         relerr  = map2(truth, value, get_relerror),
         ssdist  = map2(truth, value, symm_subspace_dist),
         cosdist = map2(truth, value, cos_dist)) %>% 
  unnest(c(relerr, ssdist, cosdist)) %>% 
  dplyr::select(-true_scores, -id, -true_patterns, -chem, -sim, -value, -truth)

save(metrics, file = "./Sep/sep_grid_metrics.RDA")

metrics %>% 
  ggplot(aes(x = model, y = relerr, fill = model)) +
  geom_boxplot() +
  facet_wrap(.~matrix, scales = "free_y")

metrics_sum = metrics %>%
  group_by(sep_num, noise_level, model, matrix) %>% 
  summarize(relerr = mean(relerr, na.rm=T),
          ssdist = mean(ssdist, na.rm=T),
          cosdist = mean(cosdist, na.rm=T))
  
metrics_sum %>%   
  filter(model == "nmfp" & matrix == "scores") %>% 
  ggplot(aes(x = sep_num, y = noise_level, fill = relerr)) +
  geom_tile() +
  geom_text(aes(label = round(relerr, 2)), size = 3.5, col = "coral") + 
  scale_x_discrete(limits = rev) +
  theme_minimal() +
  labs(x = "Number of distinct chemicals per pattern",
       y = "Noise level (as proportion of true SD)") +
  theme(legend.position = "bottom") + 
  scale_fill_distiller(palette = "YlGnBu", direction = -1)


