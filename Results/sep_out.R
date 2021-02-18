# Combine bootstrap results
# Combine VCI results
# Overall characteristics

####  Load packages ####
library(tidyverse)
library(R.matlab)
#library(openssl)

#### For all sims ####

#### Read data ####
load("./Sims/sim_dgp.RDA")
load("./Sims/sim_sep.RDA")
sim_dgp
sims

#### Read VCI data ####
load("./Sep/vci_out.RDA")
load("./Sep/sep_vci_out.RDA")
dgp_vci
sep_vci

dgp = bind_cols(sim_dgp, dgp_vci) %>% mutate(reps = 0)
sep = bind_cols(sims, sep_vci)

full = full_join(dgp, sep)

#### Normalize truth ####
full = full %>% 
  mutate(denom = map(true_patterns, function(x) apply(x,1, sum)),
         patterns_scaled = map2(true_patterns, denom, function(x,y) x/y),
         scores_scaled = map2(true_scores, denom, function(x,y) as.matrix(x) %*% diag(y)),
         id = 1:nrow(.))

#### Get Metrics ####
full = full %>% 
  mutate_at(vars(c(2:5, 7:12, 15, 16)), ~map(.,as.matrix)) %>% 
  mutate(err  = map2(scores_scaled, ewa_scaled, function(x,y) norm(x-y, "F")/norm(x, "F")),
         iqr  = map2(upperWA, lowerWA, function(x, y) mean(x-y)),
         prop = pmap(list(scores_scaled, lowerWA, upperWA),
                           function(x,y,z) sum(x >= y & x <= z)/(nrow(x)*ncol(x)) )) %>%
  unnest(c(err, prop, iqr)) 

#### Results Tables ####
full %>% 
  group_by(data, reps) %>%
  summarise(qs = quantile(prop, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(prop),
            max = max(prop),
            min = min(prop)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, reps, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(data, reps)
