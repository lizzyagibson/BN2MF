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

#### Normalize truth ####
sim_dgp = sim_dgp %>% 
  mutate(denom = map(true_patterns, function(x) apply(x,1, sum)),
         patterns_scaled = map2(true_patterns, denom, function(x,y) x/y),
         scores_scaled = map2(true_scores, denom, function(x,y) as.matrix(x) %*% diag(y)),
         id = 1:nrow(.))

#### Read VCI data ####
load("./Sep/vci_out.RDA")
load("./Sep/sep_vci_out.RDA")
dgp_vci
sep_vci

dgp = bind_cols(sim_dgp, dgp_vci) %>% mutate(reps = 0)
sep = bind_cols(sims, sep_vci)

full = full_join(dgp, sep)

#### Combine data ####
all_bs_vci = bind_cols(sim_dgp, bs_dgp, vci_dgp, bs_pred) %>% 
  mutate(vci_pred_mean = map2(vci_wa_mean, vci_h_mean, function(x,y) as.matrix(x) %*% as.matrix(y)),
         vci_pred_lower = map(vci_pred_mean, qpois, p = 0.025),
         vci_pred_upper = map(vci_pred_mean, qpois, p = 0.975)) 

bs_vci_metrics = all_bs_vci %>% 
  # For VCI predicted values, take WA*H = mean of Poisson
  # Take 2.5 and 97.5 percentiles of Poisson(WA*H) with inverse CDF
  # Scaled vs non-scaled give same predicted value
  # For bootstrapped predicted values
  # Take predicted as WA*H of each bootstrap
  # Take 2.5 and 97.5 percentiles of bootstrapped distribution (over 150 bootstraps)
  pivot_longer(c(grep("bs_", colnames(.)), grep("vci_", colnames(.))),
               names_to = c("method", "side", "matrix"),
               names_sep = "_") %>%
  mutate(value = map(value, as.matrix)) %>% 
  pivot_wider(names_from = matrix,
              values_from = value) %>% 
  mutate(truth = pmap(list(side, patterns_scaled, scores_scaled, chem),
                      function(x, y, z, a) if(x == "h") {y} else if (x == "wa") {z} else if (x == "pred") {a}),
         best = map2(median, mean, function(x,y) if(!is.null(x)) {x} else {y}),
         err  = map2(truth, best, function(x,y)
           if (ncol(x) == ncol(y)) {norm(x-y, "F")/norm(x, "F")} else {NA}),
         iqr  = map2(upper, lower, function(x, y) mean(x-y)),
         prop = pmap(list(truth, lower, upper),
                           function(x,y,z) sum(x >= y & x <= z)/(nrow(x)*ncol(x)) )) %>%
  unnest(c(err, prop, iqr)) 

#### Results Tables ####
bs_vci_metrics %>% 
  group_by(data, method, side) %>%
  summarise(qs = quantile(prop, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(prop),
            max = max(prop),
            min = min(prop)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(side, data, method, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(side, data, method)

bs_vci_metrics %>% 
  group_by(data, method, side) %>% 
  summarise(qs = quantile(err, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(err),
            max = max(err),
            min = min(err)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(side, data, method, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(side, data, method)

bs_vci_metrics %>% 
  group_by(data, method, side) %>% 
  summarise(qs = quantile(iqr, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(iqr),
            max = max(iqr),
            min = min(iqr)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(side, data, method, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(side, data, method)

#### save data ####
all_bs_vci
#save(all_bs_vci, file = "./Bootstrap/all_bs_vci.RDA")

bs_vci_metrics
#save(bs_vci_metrics, file = "./Bootstrap/bs_vci_metrics.RDA")

all_bs_vci$bs_h_lower[[201]][,1:7]
all_bs_vci$patterns_scaled[[201]][,1:7]
all_bs_vci$bs_h_upper[[201]][,1:7]
