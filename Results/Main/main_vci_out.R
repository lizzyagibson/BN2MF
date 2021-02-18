# Combine bootstrap results
# Combine VCI results
# Overall characteristics

####  Load packages ####
library(tidyverse)
library(R.matlab)
#library(openssl)

#### For all sims ####

#### Read data ####
load("./Sims/sim_old_unsep.RDA")
# Sims HAVE CHANGED!!

#### Normalize truth ####
sim_dgp = sim_dgp %>% 
  mutate(denom = map(true_patterns, function(x) apply(x,1, sum)),
         patterns_scaled = map2(true_patterns, denom, function(x,y) x/y),
         scores_scaled = map2(true_scores, denom, function(x,y) as.matrix(x) %*% diag(y)),
         id = 1:nrow(.))

#### Read VCI data ####
load("./Results/Main/vci_out.RDA")
vci_dgp <- dgp_vci %>% 
           dplyr::select(vci_h_mean = eh_scaled,
                   vci_wa_mean = ewa_scaled,
                   vci_wa_lower = lowerWA,
                   vci_wa_upper = upperWA,
                   vci_h_lower = lowerH,
                   vci_h_upper = upperH)
rm(dgp_vci)

#### Combine data ####
all_vci = bind_cols(sim_dgp, vci_dgp) %>% 
  mutate(vci_pred_mean = map2(vci_wa_mean, vci_h_mean, function(x,y) as.matrix(x) %*% as.matrix(y)),
         vci_pred_lower = map(vci_pred_mean, qpois, p = 0.025),
         vci_pred_upper = map(vci_pred_mean, qpois, p = 0.975)) 

vci_metrics = all_vci %>% 
  pivot_longer(c(grep("vci_", colnames(.))),
               names_to = c("method", "side", "matrix"),
               names_sep = "_") %>%
  mutate(value = map(value, as.matrix)) %>% 
  pivot_wider(names_from = matrix,
              values_from = value) %>% 
  mutate(truth = pmap(list(side, patterns_scaled, scores_scaled, chem),
                      function(x, y, z, a) if(x == "h") {y} else if (x == "wa") {z} else if (x == "pred") {a}),
         err  = map2(truth, mean, function(x,y)
           if (ncol(x) == ncol(y)) {norm(x-y, "F")/norm(x, "F")} else {NA}),
         iqr  = map2(upper, lower, function(x, y) mean(x-y)),
         prop = pmap(list(truth, lower, upper),
                           function(x,y,z) sum(x >= y & x <= z)/(nrow(x)*ncol(x)) )) %>%
  unnest(c(err, prop, iqr)) 

#### Results Tables ####
vci_metrics %>% 
  group_by(data, method, side) %>%
  summarise(qs = quantile(prop, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(prop),
            max = max(prop),
            min = min(prop)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(side, data, method, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(side, data, method)

vci_metrics %>% 
  group_by(data, method, side) %>% 
  summarise(qs = quantile(err, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(err),
            max = max(err),
            min = min(err)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(side, data, method, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(side, data, method)

vci_metrics %>% 
  group_by(data, method, side) %>% 
  summarise(qs = quantile(iqr, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(iqr),
            max = max(iqr),
            min = min(iqr)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(side, data, method, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(side, data, method)
