
# Combine VCI results
# Overall characteristics
# For sparse BN2MF versions

####  Load packages ####
library(tidyverse)
source("./Results/compare_functions.R")

#### Read sim data ####
load("./Sims/Main/sim_dgp.RDA")

#### Normalize truth ####
sim_dgp = sim_dgp %>% 
  mutate(denom = map(true_patterns, function(x) apply(x,1, sum)),
         patterns_scaled = map2(true_patterns, denom, function(x,y) x/y),
         scores_scaled = map2(true_scores, denom, function(x,y) as.matrix(x) %*% diag(y)),
         id = 1:nrow(.))

#### Read VCI data ####
load("./sparse_h/sparse_vci_out.RDA")
dgp_vci <- dgp_vci %>% 
  dplyr::select(vci_h_mean = eh_scaled,
                vci_wa_mean = ewa_scaled,
                vci_wa_lower = lowerWA,
                vci_wa_upper = upperWA,
                vci_h_lower = lowerH,
                vci_h_upper = upperH) %>% 
  mutate(version = "PUSH")

sparse = bind_cols(sim_dgp,dgp_vci)

#### Combine data ####
vci = sparse %>% 
      pivot_longer(grep("vci_", colnames(.)),
                   names_to = c("method", "side", "matrix"),
                   names_sep = "_") %>%
      dplyr::select(-method) %>% 
      mutate(value = map(value, as.matrix)) %>% 
      pivot_wider(names_from = matrix,
                  values_from = value) %>% 
      mutate(truth = pmap(list(side, patterns_scaled, scores_scaled),
                          function(x, y, z) if(x == "h") {y} else {z}),
             err  = map2(truth, mean, get_relerror),
             iqr  = map2(upper, lower, function(x, y) mean(x-y)),
             prop = pmap(list(truth, lower, upper), get_prop)) %>%
      unnest(c(err, prop, iqr))
vci

#### Results Tables ####
vci %>% 
  group_by(data, version, side) %>%
  summarise(qs = quantile(prop, c(0.25, 0.5, 0.75), na.rm=T), prob = c(0.25, 0.5, 0.75),
            mean = mean(prop, na.rm=T),
            max = max(prop, na.rm=T),
            min = min(prop, na.rm=T)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(side, data, version, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(side, data, version)

vci %>% 
  group_by(data, version, side) %>%
  summarise(qs = quantile(err, c(0.25, 0.5, 0.75), na.rm=T), prob = c(0.25, 0.5, 0.75),
            mean = mean(err, na.rm=T),
            max = max(err, na.rm=T),
            min = min(err, na.rm=T)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(side, data, version, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(side, data, version)

vci %>% 
  group_by(data, version, side) %>%
  summarise(qs = quantile(iqr, c(0.25, 0.5, 0.75), na.rm=T), prob = c(0.25, 0.5, 0.75),
            mean = mean(iqr, na.rm=T),
            max = max(iqr, na.rm=T),
            min = min(iqr, na.rm=T)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(side, data, version, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(side, data, version)

# OK this is obviously a rounding issue
vci$lower[[1]]

# This is all scaled
vci %>% 
  filter(side == "wa") %>% 
  select(-grep("(scal|mean|lower|upper|tru|m)", colnames(.))) %>%
  group_by(data, version, side) %>% 
  summarise(qs = quantile(err, c(0.25, 0.5, 0.75), na.rm=T), prob = c(0.25, 0.5, 0.75),
            mean = mean(err, na.rm=T),
            max = max(err, na.rm=T),
            min = min(err, na.rm=T)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, side, version, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(data, version)
  
## EX
y = sparse$vci_h_lower[[201]]
x = sparse$patterns_scaled[[201]]
z = sparse$vci_h_upper[[201]]

y[,1:5]
x[,1:5]
sparse$vci_h_mean[[201]][,1:5]
z[,1:5]

sum(x >= y & x <= z)/(nrow(x)*ncol(x))  

sort(sparse$true_patterns[[201]][1,])
sum(x[1,])
sort(x[1,])
sort(x[1,]) == 0
sort(sparse$vci_h_mean[[201]][1,]) == 0
