
## Packages
library(tidyverse)
library(R.matlab)

## Get functions
source("./Results/compare_functions.R")

sep_m = tibble() # to combine with sep_r
sep_vci = tibble() # to combine with bootsraps

how_many = 18
for (j in 1:how_many) {

  bn2mf_loadings = readMat(paste0("./Sep/sep_eh_NOTscaled", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data = everything()) %>% rename(bn2mf_loadings = data)

  bn2mf_scores = readMat(paste0("./Sep/sep_ewa_NOTscaled", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(bn2mf_scores = data)

  eh_scaled = readMat(paste0("./Sep/sep_eh_scaled", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data = everything()) %>% rename(eh_scaled = data)

  ewa_scaled = readMat(paste0("./Sep/sep_ewa_scaled", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(ewa_scaled = data)

  upperWA = readMat(paste0("./Sep/sep_upperWA_", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(upperWA = data)

  lowerWA = readMat(paste0("./Sep/sep_lowerWA_", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(lowerWA = data)

  upperH = readMat(paste0("./Sep/sep_upperH_", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(upperH = data)

  lowerH = readMat(paste0("./Sep/sep_lowerH_", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(lowerH = data)

  m_out   <- bind_cols(bn2mf_loadings, bn2mf_scores)
  vci_out <- bind_cols(eh_scaled, ewa_scaled, lowerWA, upperWA, lowerH, upperH)
  
  sep_m   <- bind_rows(sep_m, m_out)
  sep_vci <- bind_rows(sep_vci, vci_out)
  print(j)
}

sep_m <- sep_m %>%
         mutate(bn2mf_loadings = map(bn2mf_loadings, as.matrix),
                bn2mf_scores = map(bn2mf_scores, as.matrix),
                bn2mf_pred = map2(bn2mf_scores, bn2mf_loadings, function(x,y) x %*% y),
                bn2mf_rank = map(bn2mf_scores, ncol))

#load("./Sims/sim_sep.RDA")

sep = bind_cols(sims[1:how_many,], sep_m)
sep

# Calculate error metrics
sep_metrics = sep %>% 
              pivot_longer(bn2mf_loadings:bn2mf_rank,
               names_to = c("model", "object"),
               names_sep = "_") %>% 
              filter(!(object %in% c("perm", "out")))
sep_metrics

# L2 relative error -- pred vs truth
rel_err_all <- sep_metrics %>%
               filter(object == "pred") %>% 
               mutate(rel_err_all = map2(chem, value, get_relerror)) %>% 
               dplyr::select(seed, rel_err_all) %>% 
               unnest(rel_err_all)

# l2 relative error -- loadings
rel_err_loadings <- sep_metrics %>%
                filter(object == "loadings") %>% 
                mutate(rel_err_loadings = map2(true_patterns, value, get_relerror)) %>% 
                dplyr::select(seed, rel_err_loadings) %>% 
                unnest(rel_err_loadings)

# l2 relative error -- scores
rel_err_scores <- sep_metrics %>%
                  filter(object == "scores") %>% 
                  mutate(rel_err_scores = map2(true_scores, value, get_relerror)) %>% 
                  dplyr::select(seed, rel_err_scores) %>% 
                  unnest(rel_err_scores)

# cosine distance -- loadings
cos_dist_loadings <- sep_metrics %>%
                     filter(object == "loadings") %>% 
                     mutate(cos_dist_loadings = map2(true_patterns, value, cos_dist),
                            cos_dist_v_loadings = map2(true_patterns, value, cos_dist_v)) %>% 
                     dplyr::select(seed, cos_dist_loadings, cos_dist_v_loadings) %>% 
                     unnest(c(cos_dist_loadings))

# cosine distance -- scores
cos_dist_scores <- sep_metrics %>%
                   filter(object == "scores") %>% 
                   mutate(cos_dist_scores = map2(true_scores, value, cos_dist),
                         cos_dist_v_scores  = map2(true_scores, value, cos_dist_v)) %>% 
                   dplyr::select(seed, cos_dist_scores, cos_dist_v_scores) %>% 
                   unnest(c(cos_dist_scores))

# symmetric subspace distance -- loadings
ssd_loadings <- sep_metrics %>%
                filter(object == "loadings") %>% 
                mutate(ssd_loadings = map2(true_patterns, value, symm_subspace_dist)) %>% 
                dplyr::select(seed, ssd_loadings) %>% 
                unnest(ssd_loadings)

# symmetric subspace distance -- scores
ssd_scores <- sep_metrics %>%
              filter(object == "scores") %>% 
              mutate(ssd_scores = map2(true_scores, value, symm_subspace_dist)) %>% 
              dplyr::select(seed, ssd_scores) %>% 
              unnest(ssd_scores)

# rank
rank <- sep_metrics %>%
        filter(object == "rank") %>% 
        dplyr::select(seed, rank = value) %>% 
        unnest(rank)

# Combine all metrics
metrics <- rank %>% 
           full_join(., rel_err_all) %>% 
           full_join(., rel_err_loadings) %>% 
           full_join(., rel_err_scores) %>% 
           full_join(., cos_dist_loadings) %>% 
           full_join(., cos_dist_scores) %>% 
           full_join(., ssd_loadings) %>% 
           full_join(., ssd_scores)

#### Tables ####
metrics %>%
  dplyr::select(rank) %>% table()
  
##### Relative Preditive Error ####
# L2 Norm (Truth - Predicted) / L2 Norm (Truth)
metrics %>%
  dplyr::select(seed, rel_err_all) %>%
  summarise(qs = quantile(rel_err_all, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75),
            mean = mean(rel_err_all)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(`0.25`, `0.5`, mean, `0.75`) %>% 
  mutate_if(is.numeric, round, 2)

##### Relative error on loadings and scores #####
metrics %>%
  dplyr::select(seed, rel_err_loadings, rel_err_scores) %>%
  pivot_longer(c(rel_err_loadings, rel_err_scores)) %>% 
  mutate(name = str_to_title(str_remove(name, "rel_err_"))) %>% 
  group_by(name) %>% 
  summarise(qs = quantile(value, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  arrange(name) %>%
  mutate_if(is.numeric, round, 2)

##### Subspace Distance ####
# Distance between linear subspaces (orthonormal bases)
# Loading and scores
metrics %>%
  dplyr::select(seed, ssd_loadings, ssd_scores) %>%
  pivot_longer(c(ssd_loadings, ssd_scores)) %>% 
  mutate(name = str_to_title(str_remove(name, "ssd_"))) %>% 
  group_by(name) %>% 
  summarise(qs = quantile(value, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  arrange(name) %>% 
  mutate_if(is.numeric, round, 2)

#### VCI ####
sim_norm = sims %>% 
  mutate(denom = map(true_patterns, function(x) apply(x,1, sum)),
         patterns_scaled = map2(true_patterns, denom, function(x,y) x/y),
         scores_scaled = map2(true_scores, denom, function(x,y) as.matrix(x) %*% diag(y)),
         id = 1:nrow(.))

vci = bind_cols(sim_norm[1:how_many,], sep_vci) %>% 
  mutate(vci_pred_mean = map2(ewa_scaled, eh_scaled, function(x,y) as.matrix(x) %*% as.matrix(y)),
         vci_pred_lower = map(vci_pred_mean, qpois, p = 0.025),
         vci_pred_upper = map(vci_pred_mean, qpois, p = 0.975)) 

vci_metrics = vci %>% 
  mutate_all(~map(., as.matrix)) %>% 
  mutate(iqr  = map2(upperWA, lowerWA, function(x, y) mean(x-y)),
         prop = pmap(list(scores_scaled, lowerWA, upperWA),
                     function(x,y,z) sum(x >= y & x <= z)/(nrow(x)*ncol(x)) )) %>%
  unnest(c(prop, iqr)) 

#### Results ####

summary(vci_metrics$iqr)

summary(vci_metrics$prop)

