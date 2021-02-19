## Combine BN2MF with VCI, separate sim results in MATLAB
## Combine other methods, main results in R

## Take metrics for each method
  # Symmetric subspace distance
  # Cosine distance
  # Relative error
    # on sims
    # on loadings
    # on scores

# Run one time!

## Packages
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(R.matlab, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

## Get functions
source("/ifs/scratch/msph/ehs/eag2186/npbnmf/compare_functions.R")
  
# Aggregate BN2MF results from `sep_bnmf_vci.m`
sep_m = tibble() # to combine with sep_r

for (j in 1:300) {
  
  if (file.exists(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_eh_NOTscaled", j, ".mat"))) {

  bn2mf_loadings = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_eh_NOTscaled", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data = everything()) %>% rename(bn2mf_loadings = data)

  bn2mf_scores = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_ewa_NOTscaled", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(bn2mf_scores = data)

  eh_scaled = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_eh_scaled", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data = everything()) %>% rename(eh_scaled = data)

  ewa_scaled = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_ewa_scaled", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(ewa_scaled = data)

  upperWA = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_upperWA_", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(upperWA = data)

  lowerWA = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_lowerWA_", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(lowerWA = data)

  upperH = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_upperH_", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(upperH = data)

  lowerH = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/q_lowerH_", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(lowerH = data)

  m_out   <- bind_cols(bn2mf_loadings, bn2mf_scores, eh_scaled, ewa_scaled, lowerWA, upperWA, lowerH, upperH)
  
  sep_m   <- bind_rows(sep_m, m_out)
  }
  print(j)
}

sep_m <- sep_m %>%
         mutate(bn2mf_loadings = map(bn2mf_loadings, as.matrix),
                bn2mf_scores = map(bn2mf_scores, as.matrix),
                bn2mf_pred = map2(bn2mf_scores, bn2mf_loadings, function(x,y) x %*% y),
                bn2mf_rank = map(bn2mf_scores, ncol))

# Save all
save(sep_m, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/q_out.RDA")

# Calculate error metrics
sep_metrics = sep_m %>% 
              pivot_longer(pca_out:bn2mf_rank,
               names_to = c("model", "object"),
               names_sep = "_") %>% 
              filter(!(object %in% c("perm", "out")))

# L2 relative error -- pred vs truth
rel_err_all <- sep_metrics %>%
               filter(object == "pred") %>% 
               mutate(rel_err_all = map2(chem, value, get_relerror)) %>% 
               dplyr::select(seed, data, model, rel_err_all) %>% 
               unnest(rel_err_all)

# l2 relative error -- loadings
rel_err_loadings <- sep_metrics %>%
                filter(object == "loadings") %>% 
                mutate(rel_err_loadings = map2(true_patterns, value, get_relerror)) %>% 
                dplyr::select(seed, data, model, rel_err_loadings) %>% 
                unnest(rel_err_loadings)

# l2 relative error -- scores
rel_err_scores <- sep_metrics %>%
                  filter(object == "scores") %>% 
                  mutate(rel_err_scores = map2(true_scores, value, get_relerror)) %>% 
                  dplyr::select(seed, data, model, rel_err_scores) %>% 
                  unnest(rel_err_scores)

# cosine distance -- loadings
cos_dist_loadings <- sep_metrics %>%
                     filter(object == "loadings") %>% 
                     mutate(cos_dist_loadings = map2(true_patterns, value, cos_dist),
                            cos_dist_v_loadings = map2(true_patterns, value, cos_dist_v)) %>% 
                     dplyr::select(seed, data, model, cos_dist_loadings, cos_dist_v_loadings) %>% 
                     unnest(c(cos_dist_loadings))

# cosine distance -- scores
cos_dist_scores <- sep_metrics %>%
                   filter(object == "scores") %>% 
                   mutate(cos_dist_scores = map2(true_scores, value, cos_dist),
                         cos_dist_v_scores  = map2(true_scores, value, cos_dist_v)) %>% 
                   dplyr::select(seed, data, model, cos_dist_scores, cos_dist_v_scores) %>% 
                   unnest(c(cos_dist_scores))

# symmetric subspace distance -- loadings
ssd_loadings <- sep_metrics %>%
                filter(object == "loadings") %>% 
                mutate(ssd_loadings = map2(true_patterns, value, symm_subspace_dist)) %>% 
                dplyr::select(seed, data, model, ssd_loadings) %>% 
                unnest(ssd_loadings)

# symmetric subspace distance -- scores
ssd_scores <- sep_metrics %>%
              filter(object == "scores") %>% 
              mutate(ssd_scores = map2(true_scores, value, symm_subspace_dist)) %>% 
              dplyr::select(seed, data, model, ssd_scores) %>% 
              unnest(ssd_scores)

# rank
rank <- sep_metrics %>%
        filter(object == "rank") %>% 
        dplyr::select(seed, data, model, rank = value) %>% 
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

save(metrics, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/q_metrics.RDA")
