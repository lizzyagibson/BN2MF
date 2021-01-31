## Combine BN2MF with VCI, main results in MATLAB
## Combine other methods, main results in R

## Take metrics for each method
  # Symmetric subspace distance
  # Cosine distance
  # Relative error
    # on sims
    # on loadings
    # on scores

## Packages
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(R.matlab, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

## Get functions
source("/ifs/scratch/msph/ehs/eag2186/npbnmf/compare_functions.R")

## Combine R tibbles from `dgp_other_models.R`
## Model output from PCA, FA, and NMF
dgp_r = tibble()
for (i in 1:300) {
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_other_models_out/dgp_main_", i, ".RDA"))
  dgp_r = bind_rows(dgp_r, dgp_out)
  print(i)
}

dgp_r = dgp_r %>% # Replace original with reordered version if it exists
        mutate(pca_loadings   = ifelse(!is.na(pca_loadings_re),   pca_loadings_re,   pca_loadings),
               fa_loadings    = ifelse(!is.na(fa_loadings_re),    fa_loadings_re,    fa_loadings),
               nmfl2_loadings = ifelse(!is.na(nmfl2_loadings_re), nmfl2_loadings_re, nmfl2_loadings),
               nmfp_loadings  = ifelse(!is.na(nmfp_loadings_re),  nmfp_loadings_re,  nmfp_loadings),
               pca_scores   = ifelse(!is.na(pca_scores_re),   pca_scores_re,   pca_scores),
               fa_scores    = ifelse(!is.na(fa_scores_re),    fa_scores_re,    fa_scores),
               nmfl2_scores = ifelse(!is.na(nmfl2_scores_re), nmfl2_scores_re, nmfl2_scores),
               nmfp_scores  = ifelse(!is.na(nmfp_scores_re),  nmfp_scores_re,  nmfp_scores)) %>% 
        dplyr::select(-grep("_re", colnames(.)))
  
# Aggregate BN2MF results
dgp_m = tibble()
dgp_vci = tibble()
# Still have arrays for each bootstrap distribution iff needed

for (j in 1:300) {

  bn2mf_loadings = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_vci_out/dgp_eh_NOTscaled", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data = everything()) %>% rename(bn2mf_loadings = data)

  bn2mf_scores = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_vci_out/dgp_ewa_NOTscaled", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(bn2mf_scores = data)

  eh_scaled = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_vci_out/dgp_eh_scaled", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data = everything()) %>% rename(eh_scaled = data)

  ewa_scaled = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_vci_out/dgp_ewa_scaled", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(ewa_scaled = data)

  upperWA = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_vci_out/dgp_upperWA_", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(upperWA = data)

  lowerWA = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_vci_out/dgp_lowerWA_", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(lowerWA = data)

  upperH = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_vci_out/dgp_upperH_", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(upperH = data)

  lowerH = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_vci_out/dgp_lowerH_", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(lowerH = data)

  prop = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/dgp_vci_out/save_prop", j, ".mat"))[[1]] %>% 
    as_tibble(.) %>% nest(data= everything()) %>% rename(prop = data)

  m_out   <- bind_cols(bn2mf_loadings, bn2mf_scores)
  vci_out <- bind_cols(eh_scaled, ewa_scaled, lowerWA, upperWA, lowerH, upperH, prop)
  
  dgp_m   <- bind_rows(dgp_m, m_out)
  dgp_vci <- bind_rows(dgp_vci, m_out)
  print(j)
}

dgp_m <- dgp_m %>%
         mutate(bn2mf_loadings = map(bn2mf_loadings, as.matrix),
                bn2mf_scores = map(bn2mf_scores, as.matrix),
                bn2mf_pred = map2(bn2mf_scores, bn2mf_loadings, function(x,y) x %*% y),
                bn2mf_rank = map(bn2mf_scores, ncol))

dgp = bind_cols(dgp_r, dgp_m)

# Save all
save(dgp, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/main_out.RDA")

# Save VCI
save(dgp_vci, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/vci_out.RDA")

# Calculate error metrics

dgp_metrics = dgp %>% 
              pivot_longer(pca_out:bn2mf_rank,
               names_to = c("model", "object"),
               names_sep = "_") %>% 
              filter(!(object %in% c("perm", "out")))

# L2 relative error -- pred vs truth
rel_err_all <- dgp_metrics %>%
               filter(object == "pred") %>% 
               mutate(rel_err_all = map2(chem, value, get_relerror)) %>% 
               dplyr::select(seed, data, model, rel_err_all) %>% 
               unnest(rel_err_all)

# l2 relative error -- loadings
rel_err_loadings <- dgp_metrics %>%
                filter(object == "loadings") %>% 
                mutate(rel_err_loadings = map2(true_patterns, value, get_relerror)) %>% 
                dplyr::select(seed, data, model, rel_err_loadings) %>% 
                unnest(rel_err_loadings)

# l2 relative error -- scores
rel_err_scores <- dgp_metrics %>%
                  filter(object == "scores") %>% 
                  mutate(rel_err_scores = map2(true_scores, value, get_relerror)) %>% 
                  dplyr::select(seed, data, model, rel_err_scores) %>% 
                  unnest(rel_err_scores)

# cosine distance -- loadings
cos_dist_loadings <- dgp_metrics %>%
                     filter(object == "loadings") %>% 
                     mutate(cos_dist_loadings = map2(true_patterns, value, cos_dist),
                            cos_dist_v_loadings = map2(true_patterns, value, cos_dist_v)) %>% 
                     dplyr::select(seed, data, model, cos_dist_loadings, cos_dist_v_loadings) %>% 
                     unnest(c(cos_dist_loadings))

# cosine distance -- scores
cos_dist_scores <- dgp_metrics %>%
                   filter(object == "scores") %>% 
                   mutate(cos_dist_scores = map2(true_scores, value, cos_dist),
                         cos_dist_v_scores  = map2(true_scores, value, cos_dist_v)) %>% 
                   dplyr::select(seed, data, model, cos_dist_scores, cos_dist_v_scores) %>% 
                   unnest(c(cos_dist_scores))

# symmetric subspace distance -- loadings
ssd_loadings <- dgp_metrics %>%
                filter(object == "loadings") %>% 
                mutate(ssd_loadings = map2(true_patterns, value, symm_subspace_dist)) %>% 
                dplyr::select(seed, data, model, ssd_loadings) %>% 
                unnest(ssd_loadings)

# symmetric subspace distance -- scores
ssd_scores <- dgp_metrics %>%
              filter(object == "scores") %>% 
              mutate(ssd_scores = map2(true_scores, value, symm_subspace_dist)) %>% 
              dplyr::select(seed, data, model, ssd_scores) %>% 
              unnest(ssd_scores)

# rank
rank <- dgp_metrics %>%
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

save(metrics, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/main_metrics.RDA")
# Next plot metrics