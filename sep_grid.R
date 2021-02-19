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
save(sep_m, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_out_grid.RDA")
