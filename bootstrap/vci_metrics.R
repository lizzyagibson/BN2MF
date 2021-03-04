####  Load packages ####
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(R.matlab, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

#### For 600 sims ####

# 600 simulated data sets
datasets = read_csv("/ifs/scratch/msph/ehs/eag2186/Data/bs_ids.csv")

vci_bs_metrics = tibble()

for (i in 571:nrow(datasets)) {
  
  set = datasets$value[i]
  
  # VCI 95% CI
  ewa = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/sep_ewa_NOTscaled", set, ".mat"))[[1]]
  eh = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/sep_eh_NOTscaled",  set, ".mat"))[[1]]
  
  ewa_scaled = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/sep_ewa_scaled", set, ".mat"))[[1]]
  eh_scaled = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/sep_eh_scaled",  set, ".mat"))[[1]]
  
  ewa_upper = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/sep_upperWA_", set, ".mat"))[[1]]
  ewa_lower = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/sep_lowerWA_", set, ".mat"))[[1]]
  
  eh_upper = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/sep_upperH_", set, ".mat"))[[1]]
  eh_lower = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_vci_out/sep_lowerH_", set, ".mat"))[[1]]

  pred = ewa %*% eh
  pred_lower = qpois(pred, p = 0.025)
  pred_upper = qpois(pred, p = 0.975)
  
  # True Scores
  true_scores = read_csv(paste0("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_scores/scores_sep_", set, ".csv")) %>% as.matrix()
  # True patterns
  true_patterns = read_csv(paste0("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_patterns/patterns_sep_", set, ".csv")) %>% as.matrix()
  # True sims
  chem = read_csv(paste0("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_chem/chem_sep_", set, ".csv")) %>% as.matrix()
  
  #### Normalize truth ####
  denom = apply(true_patterns,1, sum)
  patterns_scaled = true_patterns/denom
  scores_scaled = as.matrix(true_scores) %*% diag(denom)
  
  prop_wa   = sum(scores_scaled >= ewa_lower & scores_scaled <= ewa_upper)/(nrow(scores_scaled)*ncol(scores_scaled))
  prop_h    = sum(patterns_scaled >= eh_lower & patterns_scaled <= eh_upper)/(nrow(patterns_scaled)*ncol(patterns_scaled))
  prop_pred = sum(chem >= pred_lower & chem <= pred_upper)/(nrow(chem)*ncol(chem))
  
  width_ewa = mean(ewa_upper - ewa_lower)
  width_eh  = mean(eh_upper - eh_lower)
  width_pred = mean(pred_upper - pred_lower)
  
  err_ewa = norm(true_scores - ewa, "F")/norm(true_scores, "F")
  err_eh = norm(true_patterns - eh, "F")/norm(true_patterns, "F")
  err_pred = norm(chem - pred, "F")/norm(chem, "F")
  
  err_ewa_scaled = norm(scores_scaled - ewa_scaled, "F")/norm(scores_scaled, "F")
  err_eh_scaled  = norm(patterns_scaled - eh_scaled, "F")/norm(patterns_scaled, "F")
  
  row = tibble(set, prop_wa, prop_h, prop_pred, 
                    err_ewa, err_ewa_scaled, err_eh, err_eh_scaled, err_pred,
                    width_ewa, width_eh, width_pred)
  
  vci_bs_metrics = bind_rows(vci_bs_metrics, row)
  print(paste0("Sim Number: ", i))
}

save(vci_bs_metrics, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/vci_bs_metrics.RDA")
