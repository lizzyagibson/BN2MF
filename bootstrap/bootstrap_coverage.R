# Get bootstrap coverage

####  Load packages ####
library(tidyverse, lib.loc = "/ifs/scratch/msph/ehs/eag2186/rpack")
library(R.matlab, lib.loc = "/ifs/scratch/msph/ehs/eag2186/rpack")

#### For 600 sims ####

# 600 simulated data sets
datasets = read_csv("/ifs/scratch/msph/ehs/eag2186/Data/bs_ids.csv")
bs_metrics = tibble()

for (i in 1:nrow(datasets)) {

  set = datasets$value[i]
  
  # BS 95% CI
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_lower_pred_", set, ".RDA"))
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_upper_pred_", set, ".RDA"))
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_mean_pred_", set, ".RDA"))
  
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_lower_wa_", set, ".RDA"))
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_upper_wa_", set, ".RDA"))
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_mean_wa_", set, ".RDA"))
  
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_lower_eh_", set, ".RDA"))
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_upper_eh_", set, ".RDA"))
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_mean_eh_", set, ".RDA"))
  
  bs_lower_pred = bs_lower_pred[,2:41]
  bs_upper_pred = bs_upper_pred[,2:41]
  bs_mean_pred  = bs_mean_pred[,2:41]
  
  bs_lower_ewa = bs_lower_wa[,2:5]
  bs_upper_ewa = bs_upper_wa[,2:5]
  bs_mean_ewa  = bs_mean_wa[,2:5]
  
  # True ewa
  true_scores = read_csv(paste0("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_scores/scores_sep_", set, ".csv")) %>% as.matrix()
  # True patterns
  true_patterns = read_csv(paste0("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_patterns/patterns_sep_", set, ".csv")) %>% as.matrix()
  # True sims
  chem = read_csv(paste0("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_chem/chem_sep_", set, ".csv")) %>% as.matrix()
  
  #### Normalize truth ####
  denom = apply(true_patterns,1, sum)
  patterns_scaled = true_patterns/denom
  scores_scaled = as.matrix(true_scores) %*% diag(denom)
  
  prop_ewa   = sum(scores_scaled >= bs_lower_ewa & scores_scaled <= bs_upper_ewa)/(nrow(scores_scaled)*ncol(scores_scaled))
  prop_eh    = sum(patterns_scaled >= bs_lower_eh & patterns_scaled <= bs_upper_eh)/(nrow(patterns_scaled)*ncol(patterns_scaled))
  prop_pred = sum(chem >= bs_lower_pred & chem <= bs_upper_pred)/(nrow(chem)*ncol(chem))
  
  width_ewa = mean(bs_upper_ewa - bs_lower_ewa)
  width_eh  = mean(bs_upper_eh - bs_lower_eh)
  width_pred = mean(bs_upper_pred - bs_lower_pred)
  
  err_pred = norm(chem - bs_mean_pred, "F")/norm(chem, "F")
  
  err_ewa_scaled = norm(scores_scaled - bs_mean_ewa, "F")/norm(scores_scaled, "F")
  err_eh_scaled  = norm(patterns_scaled - as.matrix(bs_mean_eh), "F")/norm(patterns_scaled, "F")
  
  row = tibble(set, prop_ewa, prop_eh, prop_pred, 
               err_ewa_scaled, err_eh_scaled, err_pred,
               width_ewa, width_eh, width_pred)
  
  bs_metrics = bind_rows(bs_metrics, row)
  print(paste0("Sim Number: ", i))
}

# Save prop
save(bs_metrics, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bn2mf_bs_metrics.RDA")

