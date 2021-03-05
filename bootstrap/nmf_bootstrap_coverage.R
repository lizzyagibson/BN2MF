# Get bootstrap coverage

####  Load packages ####
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

#### For 600 sims ####

# 600 simulated data sets
datasets = read_csv("/ifs/scratch/msph/ehs/eag2186/Data/bs_ids.csv")

nmf_bs_metrics = tibble()

for (i in 1:nrow(datasets)) {

  set = datasets$value[i]
  
  # BS 95% CI
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_lower_nmf_pred_", set, ".RDA"))
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_upper_nmf_pred_", set, ".RDA"))
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_mean_nmf_pred_", set, ".RDA"))
  
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_lower_nmf_score_", set, ".RDA"))
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_upper_nmf_score_", set, ".RDA"))
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_mean_nmf_score_", set, ".RDA"))
  
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_lower_nmf_loadings_", set, ".RDA"))
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_upper_nmf_loadings_", set, ".RDA"))
  load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_mean_nmf_loadings_", set, ".RDA"))

  bs_lower_pred = bs_lower_pred[,2:41]
  bs_upper_pred = bs_upper_pred[,2:41]
  bs_mean_pred  = bs_mean_pred[,2:41]
  
  bs_lower_score = bs_lower_score[,2:5]
  bs_upper_score = bs_upper_score[,2:5]
  bs_mean_score  = bs_mean_score[,2:5]
  
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
  
  prop_scores   = sum(scores_scaled >= bs_lower_score & scores_scaled <= bs_upper_score)/(nrow(scores_scaled)*ncol(scores_scaled))
  prop_load    = sum(patterns_scaled >= bs_lower_load & patterns_scaled <= bs_upper_load)/(nrow(patterns_scaled)*ncol(patterns_scaled))
  prop_pred = sum(chem >= bs_lower_pred & chem <= bs_upper_pred)/(nrow(chem)*ncol(chem))
  
  width_scores = mean(bs_upper_score - bs_lower_score)
  width_load  = mean(bs_upper_load - bs_lower_load)
  width_pred = mean(bs_upper_pred - bs_lower_pred)
  
  err_pred = norm(chem - bs_mean_pred, "F")/norm(chem, "F")
  
  err_scores_scaled = norm(scores_scaled - bs_mean_score, "F")/norm(scores_scaled, "F")
  err_load_scaled  = norm(patterns_scaled - as.matrix(bs_mean_load), "F")/norm(patterns_scaled, "F")
  
  row = tibble(set, prop_scores, prop_load, prop_pred, 
               err_scores_scaled, err_load_scaled, err_pred,
               width_scores, width_load, width_pred)
  
  nmf_bs_metrics = bind_rows(nmf_bs_metrics, row)
  print(paste0("Sim Number: ", i))
}

# Save prop
save(nmf_bs_metrics, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_bs_metrics.RDA")

