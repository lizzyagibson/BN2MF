# Combine bootstrap distributions and VCI on HPC

## Get functions
source("/ifs/scratch/msph/ehs/eag2186/npbnmf/compare_functions.R")

# First combine bootstraps

# Run 1:600
job_num = as.integer(Sys.getenv("SGE_TASK_ID"))
job_num

# 150 bootstrap samples 
bootstraps = 150

# 600 simulated data sets
datasets = read_csv("/ifs/scratch/msph/ehs/eag2186/Data/bs_ids.csv")
set = datasets$value[job_num]

# One array for all bootstraps 
bs_scores <- array(dim = c(1000, 5, bootstraps)) # 5 includes row ID
bs_loadings <- array(dim = c(4, 40, bootstraps)) # 5 includes row ID
bs_pred <- array(dim = c(1000, 41, bootstraps)) # 5 includes row ID

# Read in data output by `bootstrap_nmf.m`
for (boot in 1:bootstraps) {
  if (file.exists(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf/nmf_bs_sim_", 
                         set, "_bs_", boot, ".RDA"))) {
    load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf/nmf_bs_sim_", 
                                              set, "_bs_", boot, ".RDA"))
    nmf_boot = nmf_boot %>% rename(id = 1, bootstrap = 2, loadings_final = 3, scores_sam = 4, pred_sam = 5)
    bs_scores[,,boot] = as.matrix(nmf_boot$scores_sam[[1]])
    bs_loadings[,,boot] = as.matrix(nmf_boot$loadings_final[[1]])
    bs_pred[,,boot] = as.matrix(nmf_boot$pred_sam[[1]])
    
  }
  print(paste0("Boot Number: ", boot))
}

# Scores
# first match permutations !!!
# reorder IDs

for (boot in 1:bootstraps) {
# takes the first instance of an ID
# all instances of same ID are equal
  bs_scores[,,boot] = bs_scores[,,boot][match(1:1000, bs_scores[,,boot][,1]),]
  print(paste0("Match Score Permutations for Boot Number: ", boot))
}

# Pred
for (boot in 1:bootstraps) {
  # takes the first instance of an ID
  # all instances of same ID are equal
  bs_pred[,,boot] = bs_pred[,,boot][match(1:1000, bs_pred[,,boot][,1]),]
  print(paste0("Match Pred Permutations for Boot Number: ", boot))
}

# Create empirical confidence interval
# One matrix for each BN2MF output
# 25th, mean, mean, and 75th percentile matrices for EWA
bs_lower_score  <- apply(bs_scores, c(1,2), quantile, 0.025, na.rm = TRUE)
bs_upper_score  <- apply(bs_scores, c(1,2), quantile, 0.975, na.rm = TRUE)
bs_mean_score   <- apply(bs_scores, c(1,2), mean,          na.rm = TRUE)

bs_lower_load  <- apply(bs_loadings, c(1,2), quantile, 0.025, na.rm = TRUE)
bs_upper_load  <- apply(bs_loadings, c(1,2), quantile, 0.975, na.rm = TRUE)
bs_mean_load   <- apply(bs_loadings, c(1,2), mean,          na.rm = TRUE)

bs_lower_pred  <- apply(bs_pred, c(1,2), quantile, 0.025, na.rm = TRUE)
bs_upper_pred  <- apply(bs_pred, c(1,2), quantile, 0.975, na.rm = TRUE)
bs_mean_pred   <- apply(bs_pred, c(1,2), mean,          na.rm = TRUE)

save(bs_lower_score,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_lower_nmf_score_", set, ".RDA"))
save(bs_upper_score,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_upper_nmf_score_", set, ".RDA"))
save(bs_mean_score, file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_mean_nmf_score_", set, ".RDA"))

save(bs_lower_load,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_lower_nmf_loadings_", set, ".RDA"))
save(bs_upper_load,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_upper_nmf_loadings_", set, ".RDA"))
save(bs_mean_load, file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_mean_nmf_loadings_", set, ".RDA"))

save(bs_lower_pred,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_lower_nmf_pred_", set, ".RDA"))
save(bs_upper_pred,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_upper_nmf_pred_", set, ".RDA"))
save(bs_mean_pred, file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_mean_nmf_pred_", set, ".RDA"))

# all `bs_list...` files go into `nmf_bootstrap_out.R`
