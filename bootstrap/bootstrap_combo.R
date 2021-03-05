# Combine bootstrap distributions and VCI on HPC
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(R.matlab, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

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
bs_ewa  <- array(dim = c(1000, 5, bootstraps)) # 5 includes row ID
bs_eh   <- array(dim = c(4, 40, bootstraps)) # 5 includes row ID
bs_pred <- array(dim = c(1000, 41, bootstraps)) # 5 includes row ID

# Read in data output by `bootstrap_bn2mf.m`
for (boot in 1:bootstraps) {
  if (file.exists(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ewa/ewa_bs_sim_", 
                         set, "_bs_", boot, ".mat"))) {
  bs_ewa[,,boot] = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ewa/ewa_bs_sim_", 
                                              set, "_bs_", boot, ".mat"))[[1]]
  bs_eh[,,boot] = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_eh/eh_bs_sim_", 
                                  set, "_bs_", boot, ".mat"))[[1]]
  bs_pred[,,boot] = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_pred/pred_bs_sim_", 
                                  set, "_bs_", boot, ".mat"))[[1]]
  }
  print(paste0("Boot Number: ", boot))
}

# EWA
# first match permutations !!!
# reorder IDs

for (boot in 1:bootstraps) {
# takes the first instance of an ID
# all instances of same ID are equal
  bs_ewa[,,boot] = bs_ewa[,,boot][match(1:1000, bs_ewa[,,boot][,1]),]
  print(paste0("Match Permutations for EWA Boot Number: ", boot))
}

# pred
for (boot in 1:bootstraps) {
  # takes the first instance of an ID
  # all instances of same ID are equal
  bs_pred[,,boot] = bs_pred[,,boot][match(1:1000, bs_pred[,,boot][,1]),]
  print(paste0("Match Permutations for Pred Boot Number: ", boot))
}

# Create empirical confidence interval
# One matrix for each BN2MF output
# 25th, mean, median, and 75th percentile matrices for EWA
bs_lower_wa  <- apply(bs_ewa, c(1,2), quantile, 0.025, na.rm = TRUE)
bs_upper_wa  <- apply(bs_ewa, c(1,2), quantile, 0.975, na.rm = TRUE)
bs_mean_wa   <- apply(bs_ewa, c(1,2), mean,            na.rm = TRUE)

bs_lower_eh  <- apply(bs_eh, c(1,2), quantile, 0.025, na.rm = TRUE)
bs_upper_eh  <- apply(bs_eh, c(1,2), quantile, 0.975, na.rm = TRUE)
bs_mean_eh   <- apply(bs_eh, c(1,2), mean,          na.rm = TRUE)

bs_lower_pred  <- apply(bs_pred, c(1,2), quantile, 0.025, na.rm = TRUE)
bs_upper_pred  <- apply(bs_pred, c(1,2), quantile, 0.975, na.rm = TRUE)
bs_mean_pred   <- apply(bs_pred, c(1,2), mean,          na.rm = TRUE)

save(bs_lower_wa,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_lower_wa_", set, ".RDA"))
save(bs_upper_wa,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_upper_wa_", set, ".RDA"))
save(bs_mean_wa,   file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_mean_wa_", set, ".RDA"))

save(bs_lower_eh,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_lower_nmf_eh_", set, ".RDA"))
save(bs_upper_eh,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_upper_nmf_eh_", set, ".RDA"))
save(bs_mean_eh, file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_mean_nmf_eh_", set, ".RDA"))

save(bs_lower_pred,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_lower_nmf_pred_", set, ".RDA"))
save(bs_upper_pred,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_upper_nmf_pred_", set, ".RDA"))
save(bs_mean_pred, file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/bootstrap_ci/bs_mean_nmf_pred_", set, ".RDA"))
# all `bs_list...` files go into `bootstrap_out.R`
