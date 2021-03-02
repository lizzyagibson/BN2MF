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
bs_scores <- array(dim = c(1000, 5, bootstraps)) # 5 includes row ID

# Read in data output by `bootstrap_nmf.m`
for (boot in 1:bootstraps) {
  if (file.exists(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf/nmf_bs_sim_", 
                         set, "_bs_", boot, ".RDA"))) {
    load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf/nmf_bs_sim_", 
                                              set, "_bs_", boot, ".RDA"))
    bs_scores[,,boot] = as.matrix(nmf_boot$scores_final[[1]])
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
  print(paste0("Match Permutations for Boot Number: ", boot))
}

# Create empirical confidence interval
# One matrix for each BN2MF output
# 25th, mean, median, and 75th percentile matrices for EWA
bs_lower  <- apply(bs_scores, c(1,2), quantile, 0.025, na.rm = TRUE)
bs_upper  <- apply(bs_scores, c(1,2), quantile, 0.975, na.rm = TRUE)
#bs_mean   <- apply(bs_scores, c(1,2), mean,            na.rm = TRUE)
#bs_median <- apply(bs_scores, c(1,2), median,          na.rm = TRUE)

save(bs_lower,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_lower_nmf_", set, ".RDA"))
save(bs_upper,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_upper_nmf_", set, ".RDA"))
#save(bs_mean,   file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_mean_nmf_", set, ".RDA"))
#save(bs_median, file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/bs/nmf_ci/bs_median_nmf_", set, ".RDA"))
# all `bs_list...` files go into `nmf_bootstrap_out.R`
