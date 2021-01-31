# Combine bootstrap distributions and VCI on HPC
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(R.matlab, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

# First combine bootstraps

# 150 bootstrap samples 
bootstraps = 150

# 300 simulated data sets
datasets = 300

# List with one entry per simulation
bs_list_ewa = list(rep(NA, datasets))
bs_list_eh  = list(rep(NA, datasets))

# One array for each BN2MF output
for (run in 1:datasets) {
  bs_list_ewa[[run]] <- array(dim = c(1000, 5, bootstraps)) # 5 includes row ID
  bs_list_eh[[run]]  <- array(dim = c(4, 50, bootstraps))
}

for (run in 1:datasets) {
  
  type = ""
  if (run <= 100) {type = "dist"
  } else if (run <= 200 & run > 100) {
    type = "over"} else if (run > 200) {type = "cor"}
  
  for (boot in 1:bootstraps) {
    if (file.exists(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_", type, "_ewa/", 
                           type, "_ewa_bs_", boot, "_sim_", run, ".mat"))) {
    bs_list_ewa[[run]][,,boot] = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_", type, "_ewa/", 
                                                type, "_ewa_bs_", boot, "_sim_", run, ".mat"))[[1]]
    bs_list_eh[[run]][,,boot]  = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_", type, "_eh/", 
                                                type, "_eh_bs_", boot, "_sim_", run, ".mat"))[[1]]
    }
    }
  print(paste0("Run Number: ", run))
  }

# EH
# Create empirical confidence interval
# List of matrices, one for each BN2MF output
# 25th, mean, median, and 75th percentile matrices (4 lists) for EH
bs_list_lower_h  = list(rep(NA, datasets))
bs_list_upper_h  = list(rep(NA, datasets))
bs_list_mean_h   = list(rep(NA, datasets))
bs_list_median_h = list(rep(NA, datasets))

for (run in 1:datasets) {
  bs_list_lower_h[[run]]  <- apply(bs_list_eh[[run]], c(1,2), quantile, 0.025, na.rm = TRUE)
  bs_list_upper_h[[run]]  <- apply(bs_list_eh[[run]], c(1,2), quantile, 0.975, na.rm = TRUE)  
  bs_list_mean_h[[run]]   <- apply(bs_list_eh[[run]], c(1,2), mean, na.rm = TRUE)
  bs_list_median_h[[run]] <- apply(bs_list_eh[[run]], c(1,2), median, na.rm = TRUE)
}

save(bs_list_lower_h,  file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_lower_h.RDA")
save(bs_list_upper_h,  file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_upper_h.RDA")
save(bs_list_mean_h,   file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_mean_h.RDA")
save(bs_list_median_h, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_median_h.RDA")

# EWA
# first match permutations
# reorder IDs
for (run in 1:datasets) {
  for (boot in 1:bootstraps) {
  # takes the first instance of an ID
  # all instances of same ID are equal
    bs_list_ewa[[run]][,,boot] = bs_list_ewa[[run]][,,boot][match(1:1000, bs_list_ewa[[run]][,,boot][,1]),]
    print(paste0("Run Number: ", run, "Bootstrap: ", boot))
  }
}

# Create empirical confidence interval
# List of matrices
# One matrix for each BN2MF output
# 25th, mean, median, and 75th percentile matrices (4 lists) for EWA
bs_list_lower_wa  = list(rep(NA, datasets))
bs_list_upper_wa  = list(rep(NA, datasets))
bs_list_mean_wa   = list(rep(NA, datasets))
bs_list_median_wa = list(rep(NA, datasets))

for (run in 1:datasets) {
  bs_list_lower_wa[[run]]  <- apply(bs_list_ewa[[run]], c(1,2), quantile, 0.025, na.rm = TRUE)
  bs_list_upper_wa[[run]]  <- apply(bs_list_ewa[[run]], c(1,2), quantile, 0.975, na.rm = TRUE)
  bs_list_mean_wa[[run]]   <- apply(bs_list_ewa[[run]], c(1,2), mean,            na.rm = TRUE)
  bs_list_median_wa[[run]] <- apply(bs_list_ewa[[run]], c(1,2), median,          na.rm = TRUE)
}

save(bs_list_lower_wa,  file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_lower_wa.RDA")
save(bs_list_upper_wa,  file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_upper_wa.RDA")
save(bs_list_mean_wa,   file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_mean_wa.RDA")
save(bs_list_median_wa, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_median_wa.RDA")

