# Combine bootstrap distributions and VCI on HPC
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(R.matlab, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

# First combine bootstraps

# Run one time!

# 150 bootstrap samples 
bootstraps = 150

# 300 simulated data sets
datasets = 300

# List with one entry per simulation
bs_list_ewa  = list()
bs_list_eh   = list()
bs_list_pred = list()

# One array for each BN2MF output
for (set in 1:datasets) {
  bs_list_ewa[[set]] <- array(dim = c(1000, 5, bootstraps)) # 5 includes row ID
  bs_list_eh[[set]]  <- array(dim = c(4, 50, bootstraps))
  bs_list_pred[[set]] <- array(dim = c(1000, 50, bootstraps))
}

# Read in data output by `bootstrap_bn2mf.m`
for (set in 1:datasets) {
  
  type = ""
  if (set <= 100) {type = "dist"
  } else if (set <= 200 & set > 100) {
    type = "over"} else if (set > 200) {type = "cor"}
  
  for (boot in 1:bootstraps) {
    if (file.exists(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_", type, "_ewa/", 
                           type, "_ewa_bs_", boot, "_sim_", set, ".mat"))) {
    bs_list_ewa[[set]][,,boot] = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_", type, "_ewa/", 
                                                type, "_ewa_bs_", boot, "_sim_", set, ".mat"))[[1]]
    bs_list_eh[[set]][,,boot]  = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_", type, "_eh/", 
                                                type, "_eh_bs_", boot, "_sim_", set, ".mat"))[[1]]
    bs_list_pred[[set]][,,boot] <- bs_list_ewa[[set]][,2:5,boot] %*% bs_list_eh[[set]][,,boot]
    }
    }
  print(paste0("Reading Dataset Number: ", set))
  }

# EH
# Create empirical confidence interval
# List of matrices, one for each BN2MF output
# 25th, mean, median, and 75th percentile matrices (4 lists) for EH
bs_list_lower_h  = list()
bs_list_upper_h  = list()
bs_list_mean_h   = list()
bs_list_median_h = list()

for (set in 1:datasets) {
  bs_list_lower_h[[set]]  <- apply(bs_list_eh[[set]], c(1,2), quantile, 0.025, na.rm = TRUE)
  bs_list_upper_h[[set]]  <- apply(bs_list_eh[[set]], c(1,2), quantile, 0.975, na.rm = TRUE)  
  bs_list_mean_h[[set]]   <- apply(bs_list_eh[[set]], c(1,2), mean, na.rm = TRUE)
  bs_list_median_h[[set]] <- apply(bs_list_eh[[set]], c(1,2), median, na.rm = TRUE)
  print(paste0("EH Summary Statistics for Dataset Number: ", set))
}

# all `bs_list...` files go into `bootstrap_out.R`
save(bs_list_lower_h,  file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_lower_h.RDA")
save(bs_list_upper_h,  file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_upper_h.RDA")
save(bs_list_mean_h,   file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_mean_h.RDA")
save(bs_list_median_h, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_median_h.RDA")
rm(list = c("bs_list_lower_h", "bs_list_upper_h", "bs_list_mean_h", "bs_list_median_h"))

# Also save whole distribution for three examples
save_dist = c(2, 102, 202)
for (i in 1:length(save_dist)) {
  print(save_dist[i])
  save_eh  = bs_list_eh[[save_dist[i]]]
  save(save_eh, file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_", save_dist[i], "_eh_array.RDA"))
}
# this also goes into `bootstrap_out.R`
rm(bs_list_eh)

# EWA
# first match permutations !!!
# reorder IDs
for (set in 1:datasets) {
  for (boot in 1:bootstraps) {
  # takes the first instance of an ID
  # all instances of same ID are equal
    bs_list_ewa[[set]][,,boot] = bs_list_ewa[[set]][,,boot][match(1:1000, bs_list_ewa[[set]][,,boot][,1]),]
  }
    print(paste0("Match Permutations for Dataset Number: ", set))
}

# Create empirical confidence interval
# List of matrices
# One matrix for each BN2MF output
# 25th, mean, median, and 75th percentile matrices (4 lists) for EWA
bs_list_lower_wa  = list()
bs_list_upper_wa  = list()
bs_list_mean_wa   = list()
bs_list_median_wa = list()

for (set in 1:datasets) {
  bs_list_lower_wa[[set]]  <- apply(bs_list_ewa[[set]], c(1,2), quantile, 0.025, na.rm = TRUE)
  bs_list_upper_wa[[set]]  <- apply(bs_list_ewa[[set]], c(1,2), quantile, 0.975, na.rm = TRUE)
  bs_list_mean_wa[[set]]   <- apply(bs_list_ewa[[set]], c(1,2), mean,            na.rm = TRUE)
  bs_list_median_wa[[set]] <- apply(bs_list_ewa[[set]], c(1,2), median,          na.rm = TRUE)
  print(paste0("EWA Summary Statistics for Dataset Number: ", set))
  }

save(bs_list_lower_wa,  file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_lower_wa.RDA")
save(bs_list_upper_wa,  file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_upper_wa.RDA")
save(bs_list_mean_wa,   file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_mean_wa.RDA")
save(bs_list_median_wa, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_median_wa.RDA")
# all `bs_list...` files go into `bootstrap_out.R`
rm(list = c("bs_list_lower_wa", "bs_list_upper_wa", "bs_list_mean_wa", "bs_list_median_wa"))

# Also save whole distribution for three examples
save_dist = c(2, 102, 202)

for (i in 1:length(save_dist)) {
  print(save_dist[i])
  save_ewa = bs_list_ewa[[save_dist[i]]]
  save(save_ewa,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_", save_dist[i], "_ewa_array.RDA"))
}
# this also goes into `bootstrap_out.R`
rm(bs_list_ewa)

# Get predicted values
# List of matrices
# One matrix for each BN2MF output
# 25th, mean, median, and 75th percentile matrices (4 lists) for predicted values
bs_list_lower_pred  = list()
bs_list_upper_pred  = list()
bs_list_mean_pred   = list()
bs_list_median_pred = list()

for (set in 1:datasets) {
  bs_list_lower_pred[[set]]  <- apply(bs_list_pred[[set]], c(1,2), quantile, 0.025, na.rm = TRUE)
  bs_list_upper_pred[[set]]  <- apply(bs_list_pred[[set]], c(1,2), quantile, 0.975, na.rm = TRUE)
  bs_list_mean_pred[[set]]   <- apply(bs_list_pred[[set]], c(1,2), mean,            na.rm = TRUE)
  bs_list_median_pred[[set]] <- apply(bs_list_pred[[set]], c(1,2), median,          na.rm = TRUE)
  print(paste0("Predicted Value Summary Statistics for Dataset Number: ", set))
}

save(bs_list_lower_pred,  file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_lower_pred.RDA")
save(bs_list_upper_pred,  file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_upper_pred.RDA")
save(bs_list_mean_pred,   file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_mean_pred.RDA")
save(bs_list_median_pred, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_list_median_pred.RDA")
# all `bs_list...` files go into `bootstrap_out.R`
