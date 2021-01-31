# Run on HPC
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(R.matlab, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

# Run on array 1:3
job_num = as.integer(Sys.getenv("SGE_TASK_ID"))
job_num

# once each for distinct, overlapping, and correlated
if (job_num == 1) {
  sim_num = 2
  sim_type = "dist"
} else if (job_num == 2) {
  sim_num = 102
  sim_type = "over"
} else {
  sim_num = 202
  sim_type = "cor"}

# 150 bootstrap samples 
bootstraps = 150

# 150 bootstrapped results per array
bs_ewa <- array(dim = c(1000, 5, bootstraps)) # 5 includes row ID
bs_eh <- array(dim = c(4, 50, bootstraps))

# Read output from matlab
for (boot in 1:bootstraps) {
  if (file.exists(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_", sim_type, "_ewa/", 
                        sim_type, "_ewa_bs_", boot, "_sim_", sim_num, ".mat"))) {
    bs_ewa[,,boot] = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_", sim_type, "_ewa/", 
                                         sim_type, "_ewa_bs_", boot, "_sim_", sim_num, ".mat"))[[1]]
    bs_eh[,,boot]  = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_", sim_type, "_eh/", 
                                              sim_type, "_eh_bs_", boot, "_sim_", sim_num, ".mat"))[[1]]
  }
  print(paste0("Bootstrap Number: ", boot))
  str(bs_ewa[,,boot])
  str(bs_eh[,,boot])
}

# EH
# Create empirical confidence interval
# List of matrices, one for each BN2MF output
# 25th, mean, median, and 75th percentile matrices (4 lists) for EH
bs_lower_h  = matrix(NA, nrow = 4, ncol = 50)
bs_upper_h  = matrix(NA, nrow = 4, ncol = 50)
bs_mean_h   = matrix(NA, nrow = 4, ncol = 50)
bs_median_h = matrix(NA, nrow = 4, ncol = 50)

bs_lower_h  <- apply(bs__eh, c(1,2), quantile, 0.025, na.rm = TRUE)
bs_upper_h  <- apply(bs__eh, c(1,2), quantile, 0.975, na.rm = TRUE)  
bs_mean_h   <- apply(bs__eh, c(1,2), mean,            na.rm = TRUE)
bs_median_h <- apply(bs__eh, c(1,2), median,          na.rm = TRUE)

save(bs_lower_h,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_", sim_type, "_lower_h.RDA"))
save(bs_upper_h,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_", sim_type, "_upper_h.RDA"))
save(bs_mean_h,   file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_", sim_type, "_mean_h.RDA"))
save(bs_median_h, file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_", sim_type, "_median_h.RDA"))

# EWA
# first match permutations
for (boot in 1:bootstraps) {
  # takes the first instance of an ID
    bs_ewa[,,boot] = bs_ewa[,,boot][match(1:1000, bs_ewa[,,boot][,1]),]
    print(paste0("Bootstrap: ", boot))
  }

# Create empirical confidence interval
# List of matrices
# One matrix for each BN2MF output
# 25th, mean, median, and 75th percentile matrices (4 lists) for EWA
bs_lower_wa  = matrix(NA, nrow = 1000, ncol = 5)
bs_upper_wa  = matrix(NA, nrow = 1000, ncol = 5)
bs_mean_wa   = matrix(NA, nrow = 1000, ncol = 5)
bs_median_wa = matrix(NA, nrow = 1000, ncol = 5)

bs_lower_wa  <- apply(bs_ewa, c(1,2), quantile, 0.025, na.rm = TRUE)
bs_upper_wa  <- apply(bs_ewa, c(1,2), quantile, 0.975, na.rm = TRUE)
bs_mean_wa   <- apply(bs_ewa, c(1,2), mean,            na.rm = TRUE)
bs_median_wa <- apply(bs_ewa, c(1,2), median,          na.rm = TRUE)

save(bs_lower_wa,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_", sim_type, "_lower_wa.RDA"))
save(bs_upper_wa,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_", sim_type, "_upper_wa.RDA"))
save(bs_mean_wa,   file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_", sim_type, "_mean_wa.RDA"))
save(bs_median_wa, file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_", sim_type, "_median_wa.RDA"))

# Also save whole over for this example
save(bs_ewa,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_", sim_type, "_ewa.RDA"))
save(bs_eh,   file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_", sim_type, "_eh.RDA"))
