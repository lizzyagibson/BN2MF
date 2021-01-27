library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(R.matlab, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

# Read output from matlab
# 150 bootstrap samples 
bootstraps = 150
# for BN2MF runs on


# 150 bootstrapped results per array
bs_over_ewa <- array(dim = c(1000, 5, bootstraps)) # 5 includes row ID
bs_over_eh <- array(dim = c(4, 50, bootstraps))

for (boot in 1:bootstraps) {
  if (file.exists(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_over_ewa/", 
                        "over_ewa_bs_", boot, "_sim_102.mat"))) {
    bs_over_ewa[,,boot] = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_over_ewa/", 
                                              "over_ewa_bs_", boot, "_sim_102.mat"))[[1]]
    bs_over_eh[,,boot]  = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_over_eh/", 
                                              "over_eh_bs_", boot, "_sim_102.mat"))[[1]]
  }
  print(paste0("Bootstrap Number: ", boot))
  str(bs_over_ewa[,,boot])
  str(bs_over_eh[,,boot])
}

# EH
# Create empirical confidence interval
# List of matrices, one for each BN2MF output
# 25th, mean, median, and 75th percentile matrices (4 lists) for EH
bs_over_lower_h  = matrix(NA, nrow = 4, ncol = 50)
bs_over_upper_h  = matrix(NA, nrow = 4, ncol = 50)
bs_over_mean_h   = matrix(NA, nrow = 4, ncol = 50)
bs_over_median_h = matrix(NA, nrow = 4, ncol = 50)

bs_over_lower_h  <- apply(bs_over_eh, c(1,2), quantile, 0.025, na.rm = TRUE)
bs_over_upper_h  <- apply(bs_over_eh, c(1,2), quantile, 0.975, na.rm = TRUE)  
bs_over_mean_h   <- apply(bs_over_eh, c(1,2), mean,            na.rm = TRUE)
bs_over_median_h <- apply(bs_over_eh, c(1,2), median,          na.rm = TRUE)

save(bs_over_lower_h,  file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_over_lower_h.RDA")
save(bs_over_upper_h,  file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_over_upper_h.RDA")
save(bs_over_mean_h,   file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_over_mean_h.RDA")
save(bs_over_median_h, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_over_median_h.RDA")

# EWA
# first match permutations
for (boot in 1:bootstraps) {
  # takes the first instance of an ID
    bs_over_ewa[,,boot] = bs_over_ewa[,,boot][match(1:1000, bs_over_ewa[,,boot][,1]),]
    print(paste0("Bootstrap: ", boot))
    #print(bs_over_ewa[,,boot][,1])
  }

# Create empirical confidence interval
# List of matrices
# One matrix for each BN2MF output
# 25th, mean, median, and 75th percentile matrices (4 lists) for EWA
bs_over_lower_wa  = matrix(NA, nrow = 1000, ncol = 5)
bs_over_upper_wa  = matrix(NA, nrow = 1000, ncol = 5)
bs_over_mean_wa   = matrix(NA, nrow = 1000, ncol = 5)
bs_over_median_wa = matrix(NA, nrow = 1000, ncol = 5)

bs_over_lower_wa  <- apply(bs_over_ewa, c(1,2), quantile, 0.025, na.rm = TRUE)
bs_over_upper_wa  <- apply(bs_over_ewa, c(1,2), quantile, 0.975, na.rm = TRUE)
bs_over_mean_wa   <- apply(bs_over_ewa, c(1,2), mean,            na.rm = TRUE)
bs_over_median_wa <- apply(bs_over_ewa, c(1,2), median,          na.rm = TRUE)

save(bs_over_lower_wa,  file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_over_lower_wa.RDA")
save(bs_over_upper_wa,  file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_over_upper_wa.RDA")
save(bs_over_mean_wa,   file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_over_mean_wa.RDA")
save(bs_over_median_wa, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_over_median_wa.RDA")

# Also save whole over for this example
save(bs_over_ewa,  file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_over_ewa.RDA")
save(bs_over_eh,  file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_lists/bs_over_eh.RDA")
