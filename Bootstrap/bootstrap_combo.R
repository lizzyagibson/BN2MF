library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(R.matlab, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

# Read output from matlab
# 100 bootstrap samples 
bootstraps = 100
# for BN2MF runs on
# 100 correlated simulated data sets
datasets = 100

# List with one entry per simulation
bs_list_ewa = list(rep(NA, datasets))
bs_list_eh  = list(rep(NA, datasets))

# List of arrays
# One array for each BN2MF output
# 100 bootstrapped results per array
for (run in 1:datasets) {
  bs_list_ewa[[run]] <- array(dim = c(1000, 4, bootstraps))
  bs_list_eh[[run]]  <- array(dim = c(1000, 4, bootstraps))
}

for (run in 1:datasets) {
  sim_num = run + 200
  for (boot in 1:bootstraps) {
    if (file.exists(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_cor_ewa/cor_ewa_bs_", 
                                   boot, "sim_", sim_num, ".mat"))) {
    bs_list_ewa[[run]] = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_cor_ewa/cor_ewa_bs_", 
                                        boot, "sim_", sim_num, ".mat"))[[1]]
    bs_list_eh[[run]]  = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_cor_eh/cor_eh_bs_",  
                                        boot, "sim_", sim_num, ".mat"))[[1]]
    }
    print(paste0("Bootstrap Number: ", boot))
  }
  print(paste0("Run Number: ", run))
  }

# Create empirical confidence interval
# List of matrices
# One matrix for each BN2MF output
# 25th, mean, and 75th percentile matrices (3 lists)

bs_list_lower = list(rep(NA, datasets))
bs_list_upper = list(rep(NA, datasets))
bs_list_mean  = list(rep(NA, datasets))

for (run in 1:datasets) {
  bs_list_lower[[run]] <- apply(bs_list_ewa[[run]], c(1,2), quantile, 0.025, na.rm = TRUE)
  bs_list_upper[[run]] <- apply(bs_list_ewa[[run]], c(1,2), quantile, 0.975, na.rm = TRUE)
  bs_list_mean[[run]]  <- apply(bs_list_ewa[[run]], c(1,2), mean, na.rm = TRUE)
}

save(bs_list_lower, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_cor/bs_list_lower.RDA")
save(bs_list_upper, file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_cor/bs_list_upper.RDA")
save(bs_list_mean,  file = "/ifs/scratch/msph/ehs/eag2186/npbnmf/bootstrap_cor/bs_list_mean.RDA")



