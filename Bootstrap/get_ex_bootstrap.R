# Run on HPC
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(R.matlab, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

# Run on array 1:3
job_num = as.integer(Sys.getenv("SGE_TASK_ID"))
job_num

# this gets the bootstrapped distributions for three example outputs

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

# Read output from matlab `bootstrap_bn2mf.m`
for (boot in 1:bootstraps) {
  if (file.exists(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_", sim_type, "_ewa/", 
                        sim_type, "_ewa_bs_", boot, "_sim_", sim_num, ".mat"))) {
    bs_ewa[,,boot] = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_", sim_type, "_ewa/", 
                                         sim_type, "_ewa_bs_", boot, "_sim_", sim_num, ".mat"))[[1]]
    bs_eh[,,boot]  = readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_", sim_type, "_eh/", 
                                              sim_type, "_eh_bs_", boot, "_sim_", sim_num, ".mat"))[[1]]
  }
  print(paste0("Bootstrap Number: ", boot))
  str(bs_ewa[,,boot])
  str(bs_eh[,,boot])
}

# EWA
# first match permutations !!!
for (boot in 1:bootstraps) {
  # takes the first instance of an ID
    bs_ewa[,,boot] = bs_ewa[,,boot][match(1:1000, bs_ewa[,,boot][,1]),]
    print(paste0("Bootstrap: ", boot))
  }

# Also save whole distribution for this example
save(bs_ewa,  file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_", sim_type, "_ewa.RDA"))
save(bs_eh,   file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/main/bs/bootstrap_lists/bs_", sim_type, "_eh.RDA"))
