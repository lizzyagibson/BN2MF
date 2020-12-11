#########################################################
# Sims, npBNMF, regular NMF (L2 & Poisson), PCA, and FA #
#########################################################
# Lizzy Gibson    # 6/22/2020 ###########################
# Added BNMF      # 9/17/2020 ###########################
# Cleaned up code # 9/22/202  ###########################
#########################################################

## read job number from system environment
## This only works if run on cluster!!
job_num = as.integer(Sys.getenv("SGE_TASK_ID"))
job_num

# Packages
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(registry, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(pkgmaker, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(rngtools, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(NMF, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(GPArotation, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(psych, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

## Get functions
source("/ifs/scratch/msph/ehs/eag2186/npbnmf/compare_functions.R")

# Read in Sims
load(paste0("/ifs/scratch/msph/ehs/eag2186/Data/sim_noise.RDA"))

# Run everything

#####
# PCA
#####

output_all <- sim_noise[job_num,] %>% 
  mutate(pca_out       = map(sim, get_pca),
         pca_loadings = map(pca_out, function(x) x[[1]]),
         pca_scores    = map(pca_out, function(x) x[[2]]),
         pca_pred      = map(pca_out, function(x) x[[3]]),
         pca_rank      = map(pca_out, function(x) x[[4]]))

#####
# Factor Analysis
#####

output_all <- output_all %>% 
  mutate(fa_out       = map(sim, function(x) get_fa(x,4)),
         fa_loadings  = map(fa_out, function(x) x[[1]]),
         fa_scores    = map(fa_out, function(x) x[[2]]),
         fa_pred      = map(fa_out, function(x) x[[3]]),
         fa_rank      = map(fa_out, function(x) x[[4]]))

#####
# L2 NMF
#####

output_all <- output_all %>%
  mutate(nmfl2_out      = map(sim, function(x) get_nmf_l2(x,4)),
         nmfl2_loadings = map(nmfl2_out, function(x) x[[1]]),
         nmfl2_scores   = map(nmfl2_out, function(x) x[[2]]),
         nmfl2_pred     = map(nmfl2_out, function(x) x[[3]]),
         nmfl2_rank     = map(nmfl2_out, function(x) x[[4]]))

#####
# Poisson NMF
#####

output_all <- output_all %>% 
  mutate(nmfp_out      = map(sim, function(x) get_nmf_p(x,4)),
         nmfp_loadings = map(nmfp_out, function(x) x[[1]]),
         nmfp_scores   = map(nmfp_out, function(x) x[[2]]),
         nmfp_pred     = map(nmfp_out, function(x) x[[3]]),
         nmfp_rank     = map(nmfp_out, function(x) x[[4]]))

#####
# Symmetric Subspace Distance #
#####

output_all <- output_all %>% 
  mutate(pca_loadings_ssdist  = map2(true_patterns, pca_loadings,   symm_subspace_dist),
         pca_scores_ssdist    = map2(true_scores,   pca_scores,     symm_subspace_dist),
         fa_loadings_ssdist   = map2(true_patterns, fa_loadings,    symm_subspace_dist),
         fa_scores_ssdist     = map2(true_scores,   fa_scores,      symm_subspace_dist),
         nmfl2_loading_ssdist = map2(true_patterns, nmfl2_loadings, symm_subspace_dist),
         nmfl2_scores_ssdist  = map2(true_scores,   nmfl2_scores,   symm_subspace_dist),
         nmfp_loading_ssdist  = map2(true_patterns, nmfp_loadings,  symm_subspace_dist),
         nmfp_scores_ssdist   = map2(true_scores,   nmfp_scores,    symm_subspace_dist))

#####
# Save
#####
                           
output_all <- output_all %>% dplyr::select(-grep("_out", colnames(.)))
job_num
save(output_all, file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/noise_out/noise_out_", job_num, ".RDA"))

