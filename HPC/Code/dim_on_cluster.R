#########################################################
# Sims, npBNMF, regular NMF (L2 & Poisson), PCA, and FA #
#########################################################
# Lizzy Gibson    # 6/22/2020 ###########################
# Added BNMF      # 9/17/2020 ###########################
# Cleaned up code # 9/22/202  ###########################
#########################################################

# Packages
library(tidyverse)
# library(psych)
# library(NMF)
library(registry, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(pkgmaker, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(rngtools, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(NMF, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(GPArotation, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(psych, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

## read job number from system environment
## This only works if run on cluster!!
job_num = as.integer(Sys.getenv("SGE_TASK_ID"))
job_num

## Get functions
source("/ifs/scratch/msph/ehs/eag2186/npbnmf/compare_functions.R")

# Read in Sims
load(paste0("/ifs/scratch/msph/ehs/eag2186/Data/sim_dim", , ".RDA"))

# Run everything

#####
# PCA
#####

output_all <- to_save %>% 
  mutate(pca_out       = map(sim, get_pca),
         pca_rotations = map(pca_out, function(x) x[[1]]),
         pca_scores    = map(pca_out, function(x) x[[2]]),
         pca_pred      = map(pca_out, function(x) x[[3]]),
         pca_rank      = map(pca_out, function(x) x[[4]]))

#####
# Factor Analysis
#####

output_all <- output_all %>% 
  mutate(fa_out       = map(sim, get_fa),
         fa_rotations = map(fa_out, function(x) x[[1]]),
         fa_scores    = map(fa_out, function(x) x[[2]]),
         fa_pred      = map(fa_out, function(x) x[[3]]),
         fa_rank      = map(fa_out, function(x) x[[4]]))

#####
# L2 NMF
#####

output_all <- output_all %>%
  mutate(nmf_l2_out      = map(sim, get_nmf_l2),
         nmf_l2_loadings = map(nmf_l2_out, function(x) x[[1]]),
         nmf_l2_scores   = map(nmf_l2_out, function(x) x[[2]]),
         nmf_l2_pred     = map(nmf_l2_out, function(x) x[[3]]),
         nmf_l2_rank     = map(nmf_l2_out, function(x) x[[4]]))

#####
# Poisson NMF
#####

output_all <- output_all %>% 
  mutate(nmf_p_out      = map(sim, get_nmf_p),
         nmf_p_loadings = map(nmf_p_out, function(x) x[[1]]),
         nmf_p_scores   = map(nmf_p_out, function(x) x[[2]]),
         nmf_p_pred     = map(nmf_p_out, function(x) x[[3]]),
         nmf_p_rank     = map(nmf_p_out, function(x) x[[4]]))

#####
# Symmetric Subspace Distance #
#####

output_all <- output_all %>% 
  mutate(pca_rotation_ssdist   = map2(true_patterns, pca_rotations,   symm_subspace_dist),
         pca_scores_ssdist     = map2(true_scores,   pca_scores,      symm_subspace_dist),
         pca_uncenter_rotation_ssdist   = map2(true_patterns, pca_uncenter_rotations,   symm_subspace_dist),
         pca_uncenter_scores_ssdist     = map2(true_scores,   pca_uncenter_scores,      symm_subspace_dist),
         fa_rotations_ssdist   = map2(true_patterns, fa_rotations,    symm_subspace_dist),
         fa_scores_ssdist      = map2(true_scores,   fa_scores,       symm_subspace_dist),
         nmf_l2_loading_ssdist = map2(true_patterns, nmf_l2_loadings, symm_subspace_dist),
         nmf_l2_scores_ssdist  = map2(true_scores,   nmf_l2_scores,   symm_subspace_dist),
         nmf_p_loading_ssdist  = map2(true_patterns, nmf_p_loadings,  symm_subspace_dist),
         nmf_p_scores_ssdist   = map2(true_scores,   nmf_p_scores,    symm_subspace_dist))

#####
# Save
#####
                           
output_all <- output_all %>% dplyr::select(-grep("_out", colnames(.)))

save(output_all, file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/dim_out/dim_out_", job_num, ".RDA"))
