#########################################################
# Sims, npBNMF, regular NMF (L2 & Poisson), PCA, and FA #
#########################################################
# Lizzy Gibson    # 6/22/2020 ###########################
# Added BNMF      # 9/17/2020 ###########################
# Cleaned up code # 9/22/202  ###########################
#########################################################

# Packages
library(tidyverse)
library(psych)
library(NMF)

job_num = 784

## Get functions
source("./R/compare_functions.R")

# Read in Sims
load(paste0("./Sims/sim_dim/sim_dim", job_num, ".RDA"))

#####
# PCA
#####

output_all <- to_save %>% 
  mutate(pca_out       = map(sim, get_pca),
         pca_loadings = map(pca_out, function(x) x[[1]]),
         pca_scores    = map(pca_out, function(x) x[[2]]),
         pca_pred      = map(pca_out, function(x) x[[3]]),
         pca_rank      = map(pca_out, function(x) x[[4]]))

#####
# Factor Analysis
#####

output_all <- output_all %>% 
  mutate(fa_out       = map(sim, get_fa),
         fa_loadings = map(fa_out, function(x) x[[1]]),
         fa_scores    = map(fa_out, function(x) x[[2]]),
         fa_pred      = map(fa_out, function(x) x[[3]]),
         fa_rank      = map(fa_out, function(x) x[[4]]))

#####
# L2 NMF
#####

output_all <- output_all %>%
  mutate(nmfl2_out      = map(sim, get_nmf_l2),
         nmfl2_loadings = map(nmfl2_out, function(x) x[[1]]),
         nmfl2_scores   = map(nmfl2_out, function(x) x[[2]]),
         nmfl2_pred     = map(nmfl2_out, function(x) x[[3]]),
         nmfl2_rank     = map(nmfl2_out, function(x) x[[4]]))

#####
# Poisson NMF
#####

output_all <- output_all %>% 
  mutate(nmfp_out      = map(sim, get_nmf_p),
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

save(output_all, file = paste0("./HPC/Rout/dim_out/dim_out_", job_num, ".RDA"))
