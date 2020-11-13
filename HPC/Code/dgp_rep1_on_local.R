#########################################################
# Sims, npBNMF, regular NMF (L2 & Poisson), PCA, and FA #
#########################################################
# Lizzy Gibson    # 6/22/2020  ##########################
# Added BNMF      # 9/17/2020  ##########################
# Cleaned up code # 9/22/2020  ##########################
# Reran on local  # 10/24/2020 ##########################
#########################################################

# Packages
library(tidyverse)
library(psych)
library(NMF)
library(R.matlab)

source("./R/compare_functions.R")
source("./R/fig_set.R")

# Read in Sims
load("./Sims/sim_dgp_rep1.RDA")
check_r <- sim_dgp_rep1 %>% 
  dplyr::select(sim) %>% mutate(row = 1:nrow(.),
                                sim = map(sim, as_tibble))
check_r

# Read MATLAB sims
check_mat <- tibble()
for (i in 1:nrow(sim_dgp_rep1)) {

  check_data <- read_csv(paste0("./Sims/dgp_rep1/sim_dgp_rep1_", i, ".csv")) %>% 
    nest(sim = everything()) %>% mutate(row = i)
  check_mat <- rbind(check_mat, check_data)
}
check_mat

all.equal(check_r %>% dplyr::select(sim), check_mat %>% dplyr::select(sim))
# GOOD

#####
# Run everything
#####

#####
# PCA
#####
# dgp_rep1 <- sim_dgp_rep1 %>% 
#   mutate(pca_out       = map(sim, get_pca),
#          pca_loadings = map(pca_out, function(x) x[[1]]),
#          pca_scores    = map(pca_out, function(x) x[[2]]),
#          pca_pred      = map(pca_out, function(x) x[[3]]),
#          pca_rank      = map(pca_out, function(x) x[[4]]))

#####
# Factor Analysis
#####
# dgp_rep1 <- dgp_rep1 %>%
#   mutate(fa_out       = map(sim, get_fa),
#          fa_loadings = map(fa_out, function(x) x[[1]]),
#          fa_scores    = map(fa_out, function(x) x[[2]]),
#          fa_pred      = map(fa_out, function(x) x[[3]]),
#          fa_rank      = map(fa_out, function(x) x[[4]]))

#####
# L2 NMF
#####
# dgp_rep1 <- dgp_rep1 %>%
#   mutate(nmfl2_out      = map(sim, get_nmfl2),
#          nmfl2_loadings = map(nmfl2_out, function(x) x[[1]]),
#          nmfl2_scores   = map(nmfl2_out, function(x) x[[2]]),
#          nmfl2_pred     = map(nmfl2_out, function(x) x[[3]]),
#          nmfl2_rank     = map(nmfl2_out, function(x) x[[4]]))

#####
# Poisson NMF
#####
# dgp_rep1 <- dgp_rep1 %>% 
#   mutate(nmfp_out      = map(sim, get_nmfp),
#          nmfp_loadings = map(nmfp_out, function(x) x[[1]]),
#          nmfp_scores   = map(nmfp_out, function(x) x[[2]]),
#          nmfp_pred     = map(nmfp_out, function(x) x[[3]]),
#          nmfp_rank     = map(nmfp_out, function(x) x[[4]]))

#####
# Symmetric Subspace Distance
#####

# dgp_rep1 <- dgp_rep1 %>% 
#   mutate(pca_rotation_ssdist   = map2(true_patterns, pca_loadings,   symm_subspace_dist),
#          pca_scores_ssdist     = map2(true_scores,   pca_scores,      symm_subspace_dist),
#          fa_loadings_ssdist   = map2(true_patterns, fa_loadings,    symm_subspace_dist),
#          fa_scores_ssdist      = map2(true_scores,   fa_scores,       symm_subspace_dist),
#          nmfl2_loading_ssdist = map2(true_patterns, nmfl2_loadings, symm_subspace_dist),
#          nmfl2_scores_ssdist  = map2(true_scores,   nmfl2_scores,   symm_subspace_dist),
#          nmfp_loading_ssdist  = map2(true_patterns, nmfp_loadings,  symm_subspace_dist),
#          nmfp_scores_ssdist   = map2(true_scores,   nmfp_scores,    symm_subspace_dist))

#####
# Save
#####

# dgp_rep1 <- dgp_rep1 %>% dplyr::select(-grep("_out", colnames(.)))
# save(dgp_rep1, file = "./HPC/Rout/dgp_rep1.RDA")
load("./HPC/Rout/dgp_rep1.RDA")
dgp_rep1

#####

# Add MATLAB
dgp_bnmf_rep1 <- tibble()

for (i in 1:300) {
  eh <- readMat(here::here(paste0("./MATLAB/loop_dgp_rep1/rep1_eh_dist_", i, ".mat")))[[1]]
  eh <- as_tibble(eh) %>% nest(eh = everything()) %>% mutate(row = i)
  
  ewa <- readMat(here::here(paste0("/MATLAB/loop_dgp_rep1/rep1_ewa_dist_", i, ".mat")))[[1]]
  ewa <- as_tibble(ewa) %>% nest(ewa = everything()) %>% mutate(row = i)
  
  both <- full_join(eh, ewa, by = "row")
  dgp_bnmf_rep1 <- rbind(dgp_bnmf_rep1, both)
}

dgp_bnmf_rep1 <- dgp_bnmf_rep1 %>% 
  mutate(seed = rep(1:100, 3),
         data = c(rep("Distinct", 100), rep("Overlapping", 100), rep("Correlated", 100)))
dgp_bnmf_rep1

dgp_rep1_all <- 
  full_join(dgp_rep1, dgp_bnmf_rep1, by = c("seed", "data")) %>% 
  mutate(bnmf_pred            = map2(ewa, eh, function(x,y) as.matrix(x) %*% as.matrix(y)),
         bnmf_loading_ssdist  = map2(true_patterns, eh,  symm_subspace_dist),
         bnmf_scores_ssdist   = map2(true_scores,   ewa,    symm_subspace_dist),
         bnmf_rank = map(eh, nrow))

# dgp_rep1_all <- dgp_rep1_all %>%
#   rename(pca_loadings_ssdist = pca_rotation_ssdist)

# save(dgp_rep1_all, file = "./HPC/Rout/dgp_rep1_all.RDA")
