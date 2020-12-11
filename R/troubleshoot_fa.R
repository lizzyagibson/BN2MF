#########################################################
# Sims, npBNMF, regular NMF (L2 & Poisson), PCA, and FA #
#########################################################
# Lizzy Gibson    # 6/22/2020 ###########################
# Added BNMF      # 9/17/2020 ###########################
# Cleaned up code # 9/22/2020 ###########################
# Rerun on local  # 10/24/2020
#########################################################

# Packages
library(tidyverse)
library(psych)
library(NMF)
library(R.matlab)

source("./R/compare_functions.R")
source("./R/fig_set.R")

# Read in Sims
load("./Results/Simulations/sim_dgp_rep1.RDA")
load("./Results/Simulations/sim_dim.RDA")
sim_dgp_rep1
sim_over

# I know that sim 19 failed on hpc
sim <- sim_over$sim[[885]]
chem <- sim_over$chem[[19]]

#####
# Factor Analysis
#####
  set.seed(19888891)
  fa_3 <- fa(sim, 3, scores = "regression", rotate = "promax")
  fa_4 <- fa(sim, 4, scores = "regression", rotate = "promax")
  fa_5 <- try(fa(sim, 5, scores = "regression", rotate = "promax"))
  
  fact_3 <- factanal(sim, 3, scores = "regression", rotation = "promax", lower = 0.01)
  fact_4 <- factanal(sim, 4, scores = "regression", rotation = "promax", lower = 0.01)
  fact_5 <- factanal(sim, 5, scores = "regression", rotation = "promax", lower = 0.01)
  
  if(all(class(fa_3) == "try-error")) {fa_3$BIC = NA}
  if(all(class(fa_4) == "try-error")) {fa_4$BIC = NA}
  if(all(class(fa_5) == "try-error")) {fa_5$BIC = NA}
  
  
  # Choose the model with the LOWEST BIC
  if (min(fa_5$BIC, fa_4$BIC, fa_3$BIC, na.rm = TRUE) == fa_5$BIC & !is.na(fa_5$BIC)) {
    fa_out <- fa_5
    rank <- 5
  } else if (min(fa_5$BIC, fa_4$BIC, fa_3$BIC) == fa_4$BIC) {
    fa_out <- fa_4
    rank <- 4
  } else {
    fa_out <- fa_3
    rank <- 3
  }
  
loadings <- matrix(fa_out$loadings, ncol = ncol(fa_out$scores))
scores <- fa_out$scores
pred <- scores %*% t(loadings)

# predict function gives scores
dim(predict.psych(fa_out, sim))

### write ML into function
## repeat with ML and compare results

compare_fa <- dgp_rep1_all %>%
  mutate(fa_outML       = map(sim, function(x) get_fa(x, 4)),
         fa_loadingsML  = map(fa_outML, function(x) x[[1]]),
         fa_scoresML    = map(fa_outML, function(x) x[[2]]),
         fa_predML      = map(fa_outML, function(x) x[[3]]),
         fa_rankML      = map(fa_outML, function(x) x[[4]])) %>% 
  dplyr::select(seed, data, grep("fa", colnames(.)))

compare_fa %>% dplyr::select(grep("rank", colnames(.))) %>% 
  unnest(c(fa_rank, fa_rankML)) %>% 
  distinct()
# all rank 4

compare_fa %>% dplyr::select(grep("pred", colnames(.))) %>% 
  unnest(c(c, fa_predML))

norm(compare_fa$fa_pred[[1]], "F")
norm(compare_fa$fa_predML[[1]], "F")

# replace FA results with ML version
replace_FA <- compare_fa %>% dplyr::select(seed, data, grep("ML", colnames(.)))

dgp_no_fa <- dgp_rep1_all %>% dplyr::select(-(grep("fa_", colnames(.))))

dgp_rep1_all <- full_join(dgp_no_fa, replace_FA, by = c("seed", "data")) %>% 
  rename(fa_loadings = fa_loadingsML, fa_scores = fa_scoresML, 
         fa_pred = fa_predML, fa_rank = fa_rankML) %>% 
  dplyr::select(-fa_outML)
  
dgp_rep1_all <- dgp_rep1_all %>%
  mutate(fa_loadings_ssdist   = map2(true_patterns, fa_loadings,    symm_subspace_dist),
         fa_scores_ssdist      = map2(true_scores,   fa_scores,       symm_subspace_dist))

#save(dgp_rep1_all, file = "./Results/Solutions/dgp_rep1_all.RDA")  
  
  
