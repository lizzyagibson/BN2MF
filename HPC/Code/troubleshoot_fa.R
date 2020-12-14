#########################################################
# Factor Analysis failed on some simulations
# Solution: change default factoring method to ML (face palm)
#########################################################

# Packages
library(tidyverse)
library(NMF)

## Get functions
source("./Results/compare_functions.R")

# Read in Sims
load(paste0("./Results/Simulations/sim_dim.RDA"))

#####
# Factor Analysis
#####

dim_fa <- sim_over %>%
  mutate(fa_out      = map2(sim,   patterns, get_fa),
         fa_loadings = map(fa_out, function(x) x[[1]]),
         fa_scores   = map(fa_out, function(x) x[[2]]),
         fa_pred     = map(fa_out, function(x) x[[3]]),
         fa_rank     = map(fa_out, function(x) x[[4]]))

dim_fa %>% unnest(fa_rank) %>% 
  dplyr::select(patterns, fa_rank) %>% 
  drop_na(.) %>% 
  distinct(.)

#####
# Symmetric Subspace Distance #
#####

dim_fa <- dim_fa %>% 
  mutate(fa_loadings_ssdist = map2(true_patterns, fa_loadings, symm_subspace_dist),
         fa_scores_ssdist   = map2(true_scores,   fa_scores,   symm_subspace_dist))

dim_fa %>% slice(2019) %>% dplyr::select(9:15) %>% unnest(c(fa_loadings_ssdist, fa_scores_ssdist))

dim_fa <- dim_fa %>% 
  mutate(fa_l2er    = map2(chem, fa_pred,    get_relerror))

dim_fa <- dim_fa %>%
  mutate(fa_perm           = map2(true_patterns,  fa_loadings,    function(x,y) get_perm(x,y,nn=FALSE)),
         fa_loadings_re    = map2(fa_loadings,    fa_perm,        get_product),
         fa_scores_re      = map2(fa_scores,      fa_perm,        get_product))

dim_fa <- dim_fa %>% 
  mutate(fa_load_l2     = map2(true_patterns, fa_loadings_re,    get_relerror),
         fa_score_l2    = map2(true_scores,   fa_scores_re,      get_relerror))

dim_fa <- dim_fa %>% 
  mutate(fa_load_cos     = map2(true_patterns, fa_loadings_re,    cos_dist),
         fa_score_cos    = map2(true_scores,   fa_scores_re,      cos_dist))

dim_fa <- dim_fa %>% 
  mutate(fa_load_cos_v     = map2(true_patterns, fa_loadings_re,    cos_dist_v),
         fa_score_cos_v    = map2(true_scores,   fa_scores_re,      cos_dist_v))

dim_fa <- dim_fa %>% 
  dplyr::select(seed, participants, chemicals, patterns,
                grep("rank", colnames(.)),
                grep("_l2", colnames(.)),
                grep("ssdist", colnames(.)),
                grep("cos", colnames(.)))

save(dim_fa, file = "./Results/Sup Dim/dim_fa.RDA")

merge_fa <- dim_fa %>% dplyr::select(-true_patterns, -true_scores, -chem, -sim, -fa_out, 
                         -fa_loadings, -fa_scores, -fa_pred, -fa_perm,
                         -fa_loadings_re, -fa_scores_re) %>% 
  unnest(c(fa_rank:fa_score_cos))

metrics_to_merge <- metrics %>% dplyr::select(-grep("fa", colnames(.)))

metrics <- full_join(merge_fa, metrics_to_merge, by = c("seed", "participants", "chemicals", "patterns"))


