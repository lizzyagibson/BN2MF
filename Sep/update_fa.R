#### Packages ####
library(tidyverse)
source("./Results/compare_functions.R")

# Read in Sims
load("./Sims/sim_sep_grid.RDA")
sim_sep

# Run FA
sep_out_fa <- sim_sep %>%
  mutate(fa_out       = map(sim, function(x) get_fa(x, 4)),
         fa_loadings  = map(fa_out, function(x) x[[1]]),
         fa_scores    = map(fa_out, function(x) x[[2]]),
         fa_pred      = map(fa_out, function(x) x[[3]]),
         fa_rank      = map(fa_out, function(x) x[[4]]),
         fa_perm      = map2(true_patterns, fa_loadings, get_perm, nn = FALSE), # Rearrange
         fa_loadings  = map2(fa_loadings, fa_perm, get_perm_product),
         fa_scores    = map2(fa_scores, fa_perm, get_perm_product))
