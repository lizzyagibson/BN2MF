#########################################################
# Sims, npBNMF, regular NMF (L2 & Poisson), PCA, and FA #
#########################################################
# Lizzy Gibson    # 6/22/2020 ###########################
# Added BNMF      # 9/17/2020 ###########################
# Cleaned up code # 9/22/2020  ##########################
#########################################################

# Packages
library(tidyverse)
library(psych)
library(NMF)

# These jobs didn't finish on HPC
to_do = c(19  , 46  , 49  , 73  , 100 , 127 , 130 ,  145  , 154  , 172  , 181 ,  199  , 208  , 211 , 235
          , 253 , 262 , 280 , 289 , 292 , 307 , 316 ,  343  , 361  , 370  , 388 ,  397  , 415  , 424 , 427
          , 451 , 478 , 496 , 505 , 532 , 559 , 562 ,  577  , 586  , 589  , 604 ,  613  , 631  , 640 , 667
          , 685 , 694 , 721 , 724 , 739 , 748 , 775 ,  784  , 793  , 802  , 805 ,  829  , 856  , 859 , 883
          , 886 , 910 , 937 , 940 , 964 , 967 , 991 , 1018 , 1036 ,  1045 , 1072,  1075,  1090,  1099, 1102
          , 1117, 1126, 1129, 1153, 1180, 1198, 1207,  1210,  1234,  1252 , 1261,  1288,  1306,  1315, 1333
          , 1342, 1345, 1360, 1369, 1372, 1387, 1396,  1399,  1423,  1441 , 1450,  1468,  1477,  1495, 1504
          , 1522, 1531, 1549, 1558, 1576, 1585, 1603,  1612,  1615,  1630 , 1639,  1666,  1693,  1720, 1723
          , 1747, 1765, 1774, 1777, 1792, 1801, 1828,  1831,  1846,  1855 , 1858,  1882,  1900,  1909, 1927
          , 1936, 1963, 1981, 1990, 2008, 2017, 2020,  2035,  2044,  2062 , 2071,  2089,  2098,  2125, 2152
          , 2179, 2182, 2206, 2233, 2251, 2260, 2287,  2305,  2314,  2332 , 2341,  2359,  2368,  2386, 2395
          , 2398, 2413, 2422, 2449, 2476, 2494, 2503,  2530,  2533,  2548 , 2557,  2584,  2611,  2629, 2638
          , 2641, 2656, 2665, 2683, 2692)

## Get functions
source("./R/compare_functions.R")

# Read in Sims
load(paste0("./Sims/sim_dim.RDA"))

test_fa = sim_over %>% 
  slice(to_do)

#####
# Factor Analysis
#####

output_all <- 
  test_fa %>% 
  mutate(fa_out       = map(sim, function(x) try(get_fa(x))),
         fa_loadings  = map(fa_out, function(x) if(class(x) == "try-error") {NA} else{x[[1]]}), 
         fa_scores    = map(fa_out, function(x) if(class(x) == "try-error") {NA} else{x[[2]]}),
         fa_pred      = map(fa_out, function(x) if(class(x) == "try-error") {NA} else{x[[3]]}),
         fa_rank      = map(fa_out, function(x) if(class(x) == "try-error") {NA} else{x[[4]]}))

#####
# Symmetric Subspace Distance #
#####

output_all <- output_all %>% 
  mutate(fa_loadings_ssdist   = map2(true_patterns, fa_loadings,    
                                     function(x,y) if(is.na(y)) {NA} else{symm_subspace_dist(x,y)}),
         fa_scores_ssdist     = map2(true_scores,   fa_scores, 
                                     function(x,y) if(is.na(y)) {NA} else{symm_subspace_dist(x,y)}))

unlist(output_all$fa_rank)

output_all %>% slice(53)
to_do[53]
