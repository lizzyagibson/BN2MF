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
          , 685 , 694 , 721 , 724 , 739 , 748 , 775 
          ,  784  
          , 793  , 802  , 805 ,  829  , 856  , 859 , 883
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
source("./R/factor_correspondence.R")

# Read in Sims
load(paste0("./Sims/sim_dim.RDA"))

dim_local = sim_over %>% 
  slice(to_do) %>% 
  mutate(id = to_do)

#####
# PCA
#####

output_nofa <- dim_local %>% 
  mutate(pca_out       = map(sim, get_pca)) %>% 
  mutate(pca_loadings  = map(pca_out, function(x) x[[1]]),
         pca_scores    = map(pca_out, function(x) x[[2]]),
         pca_pred      = map(pca_out, function(x) x[[3]]),
         pca_rank      = map(pca_out, function(x) x[[4]]))

#####
# L2 NMF
#####

output_nofa <- output_nofa %>%
 mutate(nmfl2_out      = map2(sim, patterns, get_nmf_l2_dim),
        nmfl2_loadings = map(nmfl2_out, function(x) x[[1]]),
        nmfl2_scores   = map(nmfl2_out, function(x) x[[2]]),
        nmfl2_pred     = map(nmfl2_out, function(x) x[[3]]),
        nmfl2_rank     = map(nmfl2_out, function(x) x[[4]]))

#####
# Poisson NMF
#####

output_nofa <- output_nofa %>% 
 mutate(nmfp_out      = map2(sim, patterns, get_nmf_p_dim),
        nmfp_loadings = map(nmfp_out, function(x) x[[1]]),
        nmfp_scores   = map(nmfp_out, function(x) x[[2]]),
        nmfp_pred     = map(nmfp_out, function(x) x[[3]]),
        nmfp_rank     = map(nmfp_out, function(x) x[[4]]))

#####
# Symmetric Subspace Distance #
#####

output_nofa <- output_nofa %>% 
 mutate(pca_loadings_ssdist  = map2(true_patterns, pca_loadings,   symm_subspace_dist),
        pca_scores_ssdist    = map2(true_scores,   pca_scores,     symm_subspace_dist),
        nmfl2_loading_ssdist = map2(true_patterns, nmfl2_loadings, symm_subspace_dist),
        nmfl2_scores_ssdist  = map2(true_scores,   nmfl2_scores,   symm_subspace_dist),
        nmfp_loading_ssdist  = map2(true_patterns, nmfp_loadings,  symm_subspace_dist),
        nmfp_scores_ssdist   = map2(true_scores,   nmfp_scores,    symm_subspace_dist))

#####
# Save
#####
                          
output_nofa <- output_nofa %>% dplyr::select(-grep("_out", colnames(.)))

output_nofa <- output_nofa %>% 
  mutate(fa_loadings = NA,
         fa_scores = NA,
         fa_pred = NA,
         fa_rank = NA,
         fa_loadings_ssdist = NA,
         fa_scores_ssdist = NA)

# save(output_nofa, file = paste0("./HPC/Rout/dim_out_nofa.RDA"))
load("./HPC/Rout/dim_out/dim_out_nofa.RDA")
output_nofa

#####
# Combine
#####
output_nofa <- output_nofa %>% 
  mutate(pca_l2er  = map2(chem,   pca_pred,   
                          function(x,y) norm(x-y, "F")/norm(x, "F") ),
         fa_l2er   = NA,
         nmfl2_l2er = map2(chem,  nmfl2_pred, 
                           function(x,y) norm(x-y, "F")/norm(x, "F") ),
         nmfp_l2er   = map2(chem, nmfp_pred,
                            function(x,y) norm(x-y, "F")/norm(x, "F") ))

output_nofa <- output_nofa %>%
  mutate(true_patterns   = map(true_patterns, t),
         nmfl2_loadings  = map(nmfl2_loadings, t),
         nmfp_loadings   = map(nmfp_loadings, t),
         true_patterns = map(true_patterns, as.matrix),
         pca_loadings = map(pca_loadings, as.matrix)) %>%
  mutate(pca_perm = map2(true_patterns, pca_loadings,
                         function(x,y) if(ncol(y) == ncol(x))
                         {factor_correspondence(as.matrix(x),
                                                as.matrix(y), nn = FALSE)$permutation_matrix} else{NA}),
         pca_loadings_re = map2(pca_loadings, pca_perm,
                                function(x,y) if(!any(is.na(y)))
                                {x %*% y} else{NA}),
         pca_scores_re = map2(pca_scores, pca_perm,
                              function(x,y) if(!any(is.na(y)))
                              {x %*% y} else{NA}),
         fa_loadings_re = NA,
         fa_scores_re = NA,
         nmfl2_perm = map2(true_patterns, nmfl2_loadings,
                           function(x,y) if(ncol(y) == ncol(x))
                           {factor_correspondence(as.matrix(x),
                                                  as.matrix(y))$permutation_matrix} else{NA}),
         nmfl2_loadings_re = map2(nmfl2_loadings, nmfl2_perm,
                                  function(x,y) if(!any(is.na(y)))
                                  {x %*% y} else{NA}),
         nmfl2_scores_re = map2(nmfl2_scores, nmfl2_perm,
                                function(x,y) if(!any(is.na(y)))
                                {x %*% y} else{NA}),
         nmfp_perm = map2(true_patterns, nmfp_loadings,
                          function(x,y) if(ncol(y) == ncol(x))
                          {factor_correspondence(as.matrix(x),
                                                 as.matrix(y))$permutation_matrix} else{NA}),
         nmfp_loadings_re = map2(nmfp_loadings, nmfp_perm,
                                 function(x,y) if(!any(is.na(y)))
                                 {x %*% y} else{NA}),
         nmfp_scores_re = map2(nmfp_scores, nmfp_perm,
                               function(x,y) if(!any(is.na(y)))
                               {x %*% y} else{NA}))

output_nofa

output_nofa <- 
  output_nofa %>% 
  mutate(pca_load_l2 = map2(true_patterns, pca_loadings_re, 
                            function(x,y) if(!any(is.na(y))) {norm(x-y, "F")/norm(x, "F")} else{NA}),
         pca_score_l2 = map2(true_scores, pca_scores_re, 
                             function(x,y) if(!any(is.na(y))) {norm(x-y, "F")/norm(x, "F")} else{NA}),
         fa_load_l2 = NA,
         fa_score_l2 = NA,
         nmfl2_load_l2 = map2(true_patterns, nmfl2_loadings_re,
                              function(x,y) if(!any(is.na(y))) {norm(x-y, "F")/norm(x, "F")} else{NA}),
         nmfl2_score_l2 = map2(true_scores, nmfl2_scores_re,
                               function(x,y) if(!any(is.na(y))) {norm(x-y, "F")/norm(x, "F")} else{NA}),
         nmfp_load_l2 = map2(true_patterns, nmfp_loadings_re,
                             function(x,y) if(!any(is.na(y))) {norm(x-y, "F")/norm(x, "F")} else{NA}),
         nmfp_score_l2 = map2(true_scores, nmfp_scores_re,
                              function(x,y) if(!any(is.na(y))) {norm(x-y, "F")/norm(x, "F")} else{NA})
  )

output_nofa_summaries <- output_nofa %>% 
  dplyr::select(seed, participants, chemicals, patterns,
                grep("rank", colnames(.)),
                grep("_l2er", colnames(.)),
                grep("ssdist", colnames(.)),
                grep("score_l2", colnames(.)),
                grep("load_l2", colnames(.))) %>% 
  unnest(pca_rank:nmfp_load_l2)

for (i in 1:nrow(output_nofa_summaries)) {
  no_fa_summaries = output_nofa_summaries[i,]
  job_num = to_do[i]
  save(no_fa_summaries, file = paste0("./HPC/Rout/combo_dim/dim_out_", job_num, ".RDA"))
}


