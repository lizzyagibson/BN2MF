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

# These jobs have impossible ranks??
to_do_2
to_do_both = c(to_do, to_do_2)
to_do_both = sort(to_do_both)

save(to_do_both, file = "./R/to_do_fa.RDA")

## Get functions
source("./R/compare_functions.R")
source("./R/factor_correspondence.R")

# Read in Sims
load(paste0("./Sims/sim_dim.RDA"))

dim_local = sim_over %>% 
  slice(to_do_2) %>% 
  mutate(id = to_do_2)

sim = dim_local[1,]$sim[[1]]
dim(sim)

#####
# PCA
#####

output_wrongrank <- dim_local %>% 
  mutate(pca_out       = map(sim, get_pca)) %>% 
  mutate(pca_loadings  = map(pca_out, function(x) x[[1]]),
         pca_scores    = map(pca_out, function(x) x[[2]]),
         pca_pred      = map(pca_out, function(x) x[[3]]),
         pca_rank      = map(pca_out, function(x) x[[4]]))

output_wrongrank %>% 
  mutate(fa_out       = map2(sim, patterns, function(x,y) try(get_fa_dim)))

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
load("./HPC/Rout/dim_out_nofa.RDA")
output_nofa %>% 
  filter(patterns == 4)

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

# Check
check = tibble()
for (i in 1:length(to_do)) {
  job_num = to_do[i]
  load(paste0("./HPC/Rout/combo_dim/dim_out_", job_num, ".RDA"))
  no_fa_summaries = no_fa_summaries %>% mutate(id = job_num)
  check = bind_rows(check, no_fa_summaries) 
}
check %>% dplyr::select(id, everything())

