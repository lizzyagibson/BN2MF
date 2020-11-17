# Packages
library(tidyverse)
library(CVXR, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

## read job number from system environment
## This only works if run on cluster!!
job_num = as.integer(Sys.getenv("SGE_TASK_ID"))
job_num

## Get functions
source("/ifs/scratch/msph/ehs/eag2186/npbnmf/compare_functions.R")
source("/ifs/scratch/msph/ehs/eag2186/npbnmf/factor_correspondence.R")

load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/dim_out/dim_out", job_num, ".RDA")) 
load("/Users/lizzy/BN2MF/HPC/Rout/dim_out/dim_out_1_ex.RDA") 
output_all

# Calculate all error metrics

output_all <- output_all %>% 
  mutate(pca_l2er  = map2(chem,   pca_pred,   
                        function(x,y) norm(x-y, "F")/norm(x, "F") ),
         fa_l2er   = map2(chem,   fa_pred,
                        function(x,y) norm(x-y, "F")/norm(x, "F") ),
         nmfl2_l2er = map2(chem,  nmfl2_pred, 
                        function(x,y) norm(x-y, "F")/norm(x, "F") ),
         nmfp_l2er   = map2(chem, nmfp_pred,
                        function(x,y) norm(x-y, "F")/norm(x, "F") ))

output_all <- output_all %>%
  mutate(true_patterns   = map(true_patterns, t),
         nmfl2_loadings  = map(nmfl2_loadings, t),
         nmfp_loadings   = map(nmfp_loadings, t)) %>% 
  mutate(pca_loadings_re = map2(true_patterns, pca_loadings, 
                           function(x,y) if(ncol(y) == 4) 
                                         {factor_correspondence(as.matrix(x), 
                                          as.matrix(y), nn = FALSE)$rearranged} else{NA}),
         pca_scores_re = map2(true_scores, pca_scores, 
                                function(x,y) if(ncol(y) == 4) 
                                {factor_correspondence(as.matrix(x), 
                                                       as.matrix(y), nn = FALSE)$rearranged} else{NA}),
         fa_loadings_re = map2(true_patterns, fa_loadings, 
                                function(x,y) if(ncol(y) == 4) 
                                {factor_correspondence(as.matrix(x), 
                                                       as.matrix(y), nn = FALSE)$rearranged} else{NA}),
         fa_scores_re = map2(true_scores, fa_scores, 
                                function(x,y) if(ncol(y) == 4) 
                                {factor_correspondence(as.matrix(x), 
                                                       as.matrix(y), nn = FALSE)$rearranged} else{NA}),
         nmfl2_loadings_re = map2(true_patterns, nmfl2_loadings, 
                                function(x,y) if(ncol(y) == 4) 
                                {factor_correspondence(as.matrix(x), 
                                                       as.matrix(y))$rearranged} else{NA}),
         nmfl2_scores_re = map2(true_scores, nmfl2_scores, 
                                function(x,y) if(ncol(y) == 4) 
                                {factor_correspondence(as.matrix(x), 
                                                       as.matrix(y))$rearranged} else{NA}),
         nmfp_loadings_re = map2(true_patterns, nmfp_loadings, 
                                function(x,y) if(ncol(y) == 4) 
                                {factor_correspondence(as.matrix(x), 
                                                       as.matrix(y))$rearranged} else{NA}),
         nmfp_scores_re = map2(true_scores, nmfp_scores, 
                                function(x,y) if(ncol(y) == 4) 
                                {factor_correspondence(as.matrix(x), 
                                                       as.matrix(y))$rearranged} else{NA}))
         
output_all <- 
  output_all %>% 
  mutate(pca_load_l2 = map2(true_patterns, pca_loadings_re, 
                            function(x,y) if(!any(is.na(y))) {norm(x-y, "F")/norm(x, "F")} else{NA}),
         pca_score_l2 = map2(true_scores, pca_scores_re, 
                            function(x,y) if(!any(is.na(y))) {norm(x-y, "F")/norm(x, "F")} else{NA}),
        fa_load_l2 = map2(true_patterns, fa_loadings_re,
                           function(x,y) if(!any(is.na(y))) {norm(x-y, "F")/norm(x, "F")} else{NA}),
        fa_score_l2 = map2(true_scores, fa_scores_re,
                           function(x,y) if(!any(is.na(y))) {norm(x-y, "F")/norm(x, "F")} else{NA}),
        nmfl2_load_l2 = map2(true_patterns, nmfl2_loadings_re,
                           function(x,y) if(!any(is.na(y))) {norm(x-y, "F")/norm(x, "F")} else{NA}),
         nmfl2_score_l2 = map2(true_scores, nmfl2_scores_re,
                             function(x,y) if(!any(is.na(y))) {norm(x-y, "F")/norm(x, "F")} else{NA}),
         nmfp_load_l2 = map2(true_patterns, nmfp_loadings_re,
                            function(x,y) if(!any(is.na(y))) {norm(x-y, "F")/norm(x, "F")} else{NA}),
         nmfp_score_l2 = map2(true_scores, nmfp_scores_re,
                             function(x,y) if(!any(is.na(y))) {norm(x-y, "F")/norm(x, "F")} else{NA})
        )

output_all <- output_all %>% 
  dplyr::select(seed, participants, chemicals, patterns,
                grep("rank", colnames(.)),
                grep("_l2er", colnames(.)),
                grep("ssdist", colnames(.)),
                grep("score_l2", colnames(.)),
                grep("load_l2", colnames(.))) %>% 
  unnest(pca_rank:nmfp_load_l2)

save(output_all, file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/combo_dim/dim_out_", job_num, ".RDA"))


# output_all %>% unnest(pca_rank:nmfp_load_l2) %>% View()

