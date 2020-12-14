## read job number from system environment
## This only works if run on cluster!!
job_num = as.integer(Sys.getenv("SGE_TASK_ID"))
job_num

load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/noise_out/noise_out_", job_num, ".RDA")) 

# Packages
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(CVXR, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

## Get functions
source("/ifs/scratch/msph/ehs/eag2186/npbnmf/compare_functions.R")
source("/ifs/scratch/msph/ehs/eag2186/npbnmf/factor_correspondence.R")

# Calculate all error metrics

output_all <- output_all %>% 
  mutate(pca_l2er   = map2(chem, pca_pred,   get_relerror),
         fa_l2er    = map2(chem, fa_pred,    get_relerror),
         nmfl2_l2er = map2(chem, nmfl2_pred, get_relerror),
         nmfp_l2er  = map2(chem, nmfp_pred,  get_relerror))

output_all <- output_all %>%
  mutate(pca_perm          = map2(true_patterns,  pca_loadings,   function(x,y) get_perm(x,y,nn=FALSE)),
         pca_loadings_re   = map2(pca_loadings,   pca_perm,       get_product),
         pca_scores_re     = map2(pca_scores,     pca_perm,       get_product),
         fa_perm           = map2(true_patterns,  fa_loadings,    function(x,y) get_perm(x,y,nn=FALSE)),
         fa_loadings_re    = map2(fa_loadings,    fa_perm,        get_product),
         fa_scores_re      = map2(fa_scores,      fa_perm,        get_product),
         nmfl2_perm        = map2(true_patterns,  nmfl2_loadings, get_perm),
         nmfl2_loadings_re = map2(nmfl2_loadings, nmfl2_perm,     get_product),
         nmfl2_scores_re   = map2(nmfl2_scores,   nmfl2_perm,     get_product),
         nmfp_perm         = map2(true_patterns,  nmfp_loadings,  get_perm),
         nmfp_loadings_re  = map2(nmfp_loadings,  nmfp_perm,      get_product),
         nmfp_scores_re    = map2(nmfp_scores,    nmfp_perm,      get_product))
         
output_all <- 
  output_all %>% 
  mutate(pca_load_l2    = map2(true_patterns, pca_loadings_re,   get_relerror),
         pca_score_l2   = map2(true_scores,   pca_scores_re,     get_relerror),
         fa_load_l2     = map2(true_patterns, fa_loadings_re,    get_relerror),
         fa_score_l2    = map2(true_scores,   fa_scores_re,      get_relerror),
         nmfl2_load_l2  = map2(true_patterns, nmfl2_loadings_re, get_relerror),
         nmfl2_score_l2 = map2(true_scores,   nmfl2_scores_re,   get_relerror),
         nmfp_load_l2   = map2(true_patterns, nmfp_loadings_re,  get_relerror),
         nmfp_score_l2  = map2(true_scores,   nmfp_scores_re,    get_relerror))

output_all <- 
  output_all %>% 
  mutate(pca_load_cos    = map2(true_patterns, pca_loadings_re,   cos_dist),
         pca_score_cos   = map2(true_scores,   pca_scores_re,     cos_dist),
         fa_load_cos     = map2(true_patterns, fa_loadings_re,    cos_dist),
         fa_score_cos    = map2(true_scores,   fa_scores_re,      cos_dist),
         nmfl2_load_cos  = map2(true_patterns, nmfl2_loadings_re, cos_dist),
         nmfl2_score_cos = map2(true_scores,   nmfl2_scores_re,   cos_dist),
         nmfp_load_cos   = map2(true_patterns, nmfp_loadings_re,  cos_dist),
         nmfp_score_cos  = map2(true_scores,   nmfp_scores_re,    cos_dist))

output_all <- 
  output_all %>% 
  mutate(pca_load_cos_v    = map2(true_patterns, pca_loadings_re,   cos_dist_v),
         pca_score_cos_v   = map2(true_scores,   pca_scores_re,     cos_dist_v),
         fa_load_cos_v     = map2(true_patterns, fa_loadings_re,    cos_dist_v),
         fa_score_cos_v    = map2(true_scores,   fa_scores_re,      cos_dist_v),
         nmfl2_load_cos_v  = map2(true_patterns, nmfl2_loadings_re, cos_dist_v),
         nmfl2_score_cos_v = map2(true_scores,   nmfl2_scores_re,   cos_dist_v),
         nmfp_load_cos_v   = map2(true_patterns, nmfp_loadings_re,  cos_dist_v),
         nmfp_score_cos_v  = map2(true_scores,   nmfp_scores_re,    cos_dist_v))

output_all <- output_all %>% 
  dplyr::select(seed, data,
                grep("rank", colnames(.)),
                grep("_l2", colnames(.)),
                grep("ssdist", colnames(.)),
                grep("_cos", colnames(.)))

save(output_all, file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/combo_noise/noise_out_", job_num, ".RDA"))


