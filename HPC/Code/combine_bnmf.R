# Packages
library(tidyverse)
library(R.matlab, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(CVXR, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

## read job number from system environment
## This only works if run on cluster!!
job_num = as.integer(Sys.getenv("SGE_TASK_ID"))
job_num

## Get functions
source("/ifs/scratch/msph/ehs/eag2186/npbnmf/compare_functions.R")
source("/ifs/scratch/msph/ehs/eag2186/npbnmf/factor_correspondence.R")

# Aggregate BN2MF results
eh <- readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bnmf_dim_out/eh_dim_", job_num, ".mat"))[[1]]
eh <- as_tibble(eh) %>% nest(everything()) %>% rename(eh = data)

ewa <- readMat(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/bnmf_dim_out/ewa_dim_", job_num, ".mat"))[[1]]
ewa <- as_tibble(ewa) %>% nest(everything()) %>% rename(ewa = data)

bn2mf_out <- bind_cols(eh, ewa)

bn2mf_out <- bn2mf_out %>% 
  mutate(eh = map(eh, as.matrix),
         ewa = map(ewa, as.matrix),
         bn2mf_pred = map2(ewa, eh, function(x,y) as.matrix(x %*% y)),
         bn2mf_rank = map(ewa, ncol))

# Merge with Sim data
load(paste0("/ifs/scratch/msph/ehs/eag2186/Data/sim_dim/sim_dim", job_num, ".RDA"))
to_save

all = bind_cols(to_save, bn2mf_out)
all

# Calculate all error metrics

# L2 relative error
output_all <- all %>% 
  mutate(chem = map(chem, as.matrix),
    bn2mf_l2er  = map2(chem, bn2mf_pred,   
                          function(x,y) norm(x-y, "F")/norm(x, "F") ))

# factor correspondence for l2 relative error & cosine dist on loadings and scores
output_all <- output_all %>%
  mutate(true_patterns_t   = map(true_patterns, t),
         eh_t  = map(eh, t)) %>%
  mutate(bn2mf_perm = map2(true_patterns_t, eh_t, 
                         function(x,y) if(ncol(y) == ncol(x))
                         {factor_correspondence(as.matrix(x),
                                                as.matrix(y), nn = FALSE)$permutation_matrix} else{NA}),
         bn2mf_loadings_re = map2(eh_t, bn2mf_perm, 
                                function(x,y) if(!any(is.na(y))) 
                                {x %*% y} else{NA}),
         bn2mf_scores_re = map2(ewa, bn2mf_perm, 
                              function(x,y) if(!any(is.na(y))) 
                              {x %*% y} else{NA}))

# l2 relative error on loadings and scores
output_all <- output_all %>% 
  mutate(bn2mf_loadings_l2 = map2(true_patterns_t, bn2mf_loadings_re, 
                            function(x,y) if(!any(is.na(y))) {norm(x-y, "F")/norm(x, "F")} else{NA}),
         bn2mf_scores_l2 = map2(true_scores, bn2mf_scores_re, 
                             function(x,y) if(!any(is.na(y))) {norm(x-y, "F")/norm(x, "F")} else{NA}))

# cosine distance
output_all <- output_all %>%
  mutate(bn2mf_cos_dist_loadings = map2(true_patterns_t, bn2mf_loadings_re, 
                                        function(x,y) if(!any(is.na(y))) {cos_dist(x,y)} else{NA}),
         bn2mf_cos_dist_scores   = map2(true_scores, bn2mf_scores_re, 
                                        function(x,y) if(!any(is.na(y))) {cos_dist(x,y)} else{NA}))

# symmetric subspace distance
output_all <- output_all %>%
  mutate(bn2mf_loadings_ssdist  = map2(true_patterns_t, eh_t,   symm_subspace_dist),
         bn2mf_scores_ssdist    = map2(true_scores,   ewa,     symm_subspace_dist))

output_all <- output_all %>% 
  dplyr::select(seed, participants, chemicals, patterns,
                grep("rank", colnames(.)),
                grep("_l2", colnames(.)),
                grep("ssdist", colnames(.)),
                grep("cos_dist", colnames(.)))

save(output_all, file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/combo_bnmf/bnmf_dim_out_", job_num, ".RDA"))


