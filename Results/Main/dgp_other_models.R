#########################################################
# Data Generating Process Sims
# Run models: regular NMF (L2 & Poisson), PCA, and FA 
#########################################################

source("/ifs/scratch/msph/ehs/eag2186/npbnmf/compare_functions.R")

# Read in Sims
load("/ifs/scratch/msph/ehs/eag2186/Data/sim_dgp.RDA")

job_num = as.integer(Sys.getenv("SGE_TASK_ID"))
job_num

#####
# Run everything
#####

#####
# PCA
#####
dgp_out <- sim_dgp[job_num,] %>%
  mutate(pca_out       = map(sim, get_pca),
         pca_loadings  = map(pca_out, function(x) x[[1]]),
         pca_scores    = map(pca_out, function(x) x[[2]]),
         pca_pred      = map(pca_out, function(x) x[[3]]),
         pca_rank      = map(pca_out, function(x) x[[4]]))

#####
# Factor Analysis
#####
dgp_out <- dgp_out %>%
  mutate(fa_out       = map(sim, function(x) get_fa(x, 4)),
         fa_loadings  = map(fa_out, function(x) x[[1]]),
         fa_scores    = map(fa_out, function(x) x[[2]]),
         fa_pred      = map(fa_out, function(x) x[[3]]),
         fa_rank      = map(fa_out, function(x) x[[4]]))

#####
# L2 NMF
#####
dgp_out <- dgp_out %>%
  mutate(nmfl2_out      = map(sim, function(x) get_nmfl2(x, 4)),
         nmfl2_loadings = map(nmfl2_out, function(x) x[[1]]),
         nmfl2_scores   = map(nmfl2_out, function(x) x[[2]]),
         nmfl2_pred     = map(nmfl2_out, function(x) x[[3]]),
         nmfl2_rank     = map(nmfl2_out, function(x) x[[4]]))

#####
# Poisson NMF
#####
dgp_out <- dgp_out %>%
  mutate(nmfp_out      = map(sim, function(x) get_nmfp(x, 4)),
         nmfp_loadings = map(nmfp_out, function(x) x[[1]]),
         nmfp_scores   = map(nmfp_out, function(x) x[[2]]),
         nmfp_pred     = map(nmfp_out, function(x) x[[3]]),
         nmfp_rank     = map(nmfp_out, function(x) x[[4]]))

save(dgp_out, file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/dgp_main_", job_num, ".RDA"))
# Combine all