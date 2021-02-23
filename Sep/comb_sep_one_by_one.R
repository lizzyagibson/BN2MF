# Get metrics for all other models on all sims
# 1:12100

source("/ifs/scratch/msph/ehs/eag2186/npbnmf/compare_functions.R")

# Read in Sims
job_num = as.integer(Sys.getenv("SGE_TASK_ID"))
job_num

chem <- read_csv(paste0("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_chem/chem_sep_", job_num, ".csv")) %>% 
  as_tibble() %>% 
  nest(chem = c(V1:V40))

true_scores <- read_csv(paste0("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv_scores/scores_sep_", job_num, ".csv")) %>% 
  as_tibble() %>% 
  nest(true_scores = c(V1:V4))

load(paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_models_out/sep_out_", job_num, ".RDA"))

all_sep = bind_cols(chem, true_scores, sep_out)

#### Get Stats ####

#### Normalize truth ####
all_sep = all_sep %>% 
  mutate(denom = map(true_patterns, function(x) apply(x, 1, sum)),
         true_patternsscaled = map2(true_patterns, denom, function(x,y) x/y),
         true_scoresscaled = map2(true_scores, denom, function(x,y) as.matrix(x) %*% diag(y)))

#### ERROR ####
all_sep = all_sep %>% 
  dplyr::select(-grep("(out|scaled|denom|perm)", colnames(.)))
all_sep

#### Rank ####
sep_rank = all_sep %>% 
  dplyr::select(grep("rank", colnames(.))) %>% 
  pivot_longer(pca_rank:nmfp_rank,
               names_to = c("model", "drop"),
               names_sep = "_",
               values_to = "rank") %>% 
  unnest(rank) %>% 
  mutate(rank = ifelse(is.na(rank), 0, rank),
         rank_bin = ifelse(rank ==4, "right", "wrong"))

##### Metrics ####
sep_metrics = all_sep %>%
  dplyr::select(-grep("rank", colnames(.))) %>%
    pivot_longer(pca_loadings:nmfp_pred,
                 names_to = c("model", "matrix"),
                 names_sep = "_") %>%
  mutate(truth = case_when(matrix == "loadings" ~ true_patterns,
                           matrix == "scores" ~ true_scores,
                           matrix == "pred" ~ chem),
         relerr  = map2(truth, value, get_relerror),
         ssdist  = map2(truth, value, symm_subspace_dist),
         cosdist = map2(truth, value, cos_dist)) %>%
  unnest(c(relerr, ssdist, cosdist)) %>%
  dplyr::select(-true_scores, -true_patterns, -chem, -sim, -value, -truth)

#### Save ####

save(sep_rank,    file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_rank/sep_rank_", job_num, ".RDA"))
save(sep_metrics, file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/separate/sep_metrics/sep_metrics_", job_num, ".RDA"))
