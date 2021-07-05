# SCRIPT 5
# Get metrics for all other models on all sims
# 1:12100

source("./functions/compare_functions.R")

# Read in Sims
# job_num = as.integer(Sys.getenv("SGE_TASK_ID"))
job_num = 1

# Get simulated dataset ####
chem <- read_csv(paste0("./sims/csvs/chem_sep_", job_num, ".csv")) %>% 
  as_tibble() %>% nest(data = c(V1:V40)) %>% rename(chem = data)

true_scores <- read_csv(paste0("./sims/csvs/scores_sep_", job_num, ".csv")) %>% 
  as_tibble() %>% nest(data = c(V1:V4)) %>% rename(true_scores = data)

# read files made by `run_other_models.R`
load(paste0("./main/other_models/output/sep_out_", job_num, ".RDA"))

all_sep = bind_cols(chem, true_scores, sep_out)

#### Get Stats ####

#### Normalize patterns ####
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
         truth = map(truth, function(df) df[,colSums(is.na(df))<nrow(df)]),
         relerr  = map2(truth, value, get_relerror),
         ssdist  = map2(truth, value, symm_subspace_dist),
         cosdist = map2(truth, value, cos_dist)) %>%
  unnest(c(relerr, ssdist, cosdist)) %>%
  dplyr::select(-true_scores, -true_patterns, -chem, -sim, -value, -truth)

#### Save ####

save(sep_rank,    file = paste0("./main/other_models/output/other_rank_", job_num, ".RDA"))
save(sep_metrics, file = paste0("./main/other_models/output/other_metrics_", job_num, ".RDA"))
# goes into `combine_other_models.R`
