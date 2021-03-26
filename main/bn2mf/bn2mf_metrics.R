# Get metrics for all other models on all sims
# on hpc, run array 1:12100

source("./functions/compare_functions.R")

# Read in BN2MF results ####

# job_num = as.integer(Sys.getenv("SGE_TASK_ID"))

# Follow analysis with one example simulated dataset
job_num = 1

# Read in simulated dataset for comparison ####
chem <- read_csv(paste0("./sims/csvs/chem_sep_", job_num, ".csv")) %>%
  as_tibble() %>% nest(chem = c(V1:V40))

true_scores <- read_csv(paste0("./sims/csvs/scores_sep_", job_num, ".csv")) %>%
  as_tibble() %>% nest(true_scores = c(V1:V4))

sim <- read_csv(paste0("./sims/csvs/sim_sep_", job_num, ".csv")) %>%
  as_tibble() %>% nest(sim = c(V1:V40))

true_patterns <- read_csv(paste0("./sims/csvs/patterns_sep_", job_num, ".csv")) %>%
  as_tibble() %>% nest(true_patterns = c(V1:V40))

true_sep = bind_cols(chem, sim, true_patterns, true_scores)

# Get BN2MF Output ####

bn2mf_loadings = readMat(paste0("./main/bn2mf/output/eh_NOTscaled", job_num, ".mat"))[[1]] %>% 
  as_tibble(.) %>% nest(data = everything()) %>% rename(bn2mf_loadings = data)

bn2mf_scores = readMat(paste0("./main/bn2mf/output/ewa_NOTscaled", job_num, ".mat"))[[1]] %>% 
  as_tibble(.) %>% nest(data= everything()) %>% rename(bn2mf_scores = data)

ewa_scaled = readMat(paste0("./main/bn2mf/output/ewa_scaled", job_num, ".mat"))[[1]] %>% 
  as_tibble(.) %>% nest(data= everything()) %>% rename(ewa_scaled = data)

upperWA = readMat(paste0("./main/bn2mf/output/upperWA_", job_num, ".mat"))[[1]] %>% 
  as_tibble(.) %>% nest(data= everything()) %>% rename(upperWA = data)

lowerWA = readMat(paste0("./main/bn2mf/output/lowerWA_", job_num, ".mat"))[[1]] %>% 
  as_tibble(.) %>% nest(data= everything()) %>% rename(lowerWA = data)

# m has regular/unscaled loadings and scores to compare with other models
m_out   <- bind_cols(bn2mf_loadings, bn2mf_scores)

# vci has scaled scores and confidence interval to calculate coverage
vci_out <- bind_cols(ewa_scaled, lowerWA, upperWA)

sep_m <- m_out %>%
  mutate(bn2mf_loadings = map(bn2mf_loadings, as.matrix),
         bn2mf_scores = map(bn2mf_scores, as.matrix),
         bn2mf_pred = map2(bn2mf_scores, bn2mf_loadings, function(x,y) x %*% y),
         bn2mf_rank = map(bn2mf_scores, ncol))

# Get Stats ####

# Normalize truth ####
true_sep = true_sep %>%
  mutate(denom = map(true_patterns, function(x) apply(x, 1, sum)),
         true_patternsscaled = map2(true_patterns, denom, function(x,y) x/y),
         true_scoresscaled = map2(true_scores, denom, function(x,y) as.matrix(x) %*% diag(y)))

# Comb ####
all_sep = bind_cols(true_sep, sep_m)

# Rank ####
sep_rank = all_sep %>%
  dplyr::select(grep("rank", colnames(.))) %>%
  rename(rank = 1) %>% 
  unnest(rank) %>%
  mutate(model = "bn2mf",
         rank = ifelse(is.na(rank), 0, rank),
         rank_bin = ifelse(rank ==4, "right", "wrong"))

# Metrics ####
sep_metrics = all_sep %>%
  dplyr::select(-bn2mf_rank) %>%
  pivot_longer(bn2mf_loadings:bn2mf_pred,
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
  dplyr::select(-true_scores, -true_patterns, -denom, -chem, -sim, -value, -truth, -true_patternsscaled, -true_scoresscaled)

# VCI ####
vci_sep = bind_cols(true_sep, vci_out)

vci_sep = vci_sep %>% 
  mutate(denom = map(true_patterns, function(x) apply(x,1, sum)),
         patterns_scaled = map2(true_patterns, denom, function(x,y) x/y),
         scores_scaled = map2(true_scores, denom, function(x,y) as.matrix(x) %*% diag(y)))

vci_metrics = vci_sep %>% 
  mutate_all(~map(., as.matrix)) %>% 
  mutate(prop = pmap(list(scores_scaled, lowerWA, upperWA),
                     function(x,y,z) sum(x >= y & x <= z)/(nrow(x)*ncol(x)) )) %>%
  unnest(c(prop)) %>% 
  dplyr::select(prop) %>% 
  mutate(model = "bn2mf")

# Save ####
save(sep_rank,    file = paste0("./main/bn2mf/output/sep_rank_m_", job_num, ".RDA"))
save(sep_metrics, file = paste0("./main/bn2mf/output/sep_metrics_m_", job_num, ".RDA"))
save(vci_metrics, file = paste0("./main/bn2mf/output/sep_prop_m_", job_num, ".RDA"))

