library(tidyverse)

load("./HPC/Rout/dgp_rep1.RDA")
dgp_rep1

source("./R/fig_set.R")

# Relative Error
dgp_e1 <- dgp_rep1 %>% 
  dplyr::select(seed, data, sim, chem, grep("pred", colnames(.))) %>% 
  pivot_longer(c(pca_pred:nmf_p_pred)) %>% 
  mutate(l2_true = map2(chem, value, function (x,y) norm(x-y, "F")/norm(x, "F")),
         l2_sim = map2(sim, value, function (x,y) norm(x-y, "F")/norm(x, "F")),
         name = str_remove(name, "_pred")) %>% 
  dplyr::select(seed, data, name, l2_true, l2_sim) %>% 
  unnest(c(l2_sim, l2_true)) %>% 
  mutate(name = ifelse(str_detect(name, 'pca'), 'PCA', name),
         name = ifelse(str_detect(name, 'l2'), 'NMF L2', name),
         name = ifelse(str_detect(name, 'fa'), 'FA', name),
         name = ifelse(str_detect(name, '_p'), 'NMF P', name),
         name = ifelse(str_detect(name, 'bnmf'), 'BN2MF', name))

# Relative error in loadings and scores
results_rep1 <- 
  dgp_rep1 %>% 
  dplyr::select(seed, data, true_patterns, true_scores, pca_rotations, pca_scores,
                fa_rotations, fa_scores, nmf_l2_loadings, nmf_l2_scores,
                nmf_p_loadings, nmf_p_scores) %>% 
  mutate(true_patterns = map(true_patterns, t)) %>% 
  pivot_longer(true_patterns:nmf_p_scores,
               values_to = "results") %>% 
  filter(!grepl("true", name)) %>% 
  mutate(side = case_when(grepl("load|rot|eh", name) ~ "Loadings",
                          grepl("score|ewa", name) ~ "Scores"),
         model = case_when(grepl("pca", name) ~ "PCA",
                           grepl("fa", name) ~ "FA",
                           grepl("l2", name) ~ "NMF L2",
                           grepl("nmf_p", name) ~ "NMF P",
                           grepl("eh|ewa", name) ~ "BN2MF")) %>% 
  dplyr::select(-name)

truth_rep1 <- dgp_rep1 %>% 
  dplyr::select(seed, data, true_patterns, true_scores, pca_rotations, pca_scores,
                fa_rotations, fa_scores, nmf_l2_loadings, nmf_l2_scores,
                nmf_p_loadings, nmf_p_scores) %>% 
  mutate(true_patterns = map(true_patterns, t)) %>% 
  pivot_longer(true_patterns:nmf_p_scores,
               values_to = "truth") %>% 
  filter(grepl("true", name)) %>% 
  mutate(side = case_when(grepl("patterns", name) ~ "Loadings",
                          grepl("score", name) ~ "Scores")) %>% 
  dplyr::select(-name)

norm_rep1 <- 
  left_join(results_rep1, truth_rep1) %>% 
  mutate(rank = map(results, ncol)) %>% unnest(rank) %>% 
  mutate(results = case_when(grepl("NMF|2", model) & side == "Loadings" ~ map(results, t),
                             TRUE ~ map(results, as.matrix)
  )) %>% 
  mutate(rank = map(results, ncol)) %>% unnest(rank) %>%
  filter(rank == 4) %>% 
  mutate(results = map(results, as.matrix),
         l2 = map2(truth, results, function (x,y) norm(x-y, "F")/norm(x, "F")),
         l1 = map2(truth, results, function (x,y) sum(abs(x-y))/sum(abs(x)))) %>% 
  unnest(c(l1, l2)) %>% 
  dplyr::select(seed, data, side, model, results, l2, l1) 

# Investigate
dgp_e1 %>% 
  group_by(data, name) %>% 
  summarize(min = min(l2_true),
            med = median(l2_true),
            mean = mean(l2_true),
            q75 = quantile(l2_true, probs = 0.75),
            max = max(l2_true)) %>% 
  arrange(name, data)
  
norm_rep1 %>% 
  group_by(data, side, model) %>% 
  summarize(min = min(l2),
            med = median(l2),
            mean = mean(l2),
            q75 = quantile(l2, probs = 0.75),
            max = max(l2)) %>% 
  arrange(model, side, data)

# NMF Poisson
nmfp <- dgp_rep1$nmf_p_pred[[1]]
nmfp_load <- dgp_rep1$nmf_p_loadings[[1]]
nmfp_score <- dgp_rep1$nmf_p_scores[[1]]

true_load <- dgp_rep1$true_patterns[[1]]
true_score <- dgp_rep1$true_scores[[1]]
truth <- dgp_rep1$chem[[1]]

norm(truth -  nmfp, "F")/norm(truth, "F")
dgp_e1 %>% filter(name == "NMF P")

norm(true_load -  nmfp_load, "F")/norm(true_load, "F")
norm(true_score -  nmfp_score, "F")/norm(true_score, "F")
norm_rep1 %>% filter(model == "NMF P")

norm(truth -  nmfp_score%*%nmfp_load, "F")/norm(truth, "F")

head(truth)[,1:5]
head(nmfp)[,1:5]

head(true_load)[,1:5]
head(nmfp_load)[,1:5]
# Order is wrong!

head(true_score)
head(nmfp_score)
# Order is wrong!

# Need the factor correspondence code!
source("./R/factor_correspondence.R")

fact_load <- factor_correspondence(t(true_load), t(nmfp_load))
nmfp_load2 <- t(fact_load$rearranged)

fact_score <- nmfp_score %*% fact_load$permutation_matrix

norm(true_load -  nmfp_load2, "F")/norm(true_load, "F")
norm(true_score -  fact_score, "F")/norm(true_score, "F")

head(true_load)[,1:5]
head(nmfp_load2)[,1:5]
# Loadings too small

head(true_score)
head(fact_score)
# Scores too big

apply(nmfp_load, 1, sum)
apply(nmfp_score, 2, sum)

# NMF L2
nmfl2 <- dgp_rep1$nmf_l2_pred[[1]]
nml2_load <- dgp_rep1$nmf_l2_loadings[[1]]
nml2_score <- dgp_rep1$nmf_l2_scores[[1]]

true_load <- dgp_rep1$true_patterns[[1]]
true_score <- dgp_rep1$true_scores[[1]]
truth <- dgp_rep1$chem[[1]]

norm(truth -  nmfl2, "F")/norm(truth, "F")
dgp_e1 %>% filter(name == "NMF L2")

norm(true_load -  nml2_load, "F")/norm(true_load, "F")
norm(true_score -  nml2_score, "F")/norm(true_score, "F")
norm_rep1 %>% filter(model == "NMF L2")

norm(truth -  nml2_score%*%nml2_load, "F")/norm(truth, "F")

head(truth)[,1:5]
head(nmfp)[,1:5]

head(true_load)[,1:5]
head(nml2_load)[,1:5]
# Order is wrong!

head(true_score)
head(nml2_score)
# Order is wrong!

fact_load <- factor_correspondence(t(true_load), t(nml2_load))
nml2_load2 <- t(fact_load$rearranged)

fact_score <- nml2_score %*% fact_load$permutation_matrix

norm(true_load -  nml2_load2, "F")/norm(true_load, "F")
norm(true_score -  fact_score, "F")/norm(true_score, "F")

head(true_load)[,1:5]
head(nml2_load2)[,1:5]
# Loadings too big

head(true_score)
head(fact_score)
# Scores too small

apply(nml2_load, 1, sum)
apply(nml2_score, 2, sum) # Scores sum to 1 across 1000 sample size

# FA
fa <- dgp_rep1$fa_pred[[1]]
fa_load <- t(dgp_rep1$fa_rotations[[1]])
fa_score <- dgp_rep1$fa_scores[[1]]

true_load <- dgp_rep1$true_patterns[[1]]
true_score <- dgp_rep1$true_scores[[1]]
truth <- dgp_rep1$chem[[1]]

norm(truth -  fa, "F")/norm(truth, "F")
dgp_e1 %>% filter(name == "FA")

norm(true_load -  fa_load, "F")/norm(true_load, "F")
norm(true_score -  fa_score, "F")/norm(true_score, "F")
norm_rep1 %>% filter(model == "FA")

norm(truth -  fa_score%*%fa_load, "F")/norm(truth, "F")

head(truth)[,1:5]
head(fa)[,1:5]

head(true_load)[,1:5]
head(fa_load)[,1:5]
# Order is wrong!

head(true_score)
head(fa_score)
# Order is wrong!

fact_load <- factor_correspondence(t(true_load), t(fa_load))
fa_load2 <- t(fact_load$rearranged)

fact_score <- fa_score %*% fact_load$permutation_matrix

norm(true_load -  fa_load2, "F")/norm(true_load, "F")
norm(true_score -  fact_score, "F")/norm(true_score, "F")

head(true_load)[,1:5]
head(fa_load2)[,1:5]
# Loadings too small

head(true_score)
head(fact_score)
# Scores too small

apply(fa_load, 2, sum)
apply(fa_score, 2, sum)

# Conclusion: Run factor correspondence for ALL results
# Also: NMF over/underestimates 1 matrix and under/overestimates the other

