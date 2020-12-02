#####
# Vizualize Results from Simulations
# Data generating process = log-normal scores %*% piecewise uniform loadings
#####

# Load packages
library(tidyverse)

# Load data
load("./HPC/Rout/dgp_rep1_all.RDA")

# Source factor corr function to rearrange loadings/scores to match truth
source("./R/factor_correspondence.R")

# Source ggplot settings
source("./R/fig_set.R") 

#####
# Relative Preditive Error
# L2 Norm (Truth - Predicted) / L2 Norm (Truth)
#####

dgp_e1 <- dgp_rep1_all %>% 
  dplyr::select(seed, data, sim, chem, grep("pred", colnames(.))) %>% 
  pivot_longer(c(pca_pred:bnmf_pred),
               names_to = "model") %>% 
  mutate(l2_true = map2(chem, value, function (x,y) norm(x-y, "F")/norm(x, "F")),
         l2_sim  = map2(sim,  value, function (x,y) norm(x-y, "F")/norm(x, "F")),
         model = str_to_upper(str_remove(model, "_pred"))) %>% 
  unnest(c(l2_sim, l2_true))

#####
# Subspace Distance
# Distance between linear subspaces (orthonormal bases)
# Loadging and scores
#####

dgp_s1 <- dgp_rep1_all %>% dplyr::select(seed, data, grep("_ssdist", colnames(.))) %>% 
  unnest(cols = grep("_ssdist", colnames(.))) %>% 
  pivot_longer(grep("_ssdist", colnames(.)),
               names_to = c("model", "matrix", "drop"),
               values_to = "ssdist",
               names_sep = "_") %>%
  mutate(model = str_to_upper(model)) %>% 
  dplyr::select(-drop) %>% 
  mutate(matrix = case_when(matrix == "rotations" ~ "Loadings",
                            matrix == "loading" ~ "Loadings",
                            TRUE ~ matrix))

#####
# Relative error on loadings and scores
#####

# Transpose all loading matrices to have chemical rows and pattern columns
dgp_rep1_all_t <- 
  dgp_rep1_all %>%
  mutate(true_patterns   = map(true_patterns, t),
         nmfl2_loadings = map(nmfl2_loadings, t),
         nmfp_loadings  = map(nmfp_loadings, t),
         eh = map(eh, t))

dgp_rep1_all_re <- 
  dgp_rep1_all_t %>% 
  dplyr::select(seed, data, 
                true_patterns, 
                grep("scores", colnames(.)),
                grep("loadings", colnames(.)),
                bnmf_loadings = eh, bnmf_scores = ewa,
                -grep("_ssdist", colnames(.))) %>% 
  pivot_longer(pca_scores:bnmf_scores,
               names_to = c("model", "matrix"),
               names_sep = "_") %>% 
  drop_na(value) %>% 
  pivot_wider(names_from = "matrix",
              values_from = "value") %>% 
  mutate(perm = case_when(str_detect(model, "a") ~ map2(true_patterns, loadings, 
                                                   function(x,y) if(ncol(y) == 4) 
                                                   {factor_correspondence(as.matrix(x), 
                                                   as.matrix(y), nn = FALSE)$permutation_matrix} else{NA}),
                          str_detect(model, "nmf") ~ map2(true_patterns, loadings, 
                                                     function(x,y) if(ncol(y) == 4) 
                                                     {factor_correspondence(as.matrix(x), 
                                                     as.matrix(y))$permutation_matrix} else{NA})),
         loadings_re = map2(loadings, perm, 
                            function(x,y) if(ncol(x) == 4) 
                            {as.matrix(x) %*% as.matrix(y)} else{NA}),
         scores_re   = map2(scores,   perm, 
                            function(x,y) if(ncol(x) == 4) 
                            {as.matrix(x) %*% as.matrix(y)} else{NA}))
         

save(dgp_rep1_all_re, file = "./HPC/Rout/dgp_rep1_reordered.RDA")
load("./HPC/Rout/dgp_rep1_reordered.RDA")

dgp_re <- dgp_rep1_all_re %>% 
  mutate(l2_loadings = map2(true_patterns, loadings_re, function (x,y) norm(x-y, "F")/norm(x, "F")),
         l2_scores = map2(true_scores, scores_re,     function (x,y) norm(x-y, "F")/norm(x, "F"))) %>% 
  unnest(c(l2_loadings, l2_scores)) %>% 
  dplyr::select(seed, data, model, l2_loadings, l2_scores) %>% 
  pivot_longer(l2_loadings:l2_scores,
               names_prefix = "l2_",
               names_to = "matrix",
               values_to = "l2")

#####
# Cosine distance
#####

dgp_cos <- dgp_rep1_all_re %>% 
  filter(!is.na(loadings_re)) %>% 
  mutate(cos_dist_loadings = map2(true_patterns, loadings_re, cos_dist),
         cos_dist_scores = map2(true_scores, scores_re, cos_dist)) %>% 
  dplyr::select(seed, data, model, cos_dist_loadings, cos_dist_scores) %>% 
  pivot_longer(c(cos_dist_loadings, cos_dist_scores),
               names_prefix = "cos_dist_",
               names_to = "matrix",
               values_to = "cosine_dist") %>% unnest(cosine_dist)

#####
# Rank
#####

dist_rank <- dgp_rep1_all %>%
  dplyr::select(seed, data, grep("rank", colnames(.))) %>%
  unnest(grep("rank", colnames(.))) %>%
  pivot_longer(cols = pca_rank:bnmf_rank,
               names_to = c("model", "drop"),
               names_sep = "_") %>%
  mutate(value = ifelse(value > 5, "> 5", value)) %>% 
  group_by(data, model, value) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  filter(data == "Distinct") %>% 
  dplyr::select(-data) %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0)

over_rank <- dgp_rep1_all %>%
  dplyr::select(seed, data, grep("rank", colnames(.))) %>%
  unnest(grep("rank", colnames(.))) %>%
  pivot_longer(cols = pca_rank:bnmf_rank) %>%
  mutate(value = ifelse(value > 5, "> 5", value)) %>% 
  group_by(name, value, data) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  filter(data == "Overlapping") %>% 
  dplyr::select(-data) %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0)

cor_rank <- dgp_rep1_all %>%
  dplyr::select(seed, data, grep("rank", colnames(.))) %>%
  unnest(grep("rank", colnames(.))) %>%
  pivot_longer(cols = pca_rank:bnmf_rank) %>%
  mutate(value = ifelse(value > 5, "> 5", value)) %>% 
  group_by(name, value, data) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  filter(data == "Correlated") %>% 
  dplyr::select(-data) %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0)


#####
# Viz
#####

#####
# Pred error
#####

# pdf("./Figures/bnmf_error.pdf", width = 10)
dgp_e1 %>%
  mutate(data = fct_relevel(data, "Distinct", "Overlapping", "Correlated"),
         model = ifelse(model == "BNMF", "BN2MF", model)) %>% 
  ggplot(aes(x = model, y = l2_true)) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  geom_boxplot(aes(color = model, fill = model),
              alpha = 0.5, outlier.shape = NA) +
  facet_grid(. ~ data, scales = "free") + 
  scale_y_log10() +
  theme(legend.position = "none") + 
  labs(y = "Relative Predictive Error")
# dev.off()

dgp_e1 %>% 
  group_by(data, model) %>% 
  summarise(qs = quantile(l2_true, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  arrange(model)

dgp_rep1_all %>%
  dplyr::select(seed, data, grep("rank", colnames(.))) %>%
  unnest(grep("rank", colnames(.))) %>%
  pivot_longer(cols = pca_rank:bnmf_rank,
               names_to = c("model", "drop"),
               names_sep = "_",
               values_to = "rank") %>%
  mutate(rank = ifelse(rank > 5, "> 5", rank),
         model = str_to_upper(model)) %>% 
  left_join(., dgp_e1) %>% 
  group_by(data, model, rank) %>% 
  summarise(qs = quantile(l2_true, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  arrange(model)


#####
# SSD
#####

#pdf("./Figures/bnmf_ssd.pdf", width = 10, height = 10)
dgp_s1 %>%
  mutate(data = fct_relevel(data, "Distinct", "Overlapping", "Correlated"),
         model = ifelse(model == "BNMF", "BN2MF", model),
         matrix = str_to_title(matrix)) %>% 
  ggplot(aes(x = model, y = ssdist)) +
  #geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  geom_boxplot(aes(color = model, fill = model),
               alpha = 0.5, outlier.shape = NA) +
  facet_grid(matrix ~ data) + 
  labs(y = "Symmetric Subspace Distance")
#dev.off()

dgp_s1 %>% 
  group_by(data, model, matrix) %>% 
  summarise(qs = quantile(ssdist, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  arrange(matrix, model) %>% print(., n = 30)

#####
# Rank
#####

dist_rank
over_rank
cor_rank


#####
# Rel error loadings and scores
#####

dgp_re <- dgp_re %>% 
  mutate(data = as.factor(data)) 

#pdf("./Figures/bnmf_loadscore_error.pdf", width = 10, height = 10)
dgp_re %>% 
  mutate(data = fct_relevel(data, "Distinct", "Overlapping", "Correlated"),
         model = str_to_upper(model),
         model = ifelse(model == "BNMF", "BN2MF", model),
         matrix = str_to_title(matrix)) %>% 
  ggplot(aes(x = model, y = l2)) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  geom_boxplot(aes(color = model, fill = model), 
               alpha = 0.5, outlier.shape = NA, varwidth = TRUE) +
  facet_grid(matrix ~ data, scales = "free") + 
  scale_y_log10() +
  labs(y = "Relative Predictive Error")
#dev.off()

dgp_re %>% 
  group_by(data, model, matrix) %>% 
  summarise(qs = quantile(l2, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  arrange(matrix, model) %>% print(., n = 30)

#####
# Cosine distance
#####

#pdf("./Figures/bnmf_cos.pdf", width = 10, height = 10)
dgp_cos %>% 
  mutate(data = fct_relevel(data, "Distinct", "Overlapping", "Correlated"),
         model = str_to_upper(model),
         model = ifelse(model == "BNMF", "BN2MF", model),
         matrix = str_to_title(matrix)) %>% 
  ggplot(aes(x = model, y = cosine_dist)) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  geom_boxplot(aes(color = model, fill = model), 
               alpha = 0.5, outlier.shape = NA) +
  facet_grid(matrix ~ data, scales = "free") + 
  labs(y = "Cosine Similarity")
#dev.off()

dgp_cos %>% 
  #filter(grepl("bnmf", model)) %>% 
  group_by(data, model, matrix) %>% 
  summarise(qs = quantile(cosine_dist, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75),
            min = min(cosine_dist)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  arrange(matrix, model) %>% print(., n = 30)


