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
  dplyr::select(-drop)

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

dgp_rep1_all_re <- dgp_rep1_all_t %>% 
  dplyr::select(seed, data, 
                true_patterns, 
                grep("scores", colnames(.)),
                grep("loadings", colnames(.)),
                bnmf_loadings = eh, bnmf_scores = ewa,
                -grep("_ssdist", colnames(.))) %>% 
  pivot_longer(true_patterns:bnmf_scores,
               names_to = c("model", "matrix"),
               names_sep = "_") %>% 
  drop_na(value) %>% 
  mutate(matrix = ifelse(matrix == "patterns", "loadings", matrix)) %>% 
  group_by(seed, data, matrix) %>% 
  mutate(truth = value[model == "true"]) %>% 
  filter(model != "true") %>% 
  mutate(value_re = 
           case_when(str_detect(model, "a") ~ map2(truth, value, 
                                              function(x,y) if(ncol(y) == 4) 
                                              {factor_correspondence(as.matrix(x), 
                                              as.matrix(y), nn = FALSE)$rearranged} else{NA}),
                     TRUE ~ map2(truth, value, 
                           function(x,y) if(ncol(y) == 4) 
                           {factor_correspondence(as.matrix(x), 
                            as.matrix(y))$rearranged} else{NA})))

dgp_re <- dgp_rep1_all_re %>% 
  filter(!is.na(value_re)) %>%
  mutate(l2 = map2(truth, value_re, function (x,y) norm(x-y, "F")/norm(x, "F")),
         l1 = map2(truth, value_re, function (x,y) sum(abs(x-y))/sum(abs(x)))) %>% 
  unnest(c(l1, l2))

#####
# Rank
#####

dist_rank <- dgp_rep1_all %>%
  dplyr::select(seed, data, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank) %>%
  unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank)) %>%
  pivot_longer(cols = fa_rank:bnmf_rank) %>%
  mutate(value = ifelse(value > 5, "> 5", value)) %>% 
  group_by(name, value, data) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  filter(data == "Distinct") %>% 
  dplyr::select(-data) %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0)

over_rank <- dgp_rep1_all %>%
  dplyr::select(seed, data, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank) %>%
  unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank)) %>%
  pivot_longer(cols = fa_rank:bnmf_rank) %>%
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
  dplyr::select(seed, data, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank) %>%
  unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank)) %>%
  pivot_longer(cols = fa_rank:bnmf_rank) %>%
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

#pdf("./Figures/bnmf_error.pdf")
dgp_e1 %>%
  ggplot(aes(x = name, y = l2_true, color = name, fill = name)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(. ~ data, scales = "free") + 
  #geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  scale_y_log10() +
  theme(legend.position = "none") + 
  labs(y = "Relative Predictive Error",
       title = "")
#dev.off()


#####
# SSD
#####

# pdf("./Figures/bnmf_ssd.pdf", height = 10)
dgp_s1 %>% 
  # filter(!(model %in% c("PCA", "FA")))
  ggplot(aes(x = model, y = value, color = model, fill = model)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(type ~ data) + 
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  # scale_y_log10() +
  labs(y = "Symmetric Subspace Distance")
# dev.off()


#####
# Rank
#####

dist_rank
over_rank
cor_rank


#####
# Rel error loadings and scores
#####

#pdf("./Figures/bnmf_loadscore_error.pdf", height = 10)
dgp_rep1_all_re %>% 
  ggplot(aes(x = method, y = l2)) +
  geom_jitter(alpha = 0.25, size = 0.5, height = 0, width = .3) +
  geom_boxplot(aes(color = method, fill = method), 
               alpha = 0.75, outlier.shape = NA, varwidth = TRUE) +
  facet_grid(results ~ data, scales = "free") + 
  geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  scale_y_log10() +
  labs(y = "Relative Predictive Error",
       title = "")
#dev.off()

dgp_rep1_all_re %>% 
  group_by(data, results, method) %>% 
  summarize(min = min(l2),
            med = median(l2),
            mean = mean(l2),
            q75 = quantile(l2, probs = 0.75),
            max = max(l2)) %>% 
  arrange(method, data, results) %>% View()
