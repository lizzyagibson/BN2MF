#########################################################
# Sims, npBNMF, regular NMF (L2 & Poisson), PCA, and FA #
#########################################################
# Lizzy Gibson    # 6/22/2020 ###########################
# Added BNMF      # 9/17/2020 ###########################
# Cleaned up code # 9/22/2020 ###########################
# Rerun on local  # 10/24/2020
#########################################################

# Packages
library(tidyverse)
library(psych)
library(NMF)
library(R.matlab)

source("./R/compare_functions.R")

# Read in Sims
load("./Sims/sim_dgp_rep1.RDA")

# Run everything
dist <- sim_dgp_rep1$sim[[1]]
as_tibble(dist) %>% 
  select(1:5) %>% 
  pivot_longer(1:5) %>%
  ggplot(aes(x = value)) + geom_histogram() + facet_wrap(~name)

summary(apply(dist, 2, mean))
summary(apply(dist, 2, median))
summary(apply(dist, 2, sd))

#####
# PCA
#####
dgp_rep1 <- sim_dgp_rep1 %>% 
  mutate(pca_out       = map(sim, get_pca),
         pca_rotations = map(pca_out, function(x) x[[1]]),
         pca_scores    = map(pca_out, function(x) x[[2]]),
         pca_pred      = map(pca_out, function(x) x[[3]]),
         pca_rank      = map(pca_out, function(x) x[[4]]))

#####
# Factor Analysis
#####
dgp_rep1 <- dgp_rep1 %>%
  mutate(fa_out       = map(sim, get_fa),
         fa_rotations = map(fa_out, function(x) x[[1]]),
         fa_scores    = map(fa_out, function(x) x[[2]]),
         fa_pred      = map(fa_out, function(x) x[[3]]),
         fa_rank      = map(fa_out, function(x) x[[4]]))

#####
# L2 NMF
#####
dgp_rep1 <- dgp_rep1 %>%
  mutate(nmf_l2_out      = map(sim, get_nmf_l2),
         nmf_l2_loadings = map(nmf_l2_out, function(x) x[[1]]),
         nmf_l2_scores   = map(nmf_l2_out, function(x) x[[2]]),
         nmf_l2_pred     = map(nmf_l2_out, function(x) x[[3]]),
         nmf_l2_rank     = map(nmf_l2_out, function(x) x[[4]]))

#####
# Poisson NMF
#####
dgp_rep1 <- dgp_rep1 %>% 
  mutate(nmf_p_out      = map(sim, get_nmf_p),
         nmf_p_loadings = map(nmf_p_out, function(x) x[[1]]),
         nmf_p_scores   = map(nmf_p_out, function(x) x[[2]]),
         nmf_p_pred     = map(nmf_p_out, function(x) x[[3]]),
         nmf_p_rank     = map(nmf_p_out, function(x) x[[4]]))

#####
# Symmetric Subspace Distance
#####
dgp_rep1 <- dgp_rep1 %>% 
  mutate(pca_norm              = map2(sim, pca_pred,    function(x,y) norm(x-y, "F")/norm(x, "F")),
         fa_norm               = map2(sim, fa_pred,     function(x,y) norm(x-y, "F")/norm(x, "F")),
         nmf_l2_norm           = map2(sim, nmf_l2_pred, function(x,y) norm(x-y, "F")/norm(x, "F")),
         nmf_p_norm            = map2(sim, nmf_p_pred,  function(x,y) norm(x-y, "F")/norm(x, "F")),
         pca_rotation_ssdist   = map2(true_patterns, pca_rotations,   symm_subspace_dist),
         pca_scores_ssdist     = map2(true_scores,   pca_scores,      symm_subspace_dist),
         fa_rotations_ssdist   = map2(true_patterns, fa_rotations,    symm_subspace_dist),
         fa_scores_ssdist      = map2(true_scores,   fa_scores,       symm_subspace_dist),
         nmf_l2_loading_ssdist = map2(true_patterns, nmf_l2_loadings, symm_subspace_dist),
         nmf_l2_scores_ssdist  = map2(true_scores,   nmf_l2_scores,   symm_subspace_dist),
         nmf_p_loading_ssdist  = map2(true_patterns, nmf_p_loadings,  symm_subspace_dist),
         nmf_p_scores_ssdist   = map2(true_scores,   nmf_p_scores,    symm_subspace_dist))

#####
# Save
#####
dgp_rep1 <- dgp_rep1 %>% dplyr::select(-grep("_out", colnames(.)))
# save(dgp_rep1, file = "./HPC/Rout/dgp_rep1.RDA")
load("./HPC/Rout/dgp_rep1.RDA")

#####
# Add MATLAB
dgp_bnmf_rep1 <- tibble()

for (i in 1:200) {
  eh <- readMat(here::here(paste0("/MATLAB/dgp_rep1/rep1_eh_dist_", i, ".mat")))[[1]]
  eh <- as_tibble(eh) %>% nest(eh = everything()) %>% mutate(seed = i)
  
  ewa <- readMat(here::here(paste0("/MATLAB/dgp_rep1/rep1_ewa_dist_", i, ".mat")))[[1]]
  ewa <- as_tibble(ewa) %>% nest(ewa = everything()) %>% mutate(seed = i)
  both <- full_join(eh, ewa, by = "seed")
  dgp_bnmf_rep1 <- rbind(dgp_bnmf_rep1, both) %>% drop_na(seed)
}

dgp_rep1 <- left_join(dgp_rep1, dgp_bnmf_rep1, by = "seed") %>% 
  mutate(bnmf_pred            = map2(ewa, eh, function(x,y) as.matrix(x) %*% as.matrix(y)),
         bnmf_norm            = map2(sim, bnmf_pred,    function(x,y) norm(x-y, "F")/norm(x, "F")),
         bnmf_loading_ssdist  = map2(true_patterns, eh,  symm_subspace_dist),
         bnmf_scores_ssdist   = map2(true_scores,   ewa,    symm_subspace_dist),
         bnmf_rank = map(eh, nrow))

#####
# Viz
#####
# Relative Error
dgp_e1 <- dgp_rep1 %>% 
  dplyr::select(seed, data, sim, chem, grep("pred", colnames(.))) %>% 
  pivot_longer(c(pca_pred:bnmf_pred)) %>% 
  mutate(l2_true = map2(chem, value, function (x,y) norm(x-y, "F")/norm(x, "F")),
         l2_sim = map2(sim, value, function (x,y) norm(x-y, "F")/norm(x, "F")),
         name = str_remove(name, "_pred")) %>% 
  dplyr::select(seed, data, name, l2_true, l2_sim) %>% 
  unnest(c(l2_sim, l2_true))

dgp_e1 %>%
  ggplot(aes(x = name, y = l2_true, color = name, fill = name)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(. ~ data, scales = "free") + 
  geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  scale_y_log10() +
  theme(legend.position = "none") + 
  labs(y = "Relative Predictive Error",
       title = "vs PRE NOISE TRUTH")

# Subspace Distance
dgp_s1 <- dgp_rep1 %>% dplyr::select(seed, data, grep("_ssdist", colnames(.))) %>% 
  unnest(cols = grep("_ssdist", colnames(.))) %>% 
  pivot_longer(grep("_ssdist", colnames(.)),
               names_to = "type",
               values_to = "value") %>%
  mutate(model = str_sub(type, 1, 6),
         model = ifelse(str_detect(model, 'pca'), 'PCA', model),
         model = ifelse(str_detect(model, 'l2'), 'L2 NMF', model),
         model = ifelse(str_detect(model, 'fa'), 'FA', model),
         model = ifelse(str_detect(model, '_p_'), 'Poisson NMF', model),
         model = ifelse(str_detect(model, 'bnmf'), 'BNMF', model),
         type = ifelse(str_detect(type, 'scores'), 'Scores', "Loadings"))

dgp_s1 %>% 
  ggplot(aes(x = model, y = value, color = model)) +
  geom_boxplot() +
  facet_grid(data ~ type) + 
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  theme(legend.position = "none") + 
  # scale_y_log10() +
  labs(y = "Symmetric Subspace Distance")

# Rank
dgp_rep1 %>%
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

dgp_rep1 %>%
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

