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
source("./R/fig_set.R")

# Read in Sims
load("./Sims/sim_dgp_rep1.RDA")
check_r <- sim_dgp_rep1 %>% 
  dplyr::select(sim) %>% mutate(row = 1:nrow(.),
                                sim = map(sim, as_tibble))
check_r

# Read MATLAB sims
check_mat <- tibble()
for (i in 1:nrow(sim_dgp_rep1)) {

  check_data <- read_csv(paste0("./Sims/dgp_rep1/sim_dgp_rep1_", i, ".csv")) %>% 
    nest(sim = everything()) %>% mutate(row = i)
  check_mat <- rbind(check_mat, check_data)
}
check_mat

all.equal(check_r %>% dplyr::select(sim), check_mat %>% dplyr::select(sim))
# GOOD

#####
# Run everything
#####

#####
# PCA
#####
dgp_rep1_overcor <- sim_dgp_rep1[101:300,] %>% 
  mutate(pca_out       = map(sim, get_pca),
         pca_rotations = map(pca_out, function(x) x[[1]]),
         pca_scores    = map(pca_out, function(x) x[[2]]),
         pca_pred      = map(pca_out, function(x) x[[3]]),
         pca_rank      = map(pca_out, function(x) x[[4]]))

#####
# Factor Analysis
#####
dgp_rep1_overcor <- dgp_rep1_overcor %>%
  mutate(fa_out       = map(sim, get_fa),
         fa_rotations = map(fa_out, function(x) x[[1]]),
         fa_scores    = map(fa_out, function(x) x[[2]]),
         fa_pred      = map(fa_out, function(x) x[[3]]),
         fa_rank      = map(fa_out, function(x) x[[4]]))

#####
# L2 NMF
#####
dgp_rep1_overcor <- dgp_rep1_overcor %>%
  mutate(nmf_l2_out      = map(sim, get_nmf_l2),
         nmf_l2_loadings = map(nmf_l2_out, function(x) x[[1]]),
         nmf_l2_scores   = map(nmf_l2_out, function(x) x[[2]]),
         nmf_l2_pred     = map(nmf_l2_out, function(x) x[[3]]),
         nmf_l2_rank     = map(nmf_l2_out, function(x) x[[4]]))

#####
# Poisson NMF
#####
dgp_rep1_overcor <- dgp_rep1_overcor %>% 
  mutate(nmf_p_out      = map(sim, get_nmf_p),
         nmf_p_loadings = map(nmf_p_out, function(x) x[[1]]),
         nmf_p_scores   = map(nmf_p_out, function(x) x[[2]]),
         nmf_p_pred     = map(nmf_p_out, function(x) x[[3]]),
         nmf_p_rank     = map(nmf_p_out, function(x) x[[4]]))

#####
# Symmetric Subspace Distance
#####
dgp_rep1_overcor <- dgp_rep1_overcor %>% 
  mutate(pca_rotation_ssdist   = map2(true_patterns, pca_rotations,   symm_subspace_dist),
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
# dgp_rep1 <- dgp_rep1 %>% dplyr::select(-grep("_out", colnames(.)))
# save(dgp_rep1, file = "./HPC/Rout/dgp_rep1.RDA")
load("./HPC/Rout/dgp_rep1.RDA")
dgp_rep1

dgp_rep1_overcor <- dgp_rep1_overcor %>% dplyr::select(-grep("_out", colnames(.)))
save(dgp_rep1_overcor, file = "./HPC/Rout/dgp_rep1_overcor.RDA")

#####

# Add MATLAB
dgp_bnmf_rep1 <- tibble()

for (i in 1:200) {
  eh <- readMat(here::here(paste0("./MATLAB/loop_dgp_rep1/rep1_eh_dist_", i, ".mat")))[[1]]
  eh <- as_tibble(eh) %>% nest(eh = everything()) %>% mutate(row = i)
  
  ewa <- readMat(here::here(paste0("/MATLAB/loop_dgp_rep1/rep1_ewa_dist_", i, ".mat")))[[1]]
  ewa <- as_tibble(ewa) %>% nest(ewa = everything()) %>% mutate(row = i)
  both <- full_join(eh, ewa, by = "row")
  dgp_bnmf_rep1 <- rbind(dgp_bnmf_rep1, both)
}

dgp_bnmf_rep1 <- dgp_bnmf_rep1 %>% 
  mutate(seed = rep(1:100, 2),
         data = c(rep("Distinct", 100), rep("Overlapping", 100)))
dgp_bnmf_rep1

dgp_rep1 <- full_join(dgp_rep1, dgp_bnmf_rep1, by = c("seed", "data")) %>% 
  mutate(bnmf_pred            = map2(ewa, eh, function(x,y) as.matrix(x) %*% as.matrix(y)),
         bnmf_loading_ssdist  = map2(true_patterns, eh,  symm_subspace_dist),
         bnmf_scores_ssdist   = map2(true_scores,   ewa,    symm_subspace_dist),
         bnmf_rank = map(eh, nrow))

dgp_rep1

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
  unnest(c(l2_sim, l2_true)) %>% 
  mutate(name = ifelse(str_detect(name, 'pca'), 'PCA', name),
         name = ifelse(str_detect(name, 'l2'), 'NMF L2', name),
         name = ifelse(str_detect(name, 'fa'), 'FA', name),
         name = ifelse(str_detect(name, '_p'), 'NMF P', name),
         name = ifelse(str_detect(name, 'bnmf'), 'BN2MF', name))

# dat <- dgp_rep1 %>% filter(seed == 34 & data == "Overlapping")
# chem <- dat$chem[[1]]
# eh <- dat$eh[[1]]
# ewa <- dat$ewa[[1]]
# pred <- as.matrix(ewa) %*% as.matrix(eh)
# norm(chem - pred, "F")/norm(chem, "F")

# Subspace Distance
dgp_s1 <- dgp_rep1 %>% dplyr::select(seed, data, grep("_ssdist", colnames(.))) %>% 
  unnest(cols = grep("_ssdist", colnames(.))) %>% 
  pivot_longer(grep("_ssdist", colnames(.)),
               names_to = "type",
               values_to = "value") %>%
  mutate(model = str_sub(type, 1, 6),
         model = ifelse(str_detect(model, 'pca'), 'PCA', model),
         model = ifelse(str_detect(model, 'l2'), 'NMF L2', model),
         model = ifelse(str_detect(model, 'fa'), 'FA', model),
         model = ifelse(str_detect(model, '_p'), 'NMF P', model),
         model = ifelse(str_detect(model, 'bnmf'), 'BN2MF', model),
         type = ifelse(str_detect(type, 'scores'), 'Scores', "Loadings"))

results_rep1 <- 
  dgp_rep1 %>% 
  dplyr::select(seed, data, row, true_patterns, true_scores, pca_rotations, pca_scores,
                fa_rotations, fa_scores, nmf_l2_loadings, nmf_l2_scores,
                nmf_p_loadings, nmf_p_scores, eh, ewa) %>% 
  mutate(true_patterns = map(true_patterns, t)) %>% 
  pivot_longer(true_patterns:ewa,
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
  dplyr::select(seed, data, row, true_patterns, true_scores, pca_rotations, pca_scores,
                fa_rotations, fa_scores, nmf_l2_loadings, nmf_l2_scores,
                nmf_p_loadings, nmf_p_scores, eh, ewa) %>% 
  mutate(true_patterns = map(true_patterns, t)) %>% 
  pivot_longer(true_patterns:ewa,
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
  dplyr::select(seed, data, side, model, row, results, l2, l1) 

# Rank
dist_rank <- dgp_rep1 %>%
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

over_rank <- dgp_rep1 %>%
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

# Viz

# Pred error
pdf("./Figures/bnmf_error.pdf")
dgp_e1 %>%
  ggplot(aes(x = name, y = l2_true, color = name, fill = name)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(. ~ data, scales = "free") + 
  #geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  scale_y_log10() +
  theme(legend.position = "none") + 
  labs(y = "Relative Predictive Error",
       title = "")
dev.off()

# SSD
pdf("./Figures/bnmf_ssd.pdf", height = 10)
dgp_s1 %>% 
  ggplot(aes(x = model, y = value, color = model, fill = model)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(type ~ data) + 
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  # scale_y_log10() +
  labs(y = "Symmetric Subspace Distance")
dev.off()

# Rank
dist_rank
over_rank

pdf("./Figures/bnmf_loadscore_error.pdf", height = 10)
norm_rep1 %>% 
  ggplot(aes(x = model, y = l2)) +
  geom_jitter(alpha = 0.25, size = 0.5, height = 0, width = .3) +
  geom_boxplot(aes(color = model, fill = model), 
               alpha = 0.75, outlier.shape = NA, varwidth = TRUE) +
  facet_grid(side ~ data, scales = "free") + 
  geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  scale_y_log10() +
  labs(y = "Relative Predictive Error",
       title = "")
dev.off()

norm_rep1 %>% 
  group_by(data, side, model) %>% 
  summarize(min = min(l2),
            med = median(l2),
            mean = mean(l2),
            q75 = quantile(l2, probs = 0.75),
            max = max(l2)) %>% 
  arrange(model, side, data)
