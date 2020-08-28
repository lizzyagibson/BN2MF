#######################################################################################
#######################################################################################
#######################################################################################

# Aggregate HPC output

# Distinct patterns

# dist_out <- tibble()
# 
# for (i in 1:100) {
#   load(paste0("./HPC/out_", i, "_dist.RDA"))
#   dist_out[i,1:32] <- out_dist # 32 columns
# }
# load(paste0("./HPC/Rout/Standardized/out_", 1, "_dist.RDA"))
# out_dist

# bnmf_dist_out <- tibble(seed = NA, eh = NA, ewa = NA)
# 
# for (i in 1:100) {
#   eh <- readMat(here::here(paste0("./HPC/MATLABout/stand/eh_dist_", i, "_stand.mat")))[[1]]
#   eh <- as_tibble(eh) %>% nest(eh = everything()) %>% mutate(seed = i)
#   ewa <- readMat(here::here(paste0("./HPC/MATLABout/stand/ewa_dist_", i, "_stand.mat")))[[1]]
#   ewa <- as_tibble(ewa) %>% nest(ewa = everything()) %>% mutate(seed = i)
#   both <- full_join(eh, ewa, by = "seed")
#   bnmf_dist_out <- rbind(bnmf_dist_out, both) %>% drop_na(seed)
# }
# 
# dist_out <- left_join(dist_out, bnmf_dist_out, by = "seed")

# dist_out %>%
#   mutate(bnmf_rank = map(eh, nrow)) %>% 
#   dplyr::select(seed, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank) %>% 
#   unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank)) %>% 
#   pivot_longer(cols = fa_rank:bnmf_rank) %>% 
#   group_by(name, value) %>% 
#   summarise(n())

# Symmetric subspace distance
symm_subspace_dist <- function(U, V) {
  if (nrow(U) != max(nrow(U), ncol(U))) {U <- t(U)}
  if (nrow(V) != max(nrow(V), ncol(V))) {V <- t(V)}
  qrU <- qr.Q(qr(U))
  qrV <- qr.Q(qr(V))
  m <- ncol(U)
  n <- ncol(V)
  dUV <- sqrt( max(m,n) - sum((t(qrU) %*% qrV)^2) )
  ratio <- dUV/sqrt( max(m,n))
  ratio
}

# dist_out <- dist_out %>% 
#   mutate(bnmf_scores_ssdist = map2(true_scores, ewa, symm_subspace_dist),
#          bnmf_loading_ssdist = map2(true_patterns, eh, symm_subspace_dist))
# 
# dist_out %>% dplyr::select(seed, grep("_ssdist", colnames(.))) %>% 
#   unnest(cols = c(pca_rotation_ssdist, pca_scores_ssdist, fa_rotations_ssdist, 
#                   fa_scores_ssdist, nmf_l2_loading_ssdist, nmf_l2_scores_ssdist, 
#                   nmf_p_loading_ssdist, nmf_p_scores_ssdist,
#                   bnmf_scores_ssdist, bnmf_loading_ssdist)) %>% 
#   pivot_longer(grep("pca|nmf|fa", colnames(.)),
#                names_to = "type",
#                values_to = "value") %>%
#   mutate(model = str_sub(type, 1, 6),
#          model = ifelse(str_detect(model, 'pca'), 'pca', model),
#          model = ifelse(str_detect(model, 'l2'), 'nmf_l2', model),
#          model = ifelse(str_detect(model, 'fa'), 'fa', model),
#          model = ifelse(str_detect(model, '_p_'), 'nmf_p', model),
#          model = ifelse(str_detect(model, 'bnmf'), 'bnmf', model),
#          type = ifelse(str_detect(type, 'scores'), 'scores', "loadings")) %>% 
#   ggplot(aes(x = value)) +
#   geom_histogram(bins = 100, alpha = 0.75) +
#   #geom_density(color = "grey") + 
#   facet_grid(model~type, scales = "free") + 
#   geom_vline(xintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
#   theme_bw() +
#   labs(title = "Distinct Sims")

############################################################################################
# Overlapping patterns

# over_out <- tibble()
# 
# for (i in 1:100) {
#   load(paste0("./HPC/out_", i, "_over.RDA"))
#   over_out[i,1:32] <- out_over
# }

# bnmf_over_out <- tibble(seed = NA, eh = NA, ewa = NA)
# 
# for (i in 1:100) {
#   eh <- readMat(here::here(paste0("./HPC/MATLABout/stand/eh_over_", i, "_stand.mat")))[[1]]
#   eh <- as_tibble(eh) %>% nest(eh = everything()) %>% mutate(seed = i)
#   ewa <- readMat(here::here(paste0("./HPC/MATLABout/stand/ewa_over_", i, "_stand.mat")))[[1]]
#   ewa <- as_tibble(ewa) %>% nest(ewa = everything()) %>% mutate(seed = i)
#   both <- full_join(eh, ewa, by = "seed")
#   bnmf_over_out <- rbind(bnmf_over_out, both) %>% drop_na(seed)
# }
# 
# over_out <- left_join(over_out, bnmf_over_out, by = "seed")

# over_out %>%
#   mutate(bnmf_rank = map(eh, nrow)) %>% 
#   dplyr::select(seed, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank) %>% 
#   unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank)) %>% 
#   pivot_longer(cols = fa_rank:bnmf_rank) %>% 
#   group_by(name, value) %>% 
#   summarise(n())

# Symmetric subspace distance
# over_out <- over_out %>% 
#   mutate(bnmf_scores_ssdist = map2(true_scores, ewa, symm_subspace_dist),
#          bnmf_loading_ssdist = map2(true_patterns, eh, symm_subspace_dist))
# 
# over_out %>% dplyr::select(seed, grep("_ssdist", colnames(.))) %>% 
#   unnest(cols = c(pca_rotation_ssdist, pca_scores_ssdist, fa_rotations_ssdist, 
#                   fa_scores_ssdist, nmf_l2_loading_ssdist, nmf_l2_scores_ssdist, 
#                   nmf_p_loading_ssdist, nmf_p_scores_ssdist,
#                   bnmf_scores_ssdist, bnmf_loading_ssdist)) %>% 
#   pivot_longer(grep("pca|nmf|fa", colnames(.)),
#                names_to = "type",
#                values_to = "value") %>%
#   mutate(model = str_sub(type, 1, 6),
#          model = ifelse(str_detect(model, 'pca'), 'pca', model),
#          model = ifelse(str_detect(model, 'l2'), 'nmf_l2', model),
#          model = ifelse(str_detect(model, 'fa'), 'fa', model),
#          model = ifelse(str_detect(model, '_p_'), 'nmf_p', model),
#          model = ifelse(str_detect(model, 'bnmf'), 'bnmf', model),
#          type = ifelse(str_detect(type, 'scores'), 'scores', "loadings")) %>% 
#   ggplot(aes(x = value)) +
#   geom_histogram(bins = 100, alpha = 0.75) +
#   #geom_density(color = "grey") + 
#   facet_grid(model~type, scales = "free") + 
#   geom_vline(xintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
#   theme_bw() +
#   labs(title = "Overlapping Sims")

#######################################################################################

# cor_out <- tibble()
 
# for (i in 1:100) {
#   load(paste0("./HPC/Rout/Standardized/out_", i, "_cor.RDA"))
#   cor_out[i,1:32] <- out_cor
# }

# bnmf_cor_out <- tibble(seed = NA, eh = NA, ewa = NA)
 
# for (i in 1:100) {
#   eh <- readMat(here::here(paste0("./HPC/MATLABout/stand/eh_cor_", i, "_stand.mat")))[[1]]
#   eh <- as_tibble(eh) %>% nest(eh = everything()) %>% mutate(seed = i)
#   ewa <- readMat(here::here(paste0("./HPC/MATLABout/stand/ewa_cor_", i, "_stand.mat")))[[1]]
#   ewa <- as_tibble(ewa) %>% nest(ewa = everything()) %>% mutate(seed = i)
#   both <- full_join(eh, ewa, by = "seed")
#   bnmf_cor_out <- rbind(bnmf_cor_out, both) %>% drop_na(seed)
# }

# cor_out <- left_join(cor_out, bnmf_cor_out, by = "seed")

# cor_out %>% 
#   mutate(bnmf_rank = map(eh, nrow)) %>% 
#   dplyr::select(seed, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank) %>% 
#   unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank)) %>% 
#   pivot_longer(cols = fa_rank:bnmf_rank) %>% 
#   group_by(name, value) %>% 
#   summarise(n())

# Symmetric subspace distance
# cor_out <- cor_out %>% 
#   mutate(bnmf_scores_ssdist = map2(true_scores, ewa, symm_subspace_dist),
#          bnmf_loading_ssdist = map2(true_patterns, eh, symm_subspace_dist))
# 
# cor_out %>% dplyr::select(seed, grep("_ssdist", colnames(.))) %>% 
#   unnest(cols = c(pca_rotation_ssdist, pca_scores_ssdist, fa_rotations_ssdist, 
#                   fa_scores_ssdist, nmf_l2_loading_ssdist, nmf_l2_scores_ssdist, 
#                   nmf_p_loading_ssdist, nmf_p_scores_ssdist,
#                   bnmf_scores_ssdist, bnmf_loading_ssdist)) %>% 
#   pivot_longer(grep("pca|nmf|fa", colnames(.)),
#                names_to = "type",
#                values_to = "value") %>%
#   mutate(model = str_sub(type, 1, 6),
#          model = ifelse(str_detect(model, 'pca'), 'pca', model),
#          model = ifelse(str_detect(model, 'l2'), 'nmf_l2', model),
#          model = ifelse(str_detect(model, 'fa'), 'fa', model),
#          model = ifelse(str_detect(model, '_p_'), 'nmf_p', model),
#          model = ifelse(str_detect(model, 'bnmf'), 'bnmf', model),
#          type = ifelse(str_detect(type, 'scores'), 'scores', "loadings")) %>% 
#   # rbind(., npbnmf_cor) %>% 
#   # filter(type == "scores") %>% 
#   ggplot(aes(x = value)) +
#   geom_histogram(bins = 100, alpha = 0.75) +
#   #geom_density(color = "grey") + 
#   facet_grid(model~type, scales = "free") + 
#   geom_vline(xintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
#   theme_bw() +
#   labs(title = "Correlated Sims")

#######################################################################################
# Add predicted values
# dist_out <- dist_out %>% 
#   mutate(bnmf_pred = map2(ewa, eh, function(ewa, eh) as.matrix(ewa) %*% as.matrix(eh)),
#          sim_type = "distinct")
# 
# over_out <- over_out %>% 
#   mutate(bnmf_pred = map2(ewa, eh, function(ewa, eh) as.matrix(ewa) %*% as.matrix(eh)),
#          sim_type = "overlapping")
# 
# cor_out <- cor_out %>% 
#   mutate(bnmf_pred = map2(ewa, eh, function(ewa, eh) as.matrix(ewa) %*% as.matrix(eh)),
#          sim_type = "correlated")

#######################################################################################
# Aggregate distinct, overlapping, and correlated sims

# save(dist_out, file = "./HPC/Rout/dist_output_all.RDA")
load("./HPC/Rout/dist_output_all.RDA")

# save(over_out, file = "./HPC/Rout/over_output_all.RDA")
load("./HPC/Rout/over_output_all.RDA")

# save(cor_out, file = "./HPC/Rout/cor_output_all.RDA")
load("./HPC/Rout/cor_output_all.RDA")

load("./R/Sims/Iterate/sim_over.RDA")
load("./R/Sims/Iterate/sim_cor.RDA")
load("./R/Sims/Iterate/sim_dist.RDA")

dist <- left_join(dist_out, sim_dist)
corr <- left_join(cor_out, sim_cor)
over <- left_join(over_out, sim_over)

all_stand <- rbind(dist, over, corr)

# save(all_stand, file = "./HPC/Rout/all_stand.RDA")
load("./HPC/Rout/all_stand.RDA")
all_stand

####################################################################################################

# Relative Error

all_stand %>% dplyr::select(seed, sim_type, sim, grep("pred", colnames(.))) %>% 
  pivot_longer(pca_pred:bnmf_pred) %>% 
  mutate(l2 = map2(sim, value, function (x,y) norm(x-y, "F")/norm(x, "F"))) %>% 
  dplyr::select(seed, sim_type, name, l2) %>% 
  unnest(l2) %>% 
  group_by(sim_type, name) %>% 
  summarize(min = min(l2),
            mean = mean(l2),
            max = max(l2)) %>% knitr::kable()

all_stand_error <- all_stand %>% dplyr::select(seed, sim_type, sim, chem, grep("pred", colnames(.))) %>% 
  pivot_longer(pca_pred:bnmf_pred) %>% 
  mutate(l2_true = map2(chem, value, function (x,y) norm(x-y, "F")/norm(x, "F")),
         l2_sim = map2(sim, value, function (x,y) norm(x-y, "F")/norm(x, "F"))) %>% 
  dplyr::select(seed, sim_type, name, l2_true, l2_sim) %>% 
  unnest(c(l2_sim, l2_true)) %>% 
  mutate(Standardized = "Yes")
# save(all_stand_error, file = "./HPC/Rout/all_stand_error.RDA")

all_stand_error %>% 
  ggplot(aes(x = name, y = l2, color = name, fill = name)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(~sim_type, scales = "free") + 
  geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  theme_bw() +
  scale_y_log10(limits = c(0.07,1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "Relative Predictive Error")

all_stand_error %>% 
  ggplot(aes(x = l2, color = sim_type, fill = sim_type)) +
  geom_histogram(bins = 100, alpha = 0.75) +
  facet_grid(name~., scales = "free") + 
  theme_bw() + scale_x_log10() +
  labs(title = "Relative Predictive Error") +
  theme(legend.position = "bottom")

####################################################################################################
# Symmetric Subspace Distance

all_stand_wnorm <- all_stand %>% 
  mutate(true_norm = map2(true_patterns, true_scores, normalize_patterns),
         pca_norm = map2(pca_rotations, pca_scores, function(x,y) normalize_patterns(t(x), y)),
         fa_norm = map2(fa_rotations, fa_scores, function(x,y) normalize_patterns(t(x), y)),
         nmf_l2_norm = map2(nmf_l2_loadings, nmf_l2_scores, normalize_patterns),
         nmf_p_norm = map2(nmf_p_loadings, nmf_p_scores, normalize_patterns),
         npbnmf_norm = map2(eh, ewa, normalize_patterns)) %>% 
  unnest_wider(true_norm) %>% rename(true_score_norm = score_norm, true_loading_norm = loading_norm) %>% 
  unnest_wider(pca_norm) %>% rename(pca_score_norm = score_norm, pca_loading_norm = loading_norm) %>% 
  unnest_wider(fa_norm) %>% rename(fa_score_norm = score_norm, fa_loading_norm = loading_norm) %>% 
  unnest_wider(nmf_l2_norm) %>% rename(nmf_l2_score_norm = score_norm, nmf_l2_loading_norm = loading_norm) %>% 
  unnest_wider(nmf_p_norm) %>% rename(nmf_p_score_norm = score_norm, nmf_p_loading_norm = loading_norm) %>% 
  unnest_wider(npbnmf_norm) %>% rename(npbnmf_score_norm = score_norm, npbnmf_loading_norm = loading_norm) %>% 
  mutate(bnmf_norm_scores_ssdist = map2(true_score_norm, npbnmf_score_norm, symm_subspace_dist),
         bnmf_norm_loading_ssdist = map2(true_loading_norm, npbnmf_loading_norm, symm_subspace_dist),
         pca_norm_scores_ssdist = map2(true_score_norm, pca_score_norm, symm_subspace_dist),
         pca_norm_loading_ssdist = map2(true_loading_norm, pca_loading_norm, symm_subspace_dist),
         fa_norm_scores_ssdist = map2(true_score_norm, fa_score_norm, symm_subspace_dist),
         fa_norm_loading_ssdist = map2(true_loading_norm, fa_loading_norm, symm_subspace_dist),
         nmf_l2_norm_scores_ssdist = map2(true_score_norm, nmf_l2_score_norm, symm_subspace_dist),
         nmf_l2_norm_loading_ssdist = map2(true_loading_norm, nmf_l2_loading_norm, symm_subspace_dist),
         nmf_p_norm_scores_ssdist = map2(true_score_norm, nmf_p_score_norm, symm_subspace_dist),
         nmf_p_norm_loading_ssdist = map2(true_loading_norm, nmf_p_loading_norm, symm_subspace_dist))

all_stand_ssdist <- all_stand_wnorm %>% dplyr::select(seed, sim_type, grep("_ssdist", colnames(.))) %>% 
  unnest(cols = grep("_ssdist", colnames(.))) %>% 
  pivot_longer(grep("_ssdist", colnames(.)),
               names_to = "type",
               values_to = "value") %>%
  mutate(Normalized = ifelse(grepl("norm", type), "Yes", "No"),
         model = str_sub(type, 1, 6),
         model = ifelse(str_detect(model, 'pca'), 'pca', model),
         model = ifelse(str_detect(model, 'l2'), 'nmf_l2', model),
         model = ifelse(str_detect(model, 'fa'), 'fa', model),
         model = ifelse(str_detect(model, '_p_'), 'nmf_p', model),
         model = ifelse(str_detect(model, 'bnmf'), 'bnmf', model),
         type = ifelse(str_detect(type, 'scores'), 'Scores', "Loadings")) %>% 
  mutate(Standardized = "Yes")

all_stand_ssdist %>% 
  ggplot(aes(x = value)) +
  geom_histogram(aes(fill = sim_type), bins = 100, alpha = 0.75) +
  facet_grid(type ~ model) + 
  geom_vline(xintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  theme_bw() +
  labs(title = "Symmetric Subspace Distance")

all_stand_ssdist %>% 
  ggplot(aes(x = model, y = value, color = model)) +
  geom_boxplot(varwidth = TRUE) +
  facet_grid(sim_type~type) + 
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") + 
  labs(title = "Symmetric Subspace Distance")

# save(all_stand_ssdist, file = "./HPC/Rout/all_stand_ssdist.RDA")
