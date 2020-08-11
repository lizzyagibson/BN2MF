# See extract_from_HPC for full code.
# Aggregate npBNMF results
np_dist <- tibble(seed = NA, eh = NA, ewa = NA)

for (i in 1:100) {
  eh <- readMat(here::here(paste0("./HPC/MATLABout/stand/eh_dist_", i, "_stand.mat")))[[1]]
  eh <- as_tibble(eh) %>% nest(eh = everything()) %>% mutate(seed = i)
  ewa <- readMat(here::here(paste0("./HPC/MATLABout/stand/ewa_dist_", i, "_stand.mat")))[[1]]
  ewa <- as_tibble(ewa) %>% nest(ewa = everything()) %>% mutate(seed = i)
  both <- full_join(eh, ewa, by = "seed")
  np_dist <- rbind(np_dist, both) %>% drop_na(seed)
}

np_dist %>% 
  mutate(rank = map(eh, nrow)) %>% 
  unnest(rank) %>% pull(rank) %>% 
  table()

##################################

np_over <- tibble(seed = NA, eh = NA, ewa = NA)

for (i in 1:100) {
  eh <- readMat(here::here(paste0("./HPC/MATLABout/stand/eh_over_", i, "_stand.mat")))[[1]]
  eh <- as_tibble(eh) %>% nest(eh = everything()) %>% mutate(seed = i)
  ewa <- readMat(here::here(paste0("./HPC/MATLABout/stand/ewa_over_", i, "_stand.mat")))[[1]]
  ewa <- as_tibble(ewa) %>% nest(ewa = everything()) %>% mutate(seed = i)
  both <- full_join(eh, ewa, by = "seed")
  np_over <- rbind(np_over, both) %>% drop_na(seed)
}

np_over %>% 
  mutate(rank = map(eh, nrow)) %>% 
  unnest(rank) %>% pull(rank) %>% 
  table()

#############################

np_cor <- tibble(seed = NA, eh = NA, ewa = NA)

for (i in 1:100) {
  eh <- readMat(here::here(paste0("./HPC/MATLABout/stand/eh_cor_", i, "_stand.mat")))[[1]]
  eh <- as_tibble(eh) %>% nest(eh = everything()) %>% mutate(seed = i)
  ewa <- readMat(here::here(paste0("./HPC/MATLABout/stand/ewa_cor_", i, "_stand.mat")))[[1]]
  ewa <- as_tibble(ewa) %>% nest(ewa = everything()) %>% mutate(seed = i)
  both <- full_join(eh, ewa, by = "seed")
  np_cor <- rbind(np_cor, both) %>% drop_na(seed)
}

np_cor %>% 
  mutate(rank = map(eh, nrow)) %>% 
  unnest(rank) %>% pull(rank) %>% 
  table()

#############################

npbnmf_dist <- dist %>%
  mutate(number_patterns = map(eh_reordered, nrow),
         symm_ratio = map2(true_scores, ewa_reordered, symm_subspace_dist)) %>%
  dplyr::select(seed, number_patterns, symm_ratio) %>%
  unnest(cols = c(number_patterns, symm_ratio)) %>%
  rename(value = symm_ratio) %>%
  mutate(model = "npBNMF",
         type = "scores") %>% dplyr::select(-number_patterns)

# npbnmf_cor <- npb_scores_ordered_cor %>% 
#   mutate(number_patterns = map(eh_reordered, nrow),
#          symm_ratio = map2(true_scores, ewa_reordered, symm_subspace_dist)) %>% 
#   dplyr::select(seed, number_patterns, symm_ratio) %>% 
#   unnest(cols = c(number_patterns, symm_ratio)) %>% 
#   rename(value = symm_ratio) %>% 
#   mutate(model = "npBNMF",
#          type = "scores") %>% dplyr::select(-number_patterns)

# npbnmf_over <- npb_scores_ordered_over %>% 
#   mutate(number_patterns = map(eh_reordered, nrow),
#          symm_ratio = map2(true_scores, ewa_reordered, symm_subspace_dist)) %>% 
#   dplyr::select(seed, number_patterns, symm_ratio) %>% 
#   unnest(cols = c(number_patterns, symm_ratio)) %>% 
#   rename(value = symm_ratio) %>% 
#   mutate(model = "npBNMF",
#          type = "scores") %>% dplyr::select(-number_patterns)

## symmetric subspace distance
over_out %>% dplyr::select(seed, grep("_ssdist", colnames(.))) %>% 
  unnest(cols = c(pca_rotation_ssdist, pca_scores_ssdist, fa_rotations_ssdist, 
                  fa_scores_ssdist, nmf_l2_loading_ssdist, nmf_l2_scores_ssdist, 
                  nmf_p_loading_ssdist, nmf_p_scores_ssdist)) %>% 
  pivot_longer(grep("pca|nmf|fa", colnames(.)),
               names_to = "type",
               values_to = "value") %>%
  mutate(model = str_sub(type, 1, 6),
         model = ifelse(str_detect(model, 'pca'), 'pca', model),
         model = ifelse(str_detect(model, 'l2'), 'nmf_l2', model),
         model = ifelse(str_detect(model, 'fa'), 'fa', model),
         model = ifelse(str_detect(model, '_p_'), 'nmf_p', model),
         type = ifelse(str_detect(type, 'scores'), 'scores', "loadings")) %>% 
  # rbind(., npbnmf_over) %>% 
  # filter(type == "scores") %>% 
  ggplot(aes(x = value)) +
  geom_histogram(bins = 100, alpha = 0.75) +
  #geom_density(color = "grey") + 
  facet_grid(model~type, scales = "free") + 
  geom_vline(xintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  theme_bw() +
  labs(title = "Overlapping Sims")

## symmetric subspace distance
dist_out %>% dplyr::select(seed, grep("_ssdist", colnames(.))) %>% 
  unnest(cols = c(pca_rotation_ssdist, pca_scores_ssdist, fa_rotations_ssdist, 
                  fa_scores_ssdist, nmf_l2_loading_ssdist, nmf_l2_scores_ssdist, 
                  nmf_p_loading_ssdist, nmf_p_scores_ssdist)) %>% 
  pivot_longer(grep("pca|nmf|fa", colnames(.)),
               names_to = "type",
               values_to = "value") %>%
  mutate(model = str_sub(type, 1, 6),
         model = ifelse(str_detect(model, 'pca'), 'pca', model),
         model = ifelse(str_detect(model, 'l2'), 'nmf_l2', model),
         model = ifelse(str_detect(model, 'fa'), 'fa', model),
         model = ifelse(str_detect(model, '_p_'), 'nmf_p', model),
         type = ifelse(str_detect(type, 'scores'), 'scores', "loadings")) %>% 
  # rbind(., npbnmf_dist) %>% 
  # filter(type == "scores") %>% 
  ggplot(aes(x = value)) +
  geom_histogram(bins = 100, alpha = 0.75) +
  #geom_density(color = "grey") + 
  facet_grid(model~type, scales = "free") + 
  geom_vline(xintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  theme_bw() +
  labs(title = "Distinct Sims")

## symmetric subspace distance
cor_out %>% dplyr::select(seed, grep("_ssdist", colnames(.))) %>% 
  unnest(cols = c(pca_rotation_ssdist, pca_scores_ssdist, fa_rotations_ssdist, 
                  fa_scores_ssdist, nmf_l2_loading_ssdist, nmf_l2_scores_ssdist, 
                  nmf_p_loading_ssdist, nmf_p_scores_ssdist)) %>% 
  pivot_longer(grep("pca|nmf|fa", colnames(.)),
               names_to = "type",
               values_to = "value") %>%
  mutate(model = str_sub(type, 1, 6),
         model = ifelse(str_detect(model, 'pca'), 'pca', model),
         model = ifelse(str_detect(model, 'l2'), 'nmf_l2', model),
         model = ifelse(str_detect(model, 'fa'), 'fa', model),
         model = ifelse(str_detect(model, '_p_'), 'nmf_p', model),
         type = ifelse(str_detect(type, 'scores'), 'scores', "loadings")) %>% 
  # rbind(., npbnmf_cor) %>% 
  # filter(type == "scores") %>% 
  ggplot(aes(x = value)) +
  geom_histogram(bins = 100, alpha = 0.75) +
  #geom_density(color = "grey") + 
  facet_grid(model~type, scales = "free") + 
  geom_vline(xintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  theme_bw() +
  labs(title = "Correlated Sims")
