dgp_local <- dgp_local %>% mutate(row = 1:nrow(.))

results <- 
  dgp_local %>% 
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
  
truth <- dgp_local %>% 
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

new_norm <- 
  left_join(results, truth) %>% 
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

dgp_matrix_norms <- 
  dgp_local %>% 
  dplyr::select(seed, data, row, true_patterns, true_scores, pca_rotations, pca_scores,
                fa_rotations, fa_scores, nmf_l2_loadings, nmf_l2_scores,
                nmf_p_loadings, nmf_p_scores, eh, ewa) %>% 
  mutate(true_patterns = map(true_patterns, t)) %>% 
  pivot_longer(starts_with("true"),
               names_to = "matrix",
               values_to = "truth") %>% 
  pivot_longer(pca_rotations:ewa,
               names_to = "results") %>% 
  drop_na(value) %>% 
  filter( (grepl("rotations|loadings|eh", results) & matrix == "true_patterns") |
          (grepl("scores|ewa", results) & matrix == "true_scores") ) %>% 
  mutate(value = case_when(grepl("load|eh", results) ~ map(value, t),
                           TRUE ~ map(value, as.matrix))) %>% 
  mutate(rank = map(value, ncol)) %>% 
  unnest(rank) %>% 
  filter(rank == 4) %>% 
  mutate(l2 = map2(truth, value, function (x,y) norm(x-y, "F")/norm(x, "F")),
         l1 = map2(truth, value, function (x,y) sum(abs(x-y))/sum(abs(x)))) %>% 
  unnest(c(l1, l2)) %>% 
  dplyr::select(seed, data, row, matrix, results, l2, l1) 

#####
#
#####
pca_4 <- dgp_local %>% 
  dplyr::select(seed, data, row, true_patterns, true_scores, 
                pca_rotations, pca_scores, pca_rank) %>% 
  mutate(true_patterns = map(true_patterns, t),
         model = "PCA") %>% 
  rename(pred_load = 6, pred_score = 7, rank = 8)

fa_4 <- dgp_local %>% 
  dplyr::select(seed, data, row, true_patterns, true_scores, 
                fa_rotations, fa_scores, fa_rank) %>% 
  mutate(true_patterns = map(true_patterns, t),
         model = "FA") %>% 
  rename(pred_load = 6, pred_score = 7, rank = 8)

nmf_l2_4 <- dgp_local %>% 
  dplyr::select(seed, data, row, true_patterns, true_scores,
                nmf_l2_loadings, nmf_l2_scores, nmf_l2_rank) %>% 
  mutate(true_patterns = map(true_patterns, t),
         model = "NMF L2") %>% 
  rename(pred_load = 6, pred_score = 7, rank = 8)

nmf_p_4 <- dgp_local %>% 
  dplyr::select(seed, data, row, true_patterns, true_scores, 
                nmf_p_loadings, nmf_p_scores, nmf_p_rank) %>% 
  mutate(true_patterns = map(true_patterns, t),
         model = "NMF P") %>% 
  rename(pred_load = 6, pred_score = 7, rank = 8)

bnmf_4 <- dgp_local %>% 
  dplyr::select(seed, data, row, true_patterns, true_scores, 
                eh, ewa, bnmf_rank) %>% 
  mutate(true_patterns = map(true_patterns, t),
         model = "BN2MF") %>% 
  rename(pred_load = 6, pred_score = 7, rank = 8)

norm_plot <- rbind(pca_4, fa_4, nmf_l2_4, nmf_p_4, bnmf_4) %>% 
  unnest(rank) %>% 
  filter(rank == 4) %>%
  mutate(pred_load = map(pred_load, function(x) if (nrow(x) != 50) {t(x)} else{x}),
         norm_load = map2(true_patterns, pred_load, 
                          function(x,y) norm(as.matrix(x)-as.matrix(y), "F")/norm(as.matrix(x), "F")),
         pred_score = map(pred_score, function(x) if (nrow(x) != 1000) {t(x)} else{x}),
         norm_score = map2(true_scores, pred_score, 
                           function(x,y) norm(as.matrix(x)-as.matrix(y), "F")/norm(as.matrix(x), "F"))) %>% 
  unnest(c(norm_load, norm_score)) %>% 
  dplyr::select(seed, data, row, model, norm_load, norm_score) %>% 
  pivot_longer(norm_load:norm_score)

#####
# Viz
#####
norm_plot %>% 
  arrange(seed, data, name)
# Right

dgp_matrix_norms
# RIGHT

# TRUTH
fa_load <- dgp_local$fa_rotations[[1]]
fa_scor <- dgp_local$fa_scores[[1]]  
true_load <- dgp_local$true_patterns[[1]]
true_score <- dgp_local$true_scores[[1]]

norm(true_load-t(fa_load), "F")/norm(true_load, "F")
norm(true_score-fa_scor, "F")/norm(true_score, "F")
