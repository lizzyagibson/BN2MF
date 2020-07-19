#######################################################################################
#######################################################################################
#######################################################################################

# Aggregate HPC output

# Distinct patterns

# dist_out_un <- tibble()

# for (i in 1:100) {
#   load(paste0("./HPC/Rout/out_", i, "_dist_un.RDA"))
#   dist_out_un[i,1:34] <- out_dist_un # 34 columns
# }

# bnmf_dist_out_un <- tibble(seed = NA, eh = NA, ewa = NA)

# for (i in 1:100) {
#   eh <- readMat(here::here(paste0("./HPC/MATLABout/un/eh_dist_", i, "_un.mat")))[[1]]
#   eh <- as_tibble(eh) %>% nest(eh = everything()) %>% mutate(seed = i)
#   ewa <- readMat(here::here(paste0("./HPC/MATLABout/un/ewa_dist_", i, "_un.mat")))[[1]]
#   ewa <- as_tibble(ewa) %>% nest(ewa = everything()) %>% mutate(seed = i)
#   both <- full_join(eh, ewa, by = "seed")
#   bnmf_dist_out_un <- rbind(bnmf_dist_out_un, both) %>% drop_na(seed)
# }

# dist_out_un <- left_join(dist_out_un, bnmf_dist_out_un, by = "seed")

# save(dist_out_un, file = "./HPC/dist_output_all_un.RDA")
load("./HPC/Rout/dist_output_all_un.RDA")

dist_out_un %>%
  mutate(bnmf_rank = map(eh, nrow)) %>% 
  dplyr::select(seed, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank) %>% 
  unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank)) %>% 
  pivot_longer(cols = fa_rank:bnmf_rank) %>% 
  group_by(name, value) %>% 
  summarise(n())

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

dist_out_un <- dist_out_un %>% 
  mutate(bnmf_scores_ssdist = map2(scores, ewa, symm_subspace_dist),
         bnmf_loading_ssdist = map2(patterns, eh, symm_subspace_dist))

############################################################################################
# Overlapping patterns

# over_out_un <- tibble()

# for (i in 1:100) {
#   load(paste0("./HPC/Rout/out_", i, "_over_un.RDA"))
#   over_out_un[i,1:34] <- out_over_un # 34 columns
# }

# bnmf_over_out_un <- tibble(seed = NA, eh = NA, ewa = NA)

# for (i in 1:100) {
#   eh <- readMat(here::here(paste0("./HPC/MATLABout/un/eh_over_", i, "_un.mat")))[[1]]
#   eh <- as_tibble(eh) %>% nest(eh = everything()) %>% mutate(seed = i)
#   ewa <- readMat(here::here(paste0("./HPC/MATLABout/un/ewa_over_", i, "_un.mat")))[[1]]
#   ewa <- as_tibble(ewa) %>% nest(ewa = everything()) %>% mutate(seed = i)
#   both <- full_join(eh, ewa, by = "seed")
#   bnmf_over_out_un <- rbind(bnmf_over_out_un, both) %>% drop_na(seed)
# }
 
# over_out_un <- left_join(over_out_un, bnmf_over_out_un, by = "seed")

# save(over_out_un, file = "./HPC/over_output_all_un.RDA")
load("./HPC/over_output_all_un.RDA")

over_out_un %>%
  mutate(bnmf_rank = map(eh, nrow)) %>% 
  dplyr::select(seed, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank) %>% 
  unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank)) %>% 
  pivot_longer(cols = fa_rank:bnmf_rank) %>% 
  group_by(name, value) %>% 
  summarise(n())

# Symmetric subspace distance
over_out_un <- over_out_un %>% 
  mutate(bnmf_scores_ssdist = map2(scores, ewa, symm_subspace_dist),
         bnmf_loading_ssdist = map2(patterns, eh, symm_subspace_dist))

#######################################################################################

# cor_out_un <- tibble()

# for (i in 1:100) {
#   load(paste0("./HPC/Rout/out_", i, "_cor_un.RDA"))
#   cor_out_un[i,1:34] <- out_cor_un
# }

# bnmf_cor_out_un <- tibble(seed = NA, eh = NA, ewa = NA)

# for (i in 1:100) {
#   eh <- readMat(here::here(paste0("./HPC/MATLABout/un/eh_cor_", i, "_un.mat")))[[1]]
#   eh <- as_tibble(eh) %>% nest(eh = everything()) %>% mutate(seed = i)
#   ewa <- readMat(here::here(paste0("./HPC/MATLABout/un/ewa_cor_", i, "_un.mat")))[[1]]
#   ewa <- as_tibble(ewa) %>% nest(ewa = everything()) %>% mutate(seed = i)
#   both <- full_join(eh, ewa, by = "seed")
#   bnmf_cor_out_un <- rbind(bnmf_cor_out_un, both) %>% drop_na(seed)
# }
 
# cor_out_un <- left_join(cor_out_un, bnmf_cor_out_un, by = "seed")

# save(cor_out_un, file = "./HPC/cor_output_all_un.RDA")
load("./HPC/cor_output_all_un.RDA")

cor_out_un %>% 
  mutate(bnmf_rank = map(eh, nrow)) %>% 
  dplyr::select(seed, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank) %>% 
  unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank)) %>% 
  pivot_longer(cols = fa_rank:bnmf_rank) %>% 
  group_by(name, value) %>% 
  summarise(n())

# Symmetric subspace distance
cor_out_un <- cor_out_un %>% 
  mutate(bnmf_scores_ssdist = map2(scores, ewa, symm_subspace_dist),
         bnmf_loading_ssdist = map2(patterns, eh, symm_subspace_dist))

#######################################################################################
# Add predicted values
dist_out_un <- dist_out_un %>% 
  mutate(bnmf_pred = map2(ewa, eh, function(ewa, eh) as.matrix(ewa) %*% as.matrix(eh)),
         sim_type = "distinct")

over_out_un <- over_out_un %>% 
  mutate(bnmf_pred = map2(ewa, eh, function(ewa, eh) as.matrix(ewa) %*% as.matrix(eh)),
         sim_type = "overlapping")

cor_out_un <- cor_out_un %>% 
  mutate(bnmf_pred = map2(ewa, eh, function(ewa, eh) as.matrix(ewa) %*% as.matrix(eh)),
         sim_type = "correlated")

#######################################################################################
# Aggregate distinct, overlapping, and correlated sims

all_unstand <- rbind(dist_out_un, over_out_un, cor_out_un)
# save(all_unstand, file = "./HPC/all_unstand.RDA")
all_unstand

####################################################################################################
# Relative Error

all_unstand %>% dplyr::select(seed, sim_type, sim, grep("pred", colnames(.))) %>% 
  pivot_longer(pca_pred:bnmf_pred) %>% 
  mutate(l2 = map2(sim, value, function (x,y) norm(x-y, "F")/norm(x, "F"))) %>% 
  dplyr::select(seed, sim_type, name, l2) %>% 
  unnest(l2) %>% 
  group_by(sim_type, name) %>% 
  summarize(min = min(l2),
            mean = mean(l2),
            max = max(l2))

all_error <- all_unstand %>% dplyr::select(seed, sim_type, sim, grep("pred", colnames(.))) %>% 
  pivot_longer(pca_pred:bnmf_pred) %>% 
  mutate(l2 = map2(sim, value, function (x,y) norm(x-y, "F")/norm(x, "F"))) %>% 
  dplyr::select(seed, sim_type, name, l2) %>% 
  unnest(l2)

all_error %>% 
  ggplot(aes(x = name, y = l2)) +
  geom_boxplot() +
  facet_grid(~sim_type, scales = "free") + 
  geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "Relative Predictive Error")

all_error %>% 
  ggplot(aes(x = l2)) +
  geom_histogram(bins = 100, alpha = 0.75) +
  facet_grid(name~sim_type, scales = "free") + 
  geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  theme_bw() +
  labs(title = "Relative Predictive Error")

####################################################################################################
# Symmetric Subspace Distance

all_unstand %>% dplyr::select(seed, sim_type, grep("_ssdist", colnames(.))) %>% 
  unnest(cols = c(pca_rotation_ssdist, pca_scores_ssdist, fa_rotations_ssdist, 
                  fa_scores_ssdist, nmf_l2_loading_ssdist, nmf_l2_scores_ssdist, 
                  nmf_p_loading_ssdist, nmf_p_scores_ssdist,
                  bnmf_scores_ssdist, bnmf_loading_ssdist)) %>% 
  summary()

all_ssdist <- all_unstand %>% dplyr::select(seed, sim_type, grep("_ssdist", colnames(.))) %>% 
  unnest(cols = c(pca_rotation_ssdist, pca_scores_ssdist, fa_rotations_ssdist, 
                  fa_scores_ssdist, nmf_l2_loading_ssdist, nmf_l2_scores_ssdist, 
                  nmf_p_loading_ssdist, nmf_p_scores_ssdist,
                  bnmf_scores_ssdist, bnmf_loading_ssdist)) %>% 
  pivot_longer(grep("_ssdist", colnames(.)),
               names_to = "type",
               values_to = "value") %>%
  mutate(model = str_sub(type, 1, 6),
         model = ifelse(str_detect(model, 'pca'), 'pca', model),
         model = ifelse(str_detect(model, 'l2'), 'nmf_l2', model),
         model = ifelse(str_detect(model, 'fa'), 'fa', model),
         model = ifelse(str_detect(model, '_p_'), 'nmf_p', model),
         model = ifelse(str_detect(model, 'bnmf'), 'bnmf', model),
         type = ifelse(str_detect(type, 'scores'), 'scores', "loadings")) 

all_ssdist %>% 
  ggplot(aes(x = value)) +
  geom_histogram(aes(fill = sim_type), bins = 100, alpha = 0.75) +
  facet_grid(type ~ model) + 
  geom_vline(xintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  theme_bw() +
  labs(title = "Symmetric Subspace Distance")

all_ssdist %>% 
  ggplot(aes(x = model, y = value, color = model)) +
  geom_boxplot() +
  facet_grid(sim_type~type) + 
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") + 
  labs(title = "Symmetric Subspace Distance")


