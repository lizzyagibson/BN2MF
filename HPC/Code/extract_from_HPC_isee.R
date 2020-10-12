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
# load(paste0("./HPC/Rout/Unstandardized/out_", 1, "_dist_un.RDA"))
# out_dist_un

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

# dist_out_un %>%
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

# dist_out_un <- dist_out_un %>% 
#   mutate(bnmf_scores_ssdist = map2(scores, ewa, symm_subspace_dist),
#          bnmf_loading_ssdist = map2(patterns, eh, symm_subspace_dist))

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

# over_out_un %>%
#   mutate(bnmf_rank = map(eh, nrow)) %>% 
#   dplyr::select(seed, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank) %>% 
#   unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank)) %>% 
#   pivot_longer(cols = fa_rank:bnmf_rank) %>% 
#   group_by(name, value) %>% 
#   summarise(n())

# Symmetric subspace distance
# over_out_un <- over_out_un %>% 
#   mutate(bnmf_scores_ssdist = map2(scores, ewa, symm_subspace_dist),
#          bnmf_loading_ssdist = map2(patterns, eh, symm_subspace_dist))

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

# cor_out_un %>% 
#   mutate(bnmf_rank = map(eh, nrow)) %>% 
#   dplyr::select(seed, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank) %>% 
#   unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank)) %>% 
#   pivot_longer(cols = fa_rank:bnmf_rank) %>% 
#   group_by(name, value) %>% 
#   summarise(n())

# Symmetric subspace distance
# cor_out_un <- cor_out_un %>% 
#   mutate(bnmf_scores_ssdist = map2(scores, ewa, symm_subspace_dist),
#          bnmf_loading_ssdist = map2(patterns, eh, symm_subspace_dist))

#######################################################################################
# Add predicted values
# dist_out_un <- dist_out_un %>% 
#   mutate(bnmf_pred = map2(ewa, eh, function(ewa, eh) as.matrix(ewa) %*% as.matrix(eh)),
#          sim_type = "distinct")
# 
# over_out_un <- over_out_un %>% 
#   mutate(bnmf_pred = map2(ewa, eh, function(ewa, eh) as.matrix(ewa) %*% as.matrix(eh)),
#          sim_type = "overlapping")
# 
# cor_out_un <- cor_out_un %>% 
#   mutate(bnmf_pred = map2(ewa, eh, function(ewa, eh) as.matrix(ewa) %*% as.matrix(eh)),
#          sim_type = "correlated")

#######################################################################################
# Aggregate distinct, overlapping, and correlated sims
# save(dist_out_un, file = "./HPC/dist_output_all_un.RDA")
load("./HPC/Rout/dist_output_all_un.RDA")

# save(over_out_un, file = "./HPC/over_output_all_un.RDA")
load("./HPC/over_output_all_un.RDA")

# save(cor_out_un, file = "./HPC/cor_output_all_un.RDA")
load("./HPC/Rout/cor_output_all_un.RDA")
cor_out_un

# all_unstand <- rbind(dist_out_un, over_out_un, cor_out_un)
# save(all_unstand, file = "./HPC/all_unstand.RDA")
load("./Sims/all_unstand.RDA")
all_unstand

####################################################################################################
# Relative Error

all_error <- all_unstand %>% dplyr::select(seed, sim_type, chem, sim, grep("pred", colnames(.))) %>% 
  pivot_longer(pca_pred:bnmf_pred) %>% 
  mutate(l2_true = map2(chem, value, function (x,y) norm(x-y, "F")/norm(x, "F")),
         l2_sim = map2(sim, value, function (x,y) norm(x-y, "F")/norm(x, "F"))) %>% 
  dplyr::select(seed, sim_type, name, l2_sim, l2_true) %>% 
  unnest(c(l2_sim, l2_true)) %>% 
  mutate(Standardized = "No")
# save(all_unstand_error, file = "./HPC/Rout/all_unstand_error.RDA")

all_error %>% 
  ggplot(aes(x = name, y = l2_sim)) +
  geom_boxplot() +
  facet_grid(Standardized~sim_type) + 
  geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  theme_bw() + scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "Relative Predictive Error")

####################################################################################################
# Symmetric Subspace Distance

all_ssdist <- all_unstand %>% dplyr::select(seed, sim_type, grep("_ssdist", colnames(.))) %>% 
  unnest(cols = grep("_ssdist", colnames(.))) %>% 
  pivot_longer(grep("_ssdist", colnames(.)),
               names_to = "type",
               values_to = "value") %>%
  mutate(Normalized = ifelse(grepl("norm", type), "Yes", "No"),
         model = str_sub(type, 1, 6),
         model = ifelse(str_detect(model, 'pca'), 'PCA', model),
         model = ifelse(str_detect(model, 'l2'), 'NMF_l2', model),
         model = ifelse(str_detect(model, 'fa'), 'FA', model),
         model = ifelse(str_detect(model, '_p_'), 'NMF_p', model),
         model = ifelse(str_detect(model, 'bnmf'), 'npBNMF', model),
         type = ifelse(str_detect(type, 'scores'), 'Scores', "Loadings")) %>% 
  mutate(Standardized = "No")

# save(all_unstand_ssdist, file = "./HPC/Rout/all_unstand_ssdist.RDA")

all_ssdist %>% 
  ggplot(aes(x = model, y = value, color = Normalized)) +
  geom_boxplot() +
  facet_grid(sim_type~type) + 
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(title = "Symmetric Subspace Distance")

######################################################
######################################################
# Plots for ISEE
######################################################

plot_l2 <- all_error %>% 
  filter(name != "nmf_p_pred") %>% 
  mutate(sim_type = paste0(str_to_title(sim_type), " Sims")) %>% 
  mutate(name = case_when(name == "pca_pred" ~ "PCA",
                          name == "fa_pred" ~ "Factor Analysis",
                          name == "nmf_l2_pred" ~ "NMF",
                          name == "bnmf_pred" ~ "npBNMF")) %>% 
  mutate(sim_type = fct_relevel(sim_type, "Distinct Sims",
                                          "Overlapping Sims",
                                          "Correlated Sims")) %>% 
  mutate(name = fct_relevel(name, "npBNMF",
                             "PCA",
                             "Factor Analysis",
                             "NMF"))

#pdf("Figures/isee_2020_l2.pdf", width = 10)
plot_l2 %>% 
  ggplot(aes(x = name, y = l2_sim, color = name, fill = name)) +
  geom_boxplot(outlier.size = 0.25, alpha = 0.5) +
  facet_wrap(~sim_type) + 
  theme_bw(base_size = 30) + 
  scale_y_log10() +
  theme(strip.text = element_text(size = 20),
        axis.text.x = element_blank(),
        strip.background = element_rect(fill="white"),
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.1), # c(0,0) bottom left, c(1,1) top-right.
        axis.title.y = element_text(size = 30),
        legend.text = element_text(size = 25)) + 
  labs(title = "Distance from True Observations",
       y = "Relative Prediction Error", x = "",
       color = "", fill = "") + scale_color_nejm() + scale_fill_nejm()
#dev.off()

###
plot_ss <- 
  all_ssdist %>% 
  filter(Normalized == "No" & Standardized == "No") %>% 
  filter(model != "NMF_p")

#pdf("Figures/isee_2020_score.pdf", width = 10)
plot_ss %>% 
  filter(type == "Scores") %>% 
  ggplot(aes(x = model, y = value, color = model, fill = model)) +
  geom_boxplot(outlier.size = 0.25, alpha = 0.5) +
  facet_wrap(~sim_type) + 
  theme_bw(base_size = 30) +
  theme(strip.text = element_text(size = 20),
        axis.text.x = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.1), # c(0,0) bottom left, c(1,1) top-right.
        axis.title.y = element_text(size = 30),
        legend.text = element_text(size = 25)) + 
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed", size = 1.05) +
  labs(y = "Symmetric Subspace Distance",
       x = "",
       title = "Distance from True Individual Scores",
       color = "", fill = "") + scale_color_nejm() + scale_fill_nejm()
#dev.off()

#pdf("Figures/isee_2020_load.pdf", width = 10)
plot_ss %>% 
  filter(type == "Loadings") %>% 
  ggplot(aes(x = model, y = value, color = model, fill = model)) +
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed", size = 1.05) +
  geom_boxplot(outlier.size = 0.25, alpha = 0.5) +
  facet_wrap(~sim_type) + 
  theme_bw(#base_size = 30
    ) +
  theme(#strip.text = element_text(size = 20),
        axis.text.x = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.1), # c(0,0) bottom left, c(1,1) top-right.
        #axis.title.y = element_text(size = 30),
        #legend.text = element_text(size = 25)
        ) + 
  labs(y = "Symmetric Subspace Distance",
       x = "", color = "", fill = "",
       title = "Distance from True Pattern Loadings") + 
  scale_color_nejm() + scale_fill_nejm()
#dev.off()

