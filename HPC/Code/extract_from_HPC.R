#######################################################################################
#######################################################################################
#######################################################################################
library(tidyverse)

# Aggregate HPC output
comb_data <- tibble()

for (i in 1:1300) {
  if (file.exists(here::here(paste0("HPC/Rout/out_comb/out_sims", i, ".RDA")))) { 
    load(here::here(paste0("HPC/Rout/out_comb/out_sims", i, ".RDA"))) 
    }
    comb_data <- rbind(comb_data, output_all)
}

comb_data %>%
  dplyr::select(seed, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank) %>%
  unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank)) %>%
  pivot_longer(cols = fa_rank:bnmf_rank) %>%
  group_by(name, value) %>%
  summarise(n())

comb_data %>% dplyr::select(seed, grep("_ssdist", colnames(.))) %>%
  unnest(cols = c(pca_rotation_ssdist, pca_scores_ssdist, fa_rotations_ssdist,
                  fa_scores_ssdist, nmf_l2_loading_ssdist, nmf_l2_scores_ssdist,
                  nmf_p_loading_ssdist, nmf_p_scores_ssdist,
                  bnmf_scores_ssdist, bnmf_loading_ssdist)) %>%
  pivot_longer(grep("pca|nmf|fa", colnames(.)),
               names_to = "type",
               values_to = "value") %>%
  mutate(model = str_sub(type, 1, 6),
         model = ifelse(str_detect(model, 'pca'), 'pca', model),
         model = ifelse(str_detect(model, 'l2'), 'nmf_l2', model),
         model = ifelse(str_detect(model, 'fa'), 'fa', model),
         model = ifelse(str_detect(model, '_p_'), 'nmf_p', model),
         model = ifelse(str_detect(model, 'bnmf'), 'bnmf', model),
         type = ifelse(str_detect(type, 'scores'), 'scores', "loadings")) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 100, alpha = 0.75) +
  #geom_density(color = "grey") +
  facet_grid(model~type, scales = "free") +
  geom_vline(xintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  theme_bw() +
  labs(title = "Distinct Sims")

# save(comb_data, file = "./HPC/Rout/comb_data.RDA")
# load("./HPC/Rout/comb_data.RDA")
comb_data

####################################################################################################

# Relative Error
comb_data %>% dplyr::select(seed, sim_type, sim, grep("pred", colnames(.))) %>% 
  pivot_longer(pca_pred:bnmf_pred) %>% 
  mutate(l2 = map2(sim, value, function (x,y) norm(x-y, "F")/norm(x, "F"))) %>% 
  dplyr::select(seed, sim_type, name, l2) %>% 
  unnest(l2) %>% 
  group_by(sim_type, name) %>% 
  summarize(min = min(l2),
            mean = mean(l2),
            max = max(l2)) %>% knitr::kable()


comb_data_error <- comb_data %>% dplyr::select(seed, sim_type, sim, chem, grep("pred", colnames(.))) %>% 
  pivot_longer(pca_pred:bnmf_pred) %>% 
  mutate(l2_true = map2(chem, value, function (x,y) norm(x-y, "F")/norm(x, "F")),
         l2_sim = map2(sim, value, function (x,y) norm(x-y, "F")/norm(x, "F"))) %>% 
  dplyr::select(seed, sim_type, name, l2_true, l2_sim) %>% 
  unnest(c(l2_sim, l2_true)) %>% 
  mutate(Standardized = "Yes")

# save(comb_data_error, file = "./HPC/Rout/comb_data_error.RDA")

comb_data_error %>% 
  ggplot(aes(x = name, y = l2, color = name, fill = name)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(~sim_type, scales = "free") + 
  geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  theme_bw() +
  scale_y_log10(limits = c(0.07,1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(title = "Relative Predictive Error")

comb_data_error %>% 
  ggplot(aes(x = l2, color = sim_type, fill = sim_type)) +
  geom_histogram(bins = 100, alpha = 0.75) +
  facet_grid(name~., scales = "free") + 
  theme_bw() + scale_x_log10() +
  labs(title = "Relative Predictive Error") +
  theme(legend.position = "bottom")

####################################################################################################
# Symmetric Subspace Distance

comb_data_ssdist <- comb_data %>% dplyr::select(seed, sim_type, grep("_ssdist", colnames(.))) %>% 
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

comb_data_ssdist %>% 
  ggplot(aes(x = value)) +
  geom_histogram(aes(fill = sim_type), bins = 100, alpha = 0.75) +
  facet_grid(type ~ model) + 
  geom_vline(xintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  theme_bw() +
  labs(title = "Symmetric Subspace Distance")

comb_data_ssdist %>% 
  ggplot(aes(x = model, y = value, color = model)) +
  geom_boxplot(varwidth = TRUE) +
  facet_grid(sim_type~type) + 
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") + 
  labs(title = "Symmetric Subspace Distance")

# save(comb_data_ssdist, file = "./HPC/Rout/comb_data_ssdist.RDA")