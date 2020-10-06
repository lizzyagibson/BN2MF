#######################################################################################
#######################################################################################
library(tidyverse)
library(ggsci)
library(GGally)
#######################################################################################
theme_set(theme_bw(base_size = 10) + 
            theme(legend.position = "bottom",
                  strip.background =element_rect(fill="white"),
                  axis.text.x = element_text(angle = 45, hjust = 1)))
#######################################################################################
ggsci <- pal_nejm()(8)

options(
  ggplot2.discrete.colour = ggsci,
  ggplot2.discrete.fill = ggsci)
#######################################################################################

# Aggregate HPC output
load(here::here(paste0("HPC/Rout/dgp_out/dgp_sims", 1, ".RDA"))) 
dgp_data <- output_all

for (i in 2:200) {
  if (file.exists(here::here(paste0("HPC/Rout/dgp_out/dgp_sims", i, ".RDA")))) { 
    load(here::here(paste0("HPC/Rout/dgp_out/dgp_sims", i, ".RDA"))) 
    dgp_data <- full_join(dgp_data, output_all)
  }
}

dgp_data <- dgp_data %>% arrange(seed)

#####################
#####################

# What does this data look like

# Distinct
sim <- dgp_data$sim[3][[1]]

as_tibble(sim) %>% 
  pivot_longer(V1:V50) %>% 
  ggplot(aes(x = value)) + 
  geom_histogram(bins = 100) + 
  facet_wrap(~name) +
  labs(x = "Distinct Simulations",
       y = "Count")

as_tibble(sim) %>% 
  pivot_longer(V1:V50) %>% 
  ggplot(aes(x = value)) + 
  geom_histogram(bins = 100) + 
  labs(x = "Distinct Simulations",
       y = "Count")

max(sim)
summary(sim)

# 12% of sim is zero
sum(sim == 0)/(1000*50)

ggcorr(as_tibble(sim), limits = FALSE,
       hjust = 0.85, size = 2, layout.exp = 1) +
  labs(x = "Simulated Correlation Matrix") +
  theme_minimal(base_size = 10)

# Overlapping
simO <- dgp_data$sim[4][[1]]

as_tibble(simO) %>% 
  pivot_longer(V1:V50) %>% 
  ggplot(aes(x = value)) + 
  geom_histogram(bins = 100) + 
  facet_wrap(~name) +
  labs(x = "Overlapping Simulations",
       y = "Count")

as_tibble(simO) %>% 
  pivot_longer(V1:V50) %>% 
  ggplot(aes(x = value)) + 
  geom_histogram(bins = 100) + 
  labs(x = "Overlapping Simulations",
       y = "Count")

max(simO)

# 8% of sim is zero
sum(simO == 0)/(1000*50)

ggcorr(as_tibble(simO), limits = FALSE,
       hjust = 0.85, size = 2, layout.exp = 1) +
  labs(x = "Simulated Correlation Matrix") +
  theme_minimal(base_size = 10)

#####################
#####################

# Relative Error
dgp_e <- dgp_data %>% 
  dplyr::select(seed, data, sim, chem, grep("pred", colnames(.))) %>% 
  dplyr::select(-pca_uncenter_pred) %>% 
  pivot_longer(pca_pred:bnmf_pred) %>% 
  mutate(l2_true = map2(chem, value, function (x,y) norm(x-y, "F")/norm(x, "F")),
         l2_sim = map2(sim, value, function (x,y) norm(x-y, "F")/norm(x, "F")),
         name = str_remove(name, "_pred")) %>% 
  dplyr::select(seed, data, name, l2_true, l2_sim) %>% 
  unnest(c(l2_sim, l2_true))

dgp_e %>%
  ggplot(aes(x = name, y = l2_sim, color = name, fill = name)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(. ~ data, scales = "free") + 
  geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  scale_y_log10() +
  theme(legend.position = "none") + 
  labs(y = "Relative Predictive Error",
       title = "vs SIMS")

dgp_e %>%
  ggplot(aes(x = name, y = l2_true, color = name, fill = name)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(. ~ data, scales = "free") + 
  geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  scale_y_log10() +
  theme(legend.position = "none") + 
  labs(y = "Relative Predictive Error",
       title = "vs PRE NOISE TRUTH")

dgp_e %>% 
  ggplot(aes(x = l2_true, color = sim_factor, fill = sim_factor)) +
  geom_histogram(bins = 100, alpha = 0.75) +
  facet_grid(name ~ data) + 
  scale_x_log10() +
  labs(title = "Relative Predictive Error vs PRE NOISE TRUTH")

dgp_e %>% 
  ggplot(aes(x = l2_sim, color = data, fill = data)) +
  geom_histogram(bins = 100, alpha = 0.75) +
  facet_grid(name ~ .) + 
  scale_x_log10() +
  labs(title = "Relative Predictive Error vs SIMS")

#######################
#######################
# Symmetric Subspace Distance
########################
########################

dgp_s <- dgp_data %>% dplyr::select(seed, data, grep("_ssdist", colnames(.))) %>% 
  dplyr::select(-grep("uncenter", colnames(.))) %>% 
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

dgp_s %>% 
  ggplot(aes(x = value)) +
  geom_histogram(aes(fill = data), bins = 100, alpha = 0.75) +
  facet_grid(model ~ type) + 
  # scale_x_log10() +
  geom_vline(xintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  labs(title = "Symmetric Subspace Distance")

dgp_s %>% 
  ggplot(aes(x = model, y = value, color = model)) +
  geom_boxplot() +
  facet_grid(data ~ type) + 
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  theme(legend.position = "none") + 
  # scale_y_log10() +
  labs(y = "Symmetric Subspace Distance")

#####
# Rank
#####

dgp_data %>%
  dplyr::select(seed, data, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank) %>%
  unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank)) %>%
  pivot_longer(cols = fa_rank:bnmf_rank) %>%
  mutate(value = ifelse(value > 5, "> 5", as.factor(value))) %>% 
  group_by(name, value, data) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  filter(data == "Distinct") %>% 
  dplyr::select(-data) %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0) %>% 
  dplyr::select(name, `1`:`5`, `> 5`)

dgp_data %>%
  dplyr::select(seed, data, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank) %>%
  unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank)) %>%
  pivot_longer(cols = fa_rank:bnmf_rank) %>%
  mutate(value = ifelse(value > 5, "> 5", as.factor(value))) %>% 
  group_by(name, value, data) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  filter(data == "Overlapping") %>% 
  dplyr::select(-data) %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0) %>% 
  dplyr::select(name, `1`:`5`, `> 5`)