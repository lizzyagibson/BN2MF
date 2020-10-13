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
load(here::here(paste0("HPC/Rout/norm_out/norm_sims", 1, ".RDA"))) 
norm_data <- output_all

for (i in 2:600) {
  if (file.exists(here::here(paste0("HPC/Rout/norm_out/norm_sims", i, ".RDA")))) { 
             load(here::here(paste0("HPC/Rout/norm_out/norm_sims", i, ".RDA"))) 
             norm_data <- full_join(norm_data, output_all)
  }
}

norm_data <- norm_data %>% arrange(seed)

# What does this data look like

# Distinct x100
sim <- norm_data$sim[3][[1]]

as_tibble(sim) %>% 
  pivot_longer(V1:V50) %>% 
  ggplot(aes(x = value)) + 
  geom_histogram(bins = 100) + 
  labs(x = "Distinct Simulations",
       y = "Count")

max(sim)
summary(sim)

# 10% of sim is zero
sum(sim == 0)/(1000*50)

ggcorr(as_tibble(sim), limits = FALSE,
       hjust = 0.85, size = 2, layout.exp = 1) +
  labs(x = "Simulated Correlation Matrix") +
  theme_minimal(base_size = 10)

# Overlapping x100
simO <- norm_data$sim[6][[1]]

as_tibble(simO) %>% 
  pivot_longer(V1:V50) %>% 
  ggplot(aes(x = value)) + 
  geom_histogram(bins = 100) + 
  labs(x = "Overlapping Simulations",
       y = "Count")

max(simO)

# 7% of sim is zero
sum(simO == 0)/(1000*50)

ggcorr(as_tibble(simO), limits = FALSE,
       hjust = 0.85, size = 2, layout.exp = 1) +
  labs(x = "Simulated Correlation Matrix") +
  theme_minimal(base_size = 10)

#####################
#####################

# Relative Error
error <- norm_data %>% 
  dplyr::select(seed, data, sim_factor, sim, chem, grep("pred", colnames(.))) %>% 
  dplyr::select(-pca_uncenter_pred) %>% 
  pivot_longer(pca_pred:bnmf_pred) %>% 
  mutate(l2_true = map2(chem, value, function (x,y) norm(x-y, "F")/norm(x, "F")),
         l2_sim = map2(sim, value, function (x,y) norm(x-y, "F")/norm(x, "F")),    
         
         l1_true = map2(chem,value, function (x,y) sum(abs(x-y))/sum(abs(x))),
         l1_sim  = map2(sim, value, function (x,y) sum(abs(x-y))/sum(abs(x))),
         
         name = str_remove(name, "_pred")) %>% 
  dplyr::select(seed, data, sim_factor, name, l2_true, l2_sim, l1_true, l1_sim) %>% 
  unnest(c(l2_sim, l2_true, l1_true, l1_sim))

error %>%
  ggplot(aes(x = name, y = l2_sim, color = name, fill = name)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(sim_factor ~ data) + 
  geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  scale_y_log10() +
  theme(legend.position = "none") + 
  labs(y = "Relative Predictive Error (L2)",
        title = "vs SIMS")

error %>%
  ggplot(aes(x = name, y = l1_true, color = name, fill = name)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(sim_factor ~ data) + 
  geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  scale_y_log10() +
  theme(legend.position = "none") + 
  labs(y = "Relative Predictive Error (L1)",
        title = "vs PRE NOISE TRUTH")

#####
# Symmetric Subspace Distance
#####

ssdist <- norm_data %>% dplyr::select(seed, data, sim_factor, grep("_ssdist", colnames(.))) %>% 
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

ssdist %>% 
  ggplot(aes(x = model, y = value, color = model)) +
  geom_boxplot() +
  facet_grid(sim_factor ~ data + type) + 
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  theme(legend.position = "none") + 
  # scale_y_log10() +
  labs(y = "Symmetric Subspace Distance")

ssdist %>% 
  filter(sim_factor == 100) %>% 
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

norm_data %>%
  dplyr::select(seed, data, sim_factor, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank) %>%
  unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank)) %>%
  pivot_longer(cols = fa_rank:bnmf_rank) %>%
  mutate(value = ifelse(value > 5, "> 5", as.factor(value))) %>% 
  group_by(name, value, data, sim_factor) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  filter(data == "Distinct") %>% 
  dplyr::select(-data) %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0) %>% 
  arrange(sim_factor) %>% dplyr::select(name, sim_factor, `1`:`5`, `> 5`)

norm_data %>%
  dplyr::select(seed, data, sim_factor, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank) %>%
  unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank)) %>%
  pivot_longer(cols = fa_rank:bnmf_rank) %>%
  mutate(value = ifelse(value > 5, "> 5", as.factor(value))) %>% 
  group_by(name, value, data, sim_factor) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  filter(data == "Overlapping") %>% 
  dplyr::select(-data) %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0) %>% 
  arrange(sim_factor) %>% dplyr::select(name, sim_factor, `1`:`5`, `> 5`)

###
# Mean Square Error
# Median Square Error
###

med_error <- norm_data %>% 
  dplyr::select(seed, data, sim_factor, sim, chem, grep("pred", colnames(.))) %>% 
  dplyr::select(-grep("uncenter", colnames(.))) %>% 
  pivot_longer(pca_pred:bnmf_pred) %>% 
  mutate(mse_sim      = map2(sim,  value, function (x,y) mean(   (x-y)^2) ),
         med_se_sim   = map2(sim,  value, function (x,y) median( (x-y)^2) ),
         rmse_sim     = map2(sim,  value, function (x,y) sqrt(mean(   (x-y)^2)) ),
         rmed_se_sim  = map2(sim,  value, function (x,y) sqrt(median( (x-y)^2)) ),
         mse_true     = map2(chem, value, function (x,y) mean(   (x-y)^2) ),
         med_se_true  = map2(chem, value, function (x,y) median( (x-y)^2) ),
         rmse_true    = map2(chem, value, function (x,y) sqrt(mean(   (x-y)^2)) ),
         rmed_se_true = map2(chem, value, function (x,y) sqrt(median( (x-y)^2)) ),
         name = str_remove(name, "_pred")) %>% 
  dplyr::select(-sim, -chem, -value) %>% 
  unnest(c(mse_sim:rmed_se_true))

med_error

# Plot any of these metrics, similar results
# "mse_sim" "med_se_sim" "rmse_sim" "rmed_se_sim" "mse_true" "med_se_true" "rmse_true" "rmed_se_true"

med_error %>% 
  ggplot(aes(x = name, y = rmed_se_true, color = name, fill = name)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(sim_factor ~ data, scales = "free") + 
  geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  scale_y_log10() +
  theme(legend.position = "none") + 
  labs(title = "")
# could plot rmse and show a line for sd/var

#####
##### 10/6
#####

norm_data %>%
  dplyr::select(seed, data, sim_factor, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank) %>%
  unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank)) %>%
  pivot_longer(cols = fa_rank:bnmf_rank) %>%
  mutate(value = ifelse(value > 5, "> 5", as.factor(value))) %>% 
  group_by(name, value, data, sim_factor) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  filter(data == "Distinct") %>% 
  dplyr::select(-data) %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0) %>% 
  arrange(sim_factor) %>% dplyr::select(name, sim_factor, `1`:`5`, `> 5`)

error %>%
  ggplot(aes(x = name, y = l2_sim, color = name, fill = name)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(sim_factor ~ data, scales = "free") + 
  geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  scale_y_log10() +
  theme(legend.position = "none") + 
  labs(y = "Relative Predictive Error",
       title = "vs SIMS")

ssdist %>% 
  ggplot(aes(x = model, y = value, color = model)) +
  geom_boxplot() +
  facet_grid(sim_factor ~ data + type) + 
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  theme(legend.position = "none") + 
  # scale_y_log10() +
  labs(y = "Symmetric Subspace Distance")

