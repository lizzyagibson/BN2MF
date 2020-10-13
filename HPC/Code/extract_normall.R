#####
# Set-up
#####

library(tidyverse)
library(ggsci)
library(GGally)

theme_set(theme_bw(base_size = 10) + 
            theme(legend.position = "bottom",
                  strip.background =element_rect(fill="white"),
                  axis.text.x = element_text(angle = 45, hjust = 1)))

ggsci <- pal_nejm()(8)

options(
  ggplot2.discrete.colour = ggsci,
  ggplot2.discrete.fill = ggsci)

#####
# Clean
#####

normall_data <- tibble()

for (i in 1:200) {
  if (file.exists(here::here(paste0("HPC/Rout/normall/normall_sims", i, ".RDA")))) { 
    load(here::here(paste0("HPC/Rout/normall/normall_sims", i, ".RDA"))) 
    normall_data <- rbind(normall_data, output_all)
  }
}

normall <- normall_data %>% arrange(seed)

normall %>% dplyr::select(seed) %>% distinct()

#####
# What does this data look like
#####

# Distinct
normallsim <- normall$sim[3][[1]]

as_tibble(normallsim) %>% 
  pivot_longer(V1:V50) %>% 
  ggplot(aes(x = value)) + 
  geom_histogram(bins = 100) + 
  labs(x = "Distinct Simulations",
       y = "Count")

apply(normallsim, 2, mean)
apply(normallsim, 2, sd)
# Not realistic, either

# 10% of sim is zero
sum(normallsim == 0)/(1000*50)

ggcorr(as_tibble(normallsim), limits = FALSE,
       hjust = 0.85, size = 2, layout.exp = 1) +
  labs(x = "Simulated Correlation Matrix") +
  theme_minimal(base_size = 10)

# Overlapping
simO <- norm_data$sim[303][[1]]

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

# 7% of sim is zero
sum(simO == 0)/(1000*50)

ggcorr(as_tibble(simO), limits = FALSE,
       hjust = 0.85, size = 2, layout.exp = 1) +
  labs(x = "Simulated Correlation Matrix") +
  theme_minimal(base_size = 10)

#####
# Relative error
#####

normall_e <- normall %>% 
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

normall_e %>%
  ggplot(aes(x = name, y = l1_sim, color = name, fill = name)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(sim_factor ~ data) + 
  geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  scale_y_log10() +
  theme(legend.position = "none") + 
  labs(y = "Relative Predictive Error (L1)",
       title = "vs SIMS")

normall_e %>%
  ggplot(aes(x = name, y = l1_true, color = name, fill = name)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(sim_factor ~ data) + 
  geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  scale_y_log10() +
  theme(legend.position = "none") + 
  labs(y = "Relative Predictive Error (L1)",
       title = "vs PRE NOISE TRUTH")

#####
# Symmetric subspace distance
#####

normall_s <- normall %>% dplyr::select(seed, data, sim_factor, grep("_ssdist", colnames(.))) %>% 
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

normall_s %>% 
  ggplot(aes(x = model, y = value, color = model)) +
  geom_boxplot() +
  facet_grid(sim_factor ~ data + type) + 
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  theme(legend.position = "none") + 
  # scale_y_log10() +
  labs(y = "Symmetric Subspace Distance")

normall_s %>% 
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

normall %>%
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
  arrange(sim_factor) %>% dplyr::select(name, sim_factor, `1`,`2`,`3`,`4`, `> 5`)

normall %>%
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

