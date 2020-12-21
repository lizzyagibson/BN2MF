library(tidyverse)

# Source ggplot settings
source("./Results/fig_set.R") 

# Get data
noise

## Plot

# Overall relative error
# repeat for each participant sample size

#pdf("./Figures/sim_noise_relerror.pdf", width = 10, height = 10)
noise %>%
  dplyr::select(seed, data, grep("l2er", colnames(.))) %>% 
  pivot_longer(grep("l2er", colnames(.))) %>% 
  mutate(name = str_to_upper(str_remove(name, "_l2er"))) %>%
  unnest(value) %>% 
  ggplot(aes(x = name, y = value, group = name)) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  geom_boxplot(aes(color = name, fill = name),alpha = 0.5, outlier.shape = NA) +
  #geom_violin(aes(color = name, fill = name),
  #            scale = "width",alpha = 0.5, draw_quantiles = c(0.25, 0.5, 0.75)) +
  facet_grid(. ~ data, scales = "free") + 
  scale_y_log10() +
  theme(legend.position = "none") + 
  labs(y = "Relative Predictive Error")
#dev.off()

# Symmetric Subspace Distance
# repeat for each participant sample size AND chemical mixture size

#pdf("./Figures/sim_noise_ssdist.pdf", width = 10, height = 10)
noise %>%
  dplyr::select(seed, data, grep("ssdist", colnames(.))) %>% 
  pivot_longer(grep("ssdist", colnames(.)),
               names_to = c("name", "matrix", "drop"),
               names_sep = "_") %>%
  mutate(name = str_to_upper(name),
         matrix = str_to_title(ifelse(matrix == "loading", "loadings", matrix))) %>% 
  unnest(value) %>% 
  ggplot(aes(x = name, y = value, group = name)) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  geom_boxplot(aes(color = name, fill = name),alpha = 0.5, outlier.shape = NA) +
  #geom_violin(aes(color = name, fill = name),scale = "width",alpha = 0.5) +
  facet_grid(. ~ data, scales = "free") + 
  theme(legend.position = "none") + 
  labs(y = "Symmetric Subspace Distance")
#dev.off()

# relative error on loadings and scores
# repeat for each participant sample size AND chemical mixture size

#pdf("./Figures/sim_noise_ls_relerror.pdf", width = 10, height = 10)
noise %>%
  dplyr::select(seed, data,
                grep("score_l2", colnames(.)), grep("load_l2", colnames(.)),
                grep("scores_l2", colnames(.)), grep("loadings_l2", colnames(.))) %>% 
  pivot_longer(grep("_l2", colnames(.)),
               names_to = c("name", "matrix", "drop"),
               names_sep = "_") %>%
  mutate(name = str_to_upper(name),
         matrix = str_to_title(ifelse(grepl("load", matrix), "loadings", "scores"))) %>% 
  unnest(value) %>% 
  ggplot(aes(x = name, y = value, group = name)) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  geom_boxplot(aes(color = name, fill = name),alpha = 0.5, outlier.shape = NA) +
  #geom_violin(aes(color = name, fill = name),scale = "width",alpha = 0.5) +
  facet_grid(. ~ data, scales = "free") + 
  scale_y_log10() + 
  theme(legend.position = "none") + 
  labs(y = "Relative Predictive Error")
#dev.off()

# Cosine Distance
# repeat for each participant sample size AND chemical miture size
# So far only on BN2MF
noise %>%
  dplyr::select(seed, data, grep("cos_dist", colnames(.))) %>% 
  pivot_longer(grep("cos_dist", colnames(.)),
               names_to = c("name", "drop", "drop2", "matrix"),
               names_sep = "_") %>%
  mutate(name = str_to_upper(name),
         matrix = str_to_title(ifelse(grepl("load", matrix), "loadings", "scores"))) %>% 
  unnest(value) %>% 
  ggplot(aes(x = name, y = value, group = name)) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  #geom_boxplot(aes(color = name, fill = name),alpha = 0.5, outlier.shape = NA) +
  geom_violin(aes(color = name, fill = name),scale = "width",alpha = 0.5) +
  facet_grid(. ~ data, scales = "free") + 
  theme(legend.position = "none") + 
  labs(y = "Cosine Distance")

# Rank
rank_d = noise %>%
  filter(data == "Distinct") %>% 
  dplyr::select(seed, data, grep("rank", colnames(.))) %>% 
  pivot_longer(grep("rank", colnames(.)),
               names_to = c("name", "drop"),
               names_sep = "_") %>% 
  mutate(name = str_to_upper(name)) %>% 
  unnest(value) %>% 
  group_by(name, value) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0)

rank_o = noise %>%
  filter(data == "Overlapping") %>% 
  dplyr::select(seed, data, grep("rank", colnames(.))) %>% 
  pivot_longer(grep("rank", colnames(.)),
               names_to = c("name", "drop"),
               names_sep = "_") %>% 
  mutate(name = str_to_upper(name)) %>% 
  unnest(value) %>% 
  group_by(name, value) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0)

rank_c = noise %>%
  filter(data == "Correlated") %>% 
  dplyr::select(seed, data, grep("rank", colnames(.))) %>% 
  pivot_longer(grep("rank", colnames(.)),
               names_to = c("name", "drop"),
               names_sep = "_") %>% 
  mutate(name = str_to_upper(name)) %>% 
  unnest(value) %>% 
  group_by(name, value) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0)

rank_d
rank_o
rank_c

