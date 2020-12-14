library(tidyverse)

# Source ggplot settings
source("./Results/fig_set.R") 

# Get data
load("./Results/Sup Dim/sup_dim.RDA")
dim_all

## Plot

# Overall relative error
# repeat for each participant sample size
dim_all %>%
  filter(participants == "10000 Participants") %>% 
  dplyr::select(seed, participants, chemicals, patterns, grep("l2er", colnames(.))) %>% 
  pivot_longer(grep("l2er", colnames(.))) %>% 
  mutate(name = str_to_upper(str_remove(name, "_l2er"))) %>% 
  ggplot(aes(x = name, y = value, group = interaction(name,patterns))) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  geom_boxplot(aes(color = name, fill = name),alpha = 0.5, outlier.shape = NA) +
  #geom_violin(aes(color = name, fill = name),scale = "width",alpha = 0.5) +
  facet_grid(chemicals ~ patterns, scales = "free") + 
  scale_y_log10() +
  theme(legend.position = "none") + 
  labs(y = "Relative Predictive Error")

# Symmetric Subspace Distance
# repeat for each participant sample size AND chemical mixture size
dim_all %>%
  filter(participants == "10000 Participants") %>% 
  dplyr::select(seed, participants, chemicals, patterns, grep("ssdist", colnames(.))) %>% 
  pivot_longer(grep("ssdist", colnames(.)),
               names_to = c("name", "matrix", "drop"),
               names_sep = "_") %>%
  mutate(name = str_to_upper(name),
         matrix = str_to_title(ifelse(matrix == "loading", "loadings", matrix))) %>% 
  ggplot(aes(x = name, y = value, group = interaction(name,patterns))) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  #geom_boxplot(aes(color = name, fill = name),alpha = 0.5, outlier.shape = NA) +
  geom_violin(aes(color = name, fill = name),scale = "width",alpha = 0.5) +
  facet_grid(matrix ~ patterns, scales = "free") + 
  theme(legend.position = "none") + 
  labs(y = "Symmetric Subspace Distance")

# relative error on loadings and scores
# repeat for each participant sample size AND chemical mixture size
dim_all %>%
  filter(participants == "10000 Participants") %>% 
  dplyr::select(seed, participants, chemicals, patterns, 
         grep("score_l2", colnames(.)), grep("load_l2", colnames(.)),
         grep("scores_l2", colnames(.)), grep("loadings_l2", colnames(.))) %>% 
  pivot_longer(grep("_l2", colnames(.)),
               names_to = c("name", "matrix", "drop"),
               names_sep = "_") %>%
  mutate(name = str_to_upper(name),
         matrix = str_to_title(ifelse(grepl("load", matrix), "loadings", "scores"))) %>% 
  ggplot(aes(x = name, y = value, group = interaction(name,patterns))) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  #geom_boxplot(aes(color = name, fill = name),alpha = 0.5, outlier.shape = NA) +
  geom_violin(aes(color = name, fill = name),scale = "width",alpha = 0.5) +
  facet_grid(matrix ~ patterns, scales = "free") + 
  scale_y_log10() + 
  theme(legend.position = "none") + 
  labs(y = "Relative Predictive Error")

# Cosine Distance
# repeat for each participant sample size AND chemical miture size
# So far only on BN2MF
dim_all %>%
  filter(participants == "10000 Participants") %>% 
  dplyr::select(seed, participants, chemicals, patterns, grep("cos_dist", colnames(.))) %>% 
  pivot_longer(grep("cos_dist", colnames(.)),
               names_to = c("name", "drop", "drop2", "matrix"),
               names_sep = "_") %>%
mutate(name = str_to_upper(name),
       matrix = str_to_title(ifelse(grepl("load", matrix), "loadings", "scores"))) %>% 
  ggplot(aes(x = name, y = value, group = interaction(name,patterns))) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  #geom_boxplot(aes(color = name, fill = name),alpha = 0.5, outlier.shape = NA) +
  geom_violin(aes(color = name, fill = name),scale = "width",alpha = 0.5) +
  facet_grid(matrix ~ patterns, scales = "free") + 
  theme(legend.position = "none") + 
  labs(y = "Cosine Distance")

# Rank
dim_all %>%
  filter(patterns == "1 Pattern") %>% 
  dplyr::select(seed, participants, chemicals, patterns, grep("rank", colnames(.))) %>% 
  pivot_longer(grep("rank", colnames(.)),
               names_to = c("name", "drop"),
               names_sep = "_") %>% 
  mutate(name = str_to_upper(name)) %>% 
  group_by(name, value) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0)

dim_all %>%
  filter(patterns == "4 Patterns") %>% 
  dplyr::select(seed, participants, chemicals, patterns, grep("rank", colnames(.))) %>% 
  pivot_longer(grep("rank", colnames(.)),
               names_to = c("name", "drop"),
               names_sep = "_") %>% 
  mutate(name = str_to_upper(name)) %>% 
  group_by(name, value) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0)

dim_all %>%
  filter(patterns == "10 Patterns") %>% 
  dplyr::select(seed, participants, chemicals, patterns, grep("rank", colnames(.))) %>% 
  pivot_longer(grep("rank", colnames(.)),
               names_to = c("name", "drop"),
               names_sep = "_") %>% 
  mutate(name = str_to_upper(name)) %>% 
  group_by(name, value) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0)





