library(tidyverse)

# Source ggplot settings
source("./R/fig_set.R") 

# Read in metrics
n = 2700

metrics = tibble()

for (i in 1:n) {
  load(paste0("./HPC/Rout/combo_dim/dim_out_", i, ".RDA"))
  if (exists('output_all')) {
      output_all = output_all %>% mutate(id = i) %>% 
        unnest(pca_rank:nmfp_load_l2)
      metrics = bind_rows(metrics, output_all)
      rm(output_all)
      } else {
              no_fa_summaries = no_fa_summaries %>% mutate(id = i) 
              metrics = bind_rows(metrics, no_fa_summaries) 
              rm(no_fa_summaries)}
}

metrics %>%
  slice(to_do) %>% 
  dplyr::select(id, everything())

# Should have 900 in each
metrics %>% 
  filter(patterns == 1) %>% 
  #dplyr::select(seed:patterns) %>% 
  distinct()

bn2mf_m = tibble()

for (i in 1:n) {
  if (file.exists(paste0("./HPC/Rout/combo_bnmf/bnmf_dim_out_", i, ".RDA"))) {
  load(paste0("./HPC/Rout/combo_bnmf/bnmf_dim_out_", i, ".RDA"))
  output_all = output_all %>% unnest(c(bn2mf_rank:bn2mf_cos_dist_scores))
  bn2mf_m = bind_rows(bn2mf_m, output_all)
  }
}

# Should have 900 in each
bn2mf_m %>% 
  filter(patterns == 1)

all = left_join(metrics, bn2mf_m) %>% 
  mutate(patterns = paste0(patterns, " Patterns"),
         patterns = ifelse(patterns == "1 Patterns", "1 Pattern", patterns),
         participants = paste0(participants, " Participants"),
         chemicals = paste0(chemicals, " Chemicals"),
         patterns = fct_relevel(patterns, "4 Patterns", after = 1),
         chemicals = fct_inorder(chemicals))

## Plot

# Overall relative error
# repeat for each participant sample size
all %>%
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
all %>%
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
all %>%
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
all %>%
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
