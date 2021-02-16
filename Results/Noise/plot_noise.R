library(tidyverse)

# Source functions
source("./Results/compare_functions.R")

# Source ggplot settings
source("./Results/fig_set.R") 

# Load packages
library(tidyverse)
library(stargazer)

# Load data
load("./Results/Noise/noise_metrics.RDA")
metrics <- metrics %>% 
  mutate(model = str_to_upper(model))

metrics %>% 
  group_by(data, model, rank) %>% 
  summarize(n = n()) %>% 
  pivot_wider(values_from = "n",
              names_from = "rank")

##### Relative Preditive Error #####
# L2 Norm (Truth - Predicted) / L2 Norm (Truth)

# pdf("./Figures/bnmf_error.pdf", width = 10)
metrics %>%
  dplyr::select(seed, data, model, rel_err_all) %>%
  mutate(data = fct_inorder(data)) %>% 
  ggplot(aes(x = model, y = rel_err_all)) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  geom_boxplot(aes(color = model, fill = model),
               alpha = 0.5, outlier.shape = NA) +
  facet_grid(. ~ data, scales = "free") + 
  scale_y_log10() +
  theme(legend.position = "none") + 
  labs(y = "Relative Predictive Error")

##### Relative error on loadings and scores #####

#pdf("./Figures/bnmf_loadscore_error.pdf", width = 10, height = 10)
metrics %>%
  mutate(data = fct_inorder(data)) %>% 
  dplyr::select(seed, data, model, rel_err_loadings, rel_err_scores) %>%
  pivot_longer(c(rel_err_loadings, rel_err_scores)) %>% 
  mutate(name = str_to_title(str_remove(name, "rel_err_"))) %>% 
  ggplot(aes(x = model, y = value)) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  geom_boxplot(aes(color = model, fill = model), 
               alpha = 0.5, outlier.shape = NA, varwidth = TRUE) +
  facet_grid(name ~ data, scales = "free") + 
  scale_y_log10() +
  labs(y = "Relative Predictive Error")

##### Subspace Distance ####
# Distance between linear subspaces (orthonormal bases)
# Loadging and scores

#pdf("./Figures/bnmf_ssd.pdf", width = 10, height = 10)
metrics %>%
  mutate(data = fct_inorder(data)) %>% 
  dplyr::select(seed, data, model, ssd_loadings, ssd_scores) %>%
  pivot_longer(c(ssd_loadings, ssd_scores)) %>% 
  mutate(name = str_to_title(str_remove(name, "ssd_"))) %>% 
  ggplot(aes(x = model, y = value)) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  geom_boxplot(aes(color = model, fill = model),
               alpha = 0.5, outlier.shape = NA) +
  facet_grid(name ~ data) + 
  labs(y = "Symmetric Subspace Distance")

##### Cosine distance #####
metrics %>%
  mutate(data = fct_inorder(data)) %>% 
  dplyr::select(seed, data, model, cos_dist_loadings, cos_dist_scores) %>%
  pivot_longer(c(cos_dist_loadings, cos_dist_scores)) %>% 
  mutate(name = str_to_title(str_remove(name, "cos_dist_"))) %>% 
  ggplot(aes(x = model, y = value)) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  geom_boxplot(aes(color = model, fill = model), 
               alpha = 0.5, outlier.shape = NA) +
  #geom_violin(aes(color = model, fill = model), scale = "width",alpha = 0.5) +
  facet_grid(name ~ data, scales = "free") + 
  labs(y = "Cosine Distance on Matrix")

# pdf("./Figures/bnmf_cos_v.pdf", width = 10, height = 10)
metrics %>%
  dplyr::select(seed, data, model, cos_dist_v_loadings, cos_dist_v_scores) %>%
  pivot_longer(c(cos_dist_v_loadings, cos_dist_v_scores)) %>% 
  mutate(name = str_to_title(str_remove(name, "cos_dist_v_"))) %>% 
  unnest(value) %>% 
  ggplot(aes(x = model, y = value)) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  geom_boxplot(aes(color = model, fill = model), 
               alpha = 0.5, outlier.shape = NA) +
  # geom_violin(aes(color = model, fill = model), scale = "width",alpha = 0.5) +
  facet_grid(name ~ data, scales = "free") + 
  labs(y = "Cosine Distance on Vectors")
# dev.off()

##### Tables #####

##### Rank #####
metrics %>%
  dplyr::select(seed, data, model, rank) %>%
  group_by(data, model, rank) %>%
  summarise(n = n()) %>% 
  pivot_wider(names_from = rank,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0)

# 210 solutions are not rank 4
100+100+10
# this accounts for the 420 missing values (loadings and scores) in the plots above

##### Relative Preditive Error ####
# L2 Norm (Truth - Predicted) / L2 Norm (Truth)
table_l2 = metrics %>%
  dplyr::select(seed, data, model, rel_err_all) %>%
  group_by(data, model) %>% 
  summarise(qs = quantile(rel_err_all, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75),
            mean = mean(rel_err_all)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, model, `0.25`, `0.5`, mean, `0.75`) %>% 
  arrange(data) %>% 
  mutate_if(is.numeric, round, 2)

table_l2

as_data_frame(table_l2) %>% stargazer(summary = FALSE)

##### Relative error on loadings and scores #####
table_ls = metrics %>%
  dplyr::select(seed, data, model, rel_err_loadings, rel_err_scores) %>%
  pivot_longer(c(rel_err_loadings, rel_err_scores)) %>% 
  mutate(name = str_to_title(str_remove(name, "rel_err_"))) %>% 
  group_by(data, model, name) %>% 
  summarise(qs = quantile(value, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  arrange(name, data) %>%
  mutate_if(is.numeric, round, 2)

table_ls %>% print(n = 30)

as_data_frame(table_ls) %>% stargazer(summary = FALSE)

##### Subspace Distance ####
# Distance between linear subspaces (orthonormal bases)
# Loading and scores
table_ssd = metrics %>%
  dplyr::select(seed, data, model, ssd_loadings, ssd_scores) %>%
  pivot_longer(c(ssd_loadings, ssd_scores)) %>% 
  mutate(name = str_to_title(str_remove(name, "ssd_"))) %>% 
  group_by(data, model, name) %>% 
  summarise(qs = quantile(value, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  arrange(name, data) %>% 
  mutate_if(is.numeric, round, 2)
table_ssd

as_data_frame(table_ssd) %>% stargazer(summary = FALSE)

##### Cosine distance #####
table_cos = metrics %>%
  dplyr::select(seed, data, model, cos_dist_loadings, cos_dist_scores) %>%
  pivot_longer(c(cos_dist_loadings, cos_dist_scores)) %>% 
  mutate(name = str_to_title(str_remove(name, "cos_dist_"))) %>% 
  group_by(data, model, name) %>% 
  summarise(qs = quantile(value, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75),
            min = min(value)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  arrange(name, data) %>% 
  mutate_if(is.numeric, round, 2)

table_cos %>% print(n = 30)

as_data_frame(table_cos) %>% stargazer(summary = FALSE)

#### Mean (SD) ####
table_mean = metrics %>%
  dplyr::select(-rank, -grep("rel_err", colnames(.)), -cos_dist_v_scores, -cos_dist_v_loadings) %>%
  pivot_longer(c(cos_dist_loadings:ssd_scores),
               names_to = c("method", "x", "side"),
               names_sep = "_") %>% 
  mutate(side = ifelse(x %in% c("loadings", "scores"), x, side)) %>% 
  group_by(data, model, side, method) %>% 
  summarise(mean = mean(value, na.rm = T),
            sd = sd(value, na.rm = T)) %>% 
  pivot_wider(names_from = c("method", "side"),
              values_from = c("mean", "sd")) %>% 
  arrange(data) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  dplyr::select(data, model, grep("cos_score", colnames(.)), grep("ssd_score", colnames(.)),
                grep("cos_load", colnames(.)), grep("ssd_load", colnames(.)))

table_mean

as_data_frame(table_mean) %>% stargazer(summary = FALSE)

table_rel = metrics %>%
  dplyr::select(seed, data, model, grep("rel_err", colnames(.))) %>%
  pivot_longer(grep("rel_err", colnames(.)),
               names_to = c("method", "x", "side"),
               names_sep = "_") %>% 
  group_by(data, model, side, method) %>% 
  summarise(mean = mean(value, na.rm = T),
            sd = sd(value, na.rm = T)) %>% 
  pivot_wider(names_from = c("method", "side"),
              values_from = c("mean", "sd")) %>% 
  arrange(data) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  dplyr::select(data, model, grep("rel_all", colnames(.)), grep("score", colnames(.)),
                grep("load", colnames(.)))

table_rel

as_data_frame(table_rel) %>% stargazer(summary = FALSE)




#### SAVE Figures ####

# Relative Predictive Error
# L2 Norm (Truth - Predicted) / L2 Norm (Truth)
pred_error = metrics %>%
  filter(model != "NMFP") %>% 
  mutate(model = str_remove(model, "L2"),
         model = ifelse(model == "FA", "Factor Analysis", model)) %>% 
  dplyr::select(seed, data, model, rel_err_all) %>%
  mutate(data = fct_inorder(data)) %>% 
  ggplot(aes(x = model, y = rel_err_all)) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  geom_boxplot(aes(color = model, fill = model),
               alpha = 0.5, outlier.shape = NA) +
  facet_grid(. ~ data, scales = "free") + 
  scale_y_log10() +
  theme(axis.text.x = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.1), # c(0,0) bottom left, c(1,1) top-right.
        axis.title.y = element_text(size = 30),
        plot.title =element_text(size = 25, hjust = 0.5),
        legend.text = element_text(size = 25)) + 
  labs(title = "Distance from True Observations",
       y = "Relative Predictive Error",
       color = "", fill = "")

#pdf("/Users/lizzy/OneDrive/Columbia/Spring 2021/netflix/figures/pred_error.pdf", width = 10)
pred_error
#dev.off()

# Relative error on loadings and scores
score_error = metrics %>%
  filter(model != "NMFP") %>% 
  mutate(model = str_remove(model, "L2"),
         model = ifelse(model == "FA", "Factor Analysis", model)) %>% 
  dplyr::select(seed, data, model, rel_err_loadings, rel_err_scores) %>%
  mutate(data = fct_inorder(data)) %>% 
  ggplot(aes(x = model, y = rel_err_scores)) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  geom_boxplot(aes(color = model, fill = model),
               alpha = 0.5, outlier.shape = NA) +
  facet_grid(. ~ data, scales = "free") + 
  scale_y_log10() +
  theme(axis.text.x = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.1), # c(0,0) bottom left, c(1,1) top-right.
        axis.title.y = element_text(size = 30),
        plot.title =element_text(size = 25, hjust = 0.5),
        legend.text = element_text(size = 25)) + 
  labs(title = "Distance from True Scores",
       y = "Relative Predictive Error",
       color = "", fill = "")

load_error = metrics %>%
  filter(model != "NMFP") %>% 
  mutate(model = str_remove(model, "L2"),
         model = ifelse(model == "FA", "Factor Analysis", model)) %>% 
  dplyr::select(seed, data, model, rel_err_loadings, rel_err_scores) %>%
  mutate(data = fct_inorder(data)) %>% 
  ggplot(aes(x = model, y = rel_err_loadings)) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  geom_boxplot(aes(color = model, fill = model),
               alpha = 0.5, outlier.shape = NA) +
  facet_grid(. ~ data, scales = "free") + 
  scale_y_log10() +
  theme(axis.text.x = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.1), # c(0,0) bottom left, c(1,1) top-right.
        axis.title.y = element_text(size = 30),
        plot.title =element_text(size = 25, hjust = 0.5),
        legend.text = element_text(size = 25)) + 
  labs(title = "Distance from True Loadings",
       y = "Relative Predictive Error",
       color = "", fill = "")

#pdf("/Users/lizzy/OneDrive/Columbia/Spring 2021/netflix/figures/score_error.pdf", width = 10)
score_error
#dev.off()

#pdf("/Users/lizzy/OneDrive/Columbia/Spring 2021/netflix/figures/load_error.pdf", width = 10)
load_error
#dev.off()

score_error / load_error

#pdf("/Users/lizzy/OneDrive/Columbia/Spring 2021/netflix/figures/side_error.pdf", width = 10, height = 9)
metrics %>%
  filter(model != "NMFP") %>% 
  mutate(model = str_remove(model, "L2"),
         model = ifelse(model == "FA", "Factor Analysis", model)) %>% 
  dplyr::select(seed, data, model, rel_err_loadings, rel_err_scores) %>%
  pivot_longer(rel_err_loadings:rel_err_scores) %>% 
  mutate(name = str_to_title(str_remove(name, "rel_err_"))) %>% 
  mutate(data = fct_inorder(data)) %>% 
  ggplot(aes(x = model, y = value)) +
  geom_jitter(alpha = 0.15, size = 0.5, height = 0, width = .3) +
  geom_boxplot(aes(color = model, fill = model),
               alpha = 0.5, outlier.shape = NA) +
  facet_grid(name ~ data, scales = "free") + 
  scale_y_log10() +
  theme(axis.text.x = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.05), # c(0,0) bottom left, c(1,1) top-right.
        axis.title.y = element_text(size = 30),
        plot.title =element_text(size = 25, hjust = 0.5),
        legend.text = element_text(size = 25)) + 
  labs(title = "Distance from Truth",
       y = "Relative Predictive Error",
       color = "", fill = "")
#dev.off()
