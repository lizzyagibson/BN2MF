#####
# Plot diff sparse H versions
# With reg BN2MF and other methods
#####

# Load packages
library(tidyverse)

# Load data
load("./Results/Main/main_metrics.RDA")
all_metrics <- metrics %>% 
               mutate(model = str_to_upper(model),
                      model = ifelse(model == "BN2MF", "BN2MF Regular H", model)) %>% 
               filter(model == "BN2MF Regular H")

load("./sparse_h/sparse_main_metrics.RDA")
sparse_metrics <- metrics %>% 
  mutate(model = str_to_upper(model),
         model = ifelse(model == "BN2MF", "BN2MF SPARSE PUSH", model)) %>% 
  filter(model == "BN2MF SPARSE PUSH")

metrics = bind_rows(all_metrics, sparse_metrics)
metrics

# Source functions
source("./Results/compare_functions.R")

# Source ggplot settings
source("./Results/fig_set.R") 

#####
# Relative Predictive Error
# L2 Norm (Truth - Predicted) / L2 Norm (Truth)
#####

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

#####
# Relative error on loadings and scores
#####

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

#####
# Tables
#####

#####
# Rank
#####

metrics %>%
  dplyr::select(seed, data, model, rank) %>%
  group_by(data, model, rank) %>%
  summarise(n = n()) %>% 
  pivot_wider(names_from = rank,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0)

#####
# Relative Preditive Error
# L2 Norm (Truth - Predicted) / L2 Norm (Truth)
#####

metrics %>%
  dplyr::select(seed, data, model, rel_err_all) %>%
  group_by(data, model) %>% 
  summarise(qs = quantile(rel_err_all, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75),
            min = min(rel_err_all)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  arrange(model) %>% print(., n = 15)

#####
# Relative error on loadings and scores
#####

metrics %>%
  dplyr::select(seed, data, model, rel_err_loadings, rel_err_scores) %>%
  pivot_longer(c(rel_err_loadings, rel_err_scores)) %>% 
  mutate(name = str_to_title(str_remove(name, "rel_err_"))) %>% 
  group_by(data, model, name) %>% 
  summarise(qs = quantile(value, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75),
            min = min(value, na.rm=T)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  arrange(name, model) %>% print(., n = 30)
