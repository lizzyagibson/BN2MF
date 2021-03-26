# Packages ####

source("./functions/compare_functions.R")
source("./functions/fig_set.R")
library(patchwork)
options(scipen = 999)

# Load Data ####

# bootstrapped output
# bn2mf
load("./bootstrap/output/bn2mf_bs_metrics.RDA")
# nmf
load("./bootstrap/output/nmf_bs_metrics.RDA")
# vci
load("./bootstrap/output/vci_bs_metrics.RDA")

# bootstrap sub sample
load("./sims/bs_sample.RDA")
bs_ids = bs_sample %>% dplyr::select(set = id, sep_num, noise_level)

bs_metrics = bs_metrics %>% mutate(model = "BN2MF bootstrap") %>% 
  rename_all(~str_replace(., "ewa", "scores")) %>% 
  rename_all(~str_replace(., "eh", "load"))
vci_bs_metrics = vci_bs_metrics %>% mutate(model = "BN2MF variational") %>% 
  dplyr::select(-err_ewa, -err_eh) %>% 
  rename_all(~str_replace(., "(ewa|wa)", "scores")) %>% 
  rename_all(~str_replace(., "(_eh|_h)", "_load"))
nmf_bs_metrics = nmf_bs_metrics %>% mutate(model = "NMF bootstrap")

all_metrics = bind_rows(vci_bs_metrics, bs_metrics, nmf_bs_metrics) %>% full_join(., bs_ids)

# Compare ####

# Proportion in CI, SCORES
prop_table = all_metrics %>%
  group_by(sep_num, noise_level, model) %>% 
  summarise(qs = quantile(prop_scores, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  mutate_if(is.numeric, round, 2) %>% 
  rename(median = `0.5`) %>% 
  mutate_at(vars(1:2), as.factor)

blue = "#0072b5ff"

prop_table %>%
  ggplot(aes(x = sep_num, y = noise_level, fill = median)) +
  geom_tile() +
  geom_text(aes(label = median), size = 4) + 
  scale_x_discrete(limits = rev) +
  facet_wrap(~model) +
  labs(x = "Number of distinct chemicals per pattern",
       y = expression("Noise level as proportion of true "*sigma),
       fill = "Median coverage") +
  theme_test(base_size = 20) +
  theme(legend.position = "bottom") + 
  scale_fill_gradient(high = blue, low = "white") + 
  theme(legend.text = element_text(size = 10))

# Proportion in CI, LOADINGS
prop_load = all_metrics %>%
  group_by(sep_num, noise_level, model) %>% 
  summarise(qs = quantile(prop_load, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  mutate_if(is.numeric, round, 2) %>% 
  rename(median = `0.5`) %>% 
  mutate_at(vars(1:2), as.factor)

prop_load %>%
  ggplot(aes(x = sep_num, y = noise_level, fill = median)) +
  geom_tile() +
  geom_text(aes(label = median), size = 4) + 
  scale_x_discrete(limits = rev) +
  facet_wrap(~model) +
  labs(x = "Number of distinct chemicals per pattern",
       y = expression("Noise level as proportion of true "*sigma),
       fill = "Median coverage") +
  theme_test(base_size = 20) +
  theme(legend.position = "bottom") + 
  scale_fill_gradient(high = blue, low = "white") + 
  theme(legend.text = element_text(size = 10))

# Proportion in CI, SCORES
prop_pred = all_metrics %>%
  group_by(sep_num, noise_level, model) %>% 
  summarise(qs = quantile(prop_pred, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  mutate_if(is.numeric, round, 2) %>% 
  rename(median = `0.5`) %>% 
  mutate_at(vars(1:2), as.factor)

prop_pred %>%
  ggplot(aes(x = sep_num, y = noise_level, fill = median)) +
  geom_tile() +
  geom_text(aes(label = median), size = 4) + 
  scale_x_discrete(limits = rev) +
  facet_wrap(~model) +
  labs(x = "Number of distinct chemicals per pattern",
       y = expression("Noise level as proportion of true "*sigma),
       fill = "Median coverage") +
  theme_test(base_size = 20) +
  theme(legend.position = "bottom") + 
  scale_fill_gradient(high = blue, low = "white") + 
  theme(legend.text = element_text(size = 10))

# Tables ####

medians = all_metrics %>% 
  group_by(model, sep_num, noise_level) %>% 
  summarise(prop_scores       = median(prop_scores ),
            prop_load         = median(prop_load ),
            prop_pred         = median(prop_pred ),
            err_scores_scaled = median(err_scores_scaled ),
            err_load_scaled   = median(err_load_scaled ),
            err_pred          = median(err_pred ),
            width_scores      = median(width_scores ),
            width_load        = median(width_load ),
            width_pred        = median(width_pred)) %>% 
  arrange(noise_level, sep_num, model) %>% 
  mutate_at(vars(c(noise_level, sep_num, model)), as_factor) %>% 
  mutate_if(is.numeric, round, 2)

# Scores ####
score_prop = medians %>%
  ggplot(aes(x = sep_num, y = noise_level, fill = prop_scores)) +
  geom_tile() +
  geom_text(aes(label = prop_scores), size = 4) + 
  scale_x_discrete(limits = rev) +
  facet_wrap(~model) +
  labs(x = "Number of distinct chemicals per pattern",
       y = expression("Noise level as proportion of true "*sigma),
       fill = "Median coverage") +
  theme_test(base_size = 15) +
  scale_fill_gradient(high = blue, low = "white") + 
  theme(title = element_text(size = 10)) + 
  labs(subtitle = "Median coverage")
  
score_width = medians %>%
  ggplot(aes(x = sep_num, y = noise_level, fill = width_scores)) +
  geom_tile() +
  geom_text(aes(label = width_scores), size = 4) + 
  scale_x_discrete(limits = rev) +
  facet_wrap(~model) +
  labs(x = "Number of distinct chemicals per pattern",
       y = expression("Noise level as proportion of true "*sigma),
       fill = "Median width") +
  theme_test(base_size = 15) +
  theme(axis.title.y = element_text(size = 13),
        title = element_text(size = 10)) + 
  scale_fill_gradient(high = "#E18727FF", low = "white") + 
  labs(subtitle = "Median 95%CI width")

score_err = medians %>%
  ggplot(aes(x = sep_num, y = noise_level, fill = err_scores_scaled)) +
  geom_tile() +
  geom_text(aes(label = err_scores_scaled), size = 4) + 
  scale_x_discrete(limits = rev) +
  facet_wrap(~model) +
  labs(x = "Number of distinct chemicals per pattern",
       y = expression("Noise level as proportion of true "*sigma),
       fill = "Median error") +
  theme_test(base_size = 15) +
  theme(axis.title.y = element_text(size = 12),
        title = element_text(size = 10)) + 
  scale_fill_gradient(high = "#20854EFF", low = "white") + 
  labs(subtitle = "Median relative error")

pdf("./figures/bootstrap_metrics_scores.pdf", height = 10)
(score_prop +
  theme(axis.title = element_blank(),
        legend.position = "none")) / 
  (score_width +
  theme(axis.title.x = element_blank(),
        legend.position = "none")) / 
  (score_err +
     theme(axis.title.y = element_blank(),
           legend.position = "none"))
dev.off()

# Loadings ####
load_prop = medians %>%
  ggplot(aes(x = sep_num, y = noise_level, fill = prop_load)) +
  geom_tile() +
  geom_text(aes(label = prop_load), size = 4) + 
  scale_x_discrete(limits = rev) +
  facet_wrap(~model) +
  labs(x = "Number of distinct chemicals per pattern",
       y = expression("Noise level as proportion of true "*sigma),
       fill = "Median coverage") +
  theme_test(base_size = 15) +
  scale_fill_gradient(high = blue, low = "white") + 
  theme(title = element_text(size = 10)) + 
  labs(subtitle = "Median coverage")

load_width = medians %>%
  ggplot(aes(x = sep_num, y = noise_level, fill = width_load)) +
  geom_tile() +
  geom_text(aes(label = width_load), size = 4) + 
  scale_x_discrete(limits = rev) +
  facet_wrap(~model) +
  labs(x = "Number of distinct chemicals per pattern",
       y = expression("Noise level as proportion of true "*sigma),
       fill = "Median width") +
  theme_test(base_size = 15) +
  theme(axis.title.y = element_text(size = 13),
        title = element_text(size = 10)) + 
  scale_fill_gradient(high = "#E18727FF", low = "white") + 
  labs(subtitle = "Median 95%CI width")

load_err = medians %>%
  ggplot(aes(x = sep_num, y = noise_level, fill = err_load_scaled)) +
  geom_tile() +
  geom_text(aes(label = err_load_scaled), size = 4) + 
  scale_x_discrete(limits = rev) +
  facet_wrap(~model) +
  labs(x = "Number of distinct chemicals per pattern",
       y = expression("Noise level as proportion of true "*sigma),
       fill = "Median error") +
  theme_test(base_size = 15) +
  theme(axis.title.y = element_text(size = 12),
        title = element_text(size = 10)) + 
  scale_fill_gradient(high = "#20854EFF", low = "white") + 
  labs(subtitle = "Median relative error")

pdf("./figures/bootstrap_metrics_load.pdf", height = 10)
(load_prop +
    theme(axis.title = element_blank(),
          legend.position = "none")) / 
  (load_width +
     theme(axis.title.x = element_blank(),
           legend.position = "none")) / 
  (load_err +
     theme(axis.title.y = element_blank(),
           legend.position = "none"))
dev.off()

# Pred ####
pred_prop = medians %>%
  ggplot(aes(x = sep_num, y = noise_level, fill = prop_pred)) +
  geom_tile() +
  geom_text(aes(label = prop_pred), size = 4) + 
  scale_x_discrete(limits = rev) +
  facet_wrap(~model) +
  labs(x = "Number of distinct chemicals per pattern",
       y = expression("Noise level as proportion of true "*sigma),
       fill = "Median coverage") +
  theme_test(base_size = 15) +
  scale_fill_gradient(high = blue, low = "white") + 
  theme(title = element_text(size = 10)) + 
  labs(subtitle = "Median coverage")

pred_width = medians %>%
  ggplot(aes(x = sep_num, y = noise_level, fill = width_pred)) +
  geom_tile() +
  geom_text(aes(label = width_pred), size = 4) + 
  scale_x_discrete(limits = rev) +
  facet_wrap(~model) +
  labs(x = "Number of distinct chemicals per pattern",
       y = expression("Noise level as proportion of true "*sigma),
       fill = "Median width") +
  theme_test(base_size = 15) +
  theme(axis.title.y = element_text(size = 13),
        title = element_text(size = 10)) + 
  scale_fill_gradient(high = "#E18727FF", low = "white") + 
  labs(subtitle = "Median 95%CI width")

pred_err = medians %>%
  ggplot(aes(x = sep_num, y = noise_level, fill = err_pred)) +
  geom_tile() +
  geom_text(aes(label = err_pred), size = 4) + 
  scale_x_discrete(limits = rev) +
  facet_wrap(~model) +
  labs(x = "Number of distinct chemicals per pattern",
       y = expression("Noise level as proportion of true "*sigma),
       fill = "Median error") +
  theme_test(base_size = 15) +
  theme(axis.title.y = element_text(size = 12),
        title = element_text(size = 10)) + 
  scale_fill_gradient(high = "#20854EFF", low = "white") + 
  labs(subtitle = "Median relative error")

pdf("./figures/bootstrap_metrics_pred.pdf", height = 10)
(pred_prop +
    theme(axis.title = element_blank(),
          legend.position = "none")) / 
  (pred_width +
     theme(axis.title.x = element_blank(),
           legend.position = "none")) / 
  (pred_err +
     theme(axis.title.y = element_blank(),
           legend.position = "none"))
dev.off()

