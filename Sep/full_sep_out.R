#### Packages ####

library(tidyverse)
source("./Results/compare_functions.R")
source("./Results/fig_set.R")

#### Load Data ####

# BN2MF output
load("./Sep/m_rank.RDA")
load("./Sep/m_metrics.RDA")
load("./Sep/m_prop.RDA")
m_prop = m_prop %>% dplyr::select(-`NA`) %>% slice(-1) # First row repeats

m_metrics
m_rank
m_prop

# Simulations
load("./Sims/sim_full.RDA")
sim_sep

# R output of comparison models
#load("./Sep/combined_models.RDA")
#sep_r

#### Normalize truth ####
sim_sep = sim_sep %>% 
  mutate(denom = map(true_patterns, function(x) apply(x, 1, sum)),
         true_patternsscaled = map2(true_patterns, denom, function(x,y) x/y),
         true_scoresscaled = map2(true_scores, denom, function(x,y) as.matrix(x) %*% diag(y)),
         id = 1:nrow(.)) %>% 
  mutate_at(vars(1:3), as.factor)


#### VCI ####
sep_vci = bind_cols(sim_sep, m_rank, m_prop[,1]) %>% 
          dplyr::select(-true_patterns, -true_scores, -chem, -denom, -true_patternsscaled, 
                -true_scoresscaled)
sep_vci

prop_table = sep_vci %>%
  group_by(sep_num, noise_level) %>% 
  summarise(qs = quantile(prop, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  mutate_if(is.numeric, round, 2) %>% 
  rename(median = `0.5`)

#pdf("./Figures/sep_noise_heat_100.pdf")
prop_table %>% 
  ggplot(aes(x = sep_num, y = noise_level, fill = median)) +
  geom_tile() +
  geom_text(aes(label = median), size = 3.5, col = "coral") + 
  scale_x_discrete(limits = rev) +
  labs(x = "Number of distinct chemicals per pattern",
       y = "Noise level (as proportion of true SD)",
       fill = "Median coverage") +
  theme_test(base_size = 20) +
  theme(legend.position = "bottom") + 
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  theme(legend.text = element_text(size = 10))
#dev.off()

#### ERROR ####

sep_all = sep_r %>% 
  mutate_at(vars(1:3), as.factor) %>% 
  full_join(., sep_bn2mf) %>% 
  dplyr::select(!grep("(scaled|upper|lower|denom)", colnames(.))) %>% 
  dplyr::select(true_scores, id, everything())
sep_all

#### Rank ####
sep_rank = sep_all %>% 
  dplyr::select(seed, sep_num, noise_level, grep("rank", colnames(.))) %>% 
  pivot_longer(pca_rank:bn2mf_rank,
               names_to = c("model", "drop"),
               names_sep = "_",
               values_to = "rank") %>% 
  unnest(rank) %>% 
  mutate(rank = ifelse(is.na(rank), 0, rank),
         rank_bin = ifelse(rank ==4, "right", "wrong"))

sep_rank %>% 
  group_by(model, rank_bin) %>% 
  summarise(n = n())

##### Metrics ####
# metrics = sep_all %>%
#   dplyr::select(!grep("rank", colnames(.))) %>% 
#     pivot_longer(pca_pred:bn2mf_pred,
#                  names_to = c("model", "matrix"),
#                  names_sep = "_") %>% 
#   mutate(truth = case_when(matrix == "loadings" ~ true_patterns,
#                            matrix == "scores" ~ true_scores,
#                            matrix == "pred" ~ chem),
#          relerr  = map2(truth, value, get_relerror),
#          ssdist  = map2(truth, value, symm_subspace_dist),
#          cosdist = map2(truth, value, cos_dist)) %>% 
#   unnest(c(relerr, ssdist, cosdist)) %>% 
#   dplyr::select(-true_scores, -id, -true_patterns, -chem, -sim, -value, -truth)

#save(metrics, file = "./Sep/sep_grid_metrics.RDA")
load("./Sep/sep_grid_metrics.RDA")
metrics = metrics %>% 
  mutate(model = str_to_upper(model)) %>% 
  filter(sep_num %in% c(0, 10) & noise_level %in% c(0.5, 1)) %>% 
  mutate(sep = ifelse(sep_num == 10, "Distinct Patterns", "Overlapping Patterns"),
         noise = ifelse(noise_level == 0.5, "Noise +50%", "Noise +100%"))

metrics %>% 
  ggplot(aes(x = model, y = relerr, fill = model)) +
  geom_boxplot() +
  facet_wrap(.~matrix, scales = "free_y") +
  scale_y_log10() +
  theme_bw()

pdf("./Figures/sep_loadings_pred.pdf", height = 10, width = 5)
metrics %>% 
  filter(matrix == "loadings") %>% 
  ggplot(aes(x = model, y = relerr, fill = model, col = model)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.5) +
  facet_grid(sep~noise) +
  scale_y_log10() + 
  theme_bw(base_size = 20) + 
  labs(x = "", y = "Relative error", fill = "", col = "",
       title = "Loadings") + 
  theme(axis.text.x = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.direction = "horizontal",
        legend.position = c(0.365, -0.04), # c(0,0) bottom left, c(1,1) top-right.
        legend.text = element_text(size = 14))
dev.off()

pdf("./Figures/sep_scores_pred.pdf", height = 10, width = 5)
metrics %>% 
  filter(sep_num %in% c(0, 10) & noise_level %in% c(0.5, 1)) %>% 
  filter(matrix == "scores") %>% 
  ggplot(aes(x = model, y = relerr, fill = model, col = model)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.5) +
  facet_grid(sep~noise, scales = "free_y") +
  scale_y_log10() + 
  theme_bw(base_size = 20) + 
  labs(x = "", y = "Relative error", fill = "", col = "",
       title = "Scores") + 
  theme(axis.text.x = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.direction = "horizontal",
        legend.position = c(0.44, -0.04), # c(0,0) bottom left, c(1,1) top-right.
        legend.text = element_text(size = 14))
dev.off()

pdf("./Figures/sep_overall_pred.pdf", height = 10, width = 5)
metrics %>% 
  filter(sep_num %in% c(0, 10) & noise_level %in% c(0.5, 1)) %>% 
  filter(matrix == "pred") %>% 
  ggplot(aes(x = model, y = relerr, fill = model, col = model)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.5) +
  facet_grid(sep~noise, scales = "free_y") +
  scale_y_log10() + 
  theme_bw(base_size = 20) + 
  labs(x = "", y = "Relative error", fill = "", col = "",
       title = "Overall") + 
  theme(axis.text.x = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.direction = "horizontal",
        legend.position = c(0.425, -0.04), # c(0,0) bottom left, c(1,1) top-right.
        legend.text = element_text(size = 14))
dev.off()

# Summary ####
# Can make a heatmap for everything, but not ideal
metrics_sum = metrics %>%
  group_by(sep_num, noise_level, model, matrix) %>% 
  summarize(relerr = mean(relerr, na.rm=T),
            ssdist = mean(ssdist, na.rm=T),
            cosdist = mean(cosdist, na.rm=T))

metrics_sum %>%   
  filter(model == "nmfp" & matrix == "scores") %>% 
  ggplot(aes(x = sep_num, y = noise_level, fill = relerr)) +
  geom_tile() +
  geom_text(aes(label = round(relerr, 2)), size = 3.5, col = "coral") + 
  scale_x_discrete(limits = rev) +
  theme_minimal() +
  labs(x = "Number of distinct chemicals per pattern",
       y = "Noise level (as proportion of true SD)") +
  theme(legend.position = "bottom") + 
  scale_fill_distiller(palette = "YlGnBu", direction = -1)