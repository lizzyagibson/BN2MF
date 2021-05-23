# SCRIPT 7
# Packages ####

source("./functions/compare_functions_2.R")
source("./functions/fig_set.R")
options(scipen = 999)
library(RColorBrewer)

# Load Data ####

# bn2mf output
load("./main/bn2mf/output/m_rank.RDA")
load("./main/bn2mf/output/m_metrics.RDA")
load("./main/bn2mf/output/m_prop.RDA")

# Simulations
load("./sims/sim_ids.RDA")

# R output of comparison models
load("./main/other_models/output/other_metrics.RDA")
load("./main/other_models/output/other_rank.RDA")

# VCI ####
sep_vci = bind_cols(sim_ids, m_rank, m_prop[,1])
sep_vci %>% filter(sep_num %in% c(0,10) & noise_level %in% c(0.2, 0.5, 1)) %>% pull(prop) %>% summary()

# This is to plot the median coverage of BN2MF at different levels of noise and separability
prop_table = sep_vci %>%
  group_by(sep_num, noise_level) %>% 
  summarise(qs = quantile(prop, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  mutate_if(is.numeric, round, 2) %>% 
  rename(median = `0.5`) %>% 
  mutate_at(vars(1:2), as.factor) %>% 
  mutate(name = "BN2MF")

blue = "#0072b5ff"

#pdf("./figures/coverage_heat.pdf")
prop_table %>%
  ggplot(aes(x = sep_num, y = noise_level, fill = median)) +
  geom_tile() +
  geom_text(aes(label = median), size = 4) + 
  scale_x_discrete(limits = rev) +
  labs(x = "Number of distinct chemicals per pattern",
       y = expression("Noise level as proportion of true "*sigma),
       fill = "Median coverage") +
  theme_test(base_size = 20) +
  theme(legend.position = "bottom") + 
  scale_fill_gradient(high = blue, low = "white") + 
  theme(legend.text = element_text(size = 10))
#dev.off()

#pdf("./figures/coverage_heat_crop.pdf",width=10)
prop_table %>%
  filter(noise_level %in% c(0,0.1,0.2,0.3,0.4,0.5)) %>% 
  #mutate(median = as_factor(median)) %>% 
  ggplot(aes(x = sep_num, y = noise_level, fill = median)) +
  geom_tile(color="black",size=.5) +
  geom_text(aes(label = median), size = 7, color="coral") + 
  scale_x_discrete(limits = rev) +
  labs(x = "Number of distinct chemicals per pattern",
       y = expression("Noise level as proportion of true "*sigma),
       fill = "Median coverage") +
  theme_minimal(base_size = 20) +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.ticks = element_line()) + 
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  theme(legend.text = element_text(size = 10))
#dev.off()

# Rank ####
# This is to get the chosen rank of each model
bn2mf_rank = bind_cols(sim_ids, m_rank)

sim_ids_4 = bind_rows(sim_ids, sim_ids, sim_ids, sim_ids)
# ^ this gets the sim ids, noise, and sep_num to match the corresponding rows

get_rank = other_rank %>% 
              arrange(model) %>% 
              bind_cols(., sim_ids_4) %>% 
              bind_rows(bn2mf_rank) %>% 
              dplyr::select(-drop)
get_rank

get_rank %>% 
  group_by(model, rank_bin) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "rank_bin",
              values_from = "n") %>% 
  rename_all(str_to_title) %>% 
  knitr::kable()

get_rank %>% 
  mutate(rank_f = case_when(rank > 6 ~ ">6",
                            TRUE ~ as.character(rank))) %>% 
  group_by(model, rank_f) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = "rank_f",
              values_from = "n") %>% 
  mutate_if(is_integer, ~replace_na(., 0)) %>% 
  dplyr::select(model, `0`, everything()) %>% 
  knitr::kable()

get_rank %>%
  filter(rank != 4) %>% 
  group_by(model) %>% 
  summarise(n = n()) %>% 
  mutate(prop_right = 1 - n/12100)

get_rank %>%
  filter(rank != 4) %>% 
  group_by(model, noise_level) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = noise_level,
              values_from = "n") %>% 
  mutate_if(is_integer, ~replace_na(., 0)) %>% 
  dplyr::select(model:`0.1`, `0.2`:`0.3`, `0.4`:`0.6`, `0.7`, everything()) %>% 
  knitr::kable()

get_rank %>%
  filter(rank != 4) %>% 
  group_by(model, sep_num) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = sep_num,
              values_from = "n") %>% 
  mutate_if(is_integer, ~replace_na(., 0)) %>% 
  knitr::kable()

get_rank %>%
  filter(rank == 4 & model == "pca") %>% 
  group_by(model, sep_num, noise_level) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = sep_num,
              values_from = "n") %>% 
  mutate_if(is_integer, ~replace_na(., 0)) %>% 
  knitr::kable()

get_rank %>%
  filter(model == "pca") %>% 
  #mutate(rank = ifelse(rank > 5, 5, rank)) %>% 
  group_by(noise_level, rank) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = rank,
              values_from = "n") %>% 
  mutate_if(is_integer, ~replace_na(., 0)) %>% 
  knitr::kable()

# metrics ####
sim_ids_3 = bind_rows(sim_ids, sim_ids, sim_ids)
# ^ this gets the sim ids, noise, and sep_num to match the corresponding rows

bn2mf_metrics = m_metrics %>% 
                  arrange(matrix) %>% 
                  bind_cols(., sim_ids_3)

sim_ids_12 = bind_rows(sim_ids_4, sim_ids_4, sim_ids_4)
# ^ this gets the sim ids, noise, and sep_num to match the corresponding rows

get_metrics = other_metrics %>% 
                arrange(model, matrix) %>% 
                bind_cols(., sim_ids_12) %>% 
                bind_rows(bn2mf_metrics) %>% 
                mutate(model = str_to_upper(model),
                       matrix = str_to_title(matrix))
get_metrics
save(get_metrics, file = "./main/metrics.rda")

# filter to subset for plotting / tables
metrics = get_metrics %>% 
          filter(sep_num %in% c(0, 10) & noise_level %in% c(0.2, 0.5, 1)) %>% 
          mutate(sep = ifelse(sep_num == 10, "Distinct Patterns", "overlapping Patterns"),
                 noise = case_when(noise_level == 0.2 ~ "Noise +20%",
                                   noise_level == 0.5 ~ "Noise +50%", 
                                   TRUE ~ "Noise +100%"))

get_metrics %>% 
  filter(model == "NMFL2" & matrix == "Pred") %>%
  group_by(sep_num, noise_level, model) %>% 
  summarise(median = median(relerr, na.rm = T)) %>% 
  mutate_if(is.numeric, round, 2) %>% 
  mutate_at(vars(1:3), as.factor) %>% 
  #filter(model == "BN2MF") %>%
  #mutate(noise = fct_relevel(noise_level, "Noise +20%", "Noise +50%", "Noise +100%")) %>% 
  ggplot(aes(x = sep_num, y = noise_level, fill = median)) +
  geom_tile() +
  geom_text(aes(label = median), size = 4) + 
  scale_x_discrete(limits = rev) +
  labs(x = "Number of distinct chemicals per pattern",
       y = expression("Noise level as proportion of true "*sigma),
       fill = "Median coverage") +
  theme_test(base_size = 20) +
  theme(legend.position = "bottom") + 
  scale_fill_gradient(high = blue, low = "white") + 
  theme(legend.text = element_text(size = 10))

# Plots ####
# metrics %>% 
#   ggplot(aes(x = model, y = relerr, fill = model)) +
#   geom_boxplot() +
#   facet_wrap(.~matrix, scales = "free_y") +
#   scale_y_log10() +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(x = "", y = "Relative error", fill = "",
#        title = "Across all simulations (noise levels and separability)")

#pdf("./figures/loadings.pdf")
metrics %>% 
  filter(matrix == "Loadings") %>% 
  mutate(noise = fct_inorder(noise)) %>% 
  ggplot(aes(x = model, y = relerr, fill = model, col = model)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.5) +
  facet_grid(sep~noise) +
  scale_y_log10() + 
  theme_bw(base_size = 20) + 
  labs(x = "", y = "Relative Prediction Error", fill = "", col = "") + 
  theme(axis.text.x = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.05), # c(0,0) bottom left, c(1,1) top-right.
        legend.text = element_text(size = 14))
#dev.off()

#pdf("./figures/scores.pdf")
metrics %>% 
  filter(sep_num %in% c(0, 10) & noise_level %in% c(0.2, 0.5, 1)) %>% 
  filter(matrix == "Scores") %>% 
  mutate(noise = fct_inorder(noise)) %>%
  ggplot(aes(x = model, y = relerr, fill = model, col = model)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.5) +
  facet_grid(sep~noise, scales = "free_y") +
  scale_y_log10() + 
  theme_bw(base_size = 20) + 
  labs(x = "", y = "Relative Prediction Error", fill = "", col = "") + 
  theme(axis.text.x = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.05), # c(0,0) bottom left, c(1,1) top-right.
        legend.text = element_text(size = 14))
#dev.off()

#pdf("./figures/overall.pdf")
metrics %>% 
  filter(sep_num %in% c(0, 10) & noise_level %in% c(0.2, 0.5, 1)) %>% 
  filter(matrix == "Pred") %>% 
  mutate(noise = fct_inorder(noise)) %>%
  ggplot(aes(x = model, y = relerr, fill = model, col = model)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.5) +
  facet_grid(sep~noise, scales = "free_y") +
  scale_y_log10() + 
  theme_bw(base_size = 20) + 
  labs(x = "", y = "Relative Prediction Error", fill = "", col = "") + 
  theme(axis.text.x = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.05), # c(0,0) bottom left, c(1,1) top-right.
        legend.text = element_text(size = 14))
#dev.off()

# Tables ####
# these are the tables in the manuscript

metrics_sum =  metrics %>%
                group_by(sep, noise, model, matrix) %>% 
                summarize(mean_relerr = mean(relerr, na.rm=TRUE),
                          sd_relerr  = sd(relerr, na.rm = TRUE),
                          mean_ssdist = mean(ssdist, na.rm=T),
                          sd_ssdist  = sd(ssdist, na.rm = T),
                          mean_cosdist = mean(cosdist, na.rm=T),
                          sd_cosdist  = sd(cosdist, na.rm = T))
metrics_sum

# this just cleans up the output for latex
for_table = metrics_sum %>% 
            ungroup() %>% 
            mutate_if(is.numeric, ~round(., 2)) %>% 
            mutate_all(as.character) %>% 
            mutate_at(vars(5:10), ~ifelse(. == "0", "<0.01", .)) %>% # don't print zero
            mutate_at(vars(5:10), ~ifelse(grepl("(\\.\\d)$", .), str_c(., "0"), .)) %>% # if #.#, add zero to end
            mutate_at(vars(5:10), ~ifelse(!grepl("(\\.)", .), str_c(., ".00"), .)) %>%
            mutate(RelErr = str_c(mean_relerr, " (", sd_relerr, ")"),
                   CosDist = str_c(mean_cosdist, " (", sd_cosdist, ")"),
                   SSDist = str_c(mean_ssdist, " (", sd_ssdist, ")")) %>% # format values as mean (stdev) 
            dplyr::select(-grep("(mean|sd)", colnames(.))) %>% 
            mutate_at(vars(5:7), ~ifelse(is.na(.), "---", .))

# relative error table ####
error_table = for_table %>% 
              dplyr::select(-CosDist, -SSDist) %>% 
              pivot_wider(names_from = matrix,
                          values_from = RelErr) %>% 
              arrange(sep, noise, model) %>% 
              dplyr::select(sep:model, Pred, Scores, Loadings)

xtable::xtable(error_table)

# other metrics table ####
other_table = for_table %>% 
  dplyr::select(-RelErr) %>% 
  pivot_wider(names_from = matrix,
              values_from = c(CosDist, SSDist)) %>% 
  arrange(sep, noise, model) %>% 
  dplyr::select(sep:model, CosDist_Pred, CosDist_Scores, CosDist_Loadings, 
                SSDist_Pred, SSDist_Scores, SSDist_Loadings)

other_table %>% 
  print(n=30)

xtable::xtable(other_table)
