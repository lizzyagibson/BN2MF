# Packages ####

source("./Results/compare_functions.R")
source("./Results/fig_set.R")
options(scipen = 999)

# Load Data ####

# BN2MF output
load("./Sep/m_rank.RDA")
load("./Sep/m_metrics.RDA")
load("./Sep/m_prop.RDA")

m_metrics
m_rank
m_prop

# Simulations
load("./Sims/sim_full.RDA")
sim_sep = sim_sep %>% dplyr::select(1:3) %>% mutate_all(as.factor)

# R output of comparison models
load("./Sep/all_metrics.RDA")
load("./Sep/all_rank.RDA")
all_metrics
all_rank

# Bootstrapped output
load("./Sep/bs_prop.RDA")
load("./Sep/bs_sample.RDA")
bs_prop
bs_sample

# VCI ####
sep_vci = bind_cols(sim_sep, m_rank, m_prop[,1])
sep_vci

prop_table = sep_vci %>%
  group_by(sep_num, noise_level) %>% 
  summarise(qs = quantile(prop, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  mutate_if(is.numeric, round, 2) %>% 
  rename(median = `0.5`) %>% 
  mutate_at(vars(1:2), as.factor)

##pdf("./Figures/sep_noise_heat_100.pdf")
prop_table %>%
  ggplot(aes(x = sep_num, y = noise_level, fill = median)) +
  geom_tile() +
  geom_text(aes(label = median), size = 4.5, col = "coral") + 
  scale_x_discrete(limits = rev) +
  labs(x = "Number of distinct chemicals per pattern",
       y = expression("Noise level as proportion of true "*sigma),
       fill = "Median coverage") +
  theme_test(base_size = 20) +
  theme(legend.position = "bottom") + 
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  theme(legend.text = element_text(size = 30))
##dev.off()

#pdf("./Figures/sep_noise_line_100.pdf", width = 6)
prop_table %>%
  mutate(Separable = ifelse(sep_num == 0, "No", "Yes")) %>% 
  ggplot(aes(x = noise_level, y = median, group = sep_num, color = Separable)) +
  geom_errorbar(aes(ymin = `0.25`, ymax = `0.75`), width = .1) +
  geom_point() + 
  geom_line() + 
  #scale_x_discrete(limits = rev) +
  labs(fill = "Number of distinct chemicals per pattern",
       x = expression("Noise level as proportion of true "*sigma),
       y = "Variational confidence interval coverage") +
  theme_test(base_size = 20) +
  theme(legend.position = "bottom") + 
  #theme(legend.text = element_text(size = 30)) +
  scale_color_nejm()
#dev.off()

# Bootstrapped CI ####

sep_bs = bind_cols(bs_sample, bs_prop)
sep_bs %>% arrange(desc(prop))

bs_table = sep_bs %>%
  group_by(sep_num, noise_level) %>% 
  summarise(qs = quantile(prop, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, "median", 0.75)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  mutate_if(is.numeric, round, 2) %>% 
  mutate_at(vars(1:2), as.factor)

#pdf("./Figures/bs_heat_100.pdf")
bs_table %>%
  ggplot(aes(x = sep_num, y = noise_level, fill = median)) +
  geom_tile() +
  geom_text(aes(label = median), size = 7, col = "coral") + 
  scale_x_discrete(limits = rev) +
  labs(x = "Number of distinct chemicals per pattern",
       y = "Noise level (as proportion of true SD)",
       fill = "Median coverage") +
  theme_test(20) +
  theme(legend.position = "bottom") + 
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  theme(legend.text = element_text(size = 10))
#dev.off()

# Rank ####
bn2mf_rank = bind_cols(sim_sep, m_rank)

sim_sep_4 = bind_rows(sim_sep, sim_sep, sim_sep, sim_sep)

get_rank = all_rank %>% 
              arrange(model) %>% 
              bind_cols(., sim_sep_4) %>% 
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
  rename_all(str_to_title) %>% 
  mutate_if(is_integer, ~replace_na(., 0)) %>% 
  dplyr::select(Model, `0`, everything()) %>% 
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
  rename_all(str_to_title) %>% 
  mutate_if(is_integer, ~replace_na(., 0)) %>% 
  dplyr::select(Model:`0.1`, `0.2`:`0.3`, `0.4`:`0.6`, `0.7`, everything()) %>% 
  knitr::kable()

get_rank %>%
  filter(rank != 4) %>% 
  group_by(model, sep_num) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = sep_num,
              values_from = "n") %>% 
  rename_all(str_to_title) %>% 
  mutate_if(is_integer, ~replace_na(., 0)) %>% 
  knitr::kable()

get_rank %>%
  filter(rank == 4 & model == "pca") %>% 
  group_by(model, sep_num, noise_level) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = sep_num,
              values_from = "n") %>% 
  rename_all(str_to_title) %>% 
  mutate_if(is_integer, ~replace_na(., 0)) %>% 
  knitr::kable()

# Metrics ####
sim_sep_3 = bind_rows(sim_sep, sim_sep, sim_sep)

bn2mf_metrics = m_metrics %>% 
                  arrange(matrix) %>% 
                  bind_cols(., sim_sep_3)

sim_sep_12 = bind_rows(sim_sep_4, sim_sep_4, sim_sep_4)

get_metrics = all_metrics %>% 
                arrange(model, matrix) %>% 
                bind_cols(., sim_sep_12) %>% 
                bind_rows(bn2mf_metrics) %>% 
                mutate(model = str_to_upper(model),
                       matrix = str_to_title(matrix))
get_metrics

metrics = get_metrics %>% 
          filter(sep_num %in% c(0, 10) & noise_level %in% c(0.2, 0.5, 1)) %>% 
          mutate(sep = ifelse(sep_num == 10, "Distinct Patterns", "Overlapping Patterns"),
                 noise = case_when(noise_level == 0.2 ~ "Noise +20%",
                                   noise_level == 0.5 ~ "Noise +50%", 
                                   TRUE ~ "Noise +100%"))

# Plots ####
metrics %>% 
  ggplot(aes(x = model, y = relerr, fill = model)) +
  geom_boxplot() +
  facet_wrap(.~matrix, scales = "free_y") +
  scale_y_log10() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Relative error", fill = "",
       title = "Across all simulations (noise levels and separability)")

#pdf("./Figures/sep_loadings_pred_100.pdf")
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

#pdf("./Figures/sep_scores_pred_100.pdf")
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

#pdf("./Figures/sep_overall_pred100.pdf")
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
metrics_sum =  metrics %>%
                group_by(sep, noise, model, matrix) %>% 
                summarize(mean_relerr = mean(relerr, na.rm=TRUE),
                          sd_relerr  = sd(relerr, na.rm = TRUE),
                          mean_ssdist = mean(ssdist, na.rm=T),
                          sd_ssdist  = sd(ssdist, na.rm = T),
                          mean_cosdist = mean(cosdist, na.rm=T),
                          sd_cosdist  = sd(cosdist, na.rm = T))
metrics_sum

for_table = metrics_sum %>% 
            ungroup() %>% 
            mutate_if(is.numeric, ~round(., 2)) %>% 
            mutate_all(as.character) %>% 
            mutate_at(vars(5:10), ~ifelse(. == "0", "<0.01", .)) %>% 
            mutate_at(vars(5:10), ~ifelse(grepl("(\\.\\d)$", .), str_c(., "0"), .)) %>% # if #.#, add zero to end
            mutate_at(vars(5:10), ~ifelse(!grepl("(\\.)", .), str_c(., ".00"), .)) %>%
            mutate(RelErr = str_c(mean_relerr, " (", sd_relerr, ")"),
                   CosDist = str_c(mean_cosdist, " (", sd_cosdist, ")"),
                   SSDist = str_c(mean_ssdist, " (", sd_ssdist, ")")) %>% 
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

xtable::xtable(other_table)
