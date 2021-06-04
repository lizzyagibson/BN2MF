# Packages ####

source("./functions/compare_functions_2.R")
source("./functions/fig_set.R")
options(scipen = 999)

# Load Data ####

# bn2mf output
load("./misc/dim/dim_m_rank.RDA")
load("./misc/dim/dim_m_metrics.RDA")
load("./misc/dim/dim_m_prop.RDA")

# R output of comparison models
load("./misc/dim/dim_r_metrics.RDA")
dim_r_metrics

# Sim ids
load("./misc/dim/dim_ids.RDA")
dim_ids

# VCI ####
dim_vci = bind_cols(dim_ids, bn2mf_rank, bn2mf_prop[,1])
dim_vci %>% pull(prop) %>% summary()

# 6572/13500
6275/13500
# ^ didnt get the rank right

# This is to plot the median coverage of BN2MF at different levels of noise and separability
prop_table = dim_vci %>%
  group_by(overlap, noise, patterns) %>% 
  summarise(qs = quantile(prop, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  mutate_if(is.numeric, round, 2) %>% 
  rename(median = `0.5`) %>% 
  mutate_at(vars(1:2), as.factor) %>% 
  mutate(name = "BN2MF")

blue = "#0072b5ff"

#pdf("./figures/coverage_dim_k.pdf", height=5)
prop_table %>%
  mutate(patterns = str_c("K = ", patterns),
         overlap = as_factor(ifelse(overlap == 0, 10, 0)),
         patterns = fct_inorder(patterns)) %>% 
  ggplot(aes(x = overlap, y = noise, fill = median)) +
  geom_tile() +
  geom_text(aes(label = median), size = 4) + 
  facet_grid(.~patterns, scales = "free_x", space = "free") +
  labs(x = "Number of distinct chemicals per pattern",
       y = expression("Noise level as proportion of true "*sigma),
       fill = "Median coverage") +
  theme_test(base_size = 20) +
  theme(legend.position = "bottom") + 
  scale_fill_gradient(high = blue, low = "white") + 
  theme(legend.text = element_text(size = 10),
        axis.title = element_text(size = 18),
        strip.background =element_rect(fill="white"),
        legend.position = "bottom",
        legend.title = element_text(size=15))
#dev.off()

prop_part = dim_vci %>%
  #filter(patterns==10) %>% 
  group_by(overlap, noise, participants, chemicals) %>% 
  summarise(qs = quantile(prop, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  mutate_if(is.numeric, round, 2) %>% 
  rename(median = `0.5`) %>% 
  mutate_at(vars(1:2), as.factor) %>% 
  mutate(name = "BN2MF")

#pdf("./figures/coverage_dim_np.pdf")
prop_part %>%
  mutate(participants = str_c("N = ", participants),
         chemicals = str_c("P = ", chemicals),
         overlap = as_factor(ifelse(overlap == 0, 10, 0)),
         participants = fct_inorder(participants),
         chemicals = fct_inorder(chemicals)) %>% 
  ggplot(aes(x = overlap, y = noise, fill = median)) +
  geom_tile() +
  geom_text(aes(label = median), size = 4) + 
  facet_grid(chemicals~participants) +
  labs(x = "Number of distinct chemicals per pattern",
       y = expression("Noise level as proportion of true "*sigma),
       fill = "Median coverage") +
  theme_test(base_size = 20) +
  theme(legend.position = "bottom") + 
  scale_fill_gradient(high = blue, low = "white") + 
  theme(legend.text = element_text(size = 10),
        axis.title = element_text(size = 18),
        strip.background =element_rect(fill="white"),
        legend.position = "bottom",
        legend.title = element_text(size=15))
#dev.off()

# Rank ####
# This is to get the chosen rank of each model
bn2mf_rank = bind_cols(dim_ids, bn2mf_rank)

bn2mf_rank %>% 
  filter(patterns == 10) %>%
  group_by(participants, chemicals, overlap, rank) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = rank,
              values_from = n) %>% 
  mutate_all(replace_na, 0)

2301/5400
  
dim_ids_4 = bind_rows(dim_ids, dim_ids, dim_ids, dim_ids)
# ^ this gets the sim ids, noise, and overlap to match the corresponding rows

other_rank = dim_r_metrics %>% 
  dplyr::select(job_num, grep("rank", names(.))) %>% 
  pivot_longer(pca_rank:nmfp_rank,
               values_to = "rank",
               names_to = "model") %>% 
  mutate(model = str_remove(model, "_rank"))

get_rank = other_rank %>% 
              arrange(model) %>% 
              left_join(., dim_ids) %>% 
              unnest(rank) %>% 
              bind_rows(., bn2mf_rank %>% dplyr::select(-rank_bin))

get_rank %>% 
  mutate(rank = ifelse(rank > 10, 100, rank)) %>% 
  group_by(model, rank, patterns) %>% 
  summarize(n = n()) %>% 
  pivot_wider(names_from = rank,
              values_from = n)

summary(get_rank$rank)
total_models = 13500

rank4 = get_rank %>%
  filter(patterns == 4 & rank == 4) %>% 
  group_by(model, noise, overlap) %>% 
  summarise(n = n()) %>% 
  mutate(prop_right = n/((total_models/2/3)*(2/5)))

rank1 = get_rank %>%
  filter(patterns == 1 & rank == 1) %>% 
  group_by(model, noise, overlap) %>% 
  summarise(n = n()) %>% 
  mutate(prop_right = n/((total_models/2/3)*(2/5)))

rank10 = get_rank %>%
  filter(patterns == 10 & rank == 10) %>% 
  group_by(model, noise, overlap) %>% 
  summarise(n = n()) %>% 
  mutate(prop_right = n/((total_models/2/3)*(2/5)))

# rank tables ####
table1 = rank1 %>% 
  mutate(rank = round(prop_right*100,1)) %>% 
  dplyr::select(-n, -prop_right) %>% 
  pivot_wider(names_from = c(overlap), 
              values_from = rank) %>% 
  mutate(model = str_to_upper(model)) %>% 
  bind_rows(., tibble(model = "BN2MF", noise = 0.5, `0` = 0),
            tibble(model = "BN2MF", noise = 1, `0` = 0),
            tibble(model = "PCA", noise = 1, `0` = 0)) %>% 
  arrange(model, noise) %>% 
  rename(rank1 = 3)

table4 = rank4 %>% 
  mutate(rank = round(prop_right*100,1)) %>% 
  dplyr::select(-n, -prop_right) %>% 
  pivot_wider(names_from = c(overlap), 
              values_from = rank) %>% 
  mutate(model = str_to_upper(model)) %>% 
  bind_rows(., tibble(model = "PCA", noise = 1, `0` = 0)) %>% arrange(model, noise) %>% 
  rename(rank4_0 = 3, rank4_1 = 4)

table10 = rank10 %>% 
  mutate(rank = round(prop_right*100,1)) %>% 
  dplyr::select(-n, -prop_right) %>% 
  pivot_wider(names_from = c(overlap), 
              values_from = rank) %>% 
  mutate(model = str_to_upper(model)) %>% 
  bind_rows(., tibble(model = "PCA", noise = 1, `0` = 0)) %>% arrange(model, noise) %>% 
  rename(rank10_0 = 3, rank10_1 = 4)

tablerank = bind_cols(table1, table4, table10) %>% 
  dplyr::select(1,2,3,rank4_0, rank10_0, rank4_1, rank10_1) %>% 
  arrange(noise...2) %>% dplyr::select(-noise...2)

tablerank %>% stargazer::stargazer(., summary = F)

# metrics ####

# ONLY HAVE L2 REL ERROR
dim_ids_3 = bind_rows(dim_ids, dim_ids, dim_ids)
# ^ this gets the sim ids, noise, and overlap to match the corresponding rows

bn2mf_metrics = bn2mf_metrics %>% 
                  arrange(matrix) %>% 
                  bind_cols(., dim_ids_3)

get_metrics = dim_r_metrics %>% 
                left_join(., dim_ids) %>% 
    dplyr::select(-grep("rank", names(.))) %>% 
    rename(pca_pred_l2 = pca_l2er, fa_pred_l2 = fa_l2er, nmfl2_pred_l2 = nmfl2_l2er, nmfp_pred_l2 = nmfp_l2er) %>% 
    pivot_longer(pca_pred_l2:nmfp_score_l2,
                 names_to = c("model", "matrix", "metric"),
                 names_sep = "_",
                 values_to = "relerr") %>% 
    unnest(relerr) %>% 
    bind_rows(bn2mf_metrics) %>% 
                mutate(model = str_to_upper(model),
                       matrix = str_to_title(matrix))
get_metrics

# filter to subset for plotting / tables
metrics = get_metrics %>% 
          mutate(overlap = ifelse(overlap == 0, "Distinct Patterns", "Overlapping Patterns"),
                 noise = case_when(noise == 0.2 ~ "Noise +20%",
                                   noise == 0.5 ~ "Noise +50%", 
                                   TRUE ~ "Noise +100%"),
                matrix = case_when(matrix == "Load" ~ "Loadings",
                                   matrix == "Score" ~ "Scores",
                                   TRUE ~ matrix))

# Plots ####

metrics %>% 
  filter(matrix == "Loadings") %>% 
  mutate(noise = fct_inorder(noise)) %>% 
  ggplot(aes(x = model, y = relerr, fill = model, col = model)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.5) +
  facet_grid(overlap~noise) +
  scale_y_log10() + 
  theme_bw(base_size = 20) + 
  labs(x = "", y = "Relative Prediction Error", fill = "", col = "") + 
  theme(axis.text.x = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.05), # c(0,0) bottom left, c(1,1) top-right.
        legend.text = element_text(size = 14))

metrics %>% 
  filter(matrix == "Scores") %>% 
  mutate(noise = fct_inorder(noise)) %>%
  ggplot(aes(x = model, y = relerr, fill = model, col = model)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.5) +
  facet_grid(overlap~noise, scales = "free_y") +
  scale_y_log10() + 
  theme_bw(base_size = 20) + 
  labs(x = "", y = "Relative Prediction Error", fill = "", col = "") + 
  theme(axis.text.x = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.05), # c(0,0) bottom left, c(1,1) top-right.
        legend.text = element_text(size = 14))

metrics %>% 
  filter(matrix == "Pred") %>% 
  mutate(noise = fct_inorder(noise)) %>%
  ggplot(aes(x = model, y = relerr, fill = model, col = model)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.5) +
  facet_grid(overlap~noise, scales = "free_y") +
  scale_y_log10() + 
  theme_bw(base_size = 20) + 
  labs(x = "", y = "Relative Prediction Error", fill = "", col = "") + 
  theme(axis.text.x = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.direction = "horizontal",
        legend.position = c(0.5, -0.05), # c(0,0) bottom left, c(1,1) top-right.
        legend.text = element_text(size = 14))

# Tables ####
# these are the tables in the manuscript

metrics_sum =  metrics %>%
                group_by(overlap, noise, model, matrix, patterns) %>% 
                summarize(mean_relerr = mean(relerr, na.rm=TRUE),
                          sd_relerr  = sd(relerr, na.rm = TRUE))
metrics_sum

# this just cleans up the output for latex
for_table = metrics_sum %>% 
            ungroup() %>% 
            mutate_if(is.numeric, ~round(., 2)) %>% 
            mutate_all(as.character) %>% 
            mutate_at(vars(6:7), ~ifelse(. == "0", "<0.01", .)) %>% # don't print zero
            mutate_at(vars(6:7), ~ifelse(grepl("(\\.\\d)$", .), str_c(., "0"), .)) %>% # if #.#, add zero to end
            mutate_at(vars(6:7), ~ifelse(!grepl("(\\.)", .), str_c(., ".00"), .)) %>%
            mutate(RelErr = str_c(mean_relerr, " (", sd_relerr, ")")) %>% # format values as mean (stdev) 
            dplyr::select(-grep("(mean|sd)", colnames(.))) %>% 
            mutate_at(vars(6), ~ifelse(is.na(.), "---", .))

# relative error table ####
error_table = for_table %>% 
              filter(matrix == "Pred") %>% 
              pivot_wider(names_from = c(overlap, patterns),
              values_from = RelErr)

xtable::xtable(error_table)

s_table = for_table %>% 
  filter(matrix == "Scores") %>% 
  mutate(matrix = ifelse(matrix == "Loadings", "L", "S"),
         overlap = ifelse(overlap == "Distinct Patterns", "D", "O")) %>% 
  pivot_wider(names_from = c(overlap, matrix, patterns),
              values_from = RelErr)

s_table

xtable::xtable(s_table)

l_table = for_table %>% 
  filter(matrix == "Loadings") %>% 
  mutate(matrix = ifelse(matrix == "Loadings", "L", "S"),
         overlap = ifelse(overlap == "Distinct Patterns", "D", "O")) %>% 
  pivot_wider(names_from = c(overlap, matrix, patterns),
              values_from = RelErr)

l_table

xtable::xtable(l_table)

# table discussed with MAK
sup_rank4_table = metrics %>% 
  filter(patterns == 4 & noise == "Noise +50%" & matrix != "Loadings" & participants != 1000) %>% 
  group_by(participants, chemicals, overlap, model, matrix) %>% 
  summarize(mean_relerr = mean(relerr, na.rm=TRUE),
            sd_relerr  = sd(relerr, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate_if(is.numeric, ~round(., 2)) %>% 
  mutate_all(as.character) %>% 
  mutate_at(vars(6:7), ~ifelse(. == "0", "<0.01", .)) %>% # don't print zero
  mutate_at(vars(6:7), ~ifelse(grepl("(\\.\\d)$", .), str_c(., "0"), .)) %>% # if #.#, add zero to end
  mutate_at(vars(6:7), ~ifelse(!grepl("(\\.)", .), str_c(., ".00"), .)) %>%
  mutate(RelErr = str_c(mean_relerr, " (", sd_relerr, ")")) %>% # format values as mean (stdev) 
  dplyr::select(-grep("(mean|sd)", colnames(.))) %>% 
  mutate_at(vars(6), ~ifelse(is.na(.), "---", .)) %>% 
  mutate(matrix = ifelse(matrix == "Pred", "P", "S"),
         overlap = ifelse(overlap == "Distinct Patterns", "D", "O") ) %>% 
  arrange(matrix, overlap) %>% 
  pivot_wider(names_from = c(matrix, overlap, participants),
              values_from = RelErr)
              
xtable::xtable(sup_rank4_table)

sup_rank1_table = metrics %>% 
  filter(patterns == 1 & noise == "Noise +50%" & matrix != "Loadings" & participants != 1000) %>% 
  group_by(participants, chemicals, overlap, model, matrix) %>% 
  summarize(mean_relerr = mean(relerr, na.rm=TRUE),
            sd_relerr  = sd(relerr, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate_if(is.numeric, ~round(., 2)) %>% 
  mutate_all(as.character) %>% 
  mutate_at(vars(6:7), ~ifelse(. == "0", "<0.01", .)) %>% # don't print zero
  mutate_at(vars(6:7), ~ifelse(grepl("(\\.\\d)$", .), str_c(., "0"), .)) %>% # if #.#, add zero to end
  mutate_at(vars(6:7), ~ifelse(!grepl("(\\.)", .), str_c(., ".00"), .)) %>%
  mutate(RelErr = str_c(mean_relerr, " (", sd_relerr, ")")) %>% # format values as mean (stdev) 
  dplyr::select(-grep("(mean|sd)", colnames(.))) %>% 
  mutate_at(vars(6), ~ifelse(is.na(.), "---", .)) %>% 
  mutate(matrix = ifelse(matrix == "Pred", "P", "S"),
         overlap = ifelse(overlap == "Distinct Patterns", "D", "O") ) %>% 
  arrange(matrix, overlap) %>% 
  pivot_wider(names_from = c(matrix, overlap, participants),
              values_from = RelErr)

xtable::xtable(sup_rank1_table)

sup_rank10_table = metrics %>% 
  filter(patterns == 10 & noise == "Noise +50%" & matrix != "Loadings" & participants != 1000) %>% 
  group_by(participants, chemicals, overlap, model, matrix) %>% 
  summarize(mean_relerr = mean(relerr, na.rm=TRUE),
            sd_relerr  = sd(relerr, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate_if(is.numeric, ~round(., 2)) %>% 
  mutate_all(as.character) %>% 
  mutate_at(vars(6:7), ~ifelse(. == "0", "<0.01", .)) %>% # don't print zero
  mutate_at(vars(6:7), ~ifelse(grepl("(\\.\\d)$", .), str_c(., "0"), .)) %>% # if #.#, add zero to end
  mutate_at(vars(6:7), ~ifelse(!grepl("(\\.)", .), str_c(., ".00"), .)) %>%
  mutate(RelErr = str_c(mean_relerr, " (", sd_relerr, ")")) %>% # format values as mean (stdev) 
  dplyr::select(-grep("(mean|sd)", colnames(.))) %>% 
  mutate_at(vars(6), ~ifelse(is.na(.), "---", .)) %>% 
  mutate(matrix = ifelse(matrix == "Pred", "P", "S"),
         overlap = ifelse(overlap == "Distinct Patterns", "D", "O") ) %>% 
  arrange(matrix, overlap) %>% 
  pivot_wider(names_from = c(matrix, overlap, participants),
              values_from = RelErr)

xtable::xtable(sup_rank10_table)
