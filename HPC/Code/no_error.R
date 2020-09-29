#######################################################################################
#######################################################################################
library(tidyverse)
library(ggsci)
library(GGally)
#######################################################################################
theme_set(theme_bw(base_size = 15) + 
            theme(legend.position = "bottom",
                  strip.background =element_rect(fill="white"),
                  axis.text.x = element_text(angle = 45, hjust = 1)))
#######################################################################################
ggsci <- pal_nejm()(8)

options(
  ggplot2.discrete.colour = ggsci,
  ggplot2.discrete.fill = ggsci)
#######################################################################################
# Aggregate HPC output
load(here::here(paste0("HPC/Rout/chem_out/chem_sims", 1, ".RDA"))) 
chem_data <- output_all

for (i in 2:600) {
  if (file.exists(here::here(paste0("HPC/Rout/chem_out/chem_sims", i, ".RDA")))) { 
    load(here::here(paste0("HPC/Rout/chem_out/chem_sims", i, ".RDA"))) 
    chem_data <- full_join(chem_data, output_all)
  }
}

chem_data %>% arrange(seed)

#####################
#####################

# Relative Error
error_chem <- chem_data %>% 
  dplyr::select(seed, data, sim_factor, sim, chem, grep("pred", colnames(.))) %>% 
  select(-pca_uncenter_pred) %>% 
  pivot_longer(pca_pred:bnmf_pred) %>% 
  mutate(l2_true = map2(chem, value, function (x,y) norm(x-y, "F")/norm(x, "F")),
         name = str_remove(name, "_pred")) %>% 
  dplyr::select(seed, data, sim_factor, name, l2_true) %>% 
  unnest(c(l2_true))

error_chem %>% group_by(data, sim_factor, name) %>% 
  summarize(min = min(l2_true),
            med = median(l2_true),
            mean = mean(l2_true),
            max = max(l2_true)) %>% 
  filter(sim_factor == '100')

error_chem %>%
  filter(sim_factor == '100') %>% 
  ggplot(aes(x = name, y = l2_true, color = name, fill = name)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(~ data, scales = "free") + 
  geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  scale_y_log10() +
  theme(legend.position = "none") + 
  labs(title = "Relative Predictive Error vs PRE NOISE TRUTH")

error_chem %>% 
  filter(sim_factor == '100') %>% 
  ggplot(aes(x = l2_true, color = data, fill = data)) +
  geom_histogram(bins = 100, alpha = 0.75) +
  facet_grid(name~. ) + 
  scale_x_log10() +
  labs(title = "Relative Predictive Error vs PRE NOISE TRUTH")

#######################
#######################
# Symmetric Subspace Distance
########################
########################

ssdist_chem <- chem_data %>% dplyr::select(seed, data, sim_factor, grep("_ssdist", colnames(.))) %>% 
  select(-grep("uncenter", colnames(.))) %>% 
  unnest(cols = grep("_ssdist", colnames(.))) %>% 
  pivot_longer(grep("_ssdist", colnames(.)),
               names_to = "type",
               values_to = "value") %>%
  mutate(model = str_sub(type, 1, 6),
         model = ifelse(str_detect(model, 'pca'), 'PCA', model),
         model = ifelse(str_detect(model, 'l2'), 'L2 NMF', model),
         model = ifelse(str_detect(model, 'fa'), 'FA', model),
         model = ifelse(str_detect(model, '_p_'), 'Poisson NMF', model),
         model = ifelse(str_detect(model, 'bnmf'), 'BNMF', model),
         type = ifelse(str_detect(type, 'scores'), 'Scores', "Loadings"))

ssdist_chem %>% 
  filter(sim_factor == '100') %>% 
  ggplot(aes(x = value)) +
  geom_histogram(aes(fill = data), bins = 100, alpha = 0.75) +
  facet_grid(model~ type) + 
  # scale_x_log10() +
  geom_vline(xintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  labs(x = "Symmetric Subspace Distance",
       y = "Count")

ssdist_chem %>% 
  filter(sim_factor == '100') %>% 
  ggplot(aes(x = model, y = value, color = model)) +
  geom_boxplot() +
  facet_grid(data ~ type) + 
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  theme(legend.position = "none") + 
  # scale_y_log10() +
  labs(y = "Symmetric Subspace Distance")

#####
# Rank
#####

chem_data %>%
  dplyr::select(seed, data, sim_factor, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank) %>%
  unnest(c(pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank)) %>%
  pivot_longer(cols = pca_rank:bnmf_rank) %>%
  mutate(value = ifelse(value > 5, "> 5", as.factor(value))) %>% 
  group_by(name, value, data, sim_factor) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  filter(data == "Distinct" & sim_factor == "100") %>% 
  dplyr::select(-data) %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0)

chem_data %>%
  dplyr::select(seed, data, sim_factor, pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank) %>%
  unnest(c(pca_rank, nmf_l2_rank, nmf_p_rank, bnmf_rank)) %>%
  pivot_longer(cols = pca_rank:bnmf_rank) %>%
  mutate(value = ifelse(value > 5, "> 5", as.factor(value))) %>% 
  group_by(name, value, data, sim_factor) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  filter(data == "Overlapping" & sim_factor == '100') %>% 
  dplyr::select(-data) %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0)
