library(tidyverse)
library(ggsci)
library(GGally)

#####
# Load HPC output

#####
# Distinct patterns
dist_out_un <- tibble()

for (i in 1:100) {
  load(paste0("./HPC/Rout/isee_rep_1/out_", i, "_dist_un.RDA"))
  dist_out_un <- rbind(dist_out_un, out_dist_un)
}
dist_out_un <- dist_out_un %>% mutate(data = "distinct")

#####
# Overlapping patterns
over_out_un <- tibble()

for (i in 1:100) {
  load(paste0("./HPC/Rout/isee_rep_1/out_", i, "_over_un.RDA"))
  over_out_un <- rbind(over_out_un, out_over_un)
}
over_out_un <- over_out_un %>% mutate(data = "overlapping")

#####
# Correlated patterns
cor_out_un <- tibble()

for (i in 1:100) {
  load(paste0("./HPC/Rout/isee_rep_1/out_", i, "_cor_un.RDA"))
  cor_out_un <- rbind(cor_out_un, out_cor_un)
}
cor_out_un <- cor_out_un %>% mutate(data = "correlated")
  
#####
# Aggregate distinct, overlapping, and correlated sims

isee <- rbind(dist_out_un, over_out_un, cor_out_un)

#save(isee, file = "./HPC/Rout/isee.RDA")
load(file = "./HPC/Rout/isee.RDA")

isee %>% select(seed, data)
isee_hpc <- isee

#####
# Relative Error
dgp_ee <- isee %>% 
  dplyr::select(seed, data, sim, chem, grep("pred", colnames(.))) %>% 
  pivot_longer(pca_pred:nmf_p_pred) %>% 
  mutate(l2_true = map2(chem, value, function (x,y) norm(x-y, "F")/norm(x, "F")),
         l2_sim = map2(sim, value, function (x,y) norm(x-y, "F")/norm(x, "F")),
         name = str_remove(name, "_pred")) %>% 
  dplyr::select(seed, data, name, l2_true, l2_sim) %>% 
  unnest(c(l2_sim, l2_true))

dgp_ee %>%
  ggplot(aes(x = name, y = l2_true, color = name, fill = name)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(. ~ data, scales = "free") + 
  geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  scale_y_log10() +
  theme(legend.position = "none") + 
  labs(y = "Relative Predictive Error",
       title = "vs PRE NOISE TRUTH")

#####
# Symmetric Subspace Distance
dgp_ss <- isee %>% dplyr::select(seed, data, grep("_ssdist", colnames(.))) %>% 
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

dgp_ss %>% 
  ggplot(aes(x = model, y = value, color = model)) +
  geom_boxplot() +
  facet_grid(type ~ data) + 
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  theme(legend.position = "none") + 
  # scale_y_log10() +
  labs(y = "Symmetric Subspace Distance")

#####
# Rank
isee %>%
  dplyr::select(seed, data, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank) %>%
  unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank)) %>%
  pivot_longer(cols = fa_rank:nmf_p_rank) %>%
  mutate(value = ifelse(value > 5, "> 5", value)) %>% 
  group_by(name, value, data) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  filter(data == "distinct") %>% 
  dplyr::select(-data) %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0)

isee %>%
  dplyr::select(seed, data, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank) %>%
  unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank)) %>%
  pivot_longer(cols = fa_rank:nmf_p_rank) %>%
  mutate(value = ifelse(value > 5, "> 5", value)) %>% 
  group_by(name, value, data) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  filter(data == "overlapping") %>% 
  dplyr::select(-data) %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0)

all_unstand %>%
  dplyr::select(seed, data, fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank) %>%
  unnest(c(fa_rank, pca_rank, nmf_l2_rank, nmf_p_rank)) %>%
  pivot_longer(cols = fa_rank:nmf_p_rank) %>%
  mutate(value = ifelse(value > 5, "> 5", value)) %>% 
  group_by(name, value, data) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  filter(data == "correlated") %>% 
  dplyr::select(-data) %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0)
