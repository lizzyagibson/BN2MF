## Combine BN2MF with VCI, main results in MATLAB
## Combine other methods, main results in R

## Take metrics for each method
  # Symmetric subspace distance
  # Cosine distance
  # Relative error
    # on sims
    # on loadings
    # on scores

# Run one time!

## Packages
library(tidyverse)

## Get functions
source("./Results/compare_functions.R")
source("./Results/fig_set.R")

# Load all
load("./sparse_h/main_out.RDA")
sparsehH = dgp  %>% 
  rename_at(vars(7:34), ~str_c(., "_shH")) %>% 
  dplyr::select(-grep("(pca|fa|nmfl2|nmfp)", colnames(.)))

load("./sparse_h/sparse_main_out.RDA")
dgp_m = dgp_m %>% 
  rename_all(~str_c(., "_sh"))

load("./Results/main/main_out.RDA")
reg = dgp %>% 
  rename_at(vars(7:34), ~str_c(., "_reg")) %>% 
  mutate_at(grep("(a|2|p)_loadings", colnames(.)), function(x) map(x, t))

dgp_all = bind_cols(sparsehH, dgp_m) %>% 
  full_join(., reg)
dgp_all

# Calculate error metrics
dgp_metrics0 = dgp_all %>% 
              pivot_longer(bn2mf_loadings_shH:bn2mf_rank_reg,
               names_to = c("model", "object", "sparsity"),
               names_sep = "_") %>% 
              filter(!(object %in% c("perm", "out", "pred", "rank"))) %>% 
  pivot_wider(names_from = c(object),
              values_from = value)

# Normalize loadings
# Scale scores
# It really doesn't make sense to L1 scale PCA and FA components
# Just curious
dgp_metrics = dgp_metrics0 %>% 
    mutate(true_denom = map(true_patterns, function(x) apply(x,1, sum)),
           true_patterns_scaled = map2(true_patterns, true_denom, function(x,y) x/y),
           true_scores_scaled = map2(true_scores, true_denom, function(x,y) as.matrix(x) %*% diag(y)),
           denom = map(loadings, function(x) apply(x, 1, sum)),
           scaled_loadings = map2(loadings, denom, function(x,y) x/y),
           scaled_scores   = map2(scores, denom, function(x,y) x %*% diag(y))) %>% 
  rename(reg_loadings = loadings, reg_scores = scores) %>% 
  dplyr::select(-denom, -true_denom) %>% 
  pivot_longer(c(reg_loadings, reg_scores, scaled_loadings, scaled_scores),
               names_to = c("scaled", "object"),
               names_sep = "_")

# l2 relative error -- loadings
rel_err_loadings <- dgp_metrics %>%
                filter(object == "loadings") %>% 
                mutate(l2_loadings = map2(true_patterns_scaled, value, get_relerror)) %>%
                mutate(l1_loadings = map2(true_patterns_scaled, value, get_relerror_l1)) %>%
                dplyr::select(seed, data, model, sparsity, scaled, l2_loadings, l1_loadings) %>% 
                unnest(c(l2_loadings, l1_loadings))

# l2 relative error -- scores
rel_err_scores <- dgp_metrics %>%
                  filter(object == "scores") %>% 
                  mutate(l2_scores = map2(true_scores_scaled, value, get_relerror)) %>%
                  mutate(l1_scores = map2(true_scores_scaled, value, get_relerror_l1)) %>%
                  dplyr::select(seed, data, model, sparsity, scaled, l2_scores, l1_scores) %>% 
                  unnest(c(l2_scores, l1_scores))

# Combine all metrics
metrics <- full_join(rel_err_loadings, rel_err_scores)

#### Viz ####
metrics %>% 
  mutate_at(vars(2:5), as.factor) %>% 
  #filter(scaled == "scaled") %>% 
  filter(model != "pca") %>% # BC PCA error is way higher on loadings
  filter(model != "fa") %>% 
  #filter(data != "Distinct") %>% 
  mutate(model = str_c(model, "_", sparsity),
         model = str_remove(model, "_reg")) %>% 
  mutate(data = fct_inorder(data)) %>% 
  dplyr::select(seed, data, model, sparsity, scaled, l2_loadings, l2_scores) %>%
  pivot_longer(c(l2_loadings, l2_scores)) %>% 
  mutate(name = str_to_title(str_remove(name, "l2_"))) %>%
  filter(name == "Loadings") %>% 
  ggplot(aes(x = model, y = value)) +
  geom_boxplot(aes(color = model, fill = model), 
               alpha = 0.5, outlier.shape = NA) +
  facet_grid(scaled ~ data, scales = "free") + 
  scale_y_log10() +
  labs(y = "L2 Error")

metrics %>% 
  mutate_at(vars(2:5), as.factor) %>% 
  #filter(scaled == "scaled") %>% 
  #filter(model != "pca") %>% # BC PCA error is way higher on loadings
  #filter(model != "fa") %>% 
  #filter(data != "Distinct") %>% 
  mutate(model = str_c(model, "_", sparsity),
         model = str_remove(model, "_reg")) %>% 
  mutate(data = fct_inorder(data)) %>% 
  dplyr::select(seed, data, model, sparsity, scaled, l1_loadings, l1_scores) %>%
  pivot_longer(c(l1_loadings, l1_scores)) %>% 
  mutate(name = str_to_title(str_remove(name, "l1_"))) %>%
  filter(name == "Loadings") %>% 
  ggplot(aes(x = model, y = value)) +
  geom_boxplot(aes(color = model, fill = model), 
               alpha = 0.5, outlier.shape = NA) +
  facet_grid(scaled ~ data, scales = "free") + 
  scale_y_log10() +
  labs(y = "L1 Error")

#### Tables ####
metrics %>% 
  filter(model == "bn2mf") %>% 
  group_by(data, model, sparsity, scaled) %>%
  summarise(qs = quantile(rel_err_scores, c(0.25, 0.5, 0.75), na.rm=T), prob = c(0.25, 0.5, 0.75),
            mean = mean(rel_err_scores, na.rm=T),
            max = max(rel_err_scores, na.rm=T),
            min = min(rel_err_scores, na.rm=T)) %>%
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, scaled, sparsity, model, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(data, model, scaled, sparsity)

## Metrics here do not match metrics in vci_sparse
## Need to ID differences
## FIXED!
metrics

metrics %>% 
  filter(model == "bn2mf" & scaled == "scaled") %>% 
  group_by(data, model, sparsity, scaled) %>%
  summarise(qs = quantile(rel_err_scores, c(0.25, 0.5, 0.75), na.rm=T), prob = c(0.25, 0.5, 0.75),
            mean = mean(rel_err_scores, na.rm=T),
            max = max(rel_err_scores, na.rm=T),
            min = min(rel_err_scores, na.rm=T)) %>%
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, scaled, sparsity, model, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(data, model, scaled, sparsity)

metrics %>% 
  filter(model == "bn2mf" & scaled == "scaled") %>% 
  group_by(data, model, sparsity, scaled) %>%
  summarise(qs = quantile(rel_err_loadings, c(0.25, 0.5, 0.75), na.rm=T), prob = c(0.25, 0.5, 0.75),
            mean = mean(rel_err_loadings, na.rm=T),
            max = max(rel_err_loadings, na.rm=T),
            min = min(rel_err_loadings, na.rm=T)) %>%
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, scaled, sparsity, model, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(data, model, scaled, sparsity)
