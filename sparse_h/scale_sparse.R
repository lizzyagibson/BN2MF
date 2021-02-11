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
load("./sparse_h/sparse_main_out.RDA")
dgp_m = dgp_m %>% 
  rename_all(~str_c(., "_sh"))

load("./Results/main/main_out.RDA")
dgp = dgp %>% 
  rename_at(vars(7:34), ~str_c(., "_reg")) %>% 
  mutate_at(grep("(a|2|p)_loadings", colnames(.)), function(x) map(x, t))

dgp_all = bind_cols(dgp, dgp_m)
dgp_all

# Calculate error metrics
dgp_metrics0 = dgp_all %>% 
              pivot_longer(pca_out_reg:bn2mf_rank_sh,
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
  dplyr::select(-denom, -true_denom) 
  # %>% 
  # pivot_longer(c(reg_loadings, reg_scores, scaled_loadings, scaled_scores),
  #              names_to = c("scaled", "object"),
  #              names_sep = "_")

# l2 relative error -- loadings & scores
rel_err <- dgp_metrics %>%
                mutate(l2_loadings_reg    = map2(true_patterns, reg_loadings, get_relerror),
                       l2_loadings_scaled = map2(true_patterns_scaled, scaled_loadings, get_relerror),
                       l1_loadings_reg    = map2(true_patterns, reg_loadings, get_relerror),
                       l1_loadings_scaled = map2(true_patterns_scaled, scaled_loadings, get_relerror_l1),
                       l2_scores_reg    = map2(true_scores, reg_scores, get_relerror),
                       l2_scores_scaled = map2(true_scores_scaled, scaled_scores, get_relerror),
                       l1_scores_reg    = map2(true_scores, reg_scores, get_relerror),
                       l1_scores_scaled = map2(true_scores_scaled, scaled_scores, get_relerror_l1)) %>%
                dplyr::select(seed, data, model, sparsity, grep("(l1|l2)", colnames(.))) %>% 
                unnest(c(grep("(l1|l2)", colnames(.))))

# Combine all metrics
metrics <- rel_err %>% 
           pivot_longer(grep("(l1|l2)", colnames(.)),
                        names_to = c("err_type", "name", "scaled"),
                        names_sep = "_") %>% 
           mutate_at(vars(2:7), as.factor)
metrics

#### Viz ####
metrics %>% 
  #filter(scaled == "scaled") %>% 
  filter(model != "pca") %>% # BC PCA error is way higher on loadings
  filter(model != "fa") %>% 
  #filter(data != "Distinct") %>% 
  mutate(model = str_remove(str_c(model, "_", sparsity), "_reg")) %>% 
  mutate(data = fct_inorder(data)) %>% 
  filter(name == "loadings" & err_type == "l2") %>% 
  ggplot(aes(x = model, y = value)) +
  geom_boxplot(aes(color = model, fill = model), 
               alpha = 0.5, outlier.shape = NA) +
  facet_grid(scaled ~ data, scales = "free") + 
  scale_y_log10() +
  labs(y = "L2 Error")

#### Tables ####
metrics %>% 
  filter(model == "bn2mf") %>% 
  group_by(data, model, sparsity, err_type, name, scaled) %>%
  summarise(qs = quantile(value, c(0.25, 0.5, 0.75), na.rm=T), prob = c(0.25, 0.5, 0.75),
            mean = mean(value, na.rm=T),
            max = max(value, na.rm=T),
            min = min(value, na.rm=T)) %>%
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(err_type, name, data, model, scaled, sparsity, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(err_type, name, data, model, scaled, sparsity)

## Metrics here do not match metrics in vci_sparse
## Need to ID differences
## FIXED!
metrics

metrics %>% 
  filter(model == "bn2mf" & scaled == "scaled") %>% 
  group_by(data, model, sparsity, scaled) %>%
  summarise(qs = quantile(l2_scores, c(0.25, 0.5, 0.75), na.rm=T), prob = c(0.25, 0.5, 0.75),
            mean = mean(l2_scores, na.rm=T),
            max = max(l2_scores, na.rm=T),
            min = min(l2_scores, na.rm=T)) %>%
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, scaled, sparsity, model, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(data, model, scaled, sparsity)

metrics %>% 
  filter(model == "bn2mf" & scaled == "scaled") %>% 
  group_by(data, model, sparsity, scaled) %>%
  summarise(qs = quantile(l2_loadings, c(0.25, 0.5, 0.75), na.rm=T), prob = c(0.25, 0.5, 0.75),
            mean = mean(l2_loadings, na.rm=T),
            max = max(l2_loadings, na.rm=T),
            min = min(l2_loadings, na.rm=T)) %>%
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, scaled, sparsity, model, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(data, model, scaled, sparsity)
