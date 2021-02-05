## Compare scaled error wih non-scaled error

metrics
# non-scaled error

see_pred
all_dgp
# scaled error

metrics %>% 
  filter(model == "BN2MF") %>% 
  dplyr::select(1:7, -model, -rank) %>% 
  group_by(data) %>% 
  summarise(qs = quantile(rel_err_all, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(rel_err_all),
            max = max(rel_err_all),
            min = min(rel_err_all)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, min, `0.25`, `0.5`, mean, `0.75`, max)

see_pred %>% 
  filter(method == "vci") %>% 
  group_by(data, method, side) %>% 
  summarise(qs = quantile(err, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(err),
            max = max(err),
            min = min(err)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(side, data, method, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(side, data, method)
# prediction error is the same, duh

metrics %>% 
  filter(model == "BN2MF") %>% 
  dplyr::select(1:7, -model, -rank) %>% 
  group_by(data) %>% 
  summarise(qs = quantile(rel_err_loadings, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(rel_err_loadings),
            max = max(rel_err_loadings),
            min = min(rel_err_loadings)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, min, `0.25`, `0.5`, mean, `0.75`, max)

metrics %>% 
  filter(model == "BN2MF") %>% 
  dplyr::select(1:7, -model, -rank) %>% 
  group_by(data) %>% 
  summarise(qs = quantile(rel_err_scores, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(rel_err_scores),
            max = max(rel_err_scores),
            min = min(rel_err_scores)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, min, `0.25`, `0.5`, mean, `0.75`, max)

all_dgp %>% 
  filter(method == "vci") %>% 
  group_by(data, method, side) %>% 
  summarise(qs = quantile(err, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(err),
            max = max(err),
            min = min(err)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(side, data, method, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(side, data, method)

