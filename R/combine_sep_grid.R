library(tidyverse)

load("./Sep/sep_out_grid.RDA")
sep_m

load("./Sims/sim_sep_grid.RDA")
sim_sep

#### Normalize truth ####
sep = bind_cols(sim_sep, sep_m) %>% 
      mutate(denom = map(true_patterns, function(x) apply(x, 1, sum)),
         patterns_scaled = map2(true_patterns, denom, function(x,y) x/y),
         scores_scaled = map2(true_scores, denom, function(x,y) as.matrix(x) %*% diag(y)),
         id = 1:nrow(.))

#### Tables ####
sep %>%
  dplyr::select(sep_num, noise_level, bn2mf_rank) %>% unnest(bn2mf_rank) %>% table()

metrics = sep %>% 
  unnest(bn2mf_rank) %>% 
  filter(bn2mf_rank == 4) %>% 
  mutate_at(vars(10:14), ~map(., as.matrix)) %>% 
  mutate(iqr  = map2(upperWA, lowerWA, function(x, y) mean(x-y)),
         prop = pmap(list(scores_scaled, lowerWA, upperWA),
                     function(x,y,z) sum(x >= y & x <= z)/(nrow(x)*ncol(x)) )) %>%
  unnest(c(prop, iqr)) 

#### Results ####
summary(metrics$iqr)
summary(metrics$prop)

prop_table = metrics %>%
              group_by(sep_num, noise_level) %>% 
              summarise(qs = quantile(prop, c(0.25, 0.5, 0.75), na.rm = TRUE), prob = c(0.25, 0.5, 0.75)) %>% 
              pivot_wider(names_from = "prob",
                          values_from = "qs") %>% 
              mutate_if(is.numeric, round, 2) %>% 
              rename(median = `0.5`) %>% 
              mutate_at(vars(1:2), as.factor)

prop_table %>% 
  ggplot(aes(x = sep_num, y = noise_level, fill = median)) +
  geom_tile()
  



