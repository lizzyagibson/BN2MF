# PCA scree plots

library(tidyverse)
library(psych)
library(NMF)
library(R.matlab)
options(scipen = 999)
source("./R/compare_functions.R")

# Read in Sims
load("./Sims/sim_dgp_local.RDA")
sim_dgp_local

sim <- sim_dgp_local$sim[[1]]

get_scree <- function (sim) {
  # Run PCA centered, not scaled
  pca_out <- prcomp(sim)
  sv <- as_tibble(as.matrix(pca_out$sdev, ncol = 1)) %>% 
    mutate(Comp = 1:nrow(.))
  
  return(sv)
}

get_scree(sim) %>% 
  ggplot(aes(x = Comp, y = V1)) +
  geom_point() + geom_line() +
  theme_bw()

pve <- sv$V1^2/sum(sv$V1^2)

sim_dgp_local %>% 
  mutate(sv = map(sim, get_scree)) %>% 
  select(seed, data, sv) %>% 
  unnest(sv) %>% 
  filter(Comp %in% 1:6) %>% 
  mutate(seed = as.factor(seed)) %>% 
  ggplot(aes(x = Comp, y = V1, color = seed, group = interaction(seed, data))) +
  geom_point(alpha = 0.5) + geom_line(alpha = 0.5) +
  theme_bw() + 
  theme(legend.position = "none")

sim_dgp_local %>% 
  mutate(sv = map(sim, get_scree)) %>% 
  select(seed, data, sv) %>% 
  unnest(sv) %>% 
  group_by(seed, data) %>% 
  mutate(denom = sum(V1^2),
            percent = V1^2 / denom) %>%
  filter(Comp %in% 1:6) %>% 
  mutate(seed = as.factor(seed)) %>% 
  ggplot(aes(x = Comp, y = percent, color = seed, group = interaction(seed, data))) +
  geom_point(alpha = 0.5) + geom_line(alpha = 0.5) +
  theme_bw() + 
  theme(legend.position = "none")

