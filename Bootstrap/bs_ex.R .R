# Combine bootstrap results
# Combine VCI results
# Single examples for distinct, overlapping, and correlated simulations

####  Load packages ####
library(tidyverse)
library(patchwork)

#### Read bootstrap and VCI summaries ####
load("./Bootstrap/all_bs_vci.RDA")
all_bs_vci

#### Read bootstrap & VCI distribution output ####
# ID: Dist = 2; Over = 102; Cor = 202
type = c("dist", "over", "cor")
id = c(2,102,202)
bs_ex = tibble()

for (i in 1:length(type)) {
  load(paste0("./Bootstrap/example/bs_", id[i], "_ewa_array.RDA"))
  load(paste0("./Bootstrap/example/bs_", id[i], "_eh_array.RDA"))
  load(paste0("./Bootstrap/example/bs_", id[i], "_pred_array.RDA")) # save_pred
  
  vci_wa_dist = readMat(paste0("./Bootstrap/example/dgp_distWA_", id[i], ".mat"))[[1]]
  vci_h_dist  = readMat(paste0("./Bootstrap/example/dgp_distEH_", id[i], ".mat"))[[1]]
  
  vci_pred_dist = array(NA, c(1000, 50, 1000))
  for (j in 1:1000) {
    vci_pred_dist[,,j] = vci_wa_dist[,,j] %*% vci_h_dist[,,j]
  }
  
  bs_row = bind_cols(tibble(id = id[i]), 
                     tibble(bs_wa_dist = list(save_ewa)),
                     tibble(bs_h_dist = list(save_eh)),
                     tibble(vci_wa_dist = list(vci_wa_dist)),
                     tibble(vci_h_dist = list(vci_h_dist)),
                     tibble(bs_pred_dist = list(save_pred)), 
                     tibble(vci_pred_dist = list(vci_pred_dist)))
  
  bs_ex = bind_rows(bs_ex, bs_row)
}
rm(list = c("save_ewa", "save_eh", "bs_row", "vci_wa_dist", "vci_h_dist",
            "save_pred", "vci_pred_dist", "id"))
bs_ex

bs_ex = bs_ex %>% 
  mutate(bs_wa_dist = map(bs_wa_dist, function(x) x[,2:5,])) %>% 
  left_join(., all_bs_vci)

#### Example viz (Single Entry) ####

#### EWA ####
plot_wa = tibble()
prop_wa = tibble()

for (i in 1:3) {
  truth   = bs_ex[i,]$scores_scaled[[1]][735, 2]
  v_dist  = as.numeric(bs_ex[i,]$vci_wa_dist[[1]][735, 2,])
  bs_dist = bs_ex[i,]$bs_wa_dist[[1]][735, 2,]
  
  bs_ewa    = bs_ex[i,]$bs_wa_median[[1]][735, 2]
  bs_wa_25  = bs_ex[i,]$bs_wa_lower[[1]][735, 2]
  bs_wa_75  = bs_ex[i,]$bs_wa_upper[[1]][735, 2]
  
  v_ewa   = as.numeric(bs_ex[i,]$vci_wa_mean[[1]][735, 2])
  v_wa_25 = as.numeric(bs_ex[i,]$vci_wa_lower[[1]][735, 2])
  v_wa_75 = as.numeric(bs_ex[i,]$vci_wa_upper[[1]][735, 2])
  
  v_prop = sum(bs_ex[i,]$scores_scaled[[1]] <= bs_ex[i,]$vci_wa_upper[[1]] & 
                 bs_ex[i,]$scores_scaled[[1]] >= bs_ex[i,]$vci_wa_lower[[1]])/4000
  
  bs_prop = sum(bs_ex[i,]$scores_scaled[[1]] <= bs_ex[i,]$bs_wa_upper[[1]] & 
                  bs_ex[i,]$scores_scaled[[1]] >= bs_ex[i,]$bs_wa_lower[[1]])/4000
  
  add_plot = tibble(Distribution = v_dist) %>% 
    mutate(Type = "Variational") %>% 
    rbind(., tibble(Distribution = bs_dist) %>% 
            mutate(Type = "Bootstrap"))  %>% 
    drop_na(.) %>% 
    mutate(sim_type = str_to_title(type[i]),
           bs_ewa   = bs_ewa  ,
           bs_wa_25 = bs_wa_25,
           bs_wa_75 = bs_wa_75,
           v_ewa    = v_ewa   ,
           v_wa_25  = v_wa_25 ,
           v_wa_75  = v_wa_75,
           truth = truth)
  
  plot_wa = bind_rows(plot_wa, add_plot)
  prop_wa = bind_rows(prop_wa, bind_cols(bs = bs_prop, v = v_prop))
}
rm(list = c("bs_ewa", "bs_wa_75", "bs_wa_25",
            "v_ewa", "v_wa_25", "v_wa_75", "truth"))
prop_wa 
plot_wa %>% dplyr::select(sim_type, bs_ewa:truth) %>% 
  distinct( )

#### EWA Viz ####
wa_look = plot_wa %>% 
  mutate(sim_type = fct_inorder(sim_type)) %>% 
  ggplot(aes(x = Distribution)) +
  geom_rect(aes(xmin = v_wa_25,  xmax = v_wa_75,  ymin=0, ymax=Inf), fill="lightblue", alpha=0.025) +
  geom_rect(aes(xmin = bs_wa_25, xmax = bs_wa_75, ymin=0, ymax=Inf), fill="pink",      alpha=0.025) +
  geom_histogram(aes(y= after_stat(density),fill = Type), 
                 position = "identity", bins = 100, alpha = 0.5) +
  geom_density(aes(group = Type)) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() + 
  geom_vline(aes(xintercept = truth),  color = "darkgrey") +
  geom_vline(aes(xintercept = bs_ewa), linetype="dotted", color = "pink") + 
  geom_vline(aes(xintercept = v_ewa),  linetype="dotted", color = "lightblue") + 
  facet_grid(sim_type~., scales = "free") +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill="white")) +
  ylab("Density") + ggtitle("Distributions on Scores (WA)") + xlim(c(0,33))

#### EH ####
plot_h = tibble()
prop_h = tibble()

for (i in 1:3) {
  truth     = bs_ex[i,]$patterns_scaled[[1]][4, 34]
  v_dist_h  = as.numeric(bs_ex[i,]$vci_h_dist[[1]][4, 34,])
  bs_dist_h = bs_ex[i,]$bs_h_dist[[1]][4, 34,]
  
  bs_eh    = bs_ex[i,]$bs_h_median[[1]][4, 34]
  bs_h_25  = bs_ex[i,]$bs_h_lower[[1]][4, 34]
  bs_h_75  = bs_ex[i,]$bs_h_upper[[1]][4, 34]
  
  v_eh   = as.numeric(bs_ex[i,]$vci_h_mean[[1]][4, 34])
  v_h_25 = as.numeric(bs_ex[i,]$vci_h_lower[[1]][4, 34])
  v_h_75 = as.numeric(bs_ex[i,]$vci_h_upper[[1]][4, 34])
  
  v_prop = sum(bs_ex[i,]$patterns_scaled[[1]] <= bs_ex[i,]$vci_h_upper[[1]] & 
                 bs_ex[i,]$patterns_scaled[[1]] >= bs_ex[i,]$vci_h_lower[[1]])/4000
  
  bs_prop = sum(bs_ex[i,]$patterns_scaled[[1]] <= bs_ex[i,]$bs_h_upper[[1]] & 
                  bs_ex[i,]$patterns_scaled[[1]] >= bs_ex[i,]$bs_h_lower[[1]])/4000
  
  add_plot = tibble(Distribution = v_dist_h) %>% 
    mutate(Type = "Variational") %>% 
    rbind(., tibble(Distribution = bs_dist_h) %>% 
            mutate(Type = "Bootstrap"))  %>% 
    drop_na(.) %>% 
    mutate(sim_type = str_to_title(type[i]),
           bs_eh   = bs_eh  ,
           bs_h_25 = bs_h_25,
           bs_h_75 = bs_h_75,
           v_eh    = v_eh   ,
           v_h_25  = v_h_25 ,
           v_h_75  = v_h_75,
           truth = truth)
  
  plot_h = bind_rows(plot_h, add_plot)
  prop_h = bind_rows(prop_h, bind_cols(bs = bs_prop, v = v_prop))
}
rm(list = c("bs_eh", "bs_h_75", "bs_h_25",
            "v_eh", "v_h_25", "v_h_75", "truth"))
prop_h 
plot_h %>% dplyr::select(sim_type, bs_eh:truth) %>% 
  distinct( )

#### EH Viz ####
h_look = plot_h %>% 
  mutate(sim_type = fct_inorder(sim_type)) %>% 
  ggplot(aes(x = Distribution)) +
  geom_rect(aes(xmin = bs_h_25, xmax = bs_h_75, ymin=0, ymax=Inf), fill="pink",      alpha=0.025) +
  geom_rect(aes(xmin = v_h_25,  xmax = v_h_75,  ymin=0, ymax=Inf), fill="lightblue", alpha=0.025) +
  geom_histogram(aes(y=after_stat(density),fill = Type), 
                 position = "identity", bins = 100, alpha = 0.5) +
  geom_density(aes(group = Type)) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() + 
  facet_grid(sim_type~., scales = "free") +
  geom_vline(aes(xintercept = truth),  color = "darkgrey") +
  geom_vline(aes(xintercept = bs_eh), linetype="dotted", color = "pink") + 
  geom_vline(aes(xintercept = v_eh),  linetype="dotted", color = "lightblue") + 
  theme(legend.position = "bottom",
        strip.background = element_rect(fill="white")) +
  ylab("Density") + ggtitle("Distributions on Loadings (H)")

#### Predicted Values E[WAH] ####
plot_pred = tibble()
prop_pred = tibble()

for (i in 1:3) {
  truth   = bs_ex[i,]$chem[[1]][735, 2]
  v_dist  = as.numeric(bs_ex[i,]$vci_pred_dist[[1]][735, 2,])
  bs_dist = bs_ex[i,]$bs_pred_dist[[1]][735, 2,]
  
  bs_pred    = bs_ex[i,]$bs_pred_median[[1]][735, 2]
  bs_pred_25  = bs_ex[i,]$bs_pred_lower[[1]][735, 2]
  bs_pred_75  = bs_ex[i,]$bs_pred_upper[[1]][735, 2]
  
  v_pred   = as.numeric(bs_ex[i,]$vci_pred_mean[[1]][735, 2])
  v_pred_25 = as.numeric(bs_ex[i,]$vci_pred_lower[[1]][735, 2])
  v_pred_75 = as.numeric(bs_ex[i,]$vci_pred_upper[[1]][735, 2])
  
  v_prop = sum(bs_ex[i,]$chem[[1]] <= bs_ex[i,]$vci_pred_upper[[1]] & 
                 bs_ex[i,]$chem[[1]] >= bs_ex[i,]$vci_pred_lower[[1]])/50000
  
  bs_prop = sum(bs_ex[i,]$chem[[1]] <= bs_ex[i,]$bs_pred_upper[[1]] & 
                  bs_ex[i,]$chem[[1]] >= bs_ex[i,]$bs_pred_lower[[1]])/50000
  
  add_plot = tibble(Distribution = v_dist) %>% 
    mutate(Type = "Variational") %>% 
    rbind(., tibble(Distribution = bs_dist) %>% 
            mutate(Type = "Bootstrap"))  %>% 
    drop_na(.) %>% 
    mutate(sim_type = str_to_title(type[i]),
           bs_pred   = bs_pred  ,
           bs_pred_25 = bs_pred_25,
           bs_pred_75 = bs_pred_75,
           v_pred    = v_pred ,
           v_pred_25  = v_pred_25 ,
           v_pred_75  = v_pred_75,
           truth = truth)
  
  plot_pred = bind_rows(plot_pred, add_plot)
  prop_pred = bind_rows(prop_pred, bind_cols(bs = bs_prop, v = v_prop))
}
rm(list = c("bs_pred", "bs_pred_75", "bs_pred_25", "add_plot",
            "v_pred", "v_pred_25", "v_pred_75", "truth"))
prop_pred 
plot_pred %>% dplyr::select(sim_type, bs_pred:truth) %>% 
  distinct( )

#### E[WAH] Viz ####
pred_look = plot_pred %>% 
  mutate(sim_type = fct_inorder(sim_type)) %>% 
  ggplot(aes(x = Distribution)) +
  geom_rect(aes(xmin = v_pred_25,  xmax = v_pred_75,  ymin=0, ymax=Inf), fill="lightblue", alpha=0.025) +
  geom_rect(aes(xmin = bs_pred_25, xmax = bs_pred_75, ymin=0, ymax=Inf), fill="pink",      alpha=0.025) +
  geom_histogram(aes(y = after_stat(density), fill = Type), 
                 position = "identity", bins = 100, alpha = 0.5) +
  geom_density(aes(group = Type)) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() + 
  geom_vline(aes(xintercept = truth),  linetype="dotted", color = "black") +
  geom_vline(aes(xintercept = bs_pred), linetype="dotted", color = "pink") + 
  geom_vline(aes(xintercept = v_pred),  linetype="dotted", color = "lightblue") + 
  facet_grid(sim_type~., scales = "free") +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill="white")) +
  ylab("Density") + ggtitle("Distributions on Predicted Values E[WAH]") + xlim(c(0, 17))

#### Viz ####
(wa_look + theme(legend.position = "none")) / h_look



