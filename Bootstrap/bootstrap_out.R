# Combine bootstrap results
# Combine VCI results
# Single examples for distinct, overlapping, and correlated simulations
# Overall charactersistics, too

####  Load packages ####
library(tidyverse)
library(R.matlab)
library(openssl)

# create random vector for matlab on hpc
# rand_vec = rand_num(150) * ((2^31)-1)

#### For all sims ####

#### Read data ####
load("./Sims/Main/sim_dgp.RDA")

#### Normalize truth ####
sim_dgp = sim_dgp %>% 
  mutate(denom = map(true_patterns, function(x) apply(x,1, sum)),
         patterns_scaled = map2(true_patterns, denom, function(x,y) x/y),
         scores_scaled = map2(true_scores, denom, function(x,y) as.matrix(x) %*% diag(y)),
         id = 1:nrow(.))

#### Read bootstrap data ####
load("./Bootstrap/compare/bs_list_lower_wa.RDA")
load("./Bootstrap/compare/bs_list_upper_wa.RDA")
load("./Bootstrap/compare/bs_list_mean_wa.RDA")
load("./Bootstrap/compare/bs_list_median_wa.RDA")
load("./Bootstrap/compare/bs_list_lower_h.RDA")
load("./Bootstrap/compare/bs_list_upper_h.RDA")
load("./Bootstrap/compare/bs_list_mean_h.RDA")
load("./Bootstrap/compare/bs_list_median_h.RDA")
                                                       
bs_dgp = tibble(bs_h_lower = bs_list_lower_h,
                  bs_h_mean = bs_list_mean_h,
                  bs_h_upper = bs_list_upper_h,
                  bs_h_median = bs_list_median_h,
                  bs_wa_lower = bs_list_lower_wa,
                  bs_wa_mean = bs_list_mean_wa,
                  bs_wa_upper = bs_list_upper_wa,
                  bs_wa_median = bs_list_median_wa)

rm(list =  c("bs_list_lower_h",
             "bs_list_mean_h",
             "bs_list_upper_h",
             "bs_list_median_h",
             "bs_list_lower_wa",
             "bs_list_mean_wa",
             "bs_list_upper_wa",
             "bs_list_median_wa"))

#### Read VCI data ####
load("./Results/Main/vci_out.RDA")
vci_dgp <- dgp_vci %>% 
           dplyr::select(vci_h_mean = eh_scaled,
                   vci_wa_mean = ewa_scaled,
                   vci_wa_lower = lowerWA,
                   vci_wa_upper = upperWA,
                   vci_h_lower = lowerH,
                   vci_h_upper = upperH)
rm(dgp_vci)

#### Clean bootstrap data ####
bs_dgp = bs_dgp %>% 
  mutate(bs_wa_lower  = map(bs_wa_lower,  function(x) x[,2:5]),
         bs_wa_mean   = map(bs_wa_mean,   function(x) x[,2:5]),
         bs_wa_upper  = map(bs_wa_upper,  function(x) x[,2:5]),
         bs_wa_median = map(bs_wa_median, function(x) x[,2:5])) 

#### Combine data ####
all_dgp = bind_cols(sim_dgp, bs_dgp, vci_dgp) %>% 
  pivot_longer(c(grep("bs_", colnames(.)), grep("vci_", colnames(.))),
               names_to = c("method", "side", "matrix"),
               names_sep = "_") %>%
  mutate(value = map(value, as.matrix)) %>% 
  pivot_wider(names_from = matrix,
              values_from = value) %>% 
  mutate(truth = pmap(list(side, patterns_scaled, scores_scaled),
                      function(x, y, z) if(x == "h") {y} else {z}),
         best = map2(median, mean, function(x,y) if(!is.null(x)) {x} else {y}),
         err  = map2(truth, best, function(x,y)
           if (ncol(x) == ncol(y)) {norm(x-y, "F")/norm(x, "F")} else {NA}),
         iqr  = map2(upper, lower, function(x, y) mean(x-y)),
         prop = pmap(list(truth, lower, upper),
                           function(x,y,z) sum(x >= y & x <= z)/(nrow(x)*ncol(x)) )) %>%
  unnest(c(err, prop, iqr))

#### Results Tables ####
all_dgp %>% 
  group_by(data, method, side) %>%
  summarise(qs = quantile(prop, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(prop),
            max = max(prop),
            min = min(prop)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(side, data, method, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(side, data, method)

all_dgp %>% 
  group_by(data, method, side) %>% 
  summarise(qs = quantile(err, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(err),
            max = max(err),
            min = min(err)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(side, data, method, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(side, data, method)

all_dgp %>% 
  group_by(data, method, side) %>% 
  summarise(qs = quantile(iqr, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(iqr),
            max = max(iqr),
            min = min(iqr)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(side, data, method, min, `0.25`, `0.5`, mean, `0.75`, max) %>% 
  arrange(side, data, method)

#### Single Examples #####

#### Read bootstrap distribution output ####
#### Read VCI distribution output ####

# ID: Dist = 2; Over = 102; Cor = 202
type = c("dist", "over", "cor")
id = c(2,102,202)
bs_ex = tibble()

for (i in 1:length(type)) {
  load(paste0("./Bootstrap/example/bs_", id[i], "_ewa_array.RDA"))
  load(paste0("./Bootstrap/example/bs_", id[i], "_eh_array.RDA"))

  vci_wa_dist = readMat(paste0("./Bootstrap/example/dgp_distWA_", id[i], ".mat"))[[1]]
  vci_h_dist  = readMat(paste0("./Bootstrap/example/dgp_distEH_", id[i], ".mat"))[[1]]

  bs_row = bind_cols(tibble(id = id[i]), 
           tibble(bs_wa_dist = list(save_ewa)),
           tibble(bs_h_dist = list(save_eh)),
           tibble(vci_wa_dist = list(vci_wa_dist)),
           tibble(vci_h_dist = list(vci_h_dist)))
    
  bs_ex = bind_rows(bs_ex, bs_row)
}
rm(list = c("save_ewa", "save_eh", "bs_row", "vci_wa_dist", "vci_h_dist"))

bs_ex = bs_ex %>% 
  mutate(bs_wa_dist = map(bs_wa_dist, function(x) x[,2:5,])) %>% 
  left_join(., bind_cols(bs_dgp, sim_dgp, vci_dgp))

#### Example viz (Single Entry) ####
#### Loop for each kind of simulation 

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

plot_wa %>% 
  mutate(sim_type = fct_inorder(sim_type)) %>% 
  ggplot(aes(x = Distribution, fill = Type)) +
  geom_rect(aes(xmin = v_wa_25,  xmax = v_wa_75,  ymin=0, ymax=Inf), fill="lightblue", alpha=0.025) +
  geom_rect(aes(xmin = bs_wa_25, xmax = bs_wa_75, ymin=0, ymax=Inf), fill="pink",      alpha=0.025) +
  geom_histogram(aes(y=..density..,fill = Type), bins = 100) +
  #geom_density(aes(group = Type), alpha = 0.5) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() + 
  geom_vline(aes(xintercept = truth),  linetype="dotted", color = "black") +
  geom_vline(aes(xintercept = bs_ewa), linetype="dotted", color = "yellow") + 
  geom_vline(aes(xintercept = v_ewa),  linetype="dotted", color = "yellow") + 
  facet_grid(sim_type~., scales = "free") +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill="white")) +
  ylab("Density")

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

plot_h %>% 
  mutate(sim_type = fct_inorder(sim_type)) %>% 
  ggplot(aes(x = Distribution, fill = Type)) +
  geom_rect(aes(xmin = v_h_25,  xmax = v_h_75,  ymin=0, ymax=Inf), fill="lightblue", alpha=0.025) +
  geom_rect(aes(xmin = bs_h_25, xmax = bs_h_75, ymin=0, ymax=Inf), fill="pink",      alpha=0.025) +
  geom_histogram(aes(y=..density..,fill = Type), bins = 100) +
  #geom_density(aes(group = Type), alpha = 0.5) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() + 
  geom_vline(aes(xintercept = truth),  linetype="dotted", color = "black") +
  geom_vline(aes(xintercept = bs_eh), linetype="dotted", color = "yellow") + 
  geom_vline(aes(xintercept = v_eh),  linetype="dotted", color = "yellow") + 
  facet_grid(sim_type~., scales = "free") +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill="white")) +
  ylab("Density")




  
  