# Combine bootstrap results
# Combine VCI results
# Single examples for distinct, overlapping, and correlated simulations
# Overall charactersistics, too

####  Load packages ####
library(tidyverse)
library(R.matlab)
library(openssl)

# create random vector for matlab on hpc
rand_vec = rand_num(150) * ((2^31)-1)

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
load("./Bootstrap/Compare/bs_list_lower_wa.RDA")
load("./Bootstrap/Compare/bs_list_upper_wa.RDA")
load("./Bootstrap/Compare/bs_list_mean_wa.RDA")
load("./Bootstrap/Compare/bs_list_median_wa.RDA")
load("./Bootstrap/Compare/bs_list_lower_h.RDA")
load("./Bootstrap/Compare/bs_list_upper_h.RDA")
load("./Bootstrap/Compare/bs_list_mean_h.RDA")
load("./Bootstrap/Compare/bs_list_median_h.RDA")
                                                       
bs_lists = tibble(bs_h_lower = bs_list_lower_h,
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
dgp_vci <- dgp_vci %>% 
           dplyr::select(vci_h_mean = eh_scaled,
                   vci_wa_mean = ewa_scaled,
                   vci_wa_lower = lowerWA,
                   vci_wa_upper = upperWA,
                   vci_h_lower = lowerH,
                   vci_h_upper = upperH)

#### Clean bootstrap data ####
bs_dat = bind_cols(sim_dgp, bs_lists, dgp_vci) %>% 
  mutate(bs_wa_lower  = map(bs_wa_lower,  function(x) x[,2:5]),
         bs_wa_mean   = map(bs_wa_mean,   function(x) x[,2:5]),
         bs_wa_upper  = map(bs_wa_upper,  function(x) x[,2:5]),
         bs_wa_median = map(bs_wa_median, function(x) x[,2:5])) %>% 
  pivot_longer(c(grep("bs_", colnames(.)), grep("vci_", colnames(.))),
               names_to = c("method", "side", "matrix"),
               names_sep = "_") %>%
  mutate(value = map(value, as.matrix)) %>% 
  pivot_wider(names_from = matrix,
              values_from = value) %>% 
  mutate(truth = pmap(list(side, true_patterns, true_scores),
                      function(x, y, z) if(x == "h") {y} else {z}),
         best = map2(median, mean, function(x,y) if(!is.null(x)) {x} else {y}),
         err  = map2(truth, best, function(x,y)
           if (ncol(x) == ncol(y)) {norm(x-y, "F")/norm(x, "F")} else {NA}),
         iqr  = map2(upper, lower, function(x, y) mean(x-y)),
         prop = pmap(list(truth, lower, upper), 
                           function(x,y,z) sum(x >= y & x <= z)/(nrow(x)*ncol(x)))) %>% 
  unnest(c(err, prop, iqr))

bs_dat %>% 
  dplyr::select(data, method, side, id, err, iqr, prop)

#### Results Tables ####
bs_dat %>% 
  group_by(data, method, side) %>% 
  summarise(qs = quantile(err, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(err),
            max = max(err),
            min = min(err)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, min, `0.25`, `0.5`, mean, `0.75`, max)

bs_dat %>% 
  group_by(data, method, side) %>% 
  summarise(qs = quantile(prop, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(prop),
            max = max(prop),
            min = min(prop)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, min, `0.25`, `0.5`, mean, `0.75`, max)

bs_dat %>% 
  group_by(data, method, side) %>% 
  summarise(qs = quantile(iqr, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(iqr),
            max = max(iqr),
            min = min(iqr)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, min, `0.25`, `0.5`, mean, `0.75`, max)

#### Single Examples #####

#### Read bootstrap output ####
# ID: Dist = 2; Over = 102; Cor = 202
type = c("dist", "over", "cor")
matrix = c("_lower_", "_upper_", "_mean_", "_median_")
bs_ex_list = list(dist = tibble(id = 2), over = tibble(id = 102), cor = tibble(id = 202))

for (i in 1:length(type)) {
  bs_names = c("id")
  for (j in matrix) {
    load(paste0("./Bootstrap/Example Distributions/bs_", type[i], j, "wa.RDA"))
    load(paste0("./Bootstrap/Example Distributions/bs_", type[i], j, "h.RDA"))
    
    bs_ex_list[[i]] = bind_cols(bs_ex_list[[i]], tibble(list(get(paste0("bs", j, "wa"))),
                                                        list(get(paste0("bs", j, "h")))))
    bs_names = c(bs_names, paste0("bs", j, "wa"), paste0("bs", j, "h"))
  }
  load(paste0("./Bootstrap/Example Distributions/bs_", type[i], "_ewa.RDA"))
  load(paste0("./Bootstrap/Example Distributions/bs_", type[i], "_eh.RDA"))
  bs_names = c(bs_names, paste0("bs_ewa_dist"), paste0("bs_eh_dist"))
  
  bs_ex_list[[i]] = bind_cols(bs_ex_list[[i]], tibble(list(bs_ewa),
                                                      list(bs_eh)))
  names(bs_ex_list[[i]]) <- bs_names
}

bs_ex = bs_ex_list[[1]]
for (i in 2:length(bs_ex_list)) {
  bs_ex = bind_rows(bs_ex, bs_ex_list[[i]])
}

bs_ex = bs_ex %>% 
  mutate_at(vars(matches("_wa")), function(x) map(x, function(y) y[,2:5])) %>% 
  mutate(bs_ewa_dist = map(bs_ewa_dist, function(x) x[,2:5,]))
bs_ex

#### Read VCI data ####
load("./Results/Main/vci_out.RDA")
dgp_vci <- dgp_vci %>% 
  dplyr::select(vci_mean_h = eh_scaled,
                vci_mean_wa = ewa_scaled,
                vci_lower_wa = lowerWA,
                vci_upper_wa = upperWA,
                vci_lower_h = lowerH,
                vci_upper_h = upperH) %>% 
  bind_cols(., dgp_dat)

ex = bs_ex %>% left_join(., dgp_vci) %>% 
  pivot_longer(c(grep("bs_", colnames(.)), grep("vci_", colnames(.))),
               names_to = c("method", "matrix", "side"),
               names_sep = "_") %>%
  mutate(value = map(value, as.matrix)) %>% 
  pivot_wider(names_from = matrix,
              values_from = value) %>% 
  mutate(truth = pmap(list(side, true_patterns, true_scores),
                      function(x, y, z) if(x == "h") {y} else {z}),
         best = map2(median, mean, function(x,y) if(!is.null(x)) {x} else {y}),
         err  = map2(truth, best, function(x,y)
           if (ncol(x) == ncol(y)) {norm(x-y, "F")/norm(x, "F")} else {NA}),
         iqr  = map2(upper, lower, function(x, y) mean(x-y)),
         prop = pmap(list(truth, lower, upper), 
                     function(x,y,z) sum(x >= y & x <= z)/(nrow(x)*ncol(x)))) %>% 
  unnest(c(err, prop, iqr))
#### Example viz (Single Entry) ####

truth = bs_ex[1,]$scores_scaled[[1]][10,4]

bs_ewa    = bs_ex[1,]$bs_median_wa[[1]][10,4]
bs_dist   = bs_ex[1,]$bs_ewa_dist[[1]][10,4,]
bs_wa_25  = bs_ex[1,]$bs_lower_wa[[1]][10,4]
bs_wa_75  = bs_ex[1,]$bs_upper_wa[[1]][10,4]

# v_ewa = v_EWA[10,4]
# v_dist_1 = v_dist[10,4,]
# v_q25 = v_lower[10,4]
# v_q75 = v_upper[10,4]

#tibble(Distribution = v_dist_1) %>% 
#mutate(Type = "Variational") %>% 
#rbind(., 
tibble(Distribution = bs_dist) %>% 
  
  mutate(Type = "Bootstrap"
         #             )
  )  %>% 
  drop_na(.) %>% 
  ggplot(aes(x = Distribution)) +
  geom_rect(aes(xmin=bs_wa_25, xmax=bs_wa_75, ymin=0, ymax=Inf), fill="pink",      alpha=0.05) +
  #geom_rect(aes(xmin=v_q25,  xmax=v_q75,  ymin=0, ymax=Inf), fill="lightblue", alpha=0.05) +
  #geom_histogram(aes(y=..density.., group = Type, fill = Type), alpha = 0.5, bins = 100) +
  geom_density(aes(group = Type, fill = Type), alpha = 0.5) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() + 
  geom_vline(xintercept = truth,            linetype="dotted", color = "black") +
  geom_vline(xintercept = c(#v_ewa, 
    bs_ewa), linetype="dotted", color = "yellow") + 
  ylab("Density")


#sum(scores_scaled <= v_upper & scores_scaled >= v_lower)/4000

sum(bs_ex[1,]$scores_scaled[[1]] <= bs_ex[1,]$bs_upper_wa[[1]] & 
      bs_ex[1,]$scores_scaled[[1]] >= bs_ex[1,]$bs_lower_wa[[1]])/4000














  
  