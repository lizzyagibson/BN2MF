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

#### Single Examples #####

#### Read data ####
load("./Sims/Main/sim_dgp.RDA")

#### Normalize truth ####
dgp_dat = sim_dgp %>% 
  mutate(denom = map(true_patterns, function(x) apply(x,1, sum)),
         patterns_scaled = map2(true_patterns, denom, function(x,y) x/y),
         scores_scaled = map2(true_scores, denom, function(x,y) as.matrix(x) %*% diag(y)),
         id = 1:nrow(.))

# v_EWA <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_EWA.mat")[[1]]
# v_EH <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_EH.mat")[[1]]
# v_upper <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_WA_upper.mat")[[1]]
# v_lower <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_WA_lower.mat")[[1]]
# v_dist <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_var_dist_WA.mat")[[1]]

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

bs_ex = bs_ex %>% left_join(., dgp_dat, by = "id") %>% 
  mutate_at(vars(matches("_wa")), function(x) map(x, function(y) y[,2:5])) %>% 
  mutate(bs_ewa_dist = map(bs_ewa_dist, function(x) x[,2:5,]))
bs_ex

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
  
#### For all sims ####

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

#### Clean bootstrap data ####
bs_dat = bind_cols(bs_lists, dgp_dat) %>% 
         mutate(bs_wa_lower  = map(bs_wa_lower, function(x) x[,2:5]),
                 bs_wa_mean   = map(bs_wa_mean, function(x) x[,2:5]),
                 bs_wa_upper  = map(bs_wa_upper, function(x) x[,2:5]),
                 bs_wa_median = map(bs_wa_median, function(x) x[,2:5]),
                 bs_h_iqr  = map2(bs_h_upper, bs_h_lower, function(x, y) mean(x - y)),
                 bs_wa_iqr = map2(bs_wa_upper, bs_wa_lower, function(x, y) mean(x - y)),
                 bs_loading_err = map2(patterns_scaled, bs_h_median, function(x,y)
                      if (ncol(x) == ncol(y)) {norm(x-y, "F")/norm(x, "F")} else {NA}),
                 bs_score_err = map2(scores_scaled, bs_wa_median, function(x,y)
                      if (ncol(x) == ncol(y)) {norm(x-y, "F")/norm(x, "F")} else {NA}),
                 bs_wa_prop = pmap(list(scores_scaled, bs_wa_lower, bs_wa_upper), 
                      function(x,y,z) sum(x >= y & x <= z)/4000),
                 bs_h_prop = pmap(list(patterns_scaled, bs_h_lower, bs_h_upper), 
                      function(x,y,z) sum(x >= y & x <= z)/(4*50)),
                 bs_pred = map2(bs_wa_median, bs_h_median, function(x,y) x %*% y),
                 bs_err = map2(chem, bs_pred, function(x,y) norm(x-y, "F")/norm(x, "F"))) %>% 
         unnest(c(bs_loading_err, bs_score_err, bs_wa_prop, bs_h_prop, bs_err, bs_wa_iqr, bs_h_iqr))

#### Results Tables ####
bs_dat %>% 
  group_by(data) %>% 
  summarise(qs = quantile(bs_err, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(bs_err),
            max = max(bs_err),
            min = min(bs_err)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, min, `0.25`, `0.5`, mean, `0.75`, max)

bs_dat %>% 
  group_by(data) %>% 
  summarise(qs = quantile(bs_loading_err, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(bs_loading_err),
            max = max(bs_loading_err),
            min = min(bs_loading_err)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, min, `0.25`, `0.5`, mean, `0.75`, max)

bs_dat %>% 
  group_by(data) %>% 
  summarise(qs = quantile(bs_score_err, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(bs_score_err),
            max = max(bs_score_err),
            min = min(bs_score_err)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, min, `0.25`, `0.5`, mean, `0.75`, max)

bs_dat %>% 
  group_by(data) %>% 
  summarise(qs = quantile(bs_wa_prop, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(bs_wa_prop),
            max = max(bs_wa_prop),
            min = min(bs_wa_prop)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, min, `0.25`, `0.5`, mean, `0.75`, max)

bs_dat %>% 
  group_by(data) %>% 
  summarise(qs = quantile(bs_wa_iqr, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(bs_wa_iqr),
            max = max(bs_wa_iqr),
            min = min(bs_wa_iqr)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, min, `0.25`, `0.5`, mean, `0.75`, max)















  
  