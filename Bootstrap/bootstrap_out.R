#####  Load packages
library(tidyverse)
library(R.matlab)
library(openssl)

# create random vector for matlab on hpc
rand_vec = rand_num(150) * ((2^31)-1)

#####   
##### Single Sim Examples
#####

# Read data
load("./Results/Main/sim_dgp_rep1.RDA")

# Normalize truth
dgp_dat = sim_dgp_rep1 %>% 
  mutate(denom = map(true_patterns, function(x) apply(x,1, sum)),
         patterns_scaled = map2(true_patterns, denom, function(x,y) x/y),
         scores_scaled = map2(true_scores, denom, function(x,y) as.matrix(x) %*% diag(y)),
         id = 1:nrow(.))

# v_EWA <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_EWA.mat")[[1]]
# v_EH <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_EH.mat")[[1]]
# v_upper <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_WA_upper.mat")[[1]]
# v_lower <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_WA_lower.mat")[[1]]
# v_dist <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_var_dist_WA.mat")[[1]]

# Read bootstrapped BN2MF EWA output
# ID: Dist = 2; Over = 102; Cor = 202
type = c("dist", "over", "cor")
matrix = c("_lower_", "_upper_", "_mean_", "_median_")
bs_ex_list = list(dist = tibble(id = 2), over = tibble(id = 102), cor = tibble(id = 202))

for (i in 1:length(type)) {
  bs_names = c("id")
  for (j in matrix) {
  load(paste0("./Bootstrap/Compare/bs_", type[i], j, "wa.RDA"))
  load(paste0("./Bootstrap/Compare/bs_", type[i], j, "h.RDA"))
  
  bs_ex_list[[i]] = bind_cols(bs_ex_list[[i]], tibble(list(get(paste0("bs_", type[i], j, "wa"))),
                                                      list(get(paste0("bs_", type[i], j, "h")))))
  bs_names = c(bs_names, paste0("bs", j, "wa"), paste0("bs", j, "h"))
  }
  load(paste0("./Bootstrap/Compare/bs_", type[i], "_ewa.RDA"))
  load(paste0("./Bootstrap/Compare/bs_", type[i], "_eh.RDA"))
  bs_names = c(bs_names, paste0("bs_ewa_dist"), paste0("bs_eh_dist"))
  
  bs_ex_list[[i]] = bind_cols(bs_ex_list[[i]], tibble(list(get(paste0("bs_", type[i], "_ewa"))),
                                                      list(get(paste0("bs_", type[i], "_eh")))))
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

## Example viz (Single Entry)
truth = bs_ex[1,]$scores_scaled[[1]][10,4]

# v_ewa = v_EWA[10,4]
# v_dist_1 = v_dist[10,4,]
# v_q25 = v_lower[10,4]
# v_q75 = v_upper[10,4]

bs_ewa    = bs_ex[1,]$bs_median_wa[[1]][10,4]
bs_dist   = bs_ex[1,]$bs_ewa_dist[[1]][10,4,]
bs_wa_25  = bs_ex[1,]$bs_lower_wa[[1]][10,4]
bs_wa_75  = bs_ex[1,]$bs_upper_wa[[1]][10,4]

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

sum(scores_scaled <= v_upper & scores_scaled >= v_lower)/4000

sum(bs_ex[1,]$scores_scaled[[1]] <= bs_ex[1,]$bs_upper_wa[[1]] & 
      bs_ex[1,]$scores_scaled[[1]] >= bs_ex[1,]$bs_lower_wa[[1]])/4000
  
#####   
##### 300 Sims
#####

# Read data
load("./Bootstrap/Compare/bs_list_lower_wa.RDA")
load("./Bootstrap/Compare/bs_list_upper_wa.RDA")
load("./Bootstrap/Compare/bs_list_mean_wa.RDA")
load("./Bootstrap/Compare/bs_list_median_wa.RDA")
load("./Bootstrap/Compare/bs_list_lower_h.RDA")
load("./Bootstrap/Compare/bs_list_upper_h.RDA")
load("./Bootstrap/Compare/bs_list_mean_h.RDA")
load("./Bootstrap/Compare/bs_list_median_h.RDA")

bs_lists = tibble(h_lower = bs_list_lower_h,
              h_mean = bs_list_mean_h,
              h_upper = bs_list_upper_h,
              h_median = bs_list_median_h,
              wa_lower = bs_list_lower_wa,
              wa_mean = bs_list_mean_wa,
              wa_upper = bs_list_upper_wa,
              wa_median = bs_list_median_wa)

bs_dat = bind_cols(bs_lists, dgp_dat) %>% 
  mutate(wa_lower = map(wa_lower, function(x) x[,2:5]),
         wa_mean = map(wa_mean, function(x) x[,2:5]),
         wa_upper = map(wa_upper, function(x) x[,2:5]),
         wa_median = map(wa_median, function(x) x[,2:5]),
         h_iqr = map2(h_upper, h_lower, function(x, y) x - y),
         wa_iqr = map2(wa_upper, wa_lower, function(x, y) x - y),
         loading_err = map2(patterns_scaled, h_median, function(x,y)
           if (ncol(x) == ncol(y)) {norm(x-y, "F")/norm(x, "F")} else {NA}),
         score_err = map2(scores_scaled, wa_median, function(x,y)
           if (ncol(x) == ncol(y)) {norm(x-y, "F")/norm(x, "F")} else {NA}),
         wa_prop = pmap(list(scores_scaled, wa_lower, wa_upper), 
                     function(x,y,z) sum(x >= y & x <= z)/4000),
         h_prop = pmap(list(patterns_scaled, h_lower, h_upper), 
                        function(x,y,z) sum(x >= y & x <= z)/(4*50))) %>% 
  unnest(c(loading_err, score_err, wa_prop, h_prop))

bs_dat %>% 
  group_by(data) %>% 
  summarise(qs = quantile(loading_err, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(loading_err),
            max = max(loading_err),
            min = min(loading_err)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, min, `0.25`, `0.5`, mean, `0.75`, max)

bs_dat %>% 
  group_by(data) %>% 
  summarise(qs = quantile(score_err, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(score_err),
            max = max(score_err),
            min = min(score_err)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, min, `0.25`, `0.5`, mean, `0.75`, max)

bs_dat %>% 
  group_by(data) %>% 
  summarise(qs = quantile(wa_prop, c(0.25, 0.5, 0.75)), prob = c(0.25, 0.5, 0.75),
            mean = mean(wa_prop),
            max = max(wa_prop),
            min = min(wa_prop)) %>% 
  pivot_wider(names_from = "prob",
              values_from = "qs") %>% 
  dplyr::select(data, min, `0.25`, `0.5`, mean, `0.75`, max)




















  
  