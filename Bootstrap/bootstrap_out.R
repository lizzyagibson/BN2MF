#####  Load packages
library(tidyverse)
library(R.matlab)
library(openssl)
library(plotly)

# create random vector for matlab on hpc
rand_vec = rand_num(150) * ((2^31)-1)

#####   
##### Single Sim Example
#####

# Read data
chem <- read_csv("./Results/Main/Corr Ex/corr_chem.csv")
sim <- read_csv("./Results/Main/Corr Ex/corr_sim.csv")
patterns <- read_csv("./Results/Main/Corr Ex/corr_patterns.csv")
scores <- read_csv("./Results/Main/Corr Ex/corr_scores.csv")

v_EWA <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_EWA.mat")[[1]]
v_EH <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_EH.mat")[[1]]
v_upper <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_WA_upper.mat")[[1]]
v_lower <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_WA_lower.mat")[[1]]
v_dist <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_var_dist_WA.mat")[[1]]

# Normalize truth
patterns_denom      = apply(patterns, 1, sum)
patterns_scaled     = patterns/patterns_denom
patterns_denom_diag = diag(patterns_denom)
scores_scaled       = as.matrix(scores) %*% patterns_denom_diag;

# Read bootstrapped BN2MF EWA output
bs_dist <- array(dim = c(1000, 4, 500))

for (i in 1:500) {
  bs_dist[,,i] <- readMat(paste0("/Users/lizzy/BN2MF/Bootstrap/bootstrap_cor/bs_ewa_", i, ".mat"))[[1]]
  }
dim(bs_dist)

unique(bs_dist[1,1,])

# Create empirical confidence interval
bs_EWA   <- apply(bs_dist, c(1,2), mean, na.rm = TRUE)
bs_lower <- apply(bs_dist, c(1,2), quantile, 0.025, na.rm = TRUE)
bs_upper <- apply(bs_dist, c(1,2), quantile, 0.975, na.rm = TRUE)

## Example viz
truth = scores_scaled[100,4]

v_ewa = v_EWA[100,4]
v_dist_1 = v_dist[100,4,]
v_q25 = v_lower[100,4]
v_q75 = v_upper[100,4]

bs_ewa = bs_EWA[100,4]
bs_dist_1 = bs_dist[100,4,]
bs_q25 = bs_lower[100,4]
bs_q75 = bs_upper[100,4]

dist_1_1_1 = tibble(Distribution = v_dist_1) %>% 
  mutate(Type = "Variational") %>% 
  rbind(., tibble(Distribution = bs_dist_1) %>% 
          mutate(Type = "Bootstrap")) 

dist_1_1_1 %>% 
  ggplot(aes(x = Distribution)) +
  geom_rect(aes(xmin=bs_q25, xmax=bs_q75, ymin=0, ymax=Inf), fill="pink",      alpha=0.05) +
  geom_rect(aes(xmin=v_q25,  xmax=v_q75,  ymin=0, ymax=Inf), fill="lightblue", alpha=0.05) +
  #geom_histogram(aes(y=..density.., group = Type, fill = Type), alpha = 0.5, bins = 100) +
  geom_density(aes(group = Type, fill = Type), alpha = 0.5) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_bw() + 
  geom_vline(xintercept = truth,            linetype="dotted", color = "black") +
  geom_vline(xintercept = c(v_ewa, bs_ewa), linetype="dotted", color = "yellow") + 
  xlim(0, 500) + ylab("Density")

hist(dist_1_1_1[dist_1_1_1$Type == "Variational",]$Distribution)
hist(dist_1_1_1[dist_1_1_1$Type == "Bootstrap",]$Distribution)

sum(scores_scaled <= v_upper & scores_scaled >= v_lower)/4000
sum(scores_scaled <= bs_upper & scores_scaled >= bs_lower)/4000
  
#####   
##### 100 Correlated Sims
#####

# Read data
load("./Bootstrap/bs_list_lower.RDA")
load("./Bootstrap/bs_list_upper.RDA")
load("./Bootstrap/bs_list_mean.RDA")
load("./Sims/sim_dgp_rep1.RDA")

bs_dat = tibble(seed = 1:100,
              lower = bs_list_lower,
              mean = bs_list_mean,
              upper = bs_list_upper) %>% 
         left_join(., sim_dgp_rep1 %>% filter(data == "Correlated"), by = "seed")

# Normalize truth
bs_dat = bs_dat %>% 
  mutate(denom = map(true_patterns, function(x) apply(x,1, sum)),
         patterns_scaled = map2(true_patterns, denom, function(x,y) x/y),
         scores_scaled = map2(true_scores, denom, function(x,y) as.matrix(x) %*% diag(y)))
  
# Create proportion
bs_dat = bs_dat %>% 
  mutate(prop = pmap(list(scores_scaled, lower, upper), 
                     function(x,y,z) sum(x >= y & x <= z)/4000)) %>% 
  unnest(prop)

summary(bs_dat$prop)

# Check width of VCI compared with width of BS CI
bs_dat

head(bs_dat$lower[[1]])
head(bs_dat$mean[[1]])
head(bs_dat$upper[[1]])

# Are overlapping patterns orthogonal?
m = bs_dat$true_patterns[[1]]
m1 = m[1,]
m2 = m[2,]
m3 = m[3,]
m4 = m[4,]

# No, not all.
m %*% t(m)
# off diagonal are dot products



























  
  