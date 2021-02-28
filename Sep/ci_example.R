# Viz whole distributions
# VCI and bootstrap
# For 1 example sim
library(tidyverse)
library(patchwork)

# Choose sim ####

# Bootstrapped coverage
load("./Sep/bs_prop.RDA")
bs_prop

# BN2MF VCI coverage
load("./Sep/m_prop.RDA")
m_prop

# Load sims
load("./Sims/sim_full.RDA")
sim_sep = sim_sep %>% dplyr::select(1:3) %>% mutate_all(as.factor)

# Subsample for bootstrap
load("./Sep/bs_sample.RDA")
bs_sample = bs_sample %>% dplyr::select(1:3) %>% mutate_all(as.factor)

bsci = bind_cols(bs_sample, bs_prop) %>% rename(bs_prop = prop)
vci = bind_cols(sim_sep, m_prop) %>% dplyr::select(-model) %>% rename(vci_prop = prop)

left_join(bsci, vci) %>% arrange(desc(bs_prop))
# get example for simulated data set #9562
# Get full VCI and full Bootstrapped distributions

# Read in distributions ####

# Truth
chem = read_csv("./Sep/ci_ex/true_chem.csv")
patterns = read_csv("./Sep/ci_ex/true_patterns.csv")
scores = read_csv("./Sep/ci_ex/true_scores.csv")

denom = apply(patterns, 1, sum)
patterns_scaled = patterns / denom
scores_scaled = as.matrix(scores) %*% diag(denom)

# Bootstrapped distributions ####

# 150 bootstrap samples 
bootstraps = 150

set = 9562

# One array for all bootstraps 
bs_ewa <- array(dim = c(1000, 5, bootstraps)) # 5 includes row ID
bs_eh <- array(dim = c(4, 40, bootstraps))
bs_pred <- array(dim = c(1000, 40, bootstraps))

# Read in data output by `bootstrap_bn2mf.m`
for (boot in 1:bootstraps) {
  if (file.exists(paste0("./Sep/bsci_out/ewa_bs_sim_", set, "_bs_", boot, ".mat"))) {
  
    bs_ewa[,,boot]  = readMat(paste0("./Sep/bsci_out/ewa_bs_sim_", set, "_bs_", boot, ".mat"))[[1]]
    bs_eh[,,boot]    = readMat(paste0("./Sep/bsci_out/eh_bs_sim_", set, "_bs_", boot, ".mat"))[[1]]
    bs_pred[,,boot] = readMat(paste0("./Sep/bsci_out/pred_bs_sim_", set, "_bs_", boot, ".mat"))[[1]]
  }
  print(paste0("Boot Number: ", boot))
}

# EWA
# first match permutations !!!
# reorder IDs

for (boot in 1:bootstraps) {
  # takes the first instance of an ID
  # all instances of same ID are equal
  bs_ewa[,,boot] = bs_ewa[,,boot][match(1:1000, bs_ewa[,,boot][,1]),]
  print(paste0("Match Permutations for Boot Number: ", boot))
}
head(bs_ewa[,,1])

bs_ewa = bs_ewa[,2:5,]

dim(bs_ewa)
dim(bs_eh)
dim(bs_pred)

bs_ewa_median= apply(bs_ewa, c(1,2), median, na.rm = T)
bs_wa_lower  = apply(bs_ewa, c(1,2), quantile, 0.025, na.rm = T)
bs_wa_upper  = apply(bs_ewa, c(1,2), quantile, 0.975, na.rm = T)

bs_eh_median= apply(bs_eh, c(1,2), median, na.rm = T)
bs_h_lower  = apply(bs_eh, c(1,2), quantile, 0.025, na.rm = T)
bs_h_upper  = apply(bs_eh, c(1,2), quantile, 0.975, na.rm = T)

bs_pred_median= apply(bs_pred, c(1,2), median, na.rm = T)
bs_pred_lower  = apply(bs_pred, c(1,2), quantile, 0.025, na.rm = T)
bs_pred_upper  = apply(bs_pred, c(1,2), quantile, 0.975, na.rm = T)

# Variational distributions ####

vci_pred = readMat("./Sep/vci_out/q_pred_9562.mat")[[1]]
vci_eh_dist = readMat("./Sep/vci_out/sep_distEH_9562.mat")[[1]]
vci_ewa_dist = readMat("./Sep/vci_out/sep_distWA_9562.mat")[[1]]
vci_eh = readMat("./Sep/vci_out/sep_eh_scaled9562.mat")[[1]]
vci_ewa = readMat("./Sep/vci_out/sep_ewa_scaled9562.mat")[[1]]
vci_eh_lower = readMat("./Sep/vci_out/sep_lowerH_9562.mat")[[1]]
vci_ewa_lower = readMat("./Sep/vci_out/sep_lowerWA_9562.mat")[[1]]
vci_eh_upper = readMat("./Sep/vci_out/sep_upperH_9562.mat")[[1]]
vci_ewa_upper = readMat("./Sep/vci_out/sep_upperWA_9562.mat")[[1]]

vci_pred_dist = array(dim = c(1000, 40, 1000))
# for (i in 1:nrow(vci_pred)) {
#   for (j in 1:ncol(vci_pred)) {
#       vci_pred_dist[i,j,]  = as.vector(rpois(1000, vci_pred[i,j]))
# }}
# vci_pred_lower = qpois(vci_pred, p = 0.025)
# vci_pred_upper = qpois(vci_pred, p = 0.975)

# above is Poisson dist with mean = WaH (this is discrete)
# below is Wa dist * H dist (this is continuous)
for (i in 1:1000) {
  vci_pred_dist[,,i] = vci_ewa_dist[,,i] %*% vci_eh_dist[,,i]
}
str(vci_pred_dist)

vci_pred_lower = apply(vci_pred_dist, c(1,2), quantile, p = 0.025)
vci_pred_upper = apply(vci_pred_dist, c(1,2), quantile, p = 0.975)

# Get distributions ####

# EWA ####
plot_true_score  = scores_scaled[4,4]
plot_wa_v_dist  = vci_ewa_dist[4,4,]
plot_wa_bs_dist = bs_ewa[4,4,]

plot_bs_ewa    = median(plot_wa_bs_dist, na.rm = T)
plot_bs_wa_25  = quantile(plot_wa_bs_dist, 0.025, na.rm = T)
plot_bs_wa_75  = quantile(plot_wa_bs_dist, 0.975, na.rm = T)

plot_v_ewa   = vci_ewa[4,4]
plot_v_wa_25 = vci_ewa_lower[4,4]
plot_v_wa_75 = vci_ewa_upper[4,4]

ewa_plot = tibble(Distribution = plot_wa_v_dist) %>% 
  mutate(Type = "Variational") %>% 
  rbind(., tibble(Distribution = plot_wa_bs_dist) %>% 
          mutate(Type = "Bootstrap"))  %>% 
  drop_na(.)

# EH ####
plot_true_patterns  = patterns_scaled[4,7]
plot_h_v_dist  = vci_eh_dist[4,7,]
plot_h_bs_dist = bs_eh[4,7,]

plot_bs_eh    = median(plot_h_bs_dist, na.rm = T)
plot_bs_h_25  = quantile(plot_h_bs_dist, 0.025, na.rm = T)
plot_bs_h_75  = quantile(plot_h_bs_dist, 0.975, na.rm = T)

plot_v_eh   = vci_eh[4,7]
plot_v_h_25 = vci_eh_lower[4,7]
plot_v_h_75 = vci_eh_upper[4,7]

eh_plot = tibble(Distribution = plot_h_v_dist) %>% 
  mutate(Type = "Variational") %>% 
  rbind(., tibble(Distribution = plot_h_bs_dist) %>% 
          mutate(Type = "Bootstrap"))  %>% 
  drop_na(.)

# Pred ####
plot_chem  = as.numeric(chem[3,1])
plot_pred_v_dist  = vci_pred_dist[3,1,]
plot_pred_bs_dist = bs_pred[3,1,]

plot_bs_pred     = median(plot_pred_bs_dist, na.rm = T)
plot_bs_pred_25  = quantile(plot_pred_bs_dist, 0.025, na.rm = T)
plot_bs_pred_75  = quantile(plot_pred_bs_dist, 0.975, na.rm = T)

plot_v_pred    = vci_pred[3,1]
plot_v_pred_25 = vci_pred_lower[3,1]
plot_v_pred_75 = vci_pred_upper[3,1]

pred_plot = tibble(Distribution = plot_pred_v_dist) %>% 
  mutate(Type = "Variational") %>% 
  rbind(., tibble(Distribution = plot_pred_bs_dist) %>% 
          mutate(Type = "Bootstrap"))  %>% 
  drop_na(.)

# Viz ####
red = "#BC3C29FF"
blue = "#0072B5FF"

wa_look = ewa_plot %>% 
  ggplot(aes(x = Distribution)) +
  geom_histogram(aes(y= after_stat(density),fill = Type), 
                 position = "identity", bins = 100, alpha = 0.75) +
  geom_density(aes(group = Type)) +
  scale_fill_manual(values = c(red, blue)) +
  theme_bw(base_size = 15) + 
  geom_vline(aes(xintercept = plot_true_score),  linetype="dashed", color = "black") +
  geom_vline(aes(xintercept = plot_bs_ewa), linetype="dashed", color = red) + 
  geom_vline(aes(xintercept = plot_v_ewa),  linetype="dashed", color = blue) + 
  geom_vline(aes(xintercept = plot_bs_wa_25), linetype="dotted", color = red) + 
  geom_vline(aes(xintercept = plot_v_wa_25),  linetype="dotted", color = blue) + 
  geom_vline(aes(xintercept = plot_bs_wa_75), linetype="dotted", color = red) + 
  geom_vline(aes(xintercept = plot_v_wa_75),  linetype="dotted", color = blue) + 
  theme(legend.position = "bottom",
        strip.background = element_rect(fill="white")) +
  ylab("Density") + ggtitle("Distributions on Scores E[Wa]")
wa_look    

h_look = eh_plot %>% 
  ggplot(aes(x = Distribution)) +
  geom_histogram(aes(y= after_stat(density),fill = Type), 
                 position = "identity", bins = 100, alpha = 0.75) +
  geom_density(aes(group = Type)) +
  scale_fill_manual(values = c(red, blue)) +
  theme_bw(base_size = 15) + 
  geom_vline(aes(xintercept = plot_true_patterns),  linetype="dashed", color = "black") +
  geom_vline(aes(xintercept = plot_bs_eh), linetype="dashed", color = red) + 
  geom_vline(aes(xintercept = plot_v_eh),  linetype="dashed", color = blue) + 
  geom_vline(aes(xintercept = plot_bs_h_25), linetype="dotted", color = red) + 
  geom_vline(aes(xintercept = plot_v_h_25),  linetype="dotted", color = blue) + 
  geom_vline(aes(xintercept = plot_bs_h_75), linetype="dotted", color = red) + 
  geom_vline(aes(xintercept = plot_v_h_75),  linetype="dotted", color = blue) + 
  theme(legend.position = "bottom",
        strip.background = element_rect(fill="white")) +
  ylab("Density") + ggtitle("Distributions on Patterns E[H]") 
h_look    

pred_look = 
  pred_plot %>% 
  ggplot(aes(x = Distribution)) +
  geom_histogram(aes(y= after_stat(density),fill = Type), 
                 position = "identity", bins = 100, alpha = 0.75) +
  geom_density(aes(group = Type), adjust = 2) +
  scale_fill_manual(values = c(red, blue)) +
  theme_bw(base_size = 15) + 
  geom_vline(aes(xintercept = plot_chem),  linetype="dashed", color = "black") +
  geom_vline(aes(xintercept = plot_bs_pred), linetype="dashed", color = red) + 
  geom_vline(aes(xintercept = plot_v_pred),  linetype="dashed", color = blue) + 
  geom_vline(aes(xintercept = plot_bs_pred_25), linetype="dotted", color = red) + 
  geom_vline(aes(xintercept = plot_v_pred_25),  linetype="dotted", color = blue) + 
  geom_vline(aes(xintercept = plot_bs_pred_75), linetype="dotted", color = red) + 
  geom_vline(aes(xintercept = plot_v_pred_75),  linetype="dotted", color = blue) + 
  theme(legend.position = "bottom",
        strip.background = element_rect(fill="white")) +
  ylab("Density") + ggtitle("Distributions on Predicted Values E[WaH]")
pred_look    

pdf("./Figures/vci_bci.pdf", width = 10, height = 8)
(wa_look + theme(legend.position = "none")) / (h_look + theme(legend.position = "none")) / pred_look
dev.off()

# Summary

# Coverage
sum(scores_scaled >= vci_ewa_lower & scores_scaled <= vci_ewa_upper)/(1000*4)
sum(patterns_scaled >= vci_eh_lower & patterns_scaled <= vci_eh_upper)/(4*40)
sum(chem >= vci_pred_lower & chem <= vci_pred_upper)/(1000*40)

sum(scores_scaled >= bs_wa_lower & scores_scaled <= bs_wa_upper)/(1000*4)
sum(patterns_scaled >= bs_h_lower & patterns_scaled <= bs_h_upper)/(4*40)
sum(chem >= bs_pred_lower & chem <= bs_pred_upper)/(1000*40)

# Width

vci_ewa_width  = (vci_ewa_upper - vci_ewa_lower)
vci_eh_width   = (vci_eh_upper - vci_eh_lower)
vci_pred_width = (vci_pred_upper - vci_pred_lower)

bs_ewa_width  = (bs_wa_upper   - bs_wa_lower)
bs_eh_width   = (bs_h_upper    - bs_h_lower)
bs_pred_width = (bs_pred_upper - bs_pred_lower)

ewa_width  = (vci_ewa_width  - bs_ewa_width)
eh_width   = (vci_eh_width   - bs_eh_width )
pred_width = (vci_pred_width - bs_pred_width)
head(pred_width)

summary(ewa_width) # Variational usually wider (>75%)
summary(t(eh_width)) # Bootstrap always wider
summary(pred_width) # Bootstrap usually wider (>75%)
