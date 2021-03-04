# Viz whole distributions
# VCI and bootstrap
# For 1 example sim
library(tidyverse)
library(patchwork)
source("./functions/compare_functions.R")

# Choose sim ####

# Bootstrapped coverage
load("./bootstrap/output/bs_prop.RDA")
bs_prop

load("./bootstrap/output/nmf_bs_prop.RDA")
nmf_bs_prop

# BN2MF VCI coverage
load("./main/bn2mf/output/m_prop.RDA")
m_prop

# Load sims
load("./sims/sim_ids.RDA")
sim_ids

# Subsample for bootstrap
load("./sims/bs_sample.RDA")
bs_subset = bs_sample %>% dplyr::select(1:3) %>% mutate_all(as.factor)

bsci = bind_cols(bs_subset, bs_prop) %>% rename(bs_prop = prop)
vci = bind_cols(sim_ids, m_prop) %>% dplyr::select(-model) %>% rename(vci_prop = prop)
nmfci = bind_cols(bs_subset, nmf_bs_prop) %>% rename(nmf_prop = prop)

all_prop = left_join(bsci, vci) %>% left_join(., nmfci) %>% arrange(desc(bs_prop))

summary(all_prop$bs_prop - all_prop$nmf_prop)

# get example for simulated data set #9562
# Get full VCI and full Bootstrapped distributions

# Read in distributions ####

# Truth ####
chem = read_csv("./sims/ci_ex/true_chem.csv")
patterns = read_csv("./sims/ci_ex/true_patterns.csv")
scores = read_csv("./sims/ci_ex/true_scores.csv")

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
bs_pred <- array(dim = c(1000, 41, bootstraps)) # 41 includes row ID

# Read in data output by `bootstrap_bn2mf.m`
for (boot in 1:bootstraps) {

    bs_ewa[,,boot]  = readMat(paste0("./bootstrap/output/bsci_out/ewa_bs_sim_", set, "_bs_", boot, ".mat"))[[1]]
    bs_eh[,,boot]    = readMat(paste0("./bootstrap/output/bsci_out/eh_bs_sim_", set, "_bs_", boot, ".mat"))[[1]]
    bs_pred[,,boot] = cbind(bs_ewa[,1,boot], bs_ewa[,2:5,boot]%*%bs_eh[,,boot])
  
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
bs_ewa = bs_ewa[,2:5,]

for (boot in 1:bootstraps) {
  # takes the first instance of an ID
  # all instances of same ID are equal
  bs_pred[,,boot] = bs_pred[,,boot][match(1:1000, bs_pred[,,boot][,1]),]
  print(paste0("Match Permutations for Boot Number: ", boot))
}
bs_pred = bs_pred[,2:41,]

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

# NMF bootstrapped distributions ####
sim = read_csv("./sims/ci_ex/sim_sep_9562.csv")

nmf_solution = nmf(sim, 4, nrun = 100, method = "brunet")
nmf_scores <- basis(nmf_solution)
nmf_loadings <- coef(nmf_solution)
nmf_pred <- nmf_scores %*% nmf_loadings

nmf_denom = apply(nmf_loadings, 1, sum)
nmf_loadings_scaled = nmf_loadings / nmf_denom
nmf_scores_scaled = as.matrix(nmf_scores) %*% diag(nmf_denom)

fc_out = factor_correspondence(t(patterns_scaled), t(nmf_loadings_scaled))
nmf_loadings_final = t(fc_out$rearranged)
nmf_scores_final = nmf_scores_scaled %*% fc_out$permutation_matrix

# One array for all bootstraps 
nmf_bs <- tibble()

# Read in data output by `bootstrap_nmf.m`
for (boot in 1:bootstraps) {
  load(paste0("./bootstrap/output/nmf_out/nmf_bs_sim_", set, "_bs_", boot, ".RDA"))
  nmf_bs = bind_rows(nmf_bs, nmf_boot)
  print(paste0("Boot Number: ", boot))
}

# scores
# first match permutations !!!
# reorder IDs
nmf_bs_re = nmf_bs %>%
  mutate_at(vars(3:5), ~map(.,as.matrix)) %>% 
  mutate(sam = map(scores_final, function(x) x[,1]),
         pred_sam = map2(sam, pred, cbind),
         scores_perm = map(scores_final, function(x) x[match(1:1000, x[,1]),]),
         pred_perm = map(pred_sam, function(x) x[match(1:1000, x[,1]),]))

nmf_bs_score = simplify2array(nmf_bs_re$scores_perm)
dim(nmf_bs_score)
nmf_bs_score = nmf_bs_score[,2:5,]
dim(nmf_bs_score)

nmf_bs_eh = simplify2array(nmf_bs_re$loadings_final)
dim(nmf_bs_eh)

nmf_bs_pred = simplify2array(nmf_bs_re$pred_perm)
dim(nmf_bs_pred)
nmf_bs_pred = nmf_bs_pred[,2:41,]
dim(nmf_bs_score)

nmf_bs_score_median= apply(nmf_bs_score, c(1,2), median, na.rm = T)
nmf_bs_wa_lower  = apply(nmf_bs_score, c(1,2), quantile, 0.025, na.rm = T)
nmf_bs_wa_upper  = apply(nmf_bs_score, c(1,2), quantile, 0.975, na.rm = T)

nmf_bs_eh_median= apply(nmf_bs_eh, c(1,2), median, na.rm = T)
nmf_bs_h_lower  = apply(nmf_bs_eh, c(1,2), quantile, 0.025, na.rm = T)
nmf_bs_h_upper  = apply(nmf_bs_eh, c(1,2), quantile, 0.975, na.rm = T)

nmf_bs_pred_median= apply(nmf_bs_pred, c(1,2), median, na.rm = T)
nmf_bs_pred_lower  = apply(nmf_bs_pred, c(1,2), quantile, 0.025, na.rm = T)
nmf_bs_pred_upper  = apply(nmf_bs_pred, c(1,2), quantile, 0.975, na.rm = T)

# Variational distributions ####

vci_pred = readMat("./main/bn2mf/output/vci_out/q_pred_9562.mat")[[1]]
vci_eh_dist = readMat("./main/bn2mf/output/vci_out/sep_distEH_9562.mat")[[1]]
vci_ewa_dist = readMat("./main/bn2mf/output/vci_out/sep_distWA_9562.mat")[[1]]
vci_eh = readMat("./main/bn2mf/output/vci_out/sep_eh_scaled9562.mat")[[1]]
vci_ewa = readMat("./main/bn2mf/output/vci_out/sep_ewa_scaled9562.mat")[[1]]
vci_eh_lower = readMat("./main/bn2mf/output/vci_out/sep_lowerH_9562.mat")[[1]]
vci_ewa_lower = readMat("./main/bn2mf/output/vci_out/sep_lowerWA_9562.mat")[[1]]
vci_eh_upper = readMat("./main/bn2mf/output/vci_out/sep_upperH_9562.mat")[[1]]
vci_ewa_upper = readMat("./main/bn2mf/output/vci_out/sep_upperWA_9562.mat")[[1]]

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
plot_wa_v_dist   = vci_ewa_dist[4,4,]
plot_wa_bs_dist  = bs_ewa[4,4,]
plot_wa_nmf_dist = nmf_bs_score[4,4,]

plot_bs_ewa    = median(plot_wa_bs_dist, na.rm = T)
plot_bs_wa_25  = quantile(plot_wa_bs_dist, 0.025, na.rm = T)
plot_bs_wa_75  = quantile(plot_wa_bs_dist, 0.975, na.rm = T)

plot_nmf_ewa    = median(plot_wa_nmf_dist, na.rm = T)
plot_nmf_wa_25  = quantile(plot_wa_nmf_dist, 0.025, na.rm = T)
plot_nmf_wa_75  = quantile(plot_wa_nmf_dist, 0.975, na.rm = T)

plot_v_ewa   = vci_ewa[4,4]
plot_v_wa_25 = vci_ewa_lower[4,4]
plot_v_wa_75 = vci_ewa_upper[4,4]

ewa_plot = tibble(Distribution = plot_wa_v_dist) %>% 
  mutate(Type = "BN2MF Variational") %>% 
  rbind(., tibble(Distribution = plot_wa_bs_dist) %>% 
          mutate(Type = "BN2MF Bootstrap"))  %>%
  rbind(., tibble(Distribution = plot_wa_nmf_dist) %>% 
          mutate(Type = "NMF Bootstrap"))  %>% 
  drop_na(.)

# EH ####
plot_true_patterns = patterns_scaled[4,7]
plot_h_v_dist      = vci_eh_dist[4,7,]
plot_h_bs_dist     = bs_eh[4,7,]
plot_h_nmf_dist    = nmf_bs_eh[4,7,]

plot_bs_eh   = median(plot_h_bs_dist, na.rm = T)
plot_bs_h_25 = quantile(plot_h_bs_dist, 0.025, na.rm = T)
plot_bs_h_75 = quantile(plot_h_bs_dist, 0.975, na.rm = T)

plot_nmf_eh   = median(plot_h_nmf_dist, na.rm = T)
plot_nmf_h_25 = quantile(plot_h_nmf_dist, 0.025, na.rm = T)
plot_nmf_h_75 = quantile(plot_h_nmf_dist, 0.975, na.rm = T)

plot_v_eh   = vci_eh[4,7]
plot_v_h_25 = vci_eh_lower[4,7]
plot_v_h_75 = vci_eh_upper[4,7]

eh_plot = tibble(Distribution = plot_h_v_dist) %>% 
  mutate(Type = "BN2MF Variational") %>% 
  rbind(., tibble(Distribution = plot_h_bs_dist) %>% 
          mutate(Type = "BN2MF Bootstrap"))  %>% 
  rbind(., tibble(Distribution = plot_h_nmf_dist) %>% 
          mutate(Type = "NMF Bootstrap"))  %>% 
  drop_na(.)

# Pred ####
plot_chem  = as.numeric(chem[3,1])
plot_pred_v_dist   = vci_pred_dist[3,1,]
plot_pred_bs_dist  = bs_pred[3,1,]
plot_pred_nmf_dist = nmf_bs_pred[3,1,]

plot_bs_pred     = median(plot_pred_bs_dist, na.rm = T)
plot_bs_pred_25  = quantile(plot_pred_bs_dist, 0.025, na.rm = T)
plot_bs_pred_75  = quantile(plot_pred_bs_dist, 0.975, na.rm = T)

plot_nmf_pred     = median(plot_pred_nmf_dist, na.rm = T)
plot_nmf_pred_25  = quantile(plot_pred_nmf_dist, 0.025, na.rm = T)
plot_nmf_pred_75  = quantile(plot_pred_nmf_dist, 0.975, na.rm = T)

plot_v_pred    = vci_pred[3,1]
plot_v_pred_25 = vci_pred_lower[3,1]
plot_v_pred_75 = vci_pred_upper[3,1]

pred_plot = tibble(Distribution = plot_pred_v_dist) %>% 
  mutate(Type = "BN2MF Variational") %>% 
  rbind(., tibble(Distribution = plot_pred_bs_dist) %>% 
          mutate(Type = "BN2MF Bootstrap"))  %>% 
  rbind(., tibble(Distribution = plot_pred_nmf_dist) %>% 
          mutate(Type = "NMF Bootstrap"))  %>% 
  drop_na(.)

# Viz ####

# ggsci_db$"nejm"$"default" <- c(
#   "TallPoppy" = "#BC3C29", "DeepCerulean" = "#0072B5",
#   "Zest" = "#E18727", "Eucalyptus" = "#20854E",
#   "WildBlueYonder" = "#7876B1", "Gothic" = "#6F99AD",
#   "Salomie" = "#FFDC91", "FrenchRose" = "#EE4C97"
# )

green = "#20854E"
red = "#BC3C29FF"
blue = "#0072B5FF"

wa_look = 
  ewa_plot %>% 
  ggplot(aes(x = Distribution)) +
  geom_histogram(aes(y= after_stat(density),fill = Type), 
                 position = "identity", bins = 100, alpha = 0.75) +
  geom_density(aes(group = Type)) +
  scale_fill_manual(values = c(red, blue, green)) +
  theme_bw(base_size = 15) + 
  geom_vline(aes(xintercept = plot_true_score),  linetype="dashed", color = "black") +
    
  geom_vline(aes(xintercept = plot_bs_ewa), linetype="dashed", color = red) + 
  geom_vline(aes(xintercept = plot_v_ewa),  linetype="dashed", color = blue) + 
  geom_vline(aes(xintercept = plot_nmf_ewa),  linetype="dashed", color = green) +
    
  geom_vline(aes(xintercept = plot_bs_wa_25), linetype="dotted", color = red) + 
  geom_vline(aes(xintercept = plot_v_wa_25),  linetype="dotted", color = blue) + 
  geom_vline(aes(xintercept = plot_nmf_wa_25),  linetype="dotted", color = green) + 
  
  geom_vline(aes(xintercept = plot_bs_wa_75), linetype="dotted", color = red) + 
  geom_vline(aes(xintercept = plot_v_wa_75),  linetype="dotted", color = blue) + 
  geom_vline(aes(xintercept = plot_nmf_wa_75),  linetype="dotted", color = green) + 
  
  theme(legend.position = "bottom",
        strip.background = element_rect(fill="white")) +
  ylab("Density") + ggtitle("Distributions on Score E[Wa]")
wa_look    

h_look = eh_plot %>% 
  ggplot(aes(x = Distribution)) +
  geom_histogram(aes(y= after_stat(density),fill = Type), 
                 position = "identity", bins = 100, alpha = 0.75) +
  geom_density(aes(group = Type)) +
  scale_fill_manual(values = c(red, blue,green)) +
  theme_bw(base_size = 15) + 
  geom_vline(aes(xintercept = plot_true_patterns),  linetype="dashed", color = "black") +
  geom_vline(aes(xintercept = plot_bs_eh), linetype="dashed", color = red) + 
  geom_vline(aes(xintercept = plot_v_eh),  linetype="dashed", color = blue) + 
  geom_vline(aes(xintercept = plot_nmf_eh),  linetype="dashed", color = green) + 
   
  geom_vline(aes(xintercept = plot_bs_h_25), linetype="dotted", color = red) + 
  geom_vline(aes(xintercept = plot_v_h_25),  linetype="dotted", color = blue) + 
  geom_vline(aes(xintercept = plot_nmf_h_25),  linetype="dotted", color = green) + 
   
  geom_vline(aes(xintercept = plot_bs_h_75), linetype="dotted", color = red) + 
  geom_vline(aes(xintercept = plot_v_h_75),  linetype="dotted", color = blue) + 
   geom_vline(aes(xintercept = plot_nmf_h_75),  linetype="dotted", color = green) + 
   
  theme(legend.position = "bottom",
        strip.background = element_rect(fill="white")) +
  ylab("Density") + ggtitle("Distributions on Loading E[H]") 
h_look    

pred_look = 
  pred_plot %>% 
  ggplot(aes(x = Distribution)) +
  geom_histogram(aes(y= after_stat(density),fill = Type), 
                 position = "identity", bins = 100, alpha = 0.75) +
  geom_density(aes(group = Type), adjust = 2) +
  scale_fill_manual(values = c(red, blue, green)) +
  theme_bw(base_size = 15) + 
  geom_vline(aes(xintercept = plot_chem),  linetype="dashed", color = "black") +
  geom_vline(aes(xintercept = plot_bs_pred), linetype="dashed", color = red) + 
  geom_vline(aes(xintercept = plot_v_pred),  linetype="dashed", color = blue) + 
  geom_vline(aes(xintercept = plot_nmf_pred),  linetype="dashed", color = green) + 
  
  geom_vline(aes(xintercept = plot_bs_pred_25), linetype="dotted", color = red) + 
  geom_vline(aes(xintercept = plot_v_pred_25),  linetype="dotted", color = blue) + 
  geom_vline(aes(xintercept = plot_nmf_pred_25),  linetype="dotted", color = green) + 
  
  geom_vline(aes(xintercept = plot_bs_pred_75), linetype="dotted", color = red) + 
  geom_vline(aes(xintercept = plot_v_pred_75),  linetype="dotted", color = blue) + 
  geom_vline(aes(xintercept = plot_nmf_pred_75),  linetype="dotted", color = green) + 
  theme(legend.position = "bottom",
        strip.background = element_rect(fill="white")) +
  ylab("Density") + ggtitle("Distributions on Predicted Value E[Wa] * E[H]") +
  labs(fill = "")
pred_look    

#pdf("./figures/vci_bci.pdf", width = 10, height = 8)
(wa_look + theme(legend.position = "none")) / (h_look + theme(legend.position = "none")) / pred_look
#dev.off()

# Summary ####

# Coverage
sum(scores_scaled >= vci_ewa_lower & scores_scaled <= vci_ewa_upper)/(1000*4)
sum(patterns_scaled >= vci_eh_lower & patterns_scaled <= vci_eh_upper)/(4*40)
sum(chem >= vci_pred_lower & chem <= vci_pred_upper)/(1000*40)

sum(scores_scaled >= bs_wa_lower & scores_scaled <= bs_wa_upper)/(1000*4)
sum(patterns_scaled >= bs_h_lower & patterns_scaled <= bs_h_upper)/(4*40)
sum(chem >= bs_pred_lower & chem <= bs_pred_upper)/(1000*40)

sum(scores_scaled >= nmf_bs_wa_lower & scores_scaled <= nmf_bs_wa_upper)/(1000*4)
sum(patterns_scaled >= nmf_bs_h_lower & patterns_scaled <= nmf_bs_h_upper)/(4*40)
sum(chem >= nmf_bs_pred_lower & chem <= nmf_bs_pred_upper)/(1000*40)

# Width

vci_ewa_width  = (vci_ewa_upper - vci_ewa_lower)
vci_eh_width   = (vci_eh_upper - vci_eh_lower)
vci_pred_width = (vci_pred_upper - vci_pred_lower)
mean(vci_pred_width)
mean(vci_ewa_width)
mean(vci_eh_width)

bs_ewa_width  = (bs_wa_upper   - bs_wa_lower)
bs_eh_width   = (bs_h_upper    - bs_h_lower)
bs_pred_width = (bs_pred_upper - bs_pred_lower)
mean(bs_pred_width)
mean(bs_ewa_width)
mean(bs_eh_width)

nmf_ewa_width  = (nmf_bs_wa_upper   - nmf_bs_wa_lower)
nmf_eh_width   = (nmf_bs_h_upper    - nmf_bs_h_lower)
nmf_pred_width = (nmf_bs_pred_upper - nmf_bs_pred_lower)
mean(nmf_pred_width)
mean(nmf_ewa_width)
mean(nmf_eh_width)

ewa_width  = (vci_ewa_width  - bs_ewa_width)
eh_width   = (vci_eh_width   - bs_eh_width )
pred_width = (vci_pred_width - bs_pred_width)

nmf_ewa_width  = (nmf_ewa_width  - bs_ewa_width)
nmf_eh_width   = (nmf_eh_width   - bs_eh_width )
nmf_pred_width = (nmf_pred_width - bs_pred_width)

summary(ewa_width) # Variational usually wider (>75%)
summary(t(eh_width)) # Bootstrap always wider
summary(pred_width) # Bootstrap usually wider (>75%)

summary(nmf_ewa_width) # BN2MF bootstrap usually wider (>75%)
summary(t(nmf_eh_width)) # BN2MF Bootstrap always wider
summary(nmf_pred_width) # BN2MF Bootstrap always wider

norm(scores_scaled - vci_ewa, "F")/norm(scores_scaled, "F")
norm(scores_scaled - nmf_scores_final, "F")/norm(scores_scaled, "F")
norm(scores_scaled - bs_ewa_median, "F")/norm(scores_scaled, "F")
