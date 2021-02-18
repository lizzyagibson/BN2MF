# Combine uncertainty from VCI + bootstraps
library(tidyverse)
library(R.matlab)
library(patchwork)

# 1 example data set
# 150 bootstrap samples 
bootstraps = 150

# Sim data
simdata1 = read_csv("/Users/lizzy/BN2MF/Sims/Main/sim_dgp_201.csv");
patterns = read_csv("/Users/lizzy/BN2MF/Sims/Main/patterns_dgp_201.csv");
scores   = read_csv("/Users/lizzy/BN2MF/Sims/Main/scores_dgp_201.csv");

patterns_denom = apply(patterns, 1, sum)
patterns_scaled = patterns/patterns_denom
apply(patterns_scaled, 1, sum)
scores_scaled = as.matrix(scores) %*% diag(patterns_denom)

# MATLAB output
sam_dist  = readMat("./Bootstrap/vci_on_bs/sam_dist.mat")[[1]]
WA_dist_1_3 = readMat("./Bootstrap/vci_on_bs/WA_dist_1_3.mat")[[1]]

lower = readMat("./Bootstrap/vci_on_bs/lower_wa.mat")[[1]]
upper = readMat("./Bootstrap/vci_on_bs/upper_wa.mat")[[1]]
middle = readMat("./Bootstrap/vci_on_bs/ewa_scaled.mat")[[1]]

#### Viz Quantiles ####
# EWA
# first match permutations !!!
# reorder IDs
lower_bs_vci_ewa = array(dim = c(1000, 5, 150))
upper_bs_vci_ewa = array(dim = c(1000, 5, 150))
middle_bs_vci_ewa = array(dim = c(1000, 5, 150))

for (boot in 1:bootstraps) {
  # takes the first instance of an ID
  # all instances of same ID are equal
  # match returns a vector of the positions of (first) matches of its first argument in its second.
  lower_bs_vci_ewa[,,boot] = lower[,,boot][match(1:1000, sam_dist[[boot]][[1]][,1]),]
  upper_bs_vci_ewa[,,boot] = upper[,,boot][match(1:1000, sam_dist[[boot]][[1]][,1]),]
  middle_bs_vci_ewa[,,boot] = middle[,,boot][match(1:1000, sam_dist[[boot]][[1]][,1]),]                        
  print(boot)
}

check = matrix(c(rep(3, 2), rep(1, 2), rep(2, 2), rep(2,2)), nrow = 4, byrow = T)
check
match(1:4, check[,1])
check[match(1:4, check[,1]),]

# Plot
q = tibble()

for (i in 1:bootstraps) {
  q025 = lower_bs_vci_ewa[3,3,i]
  q975 = upper_bs_vci_ewa[3,3,i]
  mean = middle_bs_vci_ewa[3,3,i]
  
  tib = tibble(q025, mean, q975) %>% mutate(bootstrap = as.factor(i))
  q = bind_rows(q, tib)
}

q %>% mutate(iqr = q975-q025) %>% pull(iqr) %>% summary()

truth = scores_scaled[3,2]
truth

#### Viz Dists ####
# EWA
# first match permutations !!!
# reorder IDs
list_bs_vci_ewa = list()
 
for (boot in 1:3) {
    # takes the first instance of an ID
    # all instances of same ID are equal
    # match returns a vector of the positions of (first) matches of its first argument in its second.
  list_bs_vci_ewa[[boot]] = WA_dist_1_3[[boot]][[1]][match(1:1000, sam_dist[[boot]][[1]][,1]),,]
  print(boot)
    }

# Plot
dist_3 = tibble()

for (i in 1:3) {
  vci_dist = list_bs_vci_ewa[[i]][3,2,]
  tib = tibble(vci_dist) %>% mutate(bootstrap = as.factor(i))
  dist_3 = bind_rows(dist_3, tib)
  }

sum_dist = dist_3 %>%
  group_by(bootstrap) %>%
  summarise(mean = mean(vci_dist, na.rm = T),
            q025 = quantile(vci_dist, 0.025, na.rm = T),
            q975 = quantile(vci_dist, 0.975, na.rm = T))

#### Original Bootstrap and VCI ####
load("./Bootstrap/bs_vci_metrics.RDA")
metrics = bs_vci_metrics %>% 
  filter(seed == 1 & data == "Correlated" & side == "wa") %>% 
  dplyr::select(seed, data, scores_scaled, method, lower, upper, prop, iqr)
metrics

bs_lower = metrics$lower[[1]][3,2]
vci_lower = metrics$lower[[2]][3,2]

bs_upper = metrics$upper[[1]][3,2]
vci_upper = metrics$upper[[2]][3,2]

head(metrics$lower[[2]])
head(scores_scaled)
head(metrics$upper[[2]])

sum(scores_scaled <= metrics$upper[[2]] & scores_scaled >= metrics$lower[[2]])/4000

#### Viz All ####

quant_plot = q %>% 
  pivot_longer(q025:q975) %>% 
  ggplot(aes(x = value, group = name)) +
  geom_rect(xmin = vci_lower, xmax = vci_upper, ymin = -Inf, ymax = Inf, fill = "yellow", alpha = 0.05) +
  geom_rect(xmin = bs_lower, xmax = bs_upper, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.5) +
  geom_histogram(aes(y = ..density.., fill = name), col="black",
                 alpha = 0.75, position = "identity") +
  theme_bw() +
  geom_vline(xintercept = truth, linetype="dashed") + xlim(c(5,50))
quant_plot

ex_plot = ggplot(dist_3, aes(x = vci_dist, group = bootstrap)) +
  geom_histogram(aes(y = ..density.., fill = bootstrap, col = bootstrap),
                 alpha = 0.25, position = "identity") +
  geom_vline(data = sum_dist, aes(xintercept = mean, col=bootstrap), linetype="dashed")+
  geom_vline(data = sum_dist, aes(xintercept = q025, col=bootstrap), linetype="dotted")+
  geom_vline(data = sum_dist, aes(xintercept = q975, col=bootstrap), linetype="dotted")+
  theme_bw()+
  geom_vline(xintercept = truth, linetype="dashed") + xlim(c(5,50))
ex_plot

quant_plot / ex_plot

#### Get whole Dist ####
# MATLAB output
sam_dist  = readMat("./Bootstrap/vci_on_bs/sam_dist.mat")[[1]]
WA_dist_1 = readMat("./Bootstrap/vci_on_bs/WA_dist_1.mat")[[1]]
WA_dist_2 = readMat("./Bootstrap/vci_on_bs/WA_dist_2.mat")[[1]]
WA_dist_3 = readMat("./Bootstrap/vci_on_bs/WA_dist_3.mat")[[1]]
WA_dist_4 = readMat("./Bootstrap/vci_on_bs/WA_dist_4.mat")[[1]]
WA_dist_5 = readMat("./Bootstrap/vci_on_bs/WA_dist_5.mat")[[1]]
WA_dist_6 = readMat("./Bootstrap/vci_on_bs/WA_dist_6.mat")[[1]]

WA_dists = c(WA_dist_1, WA_dist_2, WA_dist_3, WA_dist_4, WA_dist_5, WA_dist_6)
WA_dist_1 = list()
WA_dist_2 = list()
WA_dist_3 = list()
WA_dist_4 = list()
WA_dist_5 = list()
WA_dist_6 = list()

# first match permutations !!!
# reorder IDs
reorder_WA_dists = list()
bootstraps = 150

for (boot in 1:bootstraps) {
  # takes the first instance of an ID
  # all instances of same ID are equal
  # match returns a vector of the positions of (first) matches of its first argument in its second.
  reorder_WA_dists[[boot]] = WA_dists[[boot]][[1]][match(1:1000, sam_dist[[boot]][[1]][,1]),,]
  print(boot)
}
head(reorder_WA_dists[[1]][,,1])
str(reorder_WA_dists)
WA_dists = list()
sam_dist = list()

WA_flat = array(dim = c(1000,4,150000))
dim(WA_flat)

for (i in 1:bootstraps) {
  index = ( (i-1)*1000 +1 ): (i*1000)
  WA_flat[,,index] = reorder_WA_dists[[i]]
  print(i)
}
head(WA_flat[,,1])
reorder_WA_dists = list()
dim(WA_flat)

lower_150 = apply(WA_flat, c(1,2), quantile, 0.025, na.rm=T)
upper_150 = apply(WA_flat, c(1,2), quantile, 0.975, na.rm=T)
mean_150  = apply(WA_flat, c(1,2), mean, na.rm=T)

save(lower_150, file =  "./Bootstrap/vci_on_bs/lower_150.RDA")
save(upper_150, file =  "./Bootstrap/vci_on_bs/upper_150.RDA")
save(mean_150,  file = "./Bootstrap/vci_on_bs/mean_150.RDA")

sum(scores_scaled <= upper_150 & scores_scaled >= lower_150)/4000

mean(upper_150 - lower_150)
mean(metrics$upper[[2]] - metrics$lower[[2]])

mean(lower_150)
mean(upper_150)
mean(metrics$lower[[2]])
mean(metrics$upper[[2]])

head(scores_scaled <= upper_150 & scores_scaled >= lower_150)
head(lower_150)
head(scores_scaled)
head(upper_150)



full = tibble(full_dist = WA_flat[3,2,])

full_plot = full %>% ggplot(aes(x = full_dist)) + 
  geom_histogram(aes(y = ..density..), col= "darkgrey",
                 alpha = 0.75, position = "identity", bins = 100) +
  geom_vline(xintercept = lower_150[3,2], linetype="dotted")+
  geom_vline(xintercept = upper_150[3,2], linetype="dotted")+
  theme_bw()+
  geom_vline(xintercept = truth, linetype="dashed") + xlim(c(5,50)) +
  geom_vline(xintercept = vci_lower, linetype="dotted", col="red") +
  geom_vline(xintercept = vci_upper, linetype="dotted", col="red")
full_plot

quant_plot / ex_plot / full_plot
