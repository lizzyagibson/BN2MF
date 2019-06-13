# Simulation
# 6/10/2019
options(scipen = 999)
library(MNdata)
library(reshape2)
library(compositions)
library(MASS)
library(tidyverse)
library(gridExtra)

# Simulate dataset with 1000 individuals and 35 chemicals to approx Mothers and Newborns cohort data.

#####################
## M&N Data #########
#####################

# 8 phenols, 9 phthalates, 10 pbdes, 8 pcbs

mn_data <- mn_edc %>% dplyr::select(1:52) %>%
  select_if(~n_distinct(.) > 1) %>% #drop column if all values are the same
  dplyr::select(1:18, pcb105, pcb74, pcb99, pcb118,
         pcb138_158, pcb153, pcb187, pcb180,
         grep("BDE", colnames(.))) %>%
  dplyr::select(-sid)

cormat <- round(cor(mn_data, use = "complete.obs", method = c("spearman")),2)

melted_cormat <- melt(cormat) %>% rename(Correlation = value)

labelled <- melted_cormat %>% mutate(Class1 = ifelse(str_detect(Var1, "BDE"), "PBDEs",
                                                     ifelse(str_detect(Var1, "^M"), "Phthalates",
                                                            ifelse(str_detect(Var1, "pcb"), "PCBs",
                                                                   "Phenols")))) %>%
  mutate(Class2 = ifelse(str_detect(Var2, "BDE"), "PBDEs",
                         ifelse(str_detect(Var2, "^M"), "Phthalates",
                                ifelse(str_detect(Var2, "pcb"), "PCBs",
                                       "Phenols")))) %>%
  mutate(Class2 = factor(Class2, levels = c("Phthalates", "Phenols",
                                            "PCBs", "PBDEs")))

original <- ggplot(data = labelled, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = Correlation), colour = "white") +
  scale_fill_gradient2(low = "#00BFC4", mid = "white", high = "#F8766D",
                       midpoint = 0,
                       na.value = "transparent", limits = c(-1, 1)) +
  theme_grey(base_size = 15) + labs(x = "", y = "", title = "EDC correlations in Mothers & Newborns") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside", legend.position = "bottom",
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18)) +
  facet_grid(Class2 ~ Class1, scales = "free", switch = "both")

mn_data %>% 
  gather(key = variable) %>% 
  ggplot(aes(x = value)) + geom_density() + facet_wrap(.~variable, scales = "free")

##############################
## Multivariate log normal ###
##############################

# Vector of means
mean_log <- as_vector(map_df(log(mn_data), ~mean(., na.rm = TRUE)))

# Covariance matrix
var_log <- cov(log(mn_data), use = "complete.obs")
# missing values are handled by casewise deletion (and if there are no complete cases, that gives an error)

set.seed(1988)
simulated <- rlnorm.rplus(1000, meanlog = mean_log, varlog = var_log)
colnames(simulated) <- names(mn_data)

set.seed(1988)
simulated2 <- exp(mvrnorm(n = 1000, mu = mean_log, Sigma = var_log))
colnames(simulated2) <- names(mn_data)

simulated[1:5,1:7]
simulated2[1:5,1:7]

#####################
## Heat map #########
#####################

sim_cormat <- round(cor(simulated, use = "pairwise.complete.obs", method = c("spearman")),2)

sim_melted_cormat <- melt(sim_cormat) %>% rename(Correlation = value)

sim_labelled <- sim_melted_cormat %>% mutate(Class1 = ifelse(str_detect(Var1, "BDE"), "PBDEs",
                                                     ifelse(str_detect(Var1, "^M"), "Phthalates",
                                                            ifelse(str_detect(Var1, "pcb"), "PCBs",
                                                                   "Phenols")))) %>%
  mutate(Class2 = ifelse(str_detect(Var2, "BDE"), "PBDEs",
                         ifelse(str_detect(Var2, "^M"), "Phthalates",
                                ifelse(str_detect(Var2, "pcb"), "PCBs",
                                       "Phenols")))) %>%
  mutate(Class2 = factor(Class2, levels = c("Phthalates", "Phenols",
                                            "PCBs", "PBDEs")))

sim1_plot <- ggplot(data = sim_labelled, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = Correlation), colour = "white") +
  scale_fill_gradient2(low = "#00BFC4", mid = "white", high = "#F8766D",
                       midpoint = 0,
                       na.value = "transparent", limits = c(-1, 1)) +
  theme_grey(base_size = 15) + labs(x = "", y = "", title = "Simulated multivariate log normal") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside", legend.position = "bottom",
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18)) +
  facet_grid(Class2 ~ Class1, scales = "free", switch = "both")

sim2_cormat <- round(cor(simulated2, use = "pairwise.complete.obs", method = c("spearman")),2)

sim2_melted_cormat <- melt(sim2_cormat) %>% rename(Correlation = value)

sim2_labelled <- sim2_melted_cormat %>% mutate(Class1 = ifelse(str_detect(Var1, "BDE"), "PBDEs",
                                                             ifelse(str_detect(Var1, "^M"), "Phthalates",
                                                                    ifelse(str_detect(Var1, "pcb"), "PCBs",
                                                                           "Phenols")))) %>%
  mutate(Class2 = ifelse(str_detect(Var2, "BDE"), "PBDEs",
                         ifelse(str_detect(Var2, "^M"), "Phthalates",
                                ifelse(str_detect(Var2, "pcb"), "PCBs",
                                       "Phenols")))) %>%
  mutate(Class2 = factor(Class2, levels = c("Phthalates", "Phenols",
                                            "PCBs", "PBDEs")))

sim2_plot <- ggplot(data = sim_labelled, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = Correlation), colour = "white") +
  scale_fill_gradient2(low = "#00BFC4", mid = "white", high = "#F8766D",
                       midpoint = 0,
                       na.value = "transparent", limits = c(-1, 1)) +
  theme_grey(base_size = 15) + labs(x = "", y = "", title = "Simulated exp(multivariate normal)") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside", legend.position = "bottom",
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18)) +
  facet_grid(Class2 ~ Class1, scales = "free", switch = "both")

grid.arrange(original, sim1_plot, sim2_plot, nrow = 1)
