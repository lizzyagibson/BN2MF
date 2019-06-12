# Simulation
# 6/10/2019
options(scipen = 999)
library(MNdata)
library(reshape2)
library(compositions)
library(MASS)
library(tidyverse)

# Simulate dataset with 1000 individuals and 35 chemicals to approx Mothers and Newborns cohort data.

#####################
## M&N Data #########
#####################

# 8 phenols, 9 phthalates, 10 pbdes, 8 pcbs

mn_data <- mn_edc %>% dplyr::select(1:52) %>%
  select_if(~n_distinct(.) > 1) %>% #drop column if all values are the same
  dplyr::select(1:18, ppb_lw_pcb105, ppb_lw_pcb74, ppb_lw_pcb99, ppb_lw_pcb118,
         ppb_lw_pcb138_158, ppb_lw_pcb153, ppb_lw_pcb187, ppb_lw_pcb180,
         grep("BDE", colnames(.))) %>%
  dplyr::select(-sid)

cormat <- round(cor(mn_data, use = "pairwise.complete.obs", method = c("spearman")),2)

melted_cormat <- melt(cormat) %>% rename(Correlation = value)

labelled <- melted_cormat %>% mutate(Class1 = ifelse(str_detect(Var1, "BDE"), "PBDEs",
                                                     ifelse(str_detect(Var1, "^M"), "Phthalates",
                                                            ifelse(str_detect(Var1, "^ppb"), "PCBs",
                                                                   "Phenols")))) %>%
  mutate(Class2 = ifelse(str_detect(Var2, "BDE"), "PBDEs",
                         ifelse(str_detect(Var2, "^M"), "Phthalates",
                                ifelse(str_detect(Var2, "^ppb"), "PCBs",
                                       "Phenols")))) %>%
  mutate(Class2 = factor(Class2, levels = c("Phthalates", "Phenols",
                                            "PCBs", "PBDEs")))

ggplot(data = labelled, aes(x = Var1, y = Var2)) +
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

##############################
## Correlation for Simulation
##############################

# Set up the correlation structure from M&N cohort

# Replace missing with zero for Cholesky transformation
mn_data_drop <- mn_data %>% 
  mutate_if(is.numeric, replace_na, replace = 0) # DON'T REALLY WANT TO DO THIS

# Keep missing as NA and replace zeros with 0.01 for log transformation
mn_data_log <- mn_data %>% 
  mutate(bp_3_sg_adj = ifelse(bp_3_sg_adj == 0, NA, bp_3_sg_adj),
                        tcs_sg_adj = ifelse(tcs_sg_adj == 0, NA, tcs_sg_adj),
                        b_pb_sg_adj = ifelse(b_pb_sg_adj == 0, NA, b_pb_sg_adj)) #drop zeros bc we can't log them

# Dropping these variables will also work  
# select(-BDE209, -BDE183, -BDE28, -BDE66, -BDE85, -BDE100, -BDE153, -BDE154, # drop variables with too many missing
#        -ppb_lw_pcb74, -ppb_lw_pcb105, -ppb_lw_pcb138_158, -ppb_lw_pcb153, -ppb_lw_pcb187)

##############################
## Multivariate log normal ###
##############################

mean_log <- as_vector(map_df(log(mn_data_log), ~mean(., na.rm = TRUE)))

var_log <- cor(log(mn_data_drop), use = "pairwise.complete.obs", method = c("spearman"))
#chol(var_log)

set.seed(1988)
simulated <- rlnorm.rplus(1000, meanlog = mean_log, varlog = var_log)
names(simulated) <- names(mn_data_log)

set.seed(1988)
simulated2 <- exp(mvrnorm(n = 500, mu = mean_log, Sigma = var_log))
names(simulated2) <- names(mn_data_log)

#####################
## Heat map #########
#####################

sim_cormat <- round(cor(simulated, use = "pairwise.complete.obs", method = c("spearman")),2)

sim_melted_cormat <- melt(sim_cormat) %>% rename(Correlation = value)

sim_labelled <- sim_melted_cormat %>% mutate(Class1 = ifelse(str_detect(Var1, "BDE"), "PBDEs",
                                                     ifelse(str_detect(Var1, "^M"), "Phthalates",
                                                            ifelse(str_detect(Var1, "^ppb"), "PCBs",
                                                                   "Phenols")))) %>%
  mutate(Class2 = ifelse(str_detect(Var2, "BDE"), "PBDEs",
                         ifelse(str_detect(Var2, "^M"), "Phthalates",
                                ifelse(str_detect(Var2, "^ppb"), "PCBs",
                                       "Phenols")))) %>%
  mutate(Class2 = factor(Class2, levels = c("Phthalates", "Phenols",
                                            "PCBs", "PBDEs")))

ggplot(data = sim_labelled, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = Correlation), colour = "white") +
  scale_fill_gradient2(low = "#00BFC4", mid = "white", high = "#F8766D",
                       midpoint = 0,
                       na.value = "transparent", limits = c(-1, 1)) +
  theme_grey(base_size = 15) + labs(x = "", y = "", title = "Simulated EDC correlations") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside", legend.position = "bottom",
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18)) +
  facet_grid(Class2 ~ Class1, scales = "free", switch = "both")
