# Simulation
# 6/10/2019
options(scipen = 999)
library(tidyverse)
library(MNdata)
library(reshape2)
library(compositions)

# Simulate dataset with 500 individuals and 34 chemicals to approx Mothers and Newborns cohort data.

#####################
## M&N Data #########
#####################

# 8 phenols, 9 phthalates, 9 pbdes, 8 pcbs

mn_data <- mn_edc %>% select(1:52) %>%
  select_if(~n_distinct(.) > 1) %>% #drop column if all values are the same
  select(1:18, ppb_lw_pcb105, ppb_lw_pcb74, ppb_lw_pcb99, ppb_lw_pcb118,
         ppb_lw_pcb138_158, ppb_lw_pcb153, ppb_lw_pcb187, ppb_lw_pcb180,
         grep("BDE", colnames(.))) %>%
  select(-sid)

summary(mn_data)

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
mn_data_drop <- mn_data %>% 
  mutate_if(is.numeric, replace_na, replace = 0) # DON'T REALLY WANT TO DO THIS

mn_data_log <- mn_data %>% 
  mutate(bp_3_sg_adj = ifelse(bp_3_sg_adj == 0, NA, bp_3_sg_adj),
                        tcs_sg_adj = ifelse(tcs_sg_adj == 0, NA, tcs_sg_adj),
                        b_pb_sg_adj = ifelse(b_pb_sg_adj == 0, NA, b_pb_sg_adj)) #drop zeros bc we can't log them
  
# select(-BDE209, -BDE183, -BDE28, -BDE66, -BDE85, -BDE100, -BDE153, -BDE154, # drop variables with too many missing
#        -ppb_lw_pcb74, -ppb_lw_pcb105, -ppb_lw_pcb138_158, -ppb_lw_pcb153, -ppb_lw_pcb187)

summary(mn_data_drop)

R <- cor(mn_data_drop, use = "pairwise.complete.obs", method = c("spearman"))

# Cholesky decomposition takes the original X matrix and uses the Cholesky transformation to create 
# a new, correlated, Y variable.
# The correlation structure can be specifically defined.
# Take the Cholesky decomposition of the correlation matrix.

U = t(chol(R))

# Create uncorrelated simulated variables based on their true distributions
# Log normal distributions of 48 chemicals with mean and sd of real data
numobs <- 500
nvars <- dim(R)[1]
chems <- colnames(mn_data_drop)
uncor_sim <- matrix(ncol = numobs, nrow = nvars)
mn_matrix <- as.matrix(mn_data_log)

for (i in 1:nvars) {
    set.seed(88)
    uncor_sim[i,] <- rlnorm(numobs, meanlog = mean(log(mn_matrix[,i]), na.rm = TRUE), 
                            sdlog = sd(log(mn_matrix[,i]), na.rm = TRUE))
}

# Multiply the Cholesky decomposition of the correlation matrix by the data matrix.
# The resulting matrix is a transformed dataset with the specified correlation.
X <- U %*% uncor_sim
sim_dat <- as.data.frame(t(X))

names(sim_dat) <- chems
cor(sim_dat)
head(sim_dat)

##############################
## Multivariate log normal ###
##############################

mean_log <- as_vector(map_df(log(mn_data_log), ~mean(., na.rm = TRUE)))
var_log <- cor(log(mn_data_drop), use = "pairwise.complete.obs", method = c("spearman"))

simulated <- rlnorm.rplus(500, meanlog = mean_log, varlog = var_log)

mean(mn_data$BDE47, na.rm = TRUE)
bde47 <- mn_data$BDE47 + 1
mean(bde47, na.rm = TRUE)

mean(mn_data$BDE99, na.rm = TRUE)
bde99 <- mn_data$BDE99 + 1
mean(bde99, na.rm = TRUE)

cor(mn_data$BDE47, mn_data$BDE99, use = "pairwise.complete.obs", method = c("spearman"))
cor(bde47, bde99, use = "pairwise.complete.obs", method = c("spearman"))


mn_data_one <- mn_data + 1

mean_log <- as_vector(map_df(log(mn_data_one), ~mean(., na.rm = TRUE)))
var_log <- cor(log(mn_data_one), use = "pairwise.complete.obs", method = c("spearman"))

simulated <- rlnorm.rplus(500, meanlog = mean_log, varlog = var_log)

#####################
## Heat map #########
#####################

sim_cormat <- round(cor(simulated, use = "pairwise.complete.obs", method = c("spearman")),2)

sim_melted_cormat <- melt(sim_cormat) %>% rename(Correlation = value)

ggplot(data = sim_melted_cormat, aes(x = Var1, y = Var2)) +
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
        strip.text.y = element_text(size = 18))
