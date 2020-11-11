#########################################################
# Sims, npBNMF, regular NMF (L2 & Poisson), PCA, and FA #
#########################################################
# Lizzy Gibson    # 6/22/2020 ###########################
# Added BNMF      # 9/17/2020 ###########################
# Cleaned up code # 9/22/2020 ###########################
# Rerun on local  # 10/24/2020
#########################################################

# Packages
library(tidyverse)
library(psych)
library(NMF)
library(R.matlab)

source("./R/compare_functions.R")
source("./R/fig_set.R")

# Read in Sims
load("./Sims/sim_dgp_rep1.RDA")
sim_dgp_rep1

sim <- fa_test$sim[[1]]
chem <- fa_test$chem[[1]]

#####
# Factor Analysis
#####
  set.seed(1988)
  fa_3 <- fa(sim, 3, scores = "regression", rotate = "promax")
  fa_4 <- fa(sim, 4, scores = "regression", rotate = "promax")
  fa_5 <- fa(sim, 5, scores = "regression", rotate = "promax")
  
  fact_3 <- factanal(sim, 3, scores = "regression", rotation = "promax")
  fact_4 <- factanal(sim, 4, scores = "regression", rotation = "promax")
  fact_5 <- factanal(sim, 5, scores = "regression", rotation = "promax")
  
  # Choose the model with the LOWEST BIC
  if (min(fa_5$BIC, fa_4$BIC, fa_3$BIC) == fa_5$BIC) {
    fa_out <- fa_5
    rank <- 5
  } else if (min(fa_5$BIC, fa_4$BIC, fa_3$BIC) == fa_4$BIC) {
    fa_out <- fa_4
    rank <- 4
  } else {
    fa_out <- fa_3
    rank <- 3
  }
  
loadings <- matrix(fa_out$loadings, ncol = ncol(fa_out$scores))
scores <- fa_out$scores
pred <- scores %*% t(loadings)

# predict function gives scores
dim(predict.psych(fa_out, sim))

