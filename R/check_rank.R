source("./R/compare_functions.R")
source("./R/fig_set.R")
options(scipen = 999)
library(Matrix)

# Read in Sims
load("./Sims/sim_dgp_rep1.RDA")

chem <- sim_dgp_rep1$chem[[1]]
sim <- sim_dgp_rep1$sim[[1]]
sdev <- prcomp(sim)$sdev
sdev^2/sum(sdev^2)

rankMatrix(chem)[[1]]

chem <- list()
sim <- list()
chem_rank <- c()
sim_rank <- c()
for (i in 1:200) {
  chem[[i]] <- sim_dgp_rep1$chem[[i]]
  sim[[i]] <- sim_dgp_rep1$sim[[i]]
  
  chem_rank[i] <- rankMatrix(chem[[i]])[[1]]
  sim_rank[i] <- rankMatrix(sim[[i]])[[1]]
  }

