library(R.matlab)
library(tidyverse)

mat_sim <- read_csv("./Results/Main/dgp_rep1/sim_dgp_rep1_1.csv") %>% as.matrix()
load("./Results/Main/dgp_rep1_all.RDA")
load("./Results/Main/dgp_rep1_reordered.RDA")

head(mat_sim)[,1:5]
head(dgp_rep1_all$sim[[1]])[,1:5]

all(dgp_rep1_all$sim[[1]] == mat_sim)
all_equal(as_tibble(dgp_rep1_all$sim[[1]]), as_tibble(mat_sim))
# Sim is the same

truth <- dgp_rep1_all$true_patterns[[1]]
truth_norm <- apply(truth, 2, function(x) x/norm(x, "2"))

EH <- readMat("./MATLAB/test_loadings.mat")[[1]]
eh_norm <- apply(EH, 2, function(x) x/norm(x, "2"))
upper <- readMat("./MATLAB/test_upper.mat")[[1]]
upper_norm <- apply(upper, 2, function(x) x/norm(x, "2"))
lower <- readMat("./MATLAB/test_lower.mat")[[1]]
lower_norm <- apply(lower, 2, function(x) x/norm(x, "2"))

eh_norm[,1:5]
lower_norm[,1:5]
upper_norm[,1:5]
