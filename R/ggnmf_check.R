## Comparing simulated data and results from GNMF
##
library(tidyverse)
library(R.matlab)

## Read in simulated data

sim1 <- read_csv("./Sims/sim_data_1.csv") %>% as.matrix()
sim2 <- read_csv("./Sims/sim_data_2.csv") %>% as.matrix()
sim3 <- read_csv("./Sims/sim_data_3.csv") %>% as.matrix()

## Read in MATLAB output SCORES
ewa1 <- readMat("./Data/ewa1.mat")[[1]]

ewa2 <- readMat("./Data/ewa2.mat")[[1]]

ewa3 <- readMat("./Data/ewa3.mat")[[1]]

## Read in MATLAB output LOADINGS
eh1 <- readMat("./Data/eh1.mat")[[1]]

eh2 <- readMat("./Data/eh2.mat")[[1]]

eh3 <- readMat("./Data/eh3.mat")[[1]]

## Distance between simulated data and scores*loadings

pred1 <- ewa1 %*% eh1
pred2 <- ewa2 %*% eh2
pred3 <- ewa3 %*% eh3

norm(sim1 - pred1, type = "F")/norm(sim1, type = "F")
norm(sim2 - pred2, type = "F")/norm(sim2, type = "F")
norm(sim3 - pred3, type = "F")/norm(sim3, type = "F")







