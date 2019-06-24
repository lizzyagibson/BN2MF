######################################
## Simulate Data-Generating Process ##
## 6/21/2019 #########################
######################################

options(scipen = 999)
library(MNdata)
library(reshape2)
library(compositions)
library(MASS)
library(tidyverse)
library(gridExtra)
library(bindata) # multivariate binomial dist

# Simulate dataset with 1000 individuals and 50 chemicals to approx 
# DATA GENERATING PROCESS of Mothers and Newborns cohort data.

#############
## MN Data ##
#############

# 8 phenols, 9 phthalates, 10 pbdes, 8 pcbs
# phenols and phthalates are specific gravity adjusted
# all values < LOD are replaced by participant-specific LOD/sqrt(2)

mn_data <- mn_edc %>% dplyr::select(1:52) %>%
  select_if(~n_distinct(.) > 1) %>% #drop column if all values are the same
  dplyr::select(1:18, pcb105, pcb74, pcb99, pcb118,
                pcb138_158, pcb153, pcb187, pcb180,
                grep("BDE", colnames(.))) %>%
  dplyr::select(-sid)

#######################
## Simulate Patterns ##
#######################

## Construct a binary correlation matrix for 5 TRUE exposure patterns
m <- matrix(c(1,0.25,0.25,0.25,0.25,
              0.25,1,0.25,0.25,0.25,
              0.25,0.25,1,0.25,0.25,
              0.25,0.25,0.25,1,0.25,
              0.25,0.25,0.25,0.25,1), ncol=5)   

## Simulate 50 chemical exposures, and check that they have the specified correlation structure
set.seed(1988)
patterns <- rmvbin(50, margprob = rep(0.4, 5), bincorr = m) 
sum(patterns[,5])
cor(patterns)

############################################
## Simulate Individual Scores on Patterns ##
############################################

scores <- matrix(nrow = 1000, ncol = 5)

set.seed(1988)
for (i in 1:nrow(scores)) {
  scores[i,] <- rgamma(5, shape = 0.5)
  scores[i,] <- scores[i,]/sum(scores[i,]) # Make scores sum to 1
}

head(scores)

##############################
## Multivariate log normal ###
##############################

##############################
## Multivariate log normal ###
##############################

##############################
## Multivariate log normal ###
##############################

##############################
## Multivariate log normal ###
##############################
