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

# Generate 5 non-negative vectors of length 50 “by hand” — vectors should be overlapping a bit, 
# and at this stage don’t have to look too much like pollution. 
# just make 20 (of 50) elements 1 and the rest 0, or something dumb like that

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

# generate random scores indicating the presence of each vector for each observation from a gamma distribution

hist(rgamma(5, shape = 0.5))

scores <- matrix(nrow = 1000, ncol = 5)

set.seed(1988)
for (i in 1:nrow(scores)) {
  scores[i,] <- rgamma(5, shape = 0.5)
  #scores[i,] <- scores[i,]/sum(scores[i,]) # Make scores sum to 1
}

head(scores)

#################################
## Simulate Chemical Exposures ##
#################################

# multiple scores and vectors, and add to get the mixture

chem <- scores %*% t(patterns)

cor(chem)

###########
## Noise ##
###########

# noise is harder - if your scores are big enough just use a normal distribution 
# to make sure things “work” and we can fancy it up from there

chem_n <- matrix(nrow = 1000, ncol = 50)

set.seed(1988)
for (i in 1:ncol(chem_n)) {
  chem_n[,i] <- chem[,i] + rnorm(1000)
}

head(chem_n)
cor(chem_n)

###############################
## Simuate with NMF Function ##
###############################

# internal parameters
mu.W <- 1
sd.W <- 1

# r = rank / patterns
r <- 5

# n = participants
n <- 1000

# p = chemicals
p <- 50

set.seed(1988)
g <- rmultinom(1, p, rep(1, r))			
  
# generate H
H <- matrix(0, r, p)
tmp <- 0
for( i in 1:r ){
    H[i,(tmp+1):(tmp+g[i])] <- 1
    tmp <- tmp+g[i]
  } 	

set.seed(1988)  
b <- rmultinom(1, n, rep(1, r))		

# generate W
W <- matrix(0, n, r)
tmp <- 0
for( i in 1:r ){		
    W[(tmp+1):(tmp+b[i]),i] <- abs(rnorm(b[i], mu.W, sd.W))
    tmp <- tmp + b[i]
  }	

# build the composite matrix
res <- W %*% H

# add some noise
res <- pmax(res + rmatrix(res, dist=rnorm, mean=noise$mean, sd=noise$sd), 0)	
  
# return the factors
pData <- list(Group=factor(unlist(mapply(rep, 1:r, g, SIMPLIFY=FALSE))))
fData <- list(Group=factor(unlist(mapply(rep, 1:r, b, SIMPLIFY=FALSE))))
res <- list(res, W=W, H=H, offset=offset, pData=pData, fData=fData)
  
