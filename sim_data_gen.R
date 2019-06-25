######################################
## Simulate Data-Generating Process ##
## 6/21/2019 #########################
######################################

options(scipen = 999)
library(MNdata)
library(reshape2)
library(tidyverse)
library(gridExtra)
library(bindata) # multivariate binomial dist
library(NMF) # rmatrix

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
## Look at Heatmap Structure ##
###############################

cormat <- round(cor(chem_n, use = "complete.obs"),2)

melted_cormat <- melt(cormat) %>% rename(Correlation = value)

ggplot(data = melted_cormat, aes(x = Var1, y = Var2)) +
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
        strip.text.y = element_text(size = 18))

###############################
## Simuate with NMF Function ##
###############################

# r = rank / patterns
# r specification of the factorization rank
r <- 5

# n = participants
# n number of rows of the target matrix
n <- 1000

# p = chemicals
# p number of columns of the synthetic target matrix
p <- 50

set.seed(1988)
g <- rmultinom(1, p, rep(1, r))			
# r groups of samples are generated from a draw from a multinomial 
# distribution with equal probabilities that provides their sizes

# generate H
H <- matrix(0, r, p)
tmp <- 0

set.seed(1988)
for( i in 1:r ){
    H[i,(tmp+1):(tmp+g[i])] <- 1
    tmp <- tmp+g[i]
  } 	

set.seed(1988)  
b <- rmultinom(1, n, rep(1, r))		
# r groups of samples are generated from a draw from a multinomial 
# distribution with equal probabilities that provides their sizes

# generate W
W <- matrix(0, n, r)
tmp <- 0
mu.W <- 1
sd.W <- 1

set.seed(1988)
for( i in 1:r ){		
    W[(tmp+1):(tmp+b[i]),i] <- abs(rnorm(b[i], mu.W, sd.W))
    tmp <- tmp + b[i]
  }	

# build the composite matrix
res <- W %*% H

# add some noise
noise <- list(mean=0, sd=1)
set.seed(1988)
res <- pmax(res + rmatrix(res, dist=rnorm, mean=noise$mean, sd=noise$sd), 0)	
  
# return the factors
pData <- list(Group=factor(unlist(mapply(rep, 1:r, g, SIMPLIFY=FALSE))))
fData <- list(Group=factor(unlist(mapply(rep, 1:r, b, SIMPLIFY=FALSE))))

###############################
## Look at Heatmap Structure ##
###############################

cormat2 <- round(cor(res, use = "complete.obs"),2)

melted_cormat2 <- melt(cormat2) %>% rename(Correlation = value)

ggplot(data = melted_cormat2, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = Correlation), colour = "white") +
  scale_fill_gradient2(low = "#00BFC4", mid = "white", high = "#F8766D",
                       midpoint = 0,
                       na.value = "transparent", limits = c(-1, 1)) +
  theme_grey(base_size = 15) + labs(x = "", y = "") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside", legend.position = "bottom",
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18))

