library(bindata) # multivariate binomial dist
library(NMF) # rmatrix
library(MASS) # mvrnorm
library(tidyverse)
library(R.matlab)

#####
# Distinct Patterns
### Create patterns

create_patterns_dist <- function (seed) {
  set.seed(seed)
  patterns_dist <- rbind(c(runif(22, 1, 2), rep(0, 11), rep(0, 17)),
                         c(rep(0, 22), runif(11, 1, 2), rep(0, 17)),
                         c(rep(0, 33), runif(9, 1, 2), rep(0, 8)),
                         c(rep(0, 33), rep(0, 9), runif(8, 1, 2)))
  patterns_dist
}

patterns_dist_iter <- tibble(seed = 1:100) %>% 
  mutate(patterns = map(seed, ~create_patterns_dist(seed = .x)))

### Simulate Scores

create_scores_dist <- function (seed) {
  scores_dist <- matrix(nrow = 1000, ncol = 4)
  set.seed(seed)
  scores_dist <- exp(mvrnorm(n = 1000, mu = rep(1, 4), 
                             Sigma = diag(rep(1, 4)))) # Independent log normal dist
  scores_dist
}
scores_dist_iter <- tibble(seed = 1:100) %>% 
  mutate(scores = map(seed, ~create_scores_dist(seed = .x)))

### Simulate Chemical Exposures

sim_dist0 <- full_join(scores_dist_iter, patterns_dist_iter, by = "seed")
sim_dist0 <- sim_dist0 %>% 
  mutate(chem = map2(scores, patterns, `%*%`))

### Simulate Noise

add_noise <- function (seed, chem) {
  noise <- list(mean=0, sd=1)
  set.seed(seed)
  chem_dist <- pmax(chem + rmatrix(1000, 50, dist=rnorm, mean=noise$mean, sd=noise$sd), 0)	
  chem_dist
  sim_noise <- chem_dist - chem
  list(chem_dist = chem_dist, noise = sim_noise)
}

sim_dist <- sim_dist0 %>% 
  mutate(sim = map2(seed, chem, ~add_noise(seed = .x, chem = .y)[[1]]),
         noise = map2(seed, chem, function(x, y) {add_noise(x, y)[[2]]}))

#####
# 2: Correlated Patterns

### Simulate Patterns
create_patterns_cor <- function (seed) {
  set.seed(seed)
  patterns_cor <- rbind(c(runif(22, 1, 2),   runif(11, 0.5, 1), rep(0, 17)),
                        c(rep(0, 22),        runif(11, 1, 2),   runif(9, 0.5, 1),  rep(0, 8)),
                        c(rep(0, 33),                           runif(9, 1, 2),    runif(8, 0.5, 1)),
                        c(runif(22, 0.5, 1), rep(0, 11),        rep(0, 9),         runif(8, 1, 2)))
  patterns_cor
}
patterns_cor_iter <- tibble(seed = 1:100) %>% 
  mutate(patterns = map(seed, ~create_patterns_cor(seed = .x)))

### Simulate Scores
create_scores_cor <- function (seed) {
  scores_cor <- matrix(nrow = 1000, ncol = 4)
  cov <- matrix( c(1,    0.25, 0.15, 0.15,
                   0.25, 1,    0.15, 0.15,
                   0.15, 0.15, 1,    0.25,
                   0.15, 0.15, 0.25, 1), nrow = 4)
  
  set.seed(seed)
  scores_cor <- exp(mvrnorm(n = 1000, mu = rep(1, 4), 
                            Sigma = cov)) # dists are kind of DEPENDENT
}
scores_cor_iter <- tibble(seed = 1:100) %>% 
  mutate(scores = map(seed, ~create_scores_cor(seed = .x)))

### Simulate Chemical Exposures

sim_cor0 <- full_join(scores_cor_iter, patterns_cor_iter, by = "seed")
sim_cor0 <- sim_cor0 %>% 
  mutate(chem = map2(scores, patterns, `%*%`))

### Simulate Noise
sim_cor <- sim_cor0 %>% 
  mutate(sim = map2(seed, chem, ~add_noise(seed = .x, chem = .y)[[1]]),
         noise = map2(seed, chem, function(x, y) {add_noise(x, y)[[2]]})) 

#####
# 3: Overlapping Patterns
### Simulate Chemical Exposures
sim_over0 <- full_join(scores_dist_iter, patterns_cor_iter, by = "seed")
sim_over0 <- sim_over0 %>% 
  mutate(chem = map2(scores, patterns, `%*%`))

### Simulate Noise

sim_over <- sim_over0 %>% 
  mutate(sim = map2(seed, chem, ~add_noise(seed = .x, chem = .y)[[1]]),
         noise = map2(seed, chem, function(x, y) {add_noise(x, y)[[2]]}))
### Save

isee_all <- rbind(sim_dist %>% mutate(data = "distinct"),
      sim_over %>% mutate(data = "overlapping"),
      sim_cor %>% mutate(data = "correlated"))

# save(isee_all, file = "./Sims/sim_over_old_1012.RDA")
load("./Sims/sim_over_old_1012.RDA")

#####
# Run Models

#####
# PCA
#####

## Function

get_pca <- function (sim) {
  # Run PCA centered, not scaled
  pca_out <- prcomp(sim)
  rot <- pca_out$rotation
  ex <- pca_out$x
  sv <- pca_out$sdev
  
  # Explain >=80% of var
  pve <- sv^2/sum(sv^2)
  rank <- 0
  for (i in 1:length(sv)) {
    if (sum(pve[1:i]) >= 0.8) {
      rank <- i
      break
    }}
  
  # Cut scores and patterns to rank
  rotations <- as_tibble(rot[, 1:rank])
  scores <- if (rank == 1) {matrix(ex[, 1:rank], nrow = nrow(sim))} else {ex[, 1:rank]}
  
  # Predicted values
  pred <- scores %*% t(rotations) + matrix(rep(apply(sim, 2, mean), each= nrow(scores)), nrow = nrow(scores))
  
  list(rotations = rotations, scores = scores, pred = pred, rank = rank)
}

## Run

isee_all <- isee_all %>% 
  mutate(pca_out       = map(sim, get_pca),
         pca_rotations = map(pca_out, function(x) x[[1]]),
         pca_scores    = map(pca_out, function(x) x[[2]]),
         pca_pred      = map(pca_out, function(x) x[[3]]),
         pca_rank      = map(pca_out, function(x) x[[4]]))

#####
# Factor Analysis
#####

## Function

get_fa <- function (sim) {
  # Run FA on 3, 4, and 5 factor models
  set.seed(1988)
  fa_3 <- fa(sim, 3, scores = "regression", rotate = "promax")
  fa_4 <- fa(sim, 4, scores = "regression", rotate = "promax")
  fa_5 <- fa(sim, 5, scores = "regression", rotate = "promax")
  
  # Choose the model with the lowest BIC
  if (max(fa_5$BIC, fa_4$BIC, fa_3$BIC) == fa_5$BIC) {
    fa_out <- fa_5
    rank <- 5
  } else if (max(fa_5$BIC, fa_4$BIC, fa_3$BIC) == fa_4$BIC) {
    fa_out <- fa_4
    rank <- 4
  } else {
    fa_out <- fa_3
    rank <- 3
  }
  
  loadings <- matrix(fa_out$loadings, ncol = ncol(fa_out$scores))
  fa_scores <- fa_out$scores
  pred <- fa_scores %*% t(loadings)
  list(loadings = loadings, fa_scores = fa_scores, pred = pred, rank = rank)
}

## Run

# isee_all[136,] %>% 
#   mutate(fa_out       = map(sim, get_fa))

# 114, 136, 172, 196, 197: x system is computationally singular: reciprocal condition number

isee_all <- isee_all %>% # [c(1:113, 115:135, 137:171, 173:195, 198:200),] %>% 
  mutate(fa_out       = map(sim, get_fa),
         fa_rotations = map(fa_out, function(x) x[[1]]),
         fa_scores    = map(fa_out, function(x) x[[2]]),
         fa_pred      = map(fa_out, function(x) x[[3]]),
         fa_rank      = map(fa_out, function(x) x[[4]]))

#####
# L2 NMF
#####

## Function

get_nmf_l2 <- function (sim) {
  # Run NMF with l2 method 100 times each for 3, 4, and 5 factor models
  set.seed(1988)
  nmf_3 <- nmf(sim, 3, nrun = 100, method = "lee")
  nmf_4 <- nmf(sim, 4, nrun = 100, method = "lee")
  nmf_5 <- nmf(sim, 5, nrun = 100, method = "lee")
  
  # Calculate BIC for each
  bic_3 <- sum((sim - (basis(nmf_3)%*%coef(nmf_3)))^2) + 
    (1/2)*(nrow(sim) + ncol(sim)) * 3 * log(nrow(sim) * ncol(sim))
  bic_4 <- sum((sim - (basis(nmf_4)%*%coef(nmf_4)))^2) + 
    (1/2)*(nrow(sim) + ncol(sim)) * 4 * log(nrow(sim) * ncol(sim))
  bic_5 <- sum((sim - (basis(nmf_5)%*%coef(nmf_5)))^2) +  
    (1/2)*(nrow(sim) + ncol(sim)) * 5 * log(nrow(sim) * ncol(sim))
  
  # Choose model with lowest BIC
  if (max(bic_3, bic_4, bic_5) == bic_3) {
    nmf_out <- nmf_3
    rank <- 3
  } else if (max(bic_3, bic_4, bic_5) == bic_4) {
    nmf_out <- nmf_4
    rank <- 4
  } else {
    nmf_out <- nmf_5
    rank <- 5
  }
  
  basis <- basis(nmf_out)
  coef <- coef(nmf_out)
  pred <- basis %*% coef
  
  list(coef = coef, basis = basis, pred = pred, rank = rank)
}

## Run

isee_all <- isee_all %>%
  mutate(nmf_l2_out      = map(sim, get_nmf_l2),
         nmf_l2_loadings = map(nmf_l2_out, function(x) x[[1]]),
         nmf_l2_scores   = map(nmf_l2_out, function(x) x[[2]]),
         nmf_l2_pred     = map(nmf_l2_out, function(x) x[[3]]),
         nmf_l2_rank     = map(nmf_l2_out, function(x) x[[4]]))

#####
# Poisson NMF
#####

get_nmf_p <- function (sim) {
  set.seed(1988)
  nmf_3 <- nmf(sim, 3, nrun = 100, method = "brunet")
  nmf_4 <- nmf(sim, 4, nrun = 100, method = "brunet")
  nmf_5 <- nmf(sim, 5, nrun = 100, method = "brunet")
  
  bic_3 <- -sum((sim * log(basis(nmf_3) %*% coef(nmf_3))) - (basis(nmf_3) %*% coef(nmf_3))) + 
    (1/2)*(nrow(sim) + ncol(sim)) * 3 * log(nrow(sim) * ncol(sim))
  bic_4 <- -sum((sim * log(basis(nmf_4) %*% coef(nmf_4))) - (basis(nmf_4) %*% coef(nmf_4))) + 
    (1/2)*(nrow(sim) + ncol(sim)) * 4 * log(nrow(sim) * ncol(sim))
  bic_5 <- -sum((sim * log(basis(nmf_5) %*% coef(nmf_5))) - (basis(nmf_5) %*% coef(nmf_5))) + 
    (1/2)*(nrow(sim) + ncol(sim)) * 5 * log(nrow(sim) * ncol(sim))
  
  if (max(bic_3, bic_4, bic_5) == bic_3) {
    nmf_out <- nmf_3
    rank <- 3
  } else if (max(bic_3, bic_4, bic_5) == bic_4) {
    nmf_out <- nmf_4
    rank <- 4
  } else {
    nmf_out <- nmf_5
    rank <- 5
  }
  
  basis <- basis(nmf_out)
  coef <- coef(nmf_out)
  pred <- basis %*% coef
  
  list(coef = coef, basis = basis, pred = pred, rank = rank)
}

isee_all <- isee_all %>% 
  mutate(nmf_p_out      = map(sim, get_nmf_p),
         nmf_p_loadings = map(nmf_p_out, function(x) x[[1]]),
         nmf_p_scores   = map(nmf_p_out, function(x) x[[2]]),
         nmf_p_pred     = map(nmf_p_out, function(x) x[[3]]),
         nmf_p_rank     = map(nmf_p_out, function(x) x[[4]]))

save(isee_all, file = "./HPC/Rout/isee_all.RDA")
