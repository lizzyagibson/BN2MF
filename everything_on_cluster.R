#########################################################
# Sims, npBNMF, regular NMF (L2 & Poisson), PCA, and FA #
#########################################################
# Lizzy Gibson # 6/22/2020 ##############################
#########################################################

# Packages
library(tidyverse)
library(NMF)
library(psych) # Factor Analysis

# Read in Sims
load(here::here("./R/Sims/Iterate/sim_dist.RDA"))
load(here::here("./R/Sims/Iterate/sim_cor.RDA"))
load(here::here("./R/Sims/Iterate/sim_over.RDA"))

sim_dist <- sim_dist %>% rename(true_patterns = patterns, true_scores = scores) %>% select(-chem, -noise) %>% slice(1:2)
sim_over <- sim_over %>% rename(true_patterns = patterns, true_scores = scores) %>% select(-chem, -noise) %>% slice(1:2)
sim_cor  <- sim_cor %>% rename(true_patterns = patterns, true_scores = scores) %>% select(-chem, -noise) %>% slice(1:2)

# Run everything

# PCA

## Function

get_pca <- function (sim) {
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
  rotations <- rot[, 1:rank]
  scores <- ex[, 1:rank]
  
  # Predicted values
  pred <- scores %*% t(rotations) + matrix(rep(apply(sim, 2, mean), each= nrow(scores)), nrow = nrow(scores))
  
  list(rotations = rotations, scores = scores, pred = pred, rank = rank)
}

## Run

out_dist <- sim_dist %>% 
  mutate(pca_out = map(sim, get_pca),
         pca_rotations = map(pca_out, function(x)
                          x[[1]]),
         pca_scores    = map(pca_out, function(x)
                          x[[2]]),
         pca_pred      = map (pca_out, function(x)
                          x[[3]]),
         pca_rank      = map(pca_out, function(x)
                          x[[4]]))

out_over <- sim_over %>% 
  mutate(pca_out = map(sim, get_pca),
         pca_rotations = map(pca_out, function(x)
           x[[1]]),
         pca_scores    = map(pca_out, function(x)
           x[[2]]),
         pca_pred      = map (pca_out, function(x)
           x[[3]]),
         pca_rank      = map(pca_out, function(x)
           x[[4]]))

out_cor <- sim_cor %>% 
  mutate(pca_out = map(sim, get_pca),
         pca_rotations = map(pca_out, function(x)
           x[[1]]),
         pca_scores    = map(pca_out, function(x)
           x[[2]]),
         pca_pred      = map (pca_out, function(x)
           x[[3]]),
         pca_rank      = map(pca_out, function(x)
           x[[4]]))

# Factor Analysis

## Function

# Factor Analysis

## Function

get_fa <- function (sim) {
  set.seed(1988)
  fa_3 <- fa(sim, 3, scores = "regression", rotate = "promax")
  set.seed(1988)
  fa_4 <- fa(sim, 4, scores = "regression", rotate = "promax")
  set.seed(1988)
  fa_5 <- fa(sim, 5, scores = "regression", rotate = "promax")
  
  if (fa_5$BIC < fa_4$BIC && fa_5$BIC < fa_3$BIC) {
    fa_out <- fa_5
    rank <- 5
  } else if (fa_4$BIC < fa_3$BIC && fa_4$BIC < fa_5$BIC) {
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

out_dist <- out_dist %>% 
  mutate(fa_out = map(sim, get_fa),
         fa_rotations = map(fa_out, function(x)
           x[[1]]),
         fa_scores    = map(fa_out, function(x)
           x[[2]]),
         fa_pred      = map (fa_out, function(x)
           x[[3]]),
         fa_rank      = map(fa_out, function(x)
           x[[4]]))

out_over <- out_over %>% 
  mutate(fa_out = map(sim, get_fa),
         fa_rotations = map(fa_out, function(x)
           x[[1]]),
         fa_scores    = map(fa_out, function(x)
           x[[2]]),
         fa_pred      = map (fa_out, function(x)
           x[[3]]),
         fa_rank      = map(fa_out, function(x)
           x[[4]]))

out_cor <- out_cor %>% 
  mutate(fa_out = map(sim, get_fa),
         fa_rotations = map(fa_out, function(x)
           x[[1]]),
         fa_scores    = map(fa_out, function(x)
           x[[2]]),
         fa_pred      = map (fa_out, function(x)
           x[[3]]),
         fa_rank      = map(fa_out, function(x)
           x[[4]]))

# L2 NMF

get_nmf_l2 <- function (sim) {
  set.seed(1988)
  nmf_3 <- nmf(sim, 3, nrun = 100, method = "lee")
  set.seed(1988)
  nmf_4 <- nmf(sim, 4, nrun = 100, method = "lee")
  set.seed(1988)
  nmf_5 <- nmf(sim, 5, nrun = 100, method = "lee")
  
  bic_3 <- -log(residuals(nmf_3)) + (1/2)*(nrow(sim) + ncol(sim)) * 3 * log(nrow(sim) * ncol(sim))
  bic_4 <- -log(residuals(nmf_4)) + (1/2)*(nrow(sim) + ncol(sim)) * 4 * log(nrow(sim) * ncol(sim))
  bic_5 <- -log(residuals(nmf_5)) + (1/2)*(nrow(sim) + ncol(sim)) * 5 * log(nrow(sim) * ncol(sim))
  
  if (bic_3 <= bic_4) {
    nmf_out <- nmf_3
    rank <- 3
  } else if (bic_4 < bic_3 &&
             bic_4 <= bic_5) {
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

out_dist <- out_dist %>%
  mutate(nmf_l2_out = map(sim, get_nmf_l2),
         nmf_l2_loadings = map(nmf_l2_out, function(x) 
           x[[1]]),
         nmf_l2_scores = map(nmf_l2_out, function(x) 
           x[[2]]),
         nmf_l2_pred = map(nmf_l2_out, function(x) 
           x[[3]]),
         nmf_l2_rank = map(nmf_l2_out, function(x) 
           x[[4]]))

out_over <- out_over %>%
  mutate(nmf_l2_out = map(sim, get_nmf_l2),
         nmf_l2_loadings = map(nmf_l2_out, function(x) 
           x[[1]]),
         nmf_l2_scores = map(nmf_l2_out, function(x) 
           x[[2]]),
         nmf_l2_pred = map(nmf_l2_out, function(x) 
           x[[3]]),
         nmf_l2_rank = map(nmf_l2_out, function(x) 
           x[[4]]))

out_cor <- out_cor %>%
  mutate(nmf_l2_out = map(sim, get_nmf_l2),
         nmf_l2_loadings = map(nmf_l2_out, function(x) 
           x[[1]]),
         nmf_l2_scores = map(nmf_l2_out, function(x) 
           x[[2]]),
         nmf_l2_pred = map(nmf_l2_out, function(x) 
           x[[3]]),
         nmf_l2_rank = map(nmf_l2_out, function(x) 
           x[[4]]))

# Poisson NMF

get_nmf_p <- function (sim) {
  set.seed(1988)
  nmf_3 <- nmf(sim, 3, nrun = 100, method = "brunet")
  set.seed(1988)
  nmf_4 <- nmf(sim, 4, nrun = 100, method = "brunet")
  set.seed(1988)
  nmf_5 <- nmf(sim, 5, nrun = 100, method = "brunet")
  
  bic_3 <- -log(residuals(nmf_3)) + (1/2)*(nrow(sim) + ncol(sim)) * 3 * log(nrow(sim) * ncol(sim))
  bic_4 <- -log(residuals(nmf_4)) + (1/2)*(nrow(sim) + ncol(sim)) * 4 * log(nrow(sim) * ncol(sim))
  bic_5 <- -log(residuals(nmf_5)) + (1/2)*(nrow(sim) + ncol(sim)) * 5 * log(nrow(sim) * ncol(sim))
  
  if (bic_3 <= bic_4) {
    nmf_out <- nmf_3
    rank <- 3
  } else if (bic_4 < bic_3 &&
             bic_4 <= bic_5) {
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

out_dist <- out_dist %>%
  mutate(nmf_p_out = map(sim, get_nmf_p),
         nmf_p_loadings = map(nmf_p_out, function(x) 
           x[[1]]),
         nmf_p_scores = map(nmf_p_out, function(x) 
           x[[2]]),
         nmf_p_pred = map(nmf_p_out, function(x) 
           x[[3]]),
         nmf_p_rank = map(nmf_p_out, function(x) 
           x[[4]]))

out_over <- out_over %>%
  mutate(nmf_p_out = map(sim, get_nmf_p),
         nmf_p_loadings = map(nmf_p_out, function(x) 
           x[[1]]),
         nmf_p_scores = map(nmf_p_out, function(x) 
           x[[2]]),
         nmf_p_pred = map(nmf_p_out, function(x) 
           x[[3]]),
         nmf_p_rank = map(nmf_p_out, function(x) 
           x[[4]]))

out_cor <- out_cor %>%
  mutate(nmf_p_out = map(sim, get_nmf_p),
         nmf_p_loadings = map(nmf_p_out, function(x) 
           x[[1]]),
         nmf_p_scores = map(nmf_p_out, function(x) 
           x[[2]]),
         nmf_p_pred = map(nmf_p_out, function(x) 
           x[[3]]),
         nmf_p_rank = map(nmf_p_out, function(x) 
           x[[4]]))

#####################
# Code for Matching #
#####################

# Symmetric Subspace Distance

## Function

symm_subspace_dist <- function(U, V) {
  
  if (nrow(U) != nrow(V)) stop("Matrices must have same number of participants (rows).")
  
  qrU <- qr.Q(qr(U))
  qrV <- qr.Q(qr(V))
  
  m <- ncol(U)
  n <- ncol(V)
  
  dUV <- sqrt( max(m,n) - sum((t(qrU) %*% qrV)^2) )
  
  ratio <- dUV/sqrt( max(m,n))
  
  list(symm_subspace_dist = dUV, Similarity = ratio)
  
}

## Run

### PCA

out_dist <- out_dist %>% 
  mutate(pca_score_ssd = map2(true_scores, pca_scores, symm_subspace_dist),
         pca_score_ssdist = map(pca_score_ssd, function(x)
           x[[1]]),
         pca_score_sssimilarity = map(pca_score_ssd, function(x)
           x[[2]]),
         pca_loading_ssd = map2(true_patterns, pca_rotations, function(x,y) symm_subspace_dist(t(x), y)),
         pca_loading_ssdist = map(pca_loading_ssd, function(x)
           x[[1]]),
         pca_loading_sssimilarity = map(pca_loading_ssd, function(x)
           x[[2]]))

out_over <- out_over %>% 
  mutate(pca_score_ssd = map2(true_scores, pca_scores, symm_subspace_dist),
         pca_score_ssdist = map(pca_score_ssd, function(x)
           x[[1]]),
         pca_score_sssimilarity = map(pca_score_ssd, function(x)
           x[[2]]),
         pca_loading_ssd = map2(true_patterns, pca_rotations, function(x,y) symm_subspace_dist(t(x), y)),
         pca_loading_ssdist = map(pca_loading_ssd, function(x)
           x[[1]]),
         pca_loading_sssimilarity = map(pca_loading_ssd, function(x)
           x[[2]]))

out_cor <- out_cor %>% 
  mutate(pca_score_ssd = map2(true_scores, pca_scores, symm_subspace_dist),
         pca_score_ssdist = map(pca_score_ssd, function(x)
           x[[1]]),
         pca_score_sssimilarity = map(pca_score_ssd, function(x)
           x[[2]]),
         pca_loading_ssd = map2(true_patterns, pca_rotations, function(x,y) symm_subspace_dist(t(x), y)),
         pca_loading_ssdist = map(pca_loading_ssd, function(x)
           x[[1]]),
         pca_loading_sssimilarity = map(pca_loading_ssd, function(x)
           x[[2]]))

### FA

out_dist <- out_dist %>% 
  mutate(fa_score_ssd = map2(true_scores, fa_scores, symm_subspace_dist),
         fa_score_ssdist = map(fa_score_ssd, function(x)
           x[[1]]),
         fa_score_sssimilarity = map(fa_score_ssd, function(x)
           x[[2]]),
         fa_loading_ssd = map2(true_patterns, fa_rotations, function(x,y) symm_subspace_dist(t(x), y)),
         fa_loading_ssdist = map(fa_loading_ssd, function(x)
           x[[1]]),
         fa_loading_sssimilarity = map(fa_loading_ssd, function(x)
           x[[2]]))

out_over <- out_over %>% 
  mutate(fa_score_ssd = map2(true_scores, fa_scores, symm_subspace_dist),
         fa_score_ssdist = map(fa_score_ssd, function(x)
           x[[1]]),
         fa_score_sssimilarity = map(fa_score_ssd, function(x)
           x[[2]]),
         fa_loading_ssd = map2(true_patterns, fa_rotations, function(x,y) symm_subspace_dist(t(x), y)),
         fa_loading_ssdist = map(fa_loading_ssd, function(x)
           x[[1]]),
         fa_loading_sssimilarity = map(fa_loading_ssd, function(x)
           x[[2]]))

out_cor <- out_cor %>% 
  mutate(fa_score_ssd = map2(true_scores, fa_scores, symm_subspace_dist),
         fa_score_ssdist = map(fa_score_ssd, function(x)
           x[[1]]),
         fa_score_sssimilarity = map(fa_score_ssd, function(x)
           x[[2]]),
         fa_loading_ssd = map2(true_patterns, fa_rotations, function(x,y) symm_subspace_dist(t(x), y)),
         fa_loading_ssdist = map(fa_loading_ssd, function(x)
           x[[1]]),
         fa_loading_sssimilarity = map(fa_loading_ssd, function(x)
           x[[2]]))

### NMF L2

out_dist <- out_dist %>% 
  mutate(nmf_l2_score_ssd = map2(true_scores, nmf_l2_scores, symm_subspace_dist),
         nmf_l2_score_ssdist = map(nmf_l2_score_ssd, function(x)
           x[[1]]),
         nmf_l2_score_sssimilarity = map(nmf_l2_score_ssd, function(x)
           x[[2]]),
         nmf_l2_loading_ssd = map2(true_patterns, nmf_l2_loadings, function(x,y) symm_subspace_dist(t(x), t(y))),
         nmf_l2_loading_ssdist = map(nmf_l2_loading_ssd, function(x)
           x[[1]]),
         nmf_l2_loading_sssimilarity = map(nmf_l2_loading_ssd, function(x)
           x[[2]]))

out_over <- out_over %>% 
  mutate(nmf_l2_score_ssd = map2(true_scores, nmf_l2_scores, symm_subspace_dist),
         nmf_l2_score_ssdist = map(nmf_l2_score_ssd, function(x)
           x[[1]]),
         nmf_l2_score_sssimilarity = map(nmf_l2_score_ssd, function(x)
           x[[2]]),
         nmf_l2_loading_ssd = map2(true_patterns, nmf_l2_loadings, function(x,y) symm_subspace_dist(t(x), t(y))),
         nmf_l2_loading_ssdist = map(nmf_l2_loading_ssd, function(x)
           x[[1]]),
         nmf_l2_loading_sssimilarity = map(nmf_l2_loading_ssd, function(x)
           x[[2]]))

out_cor <- out_cor %>% 
  mutate(nmf_l2_score_ssd = map2(true_scores, nmf_l2_scores, symm_subspace_dist),
         nmf_l2_score_ssdist = map(nmf_l2_score_ssd, function(x)
           x[[1]]),
         nmf_l2_score_sssimilarity = map(nmf_l2_score_ssd, function(x)
           x[[2]]),
         nmf_l2_loading_ssd = map2(true_patterns, nmf_l2_loadings, function(x,y) symm_subspace_dist(t(x), t(y))),
         nmf_l2_loading_ssdist = map(nmf_l2_loading_ssd, function(x)
           x[[1]]),
         nmf_l2_loading_sssimilarity = map(nmf_l2_loading_ssd, function(x)
           x[[2]]))

### NMF Poisson

out_dist <- out_dist %>% 
              mutate(nmf_p_score_ssd = map2(true_scores, nmf_p_scores, symm_subspace_dist),
            nmf_p_score_ssdist = map(nmf_p_score_ssd, function(x)
              x[[1]]),
            nmf_p_score_sssimilarity = map(nmf_p_score_ssd, function(x)
              x[[2]]),
            nmf_p_loading_ssd = map2(true_patterns, nmf_p_loadings, function(x,y) symm_subspace_dist(t(x), t(y))),
            nmf_p_loading_ssdist = map(nmf_p_loading_ssd, function(x)
              x[[1]]),
            nmf_p_loading_sssimilarity = map(nmf_p_loading_ssd, function(x)
              x[[2]]))

out_over <- out_over %>% 
  mutate(nmf_p_score_ssd = map2(true_scores, nmf_p_scores, symm_subspace_dist),
         nmf_p_score_ssdist = map(nmf_p_score_ssd, function(x)
           x[[1]]),
         nmf_p_score_sssimilarity = map(nmf_p_score_ssd, function(x)
           x[[2]]),
         nmf_p_loading_ssd = map2(true_patterns, nmf_p_loadings, function(x,y) symm_subspace_dist(t(x), t(y))),
         nmf_p_loading_ssdist = map(nmf_p_loading_ssd, function(x)
           x[[1]]),
         nmf_p_loading_sssimilarity = map(nmf_p_loading_ssd, function(x)
           x[[2]]))

out_cor <- out_cor %>% 
  mutate(nmf_p_score_ssd = map2(true_scores, nmf_p_scores, symm_subspace_dist),
         nmf_p_score_ssdist = map(nmf_p_score_ssd, function(x)
           x[[1]]),
         nmf_p_score_sssimilarity = map(nmf_p_score_ssd, function(x)
           x[[2]]),
         nmf_p_loading_ssd = map2(true_patterns, nmf_p_loadings, function(x,y) symm_subspace_dist(t(x), t(y))),
         nmf_p_loading_ssdist = map(nmf_p_loading_ssd, function(x)
           x[[1]]),
         nmf_p_loading_sssimilarity = map(nmf_p_loading_ssd, function(x)
           x[[2]]))

















           
           
           
