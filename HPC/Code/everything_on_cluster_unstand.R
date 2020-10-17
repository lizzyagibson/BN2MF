#########################################################
# Sims, npBNMF, regular NMF (L2 & Poisson), PCA, and FA #
#########################################################
# Lizzy Gibson # 6/22/2020 ##############################
#########################################################

# Packages
library(tidyverse)
library(registry, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(pkgmaker, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(rngtools, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(NMF, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(GPArotation, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(psych, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

# Read in Sims
load("Data/sim_dist_old_1012.RDA")
load("Data/sim_cor_old_1012.RDA")
load("Data/sim_over_old_1012.RDA")

## read job number from system environment
## This only works if run on cluster!!
job_num = as.integer(Sys.getenv("SGE_TASK_ID"))
job_num

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
  rotations <- as_tibble(rot[, 1:rank])
  scores <- ex[, 1:rank]
  
  # Predicted values
  pred <- scores %*% t(rotations) + matrix(rep(apply(sim, 2, mean), each= nrow(scores)), nrow = nrow(scores))
  
  list(rotations = rotations, scores = scores, pred = pred, rank = rank)
}

## Run

out_dist_un <- sim_dist[job_num,] %>% 
  mutate(pca_out = map(sim, get_pca),
         pca_rotations = map(pca_out, function(x)
                          x[[1]]),
         pca_scores    = map(pca_out, function(x)
                          x[[2]]),
         pca_pred      = map (pca_out, function(x)
                          x[[3]]),
         pca_rank      = map(pca_out, function(x)
                          x[[4]]))

out_over_un <- sim_over[job_num,] %>% 
  mutate(pca_out = map(sim, get_pca),
         pca_rotations = map(pca_out, function(x)
           x[[1]]),
         pca_scores    = map(pca_out, function(x)
           x[[2]]),
         pca_pred      = map (pca_out, function(x)
           x[[3]]),
         pca_rank      = map(pca_out, function(x)
           x[[4]]))

out_cor_un <- sim_cor[job_num,] %>% 
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

out_dist_un <- out_dist_un %>% 
  mutate(fa_out = map(sim, get_fa),
         fa_rotations = map(fa_out, function(x)
           x[[1]]),
         fa_scores    = map(fa_out, function(x)
           x[[2]]),
         fa_pred      = map (fa_out, function(x)
           x[[3]]),
         fa_rank      = map(fa_out, function(x)
           x[[4]]))

out_over_un <- out_over_un %>% 
  mutate(fa_out = map(sim, get_fa),
         fa_rotations = map(fa_out, function(x)
           x[[1]]),
         fa_scores    = map(fa_out, function(x)
           x[[2]]),
         fa_pred      = map (fa_out, function(x)
           x[[3]]),
         fa_rank      = map(fa_out, function(x)
           x[[4]]))

out_cor_un <- out_cor_un %>% 
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

## Function

get_nmf_l2 <- function (sim) {
  set.seed(1988)
  nmf_3 <- nmf(sim, 3, nrun = 100, method = "lee")
#  set.seed(1988)
  nmf_4 <- nmf(sim, 4, nrun = 100, method = "lee")
#  set.seed(1988)
  nmf_5 <- nmf(sim, 5, nrun = 100, method = "lee")
  
  bic_3 <- sum((sim - (basis(nmf_3)%*%coef(nmf_3)))^2) + 
    (1/2)*(nrow(sim) + ncol(sim)) * 3 * log(nrow(sim) * ncol(sim))
  bic_4 <- sum((sim - (basis(nmf_4)%*%coef(nmf_4)))^2) + 
    (1/2)*(nrow(sim) + ncol(sim)) * 4 * log(nrow(sim) * ncol(sim))
  bic_5 <- sum((sim - (basis(nmf_5)%*%coef(nmf_5)))^2) +  
    (1/2)*(nrow(sim) + ncol(sim)) * 5 * log(nrow(sim) * ncol(sim))
  
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

## Run

out_dist_un <- out_dist_un %>%
  mutate(nmf_l2_out = map(sim, get_nmf_l2),
         nmf_l2_loadings = map(nmf_l2_out, function(x) 
           x[[1]]),
         nmf_l2_scores = map(nmf_l2_out, function(x) 
           x[[2]]),
         nmf_l2_pred = map(nmf_l2_out, function(x) 
           x[[3]]),
         nmf_l2_rank = map(nmf_l2_out, function(x) 
           x[[4]]))

out_over_un <- out_over_un %>%
  mutate(nmf_l2_out = map(sim, get_nmf_l2),
         nmf_l2_loadings = map(nmf_l2_out, function(x) 
           x[[1]]),
         nmf_l2_scores = map(nmf_l2_out, function(x) 
           x[[2]]),
         nmf_l2_pred = map(nmf_l2_out, function(x) 
           x[[3]]),
         nmf_l2_rank = map(nmf_l2_out, function(x) 
           x[[4]]))

out_cor_un <- out_cor_un %>%
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
  # set.seed(1988)
  nmf_4 <- nmf(sim, 4, nrun = 100, method = "brunet")
  # set.seed(1988)
  nmf_5 <- nmf(sim, 5, nrun = 100, method = "brunet")
  
  bic_3 <- -sum((sim * log(basis(nmf_3) %*% coef(nmf_3))) - (basis(nmf_3) %*% coef(nmf_3))) + 
    (1/2)*(nrow(sim) + ncol(sim)) * 3 * log(nrow(sim) * ncol(sim))
  bic_4 <- -sum((sim * log(basis(nmf_4) %*% coef(nmf_4))) - (basis(nmf_4) %*% coef(nmf_4))) + 
    (1/2)*(nrow(sim) + ncol(sim)) * 4 * log(nrow(sim) * ncol(sim))
  bic_5 <- -sum((sim * log(basis(nmf_5) %*% coef(nmf_5))) - (basis(nmf_5) %*% coef(nmf_5))) + 
    (1/2)*(nrow(sim) + ncol(sim)) * 5 * log(nrow(sim) * ncol(sim))
  
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

out_dist_un <- out_dist_un %>%
  mutate(nmf_p_out = map(sim, get_nmf_p),
         nmf_p_loadings = map(nmf_p_out, function(x) 
           x[[1]]),
         nmf_p_scores = map(nmf_p_out, function(x) 
           x[[2]]),
         nmf_p_pred = map(nmf_p_out, function(x) 
           x[[3]]),
         nmf_p_rank = map(nmf_p_out, function(x) 
           x[[4]]))

out_over_un <- out_over_un %>%
  mutate(nmf_p_out = map(sim, get_nmf_p),
         nmf_p_loadings = map(nmf_p_out, function(x) 
           x[[1]]),
         nmf_p_scores = map(nmf_p_out, function(x) 
           x[[2]]),
         nmf_p_pred = map(nmf_p_out, function(x) 
           x[[3]]),
         nmf_p_rank = map(nmf_p_out, function(x) 
           x[[4]]))

out_cor_un <- out_cor_un %>%
  mutate(nmf_p_out = map(sim, get_nmf_p),
         nmf_p_loadings = map(nmf_p_out, function(x) 
           x[[1]]),
         nmf_p_scores = map(nmf_p_out, function(x) 
           x[[2]]),
         nmf_p_pred = map(nmf_p_out, function(x) 
           x[[3]]),
         nmf_p_rank = map(nmf_p_out, function(x) 
           x[[4]]))

###############################
# Symmetric Subspace Distance #
###############################

symm_subspace_dist <- function(U, V) {
  
  if (nrow(U) != max(nrow(U), ncol(U))) {U <- t(U)}
  if (nrow(V) != max(nrow(V), ncol(V))) {V <- t(V)}
  
  qrU <- qr.Q(qr(U))
  qrV <- qr.Q(qr(V))
  
  m <- ncol(U)
  n <- ncol(V)
  
  dUV <- sqrt( max(m,n) - sum((t(qrU) %*% qrV)^2) )
  
  ratio <- dUV/sqrt( max(m,n))
  
  ratio
  
}

out_dist_un <- out_dist_un %>% 
  mutate(pca_norm = map2(sim, pca_pred, function(x,y) norm(x-y, "F")/norm(x, "F")),
         fa_norm = map2(sim, fa_pred, function(x,y) norm(x-y, "F")/norm(x, "F")),
         nmf_l2_norm = map2(sim, nmf_l2_pred, function(x,y) norm(x-y, "F")/norm(x, "F")),
         nmf_P_norm = map2(sim, nmf_p_pred, function(x,y) norm(x-y, "F")/norm(x, "F"))) %>% 
  mutate(pca_rotation_ssdist   = map2(patterns, pca_rotations, symm_subspace_dist),
         pca_scores_ssdist     = map2(scores, pca_scores, symm_subspace_dist),
         fa_rotations_ssdist   = map2(patterns, fa_rotations, symm_subspace_dist),
         fa_scores_ssdist      = map2(scores, fa_scores, symm_subspace_dist),
         nmf_l2_loading_ssdist = map2(patterns, nmf_l2_loadings, symm_subspace_dist),
         nmf_l2_scores_ssdist  = map2(scores, nmf_l2_scores, symm_subspace_dist),
         nmf_p_loading_ssdist  = map2(patterns, nmf_p_loadings, symm_subspace_dist),
         nmf_p_scores_ssdist   = map2(scores, nmf_p_scores, symm_subspace_dist))

out_over_un <- out_over_un %>% 
  mutate(pca_norm = map2(sim, pca_pred, function(x,y) norm(x-y, "F")/norm(x, "F")),
         fa_norm = map2(sim, fa_pred, function(x,y) norm(x-y, "F")/norm(x, "F")),
         nmf_l2_norm = map2(sim, nmf_l2_pred, function(x,y) norm(x-y, "F")/norm(x, "F")),
         nmf_P_norm = map2(sim, nmf_p_pred, function(x,y) norm(x-y, "F")/norm(x, "F"))) %>% 
  mutate(pca_rotation_ssdist   = map2(patterns, pca_rotations, symm_subspace_dist),
         pca_scores_ssdist     = map2(scores, pca_scores, symm_subspace_dist),
         fa_rotations_ssdist   = map2(patterns, fa_rotations, symm_subspace_dist),
         fa_scores_ssdist      = map2(scores, fa_scores, symm_subspace_dist),
         nmf_l2_loading_ssdist = map2(patterns, nmf_l2_loadings, symm_subspace_dist),
         nmf_l2_scores_ssdist  = map2(scores, nmf_l2_scores, symm_subspace_dist),
         nmf_p_loading_ssdist  = map2(patterns, nmf_p_loadings, symm_subspace_dist),
         nmf_p_scores_ssdist   = map2(scores, nmf_p_scores, symm_subspace_dist))

out_cor_un <- out_cor_un %>% 
  mutate(pca_norm = map2(sim, pca_pred, function(x,y) norm(x-y, "F")/norm(x, "F")),
         fa_norm = map2(sim, fa_pred, function(x,y) norm(x-y, "F")/norm(x, "F")),
         nmf_l2_norm = map2(sim, nmf_l2_pred, function(x,y) norm(x-y, "F")/norm(x, "F")),
         nmf_P_norm = map2(sim, nmf_p_pred, function(x,y) norm(x-y, "F")/norm(x, "F"))) %>% 
  mutate(pca_rotation_ssdist   = map2(patterns, pca_rotations, symm_subspace_dist),
         pca_scores_ssdist     = map2(scores, pca_scores, symm_subspace_dist),
         fa_rotations_ssdist   = map2(patterns, fa_rotations, symm_subspace_dist),
         fa_scores_ssdist      = map2(scores, fa_scores, symm_subspace_dist),
         nmf_l2_loading_ssdist = map2(patterns, nmf_l2_loadings, symm_subspace_dist),
         nmf_l2_scores_ssdist  = map2(scores, nmf_l2_scores, symm_subspace_dist),
         nmf_p_loading_ssdist  = map2(patterns, nmf_p_loadings, symm_subspace_dist),
         nmf_p_scores_ssdist   = map2(scores, nmf_p_scores, symm_subspace_dist))

###
                           
out_dist_un <- out_dist_un %>% dplyr::select(-grep("_un|_out", colnames(.)))
out_over_un <- out_over_un %>% dplyr::select(-grep("_un|_out", colnames(.)))
out_cor_un <- out_cor_un %>% dplyr::select(-grep("_un|_out", colnames(.)))

save(out_dist_un, file = paste0("out_", job_num, "_dist_un.RDA"))
save(out_over_un, file = paste0("out_", job_num, "_over_un.RDA"))
save(out_cor_un, file = paste0("out_", job_num, "_cor_un.RDA"))

