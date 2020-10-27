library(bindata) # multivariate binomial dist
library(NMF) # rmatrix
library(MASS) # mvrnorm
library(tidyverse)
library(R.matlab)
library(psych)

#####
# Distinct Patterns
#####
### Create patterns

create_patterns_dist <- function (seed) {
  labels <-  c("pcb28", "pcb66", "pcb74", "pcb99", "pcb105", "pcb118", "pcb138_158", "pcb146", "pcb153", "pcb156", "pcb167",
               "pcb170", "pcb178", "pcb183", "pcb187", "pcb180", "pcb189", "pcb194", "pcb196_203", "pcb199", "pcb206", "pcb209",
               "BDE17", "BDE28", "BDE47", "BDE66", "BDE85", "BDE99", "BDE100", "BDE153", "BDE154", "BDE183", "BDE209", "MECPP",
               "MEHHP", "MEOHP", "MCPP", "MIBP", "MBP", "MBZP", "MEP", "MEHP", "dcp_24", "dcp_25", "b_pb", "bp_3", "m_pb",
               "p_pb", "tcs", "bpa")
  
  set.seed(seed)
  patterns_dist <- rbind(c(runif(22, 1, 2), rep(0, 11), rep(0, 17)),
                         c(rep(0, 22), runif(11, 1, 2), rep(0, 17)),
                         c(rep(0, 33), runif(9, 1, 2), rep(0, 8)),
                         c(rep(0, 33), rep(0, 9), runif(8, 1, 2)))
  colnames(patterns_dist) <- labels
  patterns_dist
}

patterns_dist_iter <- tibble(seed = 1:100) %>% 
  mutate(patterns = map(seed, create_patterns_dist))

### Simulate Scores

create_scores_dist <- function (seed) {
  scores_dist <- matrix(nrow = 1000, ncol = 4)
  set.seed(seed)
  scores_dist <- exp(mvrnorm(n = 1000, mu = rep(1, 4), 
                             Sigma = diag(rep(1, 4)))) # Independent log normal dist
  scores_dist
}

scores_dist_iter <- tibble(seed = 1:100) %>% 
  mutate(scores = map(seed, create_scores_dist))

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

sim_dist1 <- sim_dist0 %>% 
  mutate(sim = map2(seed, chem, ~add_noise(seed = .x, chem = .y)[[1]]),
         noise = map2(seed, chem, function(x, y) {add_noise(x, y)[[2]]}))

#####
# 2: Correlated Patterns
#####
### Simulate Patterns
create_patterns_cor <- function (seed) {
  labels <-  c("pcb28", "pcb66", "pcb74", "pcb99", "pcb105", "pcb118", "pcb138_158", "pcb146", "pcb153", "pcb156", "pcb167",
               "pcb170", "pcb178", "pcb183", "pcb187", "pcb180", "pcb189", "pcb194", "pcb196_203", "pcb199", "pcb206", "pcb209",
               "BDE17", "BDE28", "BDE47", "BDE66", "BDE85", "BDE99", "BDE100", "BDE153", "BDE154", "BDE183", "BDE209", "MECPP",
               "MEHHP", "MEOHP", "MCPP", "MIBP", "MBP", "MBZP", "MEP", "MEHP", "dcp_24", "dcp_25", "b_pb", "bp_3", "m_pb",
               "p_pb", "tcs", "bpa")
  
  set.seed(seed)
  patterns_cor <- rbind(c(runif(22, 1, 2),   runif(11, 0.5, 1), rep(0, 17)),
                        c(rep(0, 22),        runif(11, 1, 2),   runif(9, 0.5, 1),  rep(0, 8)),
                        c(rep(0, 33),                           runif(9, 1, 2),    runif(8, 0.5, 1)),
                        c(runif(22, 0.5, 1), rep(0, 11),        rep(0, 9),         runif(8, 1, 2)))
  colnames(patterns_cor) <- labels
  patterns_cor
}
patterns_cor_iter <- tibble(seed = 1:100) %>% 
  mutate(patterns = map(seed, create_patterns_cor))

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
sim_cor1 <- sim_cor0 %>% 
  mutate(sim = map2(seed, chem, ~add_noise(seed = .x, chem = .y)[[1]]),
         noise = map2(seed, chem, function(x, y) {add_noise(x, y)[[2]]})) 

#####
# 3: Overlapping Patterns
#####
### Simulate Chemical Exposures
sim_over0 <- full_join(scores_dist_iter, patterns_cor_iter, by = "seed")
sim_over0 <- sim_over0 %>% 
  mutate(chem = map2(scores, patterns, `%*%`))

### Simulate Noise

sim_over1 <- sim_over0 %>% 
  mutate(sim = map2(seed, chem, ~add_noise(seed = .x, chem = .y)[[1]]),
         noise = map2(seed, chem, function(x, y) {add_noise(x, y)[[2]]}))

#####
# Save
#####
isee_all <- rbind(sim_dist1 %>% mutate(data = "distinct"),
      sim_over1 %>% mutate(data = "overlapping"),
      sim_cor1 %>% mutate(data = "correlated"))
#####
# Load
#####
# save(isee_all, file = "./Sims/sims_isee_all.RDA")
load("./Sims/sims_isee_all.RDA")

## Compare
load("./Sims/sim_dist_old_1012.RDA")
load("./Sims/sim_over_old_1012.RDA")
load("./Sims/sim_cor_old_1012.RDA")

all_equal(sim_over, isee_all %>% filter(data == "overlapping") %>% select(-data))
all_equal(sim_dist, isee_all %>% filter(data == "distinct") %>% select(-data))
all_equal(sim_cor, isee_all %>% filter(data == "correlated") %>% select(-data))

new_over <- isee_all %>% filter(data == "overlapping") %>% select(sim)
old_over <- sim_over %>% select(sim)
all.equal(new_over[1,1][[1]], old_over[1,1][[1]])

#####
# Run Models
#####

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

isee_all <- isee_all %>% 
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

#####
# Symmetric Subspace Distance #
#####

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

isee_all <- isee_all %>% 
  mutate(pca_norm              = map2(sim, pca_pred,    function(x,y) norm(x-y, "F")/norm(x, "F")),
         fa_norm               = map2(sim, fa_pred,     function(x,y) norm(x-y, "F")/norm(x, "F")),
         nmf_l2_norm           = map2(sim, nmf_l2_pred, function(x,y) norm(x-y, "F")/norm(x, "F")),
         nmf_p_norm            = map2(sim, nmf_p_pred,  function(x,y) norm(x-y, "F")/norm(x, "F")),
         pca_rotation_ssdist   = map2(patterns, pca_rotations,   symm_subspace_dist),
         pca_scores_ssdist     = map2(scores,   pca_scores,      symm_subspace_dist),
         fa_rotations_ssdist   = map2(patterns, fa_rotations,    symm_subspace_dist),
         fa_scores_ssdist      = map2(scores,   fa_scores,       symm_subspace_dist),
         nmf_l2_loading_ssdist = map2(patterns, nmf_l2_loadings, symm_subspace_dist),
         nmf_l2_scores_ssdist  = map2(scores,   nmf_l2_scores,   symm_subspace_dist),
         nmf_p_loading_ssdist  = map2(patterns, nmf_p_loadings,  symm_subspace_dist),
         nmf_p_scores_ssdist   = map2(scores,   nmf_p_scores,    symm_subspace_dist))

#####
# Compare
#####

isee_local <- isee_all %>% select(-grep("_out", colnames(.)))

# save(isee_local, file = "./HPC/Rout/isee_local.RDA")
load("./HPC/Rout/isee_local.RDA")

load(file = "./HPC/Rout/isee.RDA")
isee_hpc <- isee

# PCA --  all same
pca_local <- isee_local %>% select(seed, data, grep("pca", colnames(.)))
pca_hpc <- isee_hpc %>% select(seed, data, grep("pca", colnames(.)))

all.equal(pca_hpc[,], pca_local[,])
all.equal(pca_local %>% unnest(pca_rank) %>% pull(pca_rank), pca_hpc %>% unnest(pca_rank) %>% pull(pca_rank))

# FA -- DIFFERENT
fa_local <- isee_local %>% select(seed, data, grep("fa", colnames(.)))
fa_hpc <- isee_hpc %>% select(seed, data, grep("fa", colnames(.)))

all.equal(fa_hpc[,], fa_local[,])
(fa_local %>% unnest(fa_rank) %>% pull(fa_rank))
(fa_hpc   %>% unnest(fa_rank) %>% pull(fa_rank))

# NMF -- DIFFERENT
nmf_l2_local <- isee_local %>% select(seed, data, grep("nmf_l2", colnames(.)))
nmf_l2_hpc <- isee_hpc %>% select(seed, data, grep("nmf_l2", colnames(.)))

all.equal(nmf_l2_hpc[,], nmf_l2_local[,])
(nmf_l2_local %>% unnest(nmf_l2_rank) %>% pull(nmf_l2_rank))
(nmf_l2_hpc   %>% unnest(nmf_l2_rank) %>% pull(nmf_l2_rank))

#####
# Viz
#####

# Relative Error
isee_e <- isee_all %>% 
  dplyr::select(seed, data, sim, chem, grep("pred", colnames(.))) %>% 
  pivot_longer(c(pca_pred:fa_pred)) %>% 
  mutate(l2_true = map2(chem, value, function (x,y) norm(x-y, "F")/norm(x, "F")),
         l2_sim = map2(sim, value, function (x,y) norm(x-y, "F")/norm(x, "F")),
         name = str_remove(name, "_pred")) %>% 
  dplyr::select(seed, data, name, l2_true, l2_sim) %>% 
  unnest(c(l2_sim, l2_true))

isee_e %>%
  ggplot(aes(x = name, y = l2_true, color = name, fill = name)) +
  geom_boxplot(alpha = 0.5) +
  facet_grid(. ~ data, scales = "free") + 
  geom_vline(xintercept = 0, color = "pink", linetype = "dashed", size = 0.5) +
  scale_y_log10() +
  theme(legend.position = "none") + 
  labs(y = "Relative Predictive Error",
       title = "vs PRE NOISE TRUTH")

isee_s <- isee_all %>% dplyr::select(seed, data, grep("_ssdist", colnames(.))) %>% 
  unnest(cols = grep("_ssdist", colnames(.))) %>% 
  pivot_longer(grep("_ssdist", colnames(.)),
               names_to = "type",
               values_to = "value") %>%
  mutate(model = str_sub(type, 1, 6),
         model = ifelse(str_detect(model, 'pca'), 'PCA', model),
         model = ifelse(str_detect(model, 'l2'), 'L2 NMF', model),
         model = ifelse(str_detect(model, 'fa'), 'FA', model),
         model = ifelse(str_detect(model, '_p_'), 'Poisson NMF', model),
         model = ifelse(str_detect(model, 'bnmf'), 'BNMF', model),
         type = ifelse(str_detect(type, 'scores'), 'Scores', "Loadings"))

isee_s %>% 
  ggplot(aes(x = model, y = value, color = model)) +
  geom_boxplot() +
  facet_grid(type~data) + 
  geom_hline(yintercept = 0.5, color = "pink", linetype = "dashed", size = 0.5) +
  theme(legend.position = "none") + 
  # scale_y_log10() +
  labs(y = "Symmetric Subspace Distance")

#####
# Rank
#####

isee_all %>%
  dplyr::select(seed, data, grep("rank", colnames(.))) %>%
  unnest(c(pca_rank:fa_rank)) %>%
  pivot_longer(cols = pca_rank:fa_rank) %>%
  mutate(value = ifelse(value > 5, "> 5", value)) %>% 
  group_by(name, value, data) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  filter(data == "distinct") %>% 
  dplyr::select(-data) %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0) %>% 
  # dplyr::select(name, `1`,`2`,`3`,`4`, `5`) %>% 
  knitr::kable(caption = "Distinct Simulations: Patterns Identified")

isee_all %>%
  dplyr::select(seed, data, grep("rank", colnames(.))) %>%
  unnest(c(pca_rank:fa_rank)) %>%
  pivot_longer(cols = pca_rank:fa_rank) %>%
  mutate(value = ifelse(value > 5, "> 5", value)) %>% 
  group_by(name, value, data) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  filter(data == "overlapping") %>% 
  dplyr::select(-data) %>% 
  pivot_wider(names_from = value,
              values_from = n) %>% 
  mutate_if(is.integer, replace_na, 0) %>% 
  # dplyr::select(name, `1`, `2`, `3`, `4`, `5`) %>% 
  knitr::kable(caption = "Overlapping Simulations: Patterns Identified")

