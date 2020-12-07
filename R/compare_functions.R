# Functions for BN2MF project

# Functions to run other models
# & to choose number of factors/components

# PCA
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
  return(list(rotations = rotations, scores = scores, pred = pred, rank = rank))
}

# Factor Analysis
get_fa <- function (sim) {
  set.seed(1988)
  fa_3 <- fa(sim, 3, scores = "regression", rotate = "promax")
  fa_4 <- fa(sim, 4, scores = "regression", rotate = "promax")
  fa_5 <- fa(sim, 5, scores = "regression", rotate = "promax")
  
  # Choose the model with the LOWEST BIC
  if (min(fa_5$BIC, fa_4$BIC, fa_3$BIC) == fa_5$BIC) {
    fa_out <- fa_5
    rank <- 5
  } else if (min(fa_5$BIC, fa_4$BIC, fa_3$BIC) == fa_4$BIC) {
    fa_out <- fa_4
    rank <- 4
  } else {
    fa_out <- fa_3
    rank <- 3
  }
  
  loadings <- matrix(fa_out$loadings, ncol = ncol(fa_out$scores))
  fa_scores <- fa_out$scores
  pred <- fa_scores %*% t(loadings)
  return(list(loadings = loadings, fa_scores = fa_scores, pred = pred, rank = rank))
}

# NMF L2
get_nmf_l2 <- function (sim) {
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
  if (min(bic_3, bic_4, bic_5) == bic_3) {
    nmf_out <- nmf_3
    rank <- 3
  } else if (min(bic_3, bic_4, bic_5) == bic_4) {
    nmf_out <- nmf_4
    rank <- 4
  } else {
    nmf_out <- nmf_5
    rank <- 5
  }
  
  basis <- basis(nmf_out)
  coef <- coef(nmf_out)
  pred <- basis %*% coef
  return(list(coef = coef, basis = basis, pred = pred, rank = rank))
}

# NMF Poisson
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
  
  if (min(bic_3, bic_4, bic_5) == bic_3) {
    nmf_out <- nmf_3
    rank <- 3
  } else if (min(bic_3, bic_4, bic_5) == bic_4) {
    nmf_out <- nmf_4
    rank <- 4
  } else {
    nmf_out <- nmf_5
    rank <- 5
  }
  
  basis <- basis(nmf_out)
  coef <- coef(nmf_out)
  pred <- basis %*% coef
  return(list(coef = coef, basis = basis, pred = pred, rank = rank))
}

# Symmetric Subspace Distance 
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

# Cosine distance
cos_dist <- function(a, b){
  # a = as.matrix(a)
  # b = as.matrix(b)
  a = as.vector(a)
  b = as.vector(b)

  #cos_dist = c()
  # for (i in 1:ncol(a)) {
  #   cos_dist[i] = (a[,i] %*% b[,i]) / (norm(a[,i], "2")*norm(b[,i], "2"))
  # }
  cos_sim = (a %*% b) / (norm(a, "2")*norm(b, "2"))
  return(1-cos_sim)
}

## For Simulations with Different Dimensions

# Factor analysis
get_fa_dim <- function (sim, patterns) {
  
  set.seed(1988)
  
  patternsm1 = ifelse(patterns == 1, 3, patterns - 1)
  patternsp1 = patterns+1
  
  fa_3 <- fa(sim, patternsm1, scores = "regression", rotate = "promax")
  fa_4 <- fa(sim, patterns, scores = "regression", rotate = "promax")
  fa_5 <- fa(sim, patternsp1, scores = "regression", rotate = "promax")
  
  # Choose the model with the LOWEST BIC
  if (min(fa_5$BIC, fa_4$BIC, fa_3$BIC) == fa_5$BIC) {
    fa_out <- fa_5
    rank <- patternsm1
  } else if (min(fa_5$BIC, fa_4$BIC, fa_3$BIC) == fa_4$BIC) {
    fa_out <- fa_4
    rank <- patterns
  } else {
    fa_out <- fa_3
    rank <- patternsp1
  }
  
  loadings <- matrix(fa_out$loadings, ncol = ncol(fa_out$scores))
  fa_scores <- fa_out$scores
  pred <- fa_scores %*% t(loadings)
  return(list(loadings = loadings, fa_scores = fa_scores, pred = pred, rank = rank))
}

# NMF L2
get_nmf_l2_dim <- function (sim, patterns) {
  
  patternsm1 = ifelse(patterns == 1, 3, patterns - 1)
  patternsp1 = patterns+1
  
  set.seed(1988)
  nmf_3 <- nmf(sim, patternsm1, nrun = 100, method = "lee")
  nmf_4 <- nmf(sim, patterns, nrun = 100, method = "lee")
  nmf_5 <- nmf(sim, patternsp1, nrun = 100, method = "lee")
  
  # Calculate BIC for each
  bic_3 <- sum((sim - (basis(nmf_3)%*%coef(nmf_3)))^2) + 
    (1/2)*(nrow(sim) + ncol(sim)) * patternsm1 * log(nrow(sim) * ncol(sim))
  bic_4 <- sum((sim - (basis(nmf_4)%*%coef(nmf_4)))^2) + 
    (1/2)*(nrow(sim) + ncol(sim)) * patterns * log(nrow(sim) * ncol(sim))
  bic_5 <- sum((sim - (basis(nmf_5)%*%coef(nmf_5)))^2) +  
    (1/2)*(nrow(sim) + ncol(sim)) * patternsp1 * log(nrow(sim) * ncol(sim))
  
  # Choose model with lowest BIC
  if (min(bic_3, bic_4, bic_5) == bic_3) {
    nmf_out <- nmf_3
    rank <- patternsm1
  } else if (min(bic_3, bic_4, bic_5) == bic_4) {
    nmf_out <- nmf_4
    rank <- patterns
  } else {
    nmf_out <- nmf_5
    rank <- patternsp1
  }
  
  basis <- basis(nmf_out)
  coef <- coef(nmf_out)
  pred <- basis %*% coef
  return(list(coef = coef, basis = basis, pred = pred, rank = rank))
}

# NMF Poisson
get_nmf_p_dim <- function (sim, patterns) {
  
  patternsm1 = ifelse(patterns == 1, 3, patterns - 1)
  patternsp1 = patterns+1
  
  set.seed(1988)
  nmf_3 <- nmf(sim, patternsm1, nrun = 100, method = "brunet")
  nmf_4 <- nmf(sim, patterns, nrun = 100, method = "brunet")
  nmf_5 <- nmf(sim, patternsp1, nrun = 100, method = "brunet")
  
  bic_3 <- -sum((sim * log(basis(nmf_3) %*% coef(nmf_3))) - (basis(nmf_3) %*% coef(nmf_3))) + 
    (1/2)*(nrow(sim) + ncol(sim)) * patternsm1 * log(nrow(sim) * ncol(sim))
  bic_4 <- -sum((sim * log(basis(nmf_4) %*% coef(nmf_4))) - (basis(nmf_4) %*% coef(nmf_4))) + 
    (1/2)*(nrow(sim) + ncol(sim)) * patterns * log(nrow(sim) * ncol(sim))
  bic_5 <- -sum((sim * log(basis(nmf_5) %*% coef(nmf_5))) - (basis(nmf_5) %*% coef(nmf_5))) + 
    (1/2)*(nrow(sim) + ncol(sim)) * patternsp1 * log(nrow(sim) * ncol(sim))
  
  if (min(bic_3, bic_4, bic_5) == bic_3) {
    nmf_out <- nmf_3
    rank <- patternsm1
  } else if (min(bic_3, bic_4, bic_5) == bic_4) {
    nmf_out <- nmf_4
    rank <- patterns
  } else {
    nmf_out <- nmf_5
    rank <- patternsp1
  }
  
  basis <- basis(nmf_out)
  coef <- coef(nmf_out)
  pred <- basis %*% coef
  return(list(coef = coef, basis = basis, pred = pred, rank = rank))
}
