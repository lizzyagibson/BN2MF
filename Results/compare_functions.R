# Functions for BN2MF project

## Functions to run other models
## & to choose number of factors/components

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

# Factor analysis
get_fa <- function (sim, patterns) {
  
  set.seed(1988)
  
  patternsm1 = ifelse(patterns == 1, 3, patterns - 1)
  patternsp1 = patterns+1
  
  fa_3 <- try(fa(sim, patternsm1, scores = "regression", rotate = "promax", fm = "ml"))
  fa_4 <- try(fa(sim, patterns,   scores = "regression", rotate = "promax", fm = "ml"))
  fa_5 <- try(fa(sim, patternsp1, scores = "regression", rotate = "promax", fm = "ml"))
  
  # Choose the model with the LOWEST BIC
  if(all(class(fa_3) == "try-error")) {fa_3$BIC = NA}
  if(all(class(fa_4) == "try-error")) {fa_4$BIC = NA}
  if(all(class(fa_5) == "try-error")) {fa_5$BIC = NA}
  
  if (min(fa_5$BIC, fa_4$BIC, fa_3$BIC, na.rm = T) == fa_5$BIC & !is.na(fa_5$BIC)) {
    fa_out <- fa_5
    rank <- patternsm1
  } else if (min(fa_5$BIC, fa_4$BIC, fa_3$BIC, na.rm = T) == fa_4$BIC & !is.na(fa_4$BIC)) {
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

# Calculate BIC
get_bic <- function(sim, patterns, scores, loadings){
  bic = sum((sim - (scores%*%loadings))^2) + 
    (1/2)*(nrow(sim) + ncol(sim)) * patterns * log(nrow(sim) * ncol(sim))
  return(bic)
}

# NMF L2
get_nmf_l2 <- function (sim, patterns) {
  
  patternsm1 = ifelse(patterns == 1, 3, patterns - 1)
  patternsp1 = patterns+1
  
  set.seed(1988)
  nmf_3 <- nmf(sim, patternsm1, nrun = 100, method = "lee")
  nmf_4 <- nmf(sim, patterns,   nrun = 100, method = "lee")
  nmf_5 <- nmf(sim, patternsp1, nrun = 100, method = "lee")
  
  # Calculate BIC for each
  bic_3 <- get_bic(sim, patternsm1, basis(nmf_3), coef(nmf_3))
  bic_4 <- get_bic(sim, patterns,   basis(nmf_4), coef(nmf_4))
  bic_5 <- get_bic(sim, patternsp1, basis(nmf_5), coef(nmf_5))
  
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
get_nmf_p <- function (sim, patterns) {
  
  patternsm1 = ifelse(patterns == 1, 3, patterns - 1)
  patternsp1 = patterns+1
  
  set.seed(1988)
  nmf_3 <- nmf(sim, patternsm1, nrun = 100, method = "brunet")
  nmf_4 <- nmf(sim, patterns,   nrun = 100, method = "brunet")
  nmf_5 <- nmf(sim, patternsp1, nrun = 100, method = "brunet")
  
  # Calculate BIC for each
  bic_3 <- get_bic(sim, patternsm1, basis(nmf_3), coef(nmf_3))
  bic_4 <- get_bic(sim, patterns,   basis(nmf_4), coef(nmf_4))
  bic_5 <- get_bic(sim, patternsp1, basis(nmf_5), coef(nmf_5))
  
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

## Metric Functions
## to compare performance

# Symmetric Subspace Distance 
symm_subspace_dist <- function(U, V) {
  
  # patterns should be columns
  if (nrow(U) < ncol(U)) {U <- t(U)}
  if (nrow(V) < ncol(V)) {V <- t(V)}
  
  qrU <- qr.Q(qr(U))
  qrV <- qr.Q(qr(V))
  
  m <- ncol(U)
  n <- ncol(V)
  
  dUV <- sqrt( max(m,n) - sum((t(qrU) %*% qrV)^2) )
  
  ratio <- dUV/sqrt( max(m,n))
  
  ratio
  
}

# Cosine distance (matrix version)
cos_dist <- function(a, b){

  if(any(is.na(a)) | any(is.na(b))) {return(NA)} else{
    
    # patterns should be columns
    if (nrow(a) < ncol(a)) {a <- t(a)}
    if (nrow(b) < ncol(b)) {b <- t(b)}
  
    if(ncol(a) != ncol(b)) { return(NA) } else{
        a = as.vector(a)
        b = as.vector(b)
        cos_sim = (a %*% b) / (norm(a, "2")*norm(b, "2"))
        return(1-cos_sim) 
    }
  }
}

# Cosine distance (vector version)
cos_dist_v <- function(a, b){
  
  if(any(is.na(a)) | any(is.na(b))) {return(NA)} else{
    
       a = as.matrix(a)
       b = as.matrix(b)

       # patterns should be columns
       if (nrow(a) < ncol(a)) {a <- t(a)}
       if (nrow(b) < ncol(b)) {b <- t(b)}
       
       if(ncol(a) != ncol(b)) { return(NA) } else{
         cos_sim = c()
         for (i in 1:ncol(a)) {
          cos_sim[i] = (a[,i] %*% b[,i]) / (norm(a[,i], "2")*norm(b[,i], "2"))
          }
         return(1-cos_sim) 
       }
  }
  }

# L2 relative error
# WITH error handling
# if the matrices are not the same size, NA
get_relerror <- function(x,y) {

  if(any(is.na(x)) | any(is.na(y))) {return(NA)} else{
    
        x = as.matrix(x)
        y = as.matrix(y)
        
        # patterns should be columns
        if (nrow(x) < ncol(x)) {x <- t(x)}
        if (nrow(y) < ncol(y)) {y <- t(y)}  
    
        if(ncol(x) != ncol(y)) { return(NA) } else{
                return(norm(x-y, "F")/norm(x, "F"))
              }
  }
  }
  
# Factor correspondence
# WITH error handling
# if matrices are not same size, NA
# If nn is false, we allow *signed permutations*
get_perm <- function(x,y, nn = TRUE) {
  
  if(any(is.na(x)) | any(is.na(y))) {return(NA)} else{
    
      x = as.matrix(x)
      y = as.matrix(y)
      
      # patterns should be columns
      if (nrow(x) < ncol(x)) {x <- t(x)}
      if (nrow(y) < ncol(y)) {y <- t(y)}  
      
      if(ncol(x) != ncol(y)) { return(NA) } else{
        return(factor_correspondence(x, y, nn)$permutation_matrix) } 
      }
}

get_product <- function(x,y) {
  if(any(is.na(x)) | any(is.na(y))) {return(NA)} else{
    
    x = as.matrix(x)
    y = as.matrix(y)
    
    # patterns should be columns
    if (nrow(x) < ncol(x)) {x <- t(x)}
    
    if(ncol(x) != ncol(y)) { return(NA) } else{ return(x %*% y) }
  }
}













