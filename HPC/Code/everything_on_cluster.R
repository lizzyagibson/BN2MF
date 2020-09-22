#########################################################
# Sims, npBNMF, regular NMF (L2 & Poisson), PCA, and FA #
#########################################################
# Lizzy Gibson    # 6/22/2020 ###########################
# Added BNMF      # 9/17/2020 ###########################
# Cleaned up code # 9/22/202  ###########################
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
load("/ifs/scratch/msph/ehs/eag2186/Data/sim_norm.RDA")

## read job number from system environment
## This only works if run on cluster!!
job_num = as.integer(Sys.getenv("SGE_TASK_ID"))
job_num

# load("/ifs/scratch/msph/ehs/eag2186/npbnmf/to_do.RDA")
# if (!(job_num %in% to_do)) {stop("Job id already ran.")}

# Run everything

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

get_pca_uncenter <- function (sim) {
  # Run PCA not centered, not scaled
  pca_out <- prcomp(sim, center = FALSE)
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

output_all <- sim_norm[job_num,] %>% 
  mutate(pca_out          = map(sim, get_pca),
         pca_rotations_un = map(pca_out, function(x) x[[1]]),
         pca_scores_un    = map(pca_out, function(x) x[[2]]),
         pca_pred         = map(pca_out, function(x) x[[3]]),
         pca_rank         = map(pca_out, function(x) x[[4]]))

output_all <- output_all %>% 
  mutate(pca_uncenter_out          = map(sim, get_pca_uncenter),
         pca_uncenter_rotations_un = map(pca_out, function(x) x[[1]]),
         pca_uncenter_scores_un    = map(pca_out, function(x) x[[2]]),
         pca_uncenter_pred         = map(pca_out, function(x) x[[3]]),
         pca_uncenter_rank         = map(pca_out, function(x) x[[4]]))

## Match Functions
## Reorder eigenvector matrix to best match truth
match_eigen <- function (truth, loading) {
  truth <- t(truth) # want patterns as columns
  loading <- as.matrix(loading) # patterns already loadings
  corr <- cor(truth, loading)
  
  # If highest correlated eigenvector is negative, flip the whole matrix
  if (corr[which.max(abs(corr))] < 0) {
    new_loading <- -loading
    corr <- cor(truth, new_loading)
    } else {new_loading <- loading}
  
  if (ncol(corr) < nrow(corr)) {corr <- corr[1:ncol(corr),]}
  
  reordered <- matrix(nrow = nrow(new_loading), ncol = ncol(new_loading))
  
  # Put eigenvector columns in order that best correlates with truth
  # Put column 1 in column it is most correlated with
  reordered[,which.max(corr[,1])] <- new_loading[,1]
  
  # Set corresponding ROW to NA so it can't be chosen again
  corr[which.max(corr[,1]),] <- NA
  
  # Put column 2 in column it is most correlated with
  if (ncol(reordered) > 1) {
    reordered[,which.max(corr[,2])] <- new_loading[,2]
    corr[which.max(corr[,2]),] <- NA
    }
  
  if (ncol(reordered) > 2) {
    reordered[,which.max(corr[,3])] <- new_loading[,3]
    corr[which.max(corr[,3]),] <- NA
  }
  
  if (ncol(reordered) > 3) {
    reordered[,which.max(corr[,4])] <- new_loading[,4]
    corr[which.max(corr[,4]),] <- NA
  }
  
  if (ncol(reordered) > 4) {
    reordered[,5:ncol(new_loading)] <- new_loading[,5:ncol(new_loading)]
  }
  
  reordered
}

match_pca_scores <- function (rotations_un, rotations, scores_un) {
  
  scores <- as_tibble(matrix(nrow = nrow(scores_un), ncol = ncol(scores_un)))
  corr <- cor(rotations, rotations_un)
  
  # Flip scores to match loading matrix
  if (corr[which(abs(corr) == 1)][1] < 0) {scores_un = -scores_un}
  
  # Match score vectors with loading vectors
  for (k in 1:ncol(rotations_un)) {
    for (j in 1:ncol(rotations)) {
            if (all(abs(rotations[, j]) == abs(rotations_un[, k]))) {
              scores[, j] = scores_un[, k]
            }
    }
    }
  scores
}

## Run

output_all <- output_all %>% 
  mutate(pca_rotations = map2(true_patterns, pca_rotations_un, match_eigen),
         pca_scores    = pmap(list(pca_rotations_un, pca_rotations, pca_scores_un), match_pca_scores))
  
output_all <- output_all %>% 
  mutate(pca_uncenter_rotations = map2(true_patterns, pca_uncenter_rotations_un, match_eigen),
         pca_uncenter_scores    = pmap(list(pca_uncenter_rotations_un, pca_uncenter_rotations, pca_uncenter_scores_un), 
                                      match_pca_scores))

#####
# Factor Analysis
#####

## Function

sim = output_all$sim[[1]]

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

output_all <- output_all %>% 
  mutate(fa_out       = map(sim, get_fa),
         fa_rotations = map(fa_out, function(x) x[[1]]),
         fa_scores    = map(fa_out, function(x) x[[2]]),
         fa_pred      = map(fa_out, function(x) x[[3]]),
         fa_rank      = map(fa_out, function(x) x[[4]]))

# FA factors are in the same order every time

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

output_all <- output_all %>%
  mutate(nmf_l2_out         = map(sim, get_nmf_l2),
         nmf_l2_loadings_un = map(nmf_l2_out, function(x) x[[1]]),
         nmf_l2_scores_un   = map(nmf_l2_out, function(x) x[[2]]),
         nmf_l2_pred        = map(nmf_l2_out, function(x) x[[3]]),
         nmf_l2_rank        = map(nmf_l2_out, function(x) x[[4]]))

## Match Functions

match_nmf_patterns <- function (truth, loading) {
  truth <- t(truth) # want patterns as columns
  
  # Flip loadings to same dimensions as Truth
  if (nrow(truth) != nrow(loading)) {loading <- t(as.matrix(loading))} else {loading <- as.matrix(loading)}
  
  corr <- cor(truth, loading)
  if (ncol(corr) < nrow(corr)) {corr <- corr[1:ncol(corr),]}
  
  if (length(corr) != 1) {
      reordered <- matrix(nrow = nrow(loading), ncol = ncol(loading))
  
      # Put column 1 in column it is most correlated with
      reordered[,which.max(corr[,1])] <- loading[,1]
  
      # Set corresponding ROW to NA so it can't be chosen again
      corr[which.max(corr[,1]),] <- NA
   
      # Put column 2 in column it is most correlated with
      if (ncol(reordered) > 1) {
         reordered[,which.max(corr[,2])] <- loading[,2]
         corr[which.max(corr[,2]),] <- NA
       }
       
      if (ncol(reordered) > 2) {
        reordered[,which.max(corr[,3])] <- loading[,3]
        corr[which.max(corr[,3]),] <- NA
        }
      
      if (ncol(reordered) > 3) {
        reordered[,which.max(corr[,4])] <- loading[,4]
        corr[which.max(corr[,4]),] <- NA
      }
      
      if (ncol(reordered) > 4) {reordered[,5:ncol(reordered)] <- loading[,5:ncol(reordered)]}
  } else {reordered = loading}
      
      reordered
}

match_nmf_scores <- function (rotations_un, rotations, scores_un) {
  
  if (nrow(scores_un) != 1000) {scores_un = t(scores_un)} # want patterns as columns
  scores <- as_tibble(matrix(nrow = nrow(scores_un), ncol = ncol(scores_un)))
  
  for (k in 1:nrow(rotations_un)) {
    for (j in 1:ncol(rotations)) {
        if (all(rotations_un[k, ] == rotations[, j])) {
          scores[, j] = scores_un[, k]
        }}
      }
  scores
}

## Run

output_all <- output_all %>%
  mutate(nmf_l2_loadings = map2(true_patterns, nmf_l2_loadings_un, match_nmf_patterns),
         nmf_l2_scores   = pmap(list(nmf_l2_loadings_un, nmf_l2_loadings, nmf_l2_scores_un), 
                              match_nmf_scores))

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

output_all <- output_all %>% 
  mutate(nmf_p_out         = map(sim, get_nmf_p),
         nmf_p_loadings_un = map(nmf_p_out, function(x) x[[1]]),
         nmf_p_scores_un = map(nmf_p_out, function(x) x[[2]]),
         nmf_p_pred = map(nmf_p_out, function(x) x[[3]]),
         nmf_p_rank = map(nmf_p_out, function(x) x[[4]]))

output_all <- output_all %>%
  mutate(nmf_p_loadings = map2(true_patterns, nmf_p_loadings_un, match_nmf_patterns),
         nmf_p_scores = pmap(list(nmf_p_loadings_un, nmf_p_loadings, nmf_p_scores_un), 
                              match_nmf_scores))

#####
# np BNMF
#####

##### Read function
NPBayesNMF <- function(X) {
  X = as.matrix(X)
  bnp_switch = 1
  dim = nrow(X)
  N = ncol(X)
  Kinit = ncol(X)
  
  nruns = 100
  end_score = matrix(rep(0, times = nruns))
  
  EA = matrix()
  EWA = matrix()
  EH = matrix()
  EW = matrix()
  
  for (i in 1:nruns) {
    
    K = Kinit
    
    w01 = 1
    w02 = 1
    
    W1 = matrix(rgamma(dim*Kinit, shape = dim, scale = 1/dim), nrow = dim, ncol = Kinit)
    W2 = dim*matrix(1, nrow = dim, ncol = Kinit)
    
    a01 = bnp_switch*1/Kinit + (1-bnp_switch)
    a02 = 1
    
    A1 = a01 + bnp_switch*1000*matrix(1, 1, Kinit)/Kinit
    A2 = a02 + bnp_switch*1000*matrix(1, 1, Kinit)
    
    h01 = 1 # Non-sparse prior
    # h01 = 1/Kinit # Sparse prior
    h02 = 1
    
    H1 = matrix(1, Kinit, N)
    H2 = matrix(1, Kinit, N)
    
    num_iter = 100000
    
    score = vector("numeric", length = num_iter)
    
    for (iter in 1:num_iter) {
      
      EW = W1 / W2
      
      # The equivalent of Matlab's repmat(a,2,3) in base R is kronecker(matrix(1,2,3),a)
      X_reshape = kronecker(array(1, dim = c(K,1,1)), array(t(X), c(1, N, dim)))
      ElnWA = digamma(W1) - log(W2) + kronecker(array(1, c(dim,1)), digamma(A1)-log(A2))
      ElnWA_reshape = kronecker(array(1, dim = c(1,N,1)), array(t(ElnWA), c(K, 1, dim)))
      
      t1 = array(apply(ElnWA_reshape,c(2,3), max), c(1, Kinit, dim))
      
      # % Expected value of log W * expected value of log A
      ElnWA_reshape = ElnWA_reshape - kronecker(array(1, dim = c(K,1,1)), t1)
      # % expected value of log H
      ElnH = digamma(H1) - log(H2)                               
      
      P = ElnWA_reshape + kronecker(array(1, dim=c(1,1,dim)), ElnH)
      P = exp(P)
      P = P / kronecker(array(1, dim = c(K,1,1)), array(apply(P, 3, colSums), dim = c(1, Kinit, dim)))
      
      # % P is a probability to put a lower bound on the ELBO
      # % expected value of log WAH is concave
      # % include P to make expectation tractable
      # % These are update steps from optimizing the ELBO 
      # % take the gradient with respect to parameter, set to zero, solve
      
      H1 = h01 + apply(P * X_reshape, c(1,2), sum)
      H2 = h02 + t(kronecker(matrix(1, nrow = N, ncol = 1), 
                             matrix(colSums(EW * kronecker(matrix(1, nrow = dim, ncol = 1), A1/A2)), nrow=1)))
      
      W1 = w01 + t(array(apply(X_reshape*P, c(1,3), sum), dim = c(K, dim)))
      W2 = w02 + kronecker(matrix(1, dim, 1), t(rowSums((H1/H2) * kronecker(matrix(1, 1, N), t(A1/A2)))))
      
      A1 = a01 + bnp_switch * t(rowSums(apply(X_reshape*P, c(1,2), sum)))
      A2 = a02 + bnp_switch * (colSums(W1/W2) * t(rowSums(H1/H2)))
      
      # % This is the sparse prior on A, pushing A to zero
      # If all patterns are zero, don't push anything
      # prohibits pushing all patterns to zero
      if ( all(A1/A2 < 10^-3) ) {idx_prune = integer(0)} else {idx_prune = which(A1/A2 < 10^-3)} 
      
      if (length(idx_prune) >= 1) {
        W1 = matrix(W1[,-idx_prune], nrow = dim)
        W2 = matrix(W2[,-idx_prune], nrow = dim)
        A1 = matrix(A1[,-idx_prune], nrow = 1)
        A2 = matrix(A2[,-idx_prune], nrow = 1)
        H1 = matrix(H1[-idx_prune,], ncol = N)
        H2 = matrix(H2[-idx_prune,], ncol = N)
      }
      
      K = length(A1)
      
      score[iter] = if (ncol(A1/A2) > 1) {
        sum( abs( X - (W1/W2) %*% diag(as.vector(A1/A2)) %*% (H1/H2) ) )
      } else if (ncol(A1/A2) == 1) {
        sum( abs( X - (W1/W2) %*% as.vector(A1/A2) %*% (H1/H2) ) )
      }
      
      if (iter %% 100 == 0) {print(paste0("Run Number: ", i, "; Iter Number: ", iter, "; Iter Score: ", round(score[iter], 4)))}
      if (iter > 1 && abs(score[iter-1] - score[iter]) < 1e-5) {
        print(paste0('Converged in ', iter,' iterations.'))
        break} # Convergence criteria!
    }
    
    if( iter == num_iter ) warning('Maximum iterations reached. BNMF did not converge.')
    end_score[i] = score[tail(which(score != 0),1)]
    
    #print(paste0("Run Number: ", i, "; Iter Number: ", iter, "; Final Score: ", round(end_score[i], 4)))
    
    # % Among the results, use the fitted variational parameters that achieve the highest ELBO
    if (i == 1 | (i > 1 && (end_score[i] >= max(end_score)))) {
      EA = A1/A2
      EWA = (W1/W2)*diag(A1/A2)
      EH = H1/H2
      EW = (W1/W2)
      varA = A1/(A2^2)
      varW = W1/(W2^2)
      alphaH = H1
      betaH = H2
    }
    H_CI_low  <- qgamma(0.025, shape = alphaH, rate = betaH)
    H_CI_high <- qgamma(0.975, shape = alphaH, rate = betaH)
  }
  
  list(EWA = EWA, EH = EH, H_CI_low = H_CI_low, H_CI_high = H_CI_high)
}

get_bnmf <- function (sim) {
  set.seed(1988)
  bnmf_out <- NPBayesNMF(sim)
  
  ewa <- bnmf_out$EWA
  eh <- bnmf_out$EH
  pred <- ewa %*% eh
  rank = ncol(ewa)
  
  list(ewa = ewa, eh = eh, pred = pred, rank = rank)
}

output_all <- output_all %>% 
              mutate(bnmf_out         = map(sim, get_bnmf),
                     bnmf_scores_un   = map(bnmf_out, function(x) x[[1]]),
                     bnmf_loadings_un = map(bnmf_out, function(x) x[[2]]),
                     bnmf_pred        = map(bnmf_out, function(x) x[[3]]),
                     bnmf_rank        = map(bnmf_out, function(x) x[[4]]))
            
output_all <- output_all %>%
  mutate(bnmf_loadings = map2(true_patterns, bnmf_loadings_un, match_nmf_patterns),
         bnmf_scores   = pmap(list(bnmf_loadings_un, bnmf_loadings, bnmf_scores_un), 
                             match_nmf_scores))


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

output_all <- output_all %>% 
  mutate(pca_norm              = map2(sim, pca_pred,    function(x,y) norm(x-y, "F")/norm(x, "F")),
         fa_norm               = map2(sim, fa_pred,     function(x,y) norm(x-y, "F")/norm(x, "F")),
         nmf_l2_norm           = map2(sim, nmf_l2_pred, function(x,y) norm(x-y, "F")/norm(x, "F")),
         nmf_p_norm            = map2(sim, nmf_p_pred,  function(x,y) norm(x-y, "F")/norm(x, "F")),
         bnmf_norm             = map2(sim, bnmf_pred,   function(x,y) norm(x-y, "F")/norm(x, "F")),
         pca_rotation_ssdist   = map2(true_patterns, pca_rotations,   symm_subspace_dist),
         pca_scores_ssdist     = map2(true_scores,   pca_scores,      symm_subspace_dist),
         fa_rotations_ssdist   = map2(true_patterns, fa_rotations,    symm_subspace_dist),
         fa_scores_ssdist      = map2(true_scores,   fa_scores,       symm_subspace_dist),
         nmf_l2_loading_ssdist = map2(true_patterns, nmf_l2_loadings, symm_subspace_dist),
         nmf_l2_scores_ssdist  = map2(true_scores,   nmf_l2_scores,   symm_subspace_dist),
         nmf_p_loading_ssdist  = map2(true_patterns, nmf_p_loadings,  symm_subspace_dist),
         nmf_p_scores_ssdist   = map2(true_scores,   nmf_p_scores,    symm_subspace_dist),
         bnmf_loading_ssdist   = map2(true_patterns, bnmf_loadings,   symm_subspace_dist),
         bnmf_scores_ssdist    = map2(true_scores,   bnmf_scores,     symm_subspace_dist))

#####
# Save
#####
                           
output_all <- output_all %>% dplyr::select(-grep("_un|_out", colnames(.)))

save(output_all, file = paste0("/ifs/scratch/msph/ehs/eag2186/npbnmf/norm_out/norm_sims", job_num, ".RDA"))

