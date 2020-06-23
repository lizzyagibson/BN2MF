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

save(out_dist, file = paste0("out_", job_num, "_dist.RDA"))
save(out_over, file = paste0("out_", job_num, "_over.RDA"))
save(out_cor, file = paste0("out_", job_num, "_cor.RDA"))


















