dgp <- 
  dgp %>%
  mutate(true_patterns = map(true_patterns, t),
         nmf_l2_loadings  = map(nmf_l2_loadings, t),
         nmf_p_loadings  = map(nmf_p_loadings, t),
         eh = map(eh, t))

dgp_re <- dgp %>% 
  dplyr::select(seed, data, true_patterns, true_scores, pca_rotations, pca_scores,
                fa_rotations, fa_scores, nmfl2_loadings = nmf_l2_loadings, nmfl2_scores = nmf_l2_scores,
                nmfp_loadings = nmf_p_loadings, nmfp_scores = nmf_p_scores, bnmf_loadings = eh, bnmf_scores = ewa) %>% 
  pivot_longer(true_patterns:bnmf_scores,
               names_to = c("method", "results"),
               names_sep = "_") %>% 
  drop_na(value) %>% 
  mutate(results = ifelse(results == "rotations" | results == "patterns", "loadings", results)) %>% 
  group_by(seed, data, results) %>% 
  mutate(truth = value[method == "true"]) %>% 
  filter(method != "true") %>% 
  mutate(value_re = map2(truth, value, 
                        function(x,y) if(ncol(y) == 4) 
                        {factor_correspondence(as.matrix(x), as.matrix(y))$rearranged} 
                        else{NA} ))

save(dgp_re, file = "./HPC/Rout/dgp_re.RDA")
