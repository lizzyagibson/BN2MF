  library(tidyverse)
  library(R.matlab)
  library(openssl)
  library(plotly)
  
  # create random vector for matlab on hpc
  rand_vec = rand_num(1000) * ((2^31)-1)
   
  # Read data
  chem <- read_csv("./Results/Main/Corr Ex/corr_chem.csv")
  sim <- read_csv("./Results/Main/Corr Ex/corr_sim.csv")
  patterns <- read_csv("./Results/Main/Corr Ex/corr_patterns.csv")
  scores <- read_csv("./Results/Main/Corr Ex/corr_scores.csv")
  
  v_EWA <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_EWA.mat")[[1]]
  v_EH <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_EH.mat")[[1]]
  v_upper <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_WA_upper.mat")[[1]]
  v_lower <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_WA_lower.mat")[[1]]
  v_dist <- readMat("/Users/lizzy/BN2MF/Bootstrap/compare/cor_var_dist_WA.mat")[[1]]
  
  # Normalize truth
  patterns_denom      = apply(patterns, 1, sum)
  patterns_scaled     = patterns/patterns_denom
  patterns_denom_diag = diag(patterns_denom)
  scores_scaled       = as.matrix(scores) %*% patterns_denom_diag;
  
  # Read bootstrapped BN2MF EWA output
  bs_dist <- array(dim = c(1000, 4, 500))
  
  for (i in 1:500) {
    bs_dist[,,i] <- readMat(paste0("/Users/lizzy/BN2MF/Bootstrap/bootstrap_cor/bs_ewa_", i, ".mat"))[[1]]
    }
  
  head(bs_dist[,,14])
  
  # Create empirical confidence interval
  
  bs_EWA   <- apply(bs_dist, c(1,2), mean, na.rm = TRUE)
  bs_lower <- apply(bs_dist, c(1,2), quantile, 0.025, na.rm = TRUE)
  bs_upper <- apply(bs_dist, c(1,2), quantile, 0.975, na.rm = TRUE)
  
  head(bs_lower)
  head(bs_EWA)
  head(bs_upper)
  head(scores_scaled)
  head(v_EWA)
  
  ## Example viz
  truth = scores_scaled[1,1]
  
  v_ewa = v_EWA[1,1]
  v_dist_1 = v_dist[1,1,]
  v_q25 = v_lower[1,1]
  v_q75 = v_upper[1,1]
  
  bs_ewa = bs_EWA[1,1]
  bs_dist_1 = bs_dist[1,1,]
  bs_q25 = bs_lower[1,1]
  bs_q75 = bs_upper[1,1]
  
  histo <- tibble(Distribution = v_dist_1) %>% 
    mutate(Type = "Variational") %>% 
    rbind(., tibble(Distribution = bs_dist_1) %>% 
            mutate(Type = "Bootstrap")) %>% 
    ggplot(aes(x = Distribution)) +
    geom_rect(aes(xmin=bs_q25, xmax=bs_q75, ymin=0, ymax=Inf), fill="pink",      alpha=0.05) +
    geom_rect(aes(xmin=v_q25,  xmax=v_q75,  ymin=0, ymax=Inf), fill="lightblue", alpha=0.05) +
    geom_histogram(aes(y=..density.., group = Type, fill = Type), alpha = 0.5, bins = 100) +
    scale_fill_manual(values = c("red", "blue")) +
    theme_bw() + 
    geom_vline(xintercept = truth,            linetype="dotted", color = "black") +
    geom_vline(xintercept = c(v_ewa, bs_ewa), linetype="dotted", color = "yellow") + 
    xlim(0, 450) + ylab("Density")
  
  histo
  
  sum(scores_scaled <= v_upper & scores_scaled >= v_lower)/4000
  sum(scores_scaled <= bs_upper & scores_scaled >= bs_lower)/4000
  
  
  