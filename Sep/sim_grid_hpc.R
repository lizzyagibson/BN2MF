
# Packages
library(tidyverse, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")
library(LearnBayes, lib.loc = "/ifs/home/msph/ehs/eag2186/local/hpc/")

## 1: Distinct Patterns
### Simulate Patterns
create_patterns <- function (seed, reps) {
  set.seed(seed)
  pat1=c(0,0,5,10)
  pat2=c(0,5,10,0)
  pat3=c(5,10,0,0)
  pat4=c(10,0,0,5)
  
  start = c()
  
  if (reps != 0) { 
    for (i in 1:reps) {
      start = cbind(start,diag(1,4))
    }}
  
  patterns = cbind(start, t(rdirichlet(10-reps, pat1)), t(rdirichlet(10-reps, pat2)),
                   t(rdirichlet(10-reps, pat3)), t(rdirichlet(10-reps, pat4)))
  return(patterns)
}

seed = 1:100
sep_num = 0:10
noise_level = seq(0,1, 0.1)
all_comb = expand_grid(seed, sep_num, noise_level)

patterns_iter <- all_comb %>% 
  mutate(true_patterns = map2(seed, sep_num, create_patterns))

### Simulate Scores
create_scores <- function (seed) {
  set.seed(seed)
  scores <- matrix(exp(rnorm(1000*4)), ncol = 4)
  return(as.matrix(scores))
}

scores_iter <- patterns_iter %>% 
  mutate(true_scores = map(seed, create_scores))

### Simulate Chemical Exposures
sim_sep <- scores_iter %>% 
  mutate(chem = map2(true_scores, true_patterns, `%*%`))

### Simulate Noise
add_noise <- function (seed, chem, noise_level) {
  n = nrow(chem)
  p = ncol(chem)
  noise <- matrix(NA, nrow = n, ncol = p)
  stdev = apply(chem, 2, sd)
  
  for (i in 1:p) {
    noise[,i] <- (rnorm(n, mean = 0, sd = (stdev[i]*noise_level)))
  }
  
  sim = pmax(chem + noise, 0)
  sim
}

sim_sep <- sim_sep %>% 
  mutate(sim = pmap(list(seed, chem, noise_level), add_noise))

## Save
sim_sep

# for (i in 1:nrow(sim_sep)) {
#   write_csv(as_tibble(sim_sep$sim[[i]]),           paste0("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv/sim_sep_", i, ".csv"))
#   write_csv(as_tibble(sim_sep$chem[[i]]),          paste0("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv/chem_sep_", i, ".csv"))
#   write_csv(as_tibble(sim_sep$true_patterns[[i]]), paste0("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv/patterns_sep_", i, ".csv"))
#   write_csv(as_tibble(sim_sep$true_scores[[i]]),   paste0("/ifs/scratch/msph/ehs/eag2186/Data/sep_csv/scores_sep_", i, ".csv"))
# }
# 
# save(sim_sep, file = "/ifs/scratch/msph/ehs/eag2186/Data/sim_sep_grid.RDA")
save(sim_sep, file = "./Sims/sim_full.RDA")

# Bootstrap subsample ####
bs_sample = sim_sep %>% 
              mutate(id = 1:nrow(.)) %>% 
              filter(sep_num %in% c(0,10) & noise_level %in% c(0.2, 0.5, 1))

bs_ids = as_tibble(bs_sample$id)
#write_csv(bs_ids, "./Sep/bs_ids.csv")
#save(bs_sample, file = "./Sep/bs_sample.RDA")

boots = 1:150
expand_grid(bs_ids, boots)
