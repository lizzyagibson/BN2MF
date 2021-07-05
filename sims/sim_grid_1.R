# SCRIPT 1
# Create simulated datasets

# Get functions
source("./functions/compare_functions_2.R")

# Simulate Patterns
create_patterns <- function (seed, reps) {
  set.seed(seed)
  
  # 4 patterns
  # each chemical loads high on one, medium on another
  # zero on lat two patterns
  pat1=c(0,0,5,10)
  pat2=c(0,5,10,0)
  pat3=c(5,10,0,0)
  pat4=c(10,0,0,5)
  
  # `inseparable` is how overlapping the patterns are
  # inseparable = 0 -- completely overlapping
  # inseparable = 10 completely distinct
  separable = c()
  if (inseparable != 0) { 
    for (i in 1:inseparable) {
      separable = cbind(separable,diag(1,4))
    }}
  
  patterns = cbind(separable, t(rdirichlet(10-inseparable, pat1)), t(rdirichlet(10-inseparable, pat2)),
                   t(rdirichlet(10-inseparable, pat3)), t(rdirichlet(10-inseparable, pat4)))
  return(patterns)
}

# 100 random samples from each data generating process
seed = 1:100
# gradient from overlapping to distinct
sep_num = 0:10
# gradient from no noise to high noise
noise_level = seq(0,1, 0.1)
# get every combination
all_comb = expand_grid(seed, sep_num, noise_level)

# simulate patterns
patterns_iter <- all_comb %>% 
  mutate(true_patterns = map2(seed, sep_num, create_patterns))

# Simulate Scores
create_scores <- function (seed) {
  set.seed(seed)
  # independent draws from standard log-normal
  scores <- matrix(exp(rnorm(1000*4)), ncol = 4)
  return(as.matrix(scores))
}

scores_iter <- patterns_iter %>% 
  mutate(true_scores = map(seed, create_scores))

# Simulate Chemical Exposures
# matrix multiply scores times loadings
sim_sep <- scores_iter %>% 
  mutate(chem = map2(true_scores, true_patterns, `%*%`))

# Simulate Noise
add_noise <- function (seed, chem, noise_level) {
  n = nrow(chem)
  p = ncol(chem)
  noise <- matrix(NA, nrow = n, ncol = p)
  
  # stdev of each column in the true simulation before noise was added
  stdev = apply(chem, 2, sd)
  
  # add noise from normal dist, mean = 0
  # sd = proprtion of true noise
  for (i in 1:p) {
    noise[,i] <- (rnorm(n, mean = 0, sd = (stdev[i]*noise_level)))
  }
  
  # if negative, push to zero
  sim = pmax(chem + noise, 0)
  sim
}

# add noise
sim_sep <- sim_sep %>% 
  mutate(sim = pmap(list(seed, chem, noise_level), add_noise))
# ^ this file is 5G so I did not push it to github
# save it to run other scripts

# save nested dataframe 
# save(sim_sep, file = "./sims/sim_sep.RDA")

sim_ids = sim_sep %>% dplyr::select(1:3) %>% mutate_all(as.factor)

#save csv files
# for (i in 1:nrow(sim_sep)) {
#   write_csv(as_tibble(sim_sep$sim[[i]]),           paste0("./sims/csvs/sim_sep_", i, ".csv"))
#   write_csv(as_tibble(sim_sep$chem[[i]]),          paste0("./sims/csvs/chem_sep_", i, ".csv"))
#   write_csv(as_tibble(sim_sep$true_patterns[[i]]), paste0("./sims/csvs/patterns_sep_", i, ".csv"))
#   write_csv(as_tibble(sim_sep$true_scores[[i]]),   paste0("./sims/csvs/scores_sep_", i, ".csv"))
# }

# save ids
# save(sim_ids, file = "./sims/sim_ids.RDA")

# Bootstrap subsample
# distinct and overlapping
# noise = 0.2, 0.5, and 1
bs_sample = sim_sep %>% 
              mutate(id = 1:nrow(.)) %>% 
              filter(sep_num %in% c(0,10) & noise_level %in% c(0.2, 0.5, 1))
# ^ this file is also big so I did not push it to gethub
# save it to run other scripts
# save(bs_sample, file = "./sims/bs_sample.RDA")

# save bootstrap ids
bs_ids = as_tibble(bs_sample$id)
# write_csv(bs_ids, "./sims/bs_ids.csv")


