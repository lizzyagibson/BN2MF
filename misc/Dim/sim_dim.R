## Get functions
source("/ifs/scratch/msph/ehs/eag2186/npbnmf/compare_functions.R")

# Create simulated datasets

patterns <- c(1, 4, 10)
participants <- c(200, 1000, 10000)
chemicals <- c(20, 40, 100)
seed <- c(1:100)
overlap = c(0,1)

loop_scores   <- expand_grid(seed, patterns, participants, chemicals, overlap) %>% filter(!(patterns == 1 & overlap == 1))
loop_patterns <- expand_grid(seed, participants, chemicals)

# Simulate Patterns ####

## Construct a matrix for 1,4,10 TRUE exposure patterns
create_1pattern <- function (seed, chemicals) { # nothing to overlapp
  set.seed(seed)
  patterns <- matrix(c(rep(1, chemicals/2), rep(0, chemicals/2)), nrow = 1)
  patterns
}

create_4patterns <- function (seed, chemicals, overlap) { # overlapping or distinct
  set.seed(seed)

  group <- chemicals / 4
  
  pat1=c(0,0,5,10)
  pat2=c(0,5,10,0)
  pat3=c(5,10,0,0)
  pat4=c(10,0,0,5)

  if (overlap == 1) {
    patterns = cbind(t(rdirichlet(group, pat1)), t(rdirichlet(group, pat2)),
                     t(rdirichlet(group, pat3)), t(rdirichlet(group, pat4))) # completely overlapping
  } else if (overlap == 0) {
    patterns = rbind( c(rep(1, group), rep(0, chemicals-group)),
             c(rep(0, group), rep(1, group), rep(0, 2*group)),
             c(rep(0, 2*group), rep(1, group), rep(0, group)),
             c(rep(0, chemicals-group), rep(1, group))) # completely distinct
  }
  return(patterns)
}

create_10patterns <- function (seed, chemicals, overlap) {
  set.seed(seed)
  group <- chemicals / 10

  pat4 =c(10,0,0,0,0,0,0,0,0,5)
  pat3 =c(5,10,0,0,0,0,0,0,0,0)
  pat2 =c(0,5,10,0,0,0,0,0,0,0)
  pat1 =c(0,0,5,10,0,0,0,0,0,0)
  pat5 =c(0,0,0,5,10,0,0,0,0,0)
  pat6 =c(0,0,0,0,5,10,0,0,0,0)
  pat7 =c(0,0,0,0,0,5,10,0,0,0)
  pat8 =c(0,0,0,0,0,0,5,10,0,0)
  pat9 =c(0,0,0,0,0,0,0,5,10,0)
  pat10=c(0,0,0,0,0,0,0,0,5,10)

  if (overlap == 1) {
    patterns = cbind(t(rdirichlet(group, pat1)), t(rdirichlet(group, pat2)),
                     t(rdirichlet(group, pat3)), t(rdirichlet(group, pat4)),
                     t(rdirichlet(group, pat5)), t(rdirichlet(group, pat6)),
                     t(rdirichlet(group, pat7)), t(rdirichlet(group, pat8)),
                     t(rdirichlet(group, pat9)), t(rdirichlet(group, pat10))) # completely overlapping
  } else if (overlap == 0) {
    patterns = rbind( c(rep(1, group), rep(0, 9*group)),
                      c(rep(0, group), rep(1, group), rep(0, 8*group)),
                      c(rep(0, 2*group), rep(1, group), rep(0, 7*group)),
                      c(rep(0, 3*group), rep(1, group), rep(0, 6*group)),
                      c(rep(0, 4*group), rep(1, group), rep(0, 5*group)),
                      c(rep(0, 5*group), rep(1, group), rep(0, 4*group)),
                      c(rep(0, 6*group), rep(1, group), rep(0, 3*group)),
                      c(rep(0, 7*group), rep(1, group), rep(0, 2*group)),
                      c(rep(0, 8*group), rep(1, group), rep(0, group)),
                      c(rep(0, 9*group), rep(1, group))) # completely distinct
  }
  return(patterns)
}

patterns_iter <- loop_patterns %>% 
  mutate(four_1 = map2(seed, chemicals, create_4patterns, 1),
         four_0 = map2(seed, chemicals, create_4patterns, 0),
         one_0   = map2(seed, chemicals, create_1pattern),
         ten_1  = map2(seed, chemicals, create_10patterns, 1),
         ten_0  = map2(seed, chemicals, create_10patterns, 0)) %>% 
  pivot_longer(c(four_1:ten_0),
               values_to = "true_patterns",
               names_to = c("patterns", "overlap"),
               names_sep = "_") %>% 
  mutate(patterns = case_when(patterns == "four" ~ 4,
                              patterns == "one" ~ 1,
                              patterns == "ten" ~ 10),
         overlap = as.numeric(overlap))

# Simulate Scores ####
create_scores <- function (seed, participants, patterns) {
  set.seed(seed)
  # independent draws from standard log-normal
  scores <- matrix(exp(rnorm(participants*patterns)), ncol = patterns)
  return(as.matrix(scores))
}

scores_iter <- loop_scores %>% 
  mutate(true_scores = pmap(list(seed, participants, patterns), create_scores))

# Simulate Chemical Exposures
# matrix multiply scores times loadings
sim_dim <- left_join(patterns_iter, scores_iter) %>% 
  mutate(chem = map2(true_scores, true_patterns, `%*%`))

# Simulate Noise ####
add_noise <- function (chem, noise_level) {
  n = nrow(chem)
  p = ncol(chem)
  noise <- matrix(NA, nrow = n, ncol = p)
  
  # stdev of each column in the true simulation before noise was added
  stdev = apply(chem, 2, sd)
  
  if(length(stdev) != sum(stdev != 0)) {stdev = rep(stdev[stdev != 0], 2)}
  
  # add noise from normal dist, mean = 0
  # sd = proprtion of true noise
  
  # for 1 pattern, some rows = all zeros, sd = 0
  # replace
  for (i in 1:p) {
    noise[,i] <- (rnorm(n, mean = 0, sd = (stdev[i]*noise_level)))
  }
  
  # if negative, push to zero
  sim = pmax(chem + noise, 0)
  sim
}

# Add Noise ####
sim_dim <- sim_dim %>% 
  mutate(sim_0.2 = map(chem, add_noise, 0.2),
         sim_0.5 = map(chem, add_noise, 0.5),
         sim_1   = map(chem, add_noise, 1)) %>% 
  pivot_longer(sim_0.2:sim_1,
               values_to = "sim",
               names_to  = "noise") %>% 
  mutate(noise = as.numeric(str_remove(noise, "sim_")))

# save nested dataframe 
# save(sim_dim, file = "/ifs/scratch/msph/ehs/eag2186/Data/sim_dim.RDA")

# save ids, too
# dim_ids = sim_dim %>% 
#          dplyr::select(seed, participants, chemicals, patterns, overlap, noise) %>% 
#          mutate(job_num = 1:nrow(.))  
# save(dim_ids, file = "./misc/dim/dim_ids.RDA")

# dim_10 = dim_ids %>% filter(patterns == 10) %>% pull(job_num)
# write.csv(dim_10, file = "./misc/dim/dim_10.csv", row.names = F)

#save csv files
for (i in 1:nrow(sim_dim)) {
  if (sim_dim$patterns[[i]] == 10) {
  write_csv(as_tibble(sim_dim$sim[[i]]),           paste0("/ifs/scratch/msph/ehs/eag2186/Data/sim_dim/sim_dim_", i, ".csv"))
  write_csv(as_tibble(sim_dim$chem[[i]]),          paste0("/ifs/scratch/msph/ehs/eag2186/Data/chem_dim/chem_dim_", i, ".csv"))
  write_csv(as_tibble(sim_dim$true_patterns[[i]]), paste0("/ifs/scratch/msph/ehs/eag2186/Data/patterns_dim/patterns_dim_", i, ".csv"))
  write_csv(as_tibble(sim_dim$true_scores[[i]]),   paste0("/ifs/scratch/msph/ehs/eag2186/Data/scores_dim/scores_dim_", i, ".csv"))
  }
}

