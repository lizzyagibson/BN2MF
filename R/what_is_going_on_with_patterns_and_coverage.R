## Compare single examples of overlapping sims

add_noise_sd1 <- function (seed, chemx) {
  n = nrow(chemx)
  p = ncol(chemx)
  noise <- matrix(NA, nrow = n, ncol = p)
  
  for (i in 1:p) {
    noise[,i] <- (rnorm(n, mean = 0, sd = 1))
  }
  sim = pmax(chemx + noise, 0)
  sim
}

add_noise_rel <- function (seed, chem, noise_level) {
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

load("./Sims/sim_old_unsep.RDA")
sim_dgp
load("./Sims/sim_sep_grid.RDA")
out.sep = sim_sep %>% 
          filter(sep_num == 0 & seed == 1) %>% 
          filter(noise_level == "0.6")
metrics %>% 
  filter(sep_num == 0 & seed == 1) %>% 
  filter(noise_level == "0.6") %>% pull(prop)

new_patterns = out.sep$true_patterns[[1]]

pat1=c(0,0,5,10)
pat2=c(0,5,10,0)
pat3=c(5,10,0,0)
pat4=c(10,0,0,5)
new_patterns50 = cbind(new_patterns, t(rdirichlet(3, pat1)), t(rdirichlet(3, pat2)),
                 t(rdirichlet(2, pat3)), t(rdirichlet(2, pat4)))

new_scores   = out.sep$true_scores[[1]]
old_patterns = sim_dgp$true_patterns[[101]]
chem_old_true = sim_dgp$chem[[101]]
old_patterns_std = apply(old_patterns, 2, function(x) x/sum(x))
old_patterns_less = old_patterns[,c(1:10, 13:22, 25:34,38:47)]

old_scores = sim_dgp$true_scores[[101]]

chem_new     = out.sep$chem[[1]]
sim_new      = out.sep$sim[[1]]

chem_old = old_scores%*%new_patterns50
sim_old = add_noise_sd1(1, chem_old)

old_patterns %*% t(old_patterns)
old_patterns_std %*% t(old_patterns_std)
new_patterns %*% t(new_patterns)

sort(apply(old_patterns, 2, sd))
sort(apply(chem_old_true, 2, sd))

sort(apply(new_patterns50, 2, sd))
sort(apply(chem_new, 2, sd))

write_csv(as_tibble(sim_new),      "./Sims/Test/test_sim_101.csv")
write_csv(as_tibble(chem_new),     "./Sims/Test/test_chem_101.csv")
write_csv(as_tibble(new_patterns), "./Sims/Test/test_patterns_101.csv")
write_csv(as_tibble(new_scores),   "./Sims/Test/test_scores_101.csv")

mean(old_patterns_std %*% t(old_patterns_std))
mean(new_patterns50 %*% t(new_patterns50))

sort(apply(old_patterns, 2, var))
sort(apply(old_patterns_std, 2, var))
sort(apply(new_patterns50, 2, sd))

sort(apply(chem_old_std, 2, var))
0.7*sort(apply(chem_new_50, 2, var))

old_patterns_std
new_patterns50

norm(old_patterns_std, "F")
norm(new_patterns50, "F")

chem_old_std = old_scores%*%old_patterns_std
sim_old_std = add_noise_sd1(1, chem_old_std)

chem_new_50 = new_scores%*%new_patterns50
sim_new_50 = add_noise_rel(1, chem_new_50, 0.7)

norm(sim_old_std, "F")
norm(sim_new_50, "F")
