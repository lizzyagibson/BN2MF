normed <- c()

for (i in 1:200) {
  eh <- readMat(here::here(paste0("./MATLAB/dgp_hpc_rep1_100/rep1_eh_dist_", i, ".mat")))[[1]]
  ewa <- readMat(here::here(paste0("/MATLAB/dgp_hpc_rep1_100/rep1_ewa_dist_", i, ".mat")))[[1]]
  pred = ewa %*% eh
  sim <- read_csv(paste0("./Sims/dgp_rep1/sim_dgp_rep1_", i, ".csv")) %>% as.matrix()
  normed[i] <- norm(sim - pred, "F")/norm(sim, "F")
}
summary(normed)

#####

dgp_bnmf_rep1 <- tibble()

for (i in 1:200) {
  eh <- readMat(here::here(paste0("./MATLAB/dgp_hpc_rep1_100/rep1_eh_dist_", i, ".mat")))[[1]]
  eh <- as_tibble(eh) %>% nest(eh = everything()) %>% mutate(row = i)
  
  ewa <- readMat(here::here(paste0("/MATLAB/dgp_hpc_rep1_100/rep1_ewa_dist_", i, ".mat")))[[1]]
  ewa <- as_tibble(ewa) %>% nest(ewa = everything()) %>% mutate(row = i)
  
  sim <- read_csv(paste0("./Sims/dgp_rep1/sim_dgp_rep1_", i, ".csv")) %>% 
    nest(sim = everything()) %>% mutate(row = i)
  
  both <- full_join(eh, ewa, by = "row") %>% full_join(., sim, by = "row")
  dgp_bnmf_rep1 <- rbind(dgp_bnmf_rep1, both)
}

dgp_bnmf_rep1 <- dgp_bnmf_rep1 %>% 
  mutate(seed = rep(1:100, 2),
         data = c(rep("Distinct", 100), rep("Overlapping", 100)))

dgp_bnmf_rep1

load("./Sims/sim_dgp_rep1_test.RDA")

sim_test <- sim_dgp_rep1_test$sim[[1]]
sim_dgb <- sim_dgp_rep1$sim[[1]]
sim <- read_csv(paste0("./Sims/dgp_rep1/sim_dgp_rep1_1.csv"))
head(sim)[,1:5]
head(sim_dgb)[,1:5]

test_sim <- sim_dgp_rep1 %>% dplyr::select(seed, data, sim) %>% mutate(sim = map(sim, as_tibble))
bnmf_sim <- dgp_bnmf_rep1 %>% dplyr::select(seed, data, sim)

all_equal(test_sim$sim[[100]], bnmf_sim$sim[[100]])
# NO
