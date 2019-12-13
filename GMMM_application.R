# Simulated Mothers and Newborns data
dat <- read_csv('./Sims/sim_data_1.csv') %>% as.matrix()

# Model
out_4 <- GMMM(x = dat, k = 4, sigma_2 = 1, m = 0, lambda_2 = 1)
out_4$FinalState

# Data Cleaning
out_mu <- out_4$TheStateRecord[grep("mu", names(out_4$TheStateRecord))]
out_mu <- array(as.numeric(unlist(out_mu)), dim=c(4, ncol(dat), out_4$StepCount))

out_theta <- out_4$TheStateRecord[grep("theta", names(out_4$TheStateRecord))]
out_theta <- array(as.numeric(unlist(out_theta)), dim=c(nrow(dat), 4, out_4$StepCount))

out_z <- out_4$TheStateRecord[grep("z", names(out_4$TheStateRecord))]
out_z <- array(as.numeric(unlist(out_z)), dim=c(nrow(dat), ncol(dat), out_4$StepCount))

out_lp <- out_4$TheStateRecord[grep("dev", names(out_4$TheStateRecord))] %>% as_vector() %>% 
  cbind(., 1:out_4$StepCount) %>% 
  as_tibble() %>% rename(deviance = 1, iteration = 2)

# Results
out_lp %>% 
  ggplot(aes(x = iteration, y = deviance)) +
  geom_line() +
  theme_minimal() +
  labs(x = "Iterations", y = "Deviance")

# Need to drop burn in !!!
mu_mc <- apply(out_mu, c(1,2), mean)
theta_mc <- apply(out_theta, c(1,2), mean)       
z_mc <- apply(out_z, c(1,2), mean)

# Need to drop burn in !!!
mu_mc_sd <- apply(out_mu, c(1,2), function(x) sd(as.vector(x)))
theta_mc_sd <- apply(out_theta, c(1,2), function(x) sd(as.vector(x)))
z_mc_sd <- apply(out_z, c(1,2), function(x) sd(as.vector(x)))

# Mu
t(mu_mc) %>% as_tibble() %>% 
  rename(ave1, ave2, ave3, ave4) %>% 
  cbind(., t(mu_mc_sd))







