# Simulated Mothers and Newborns data
dat <- read_csv('./Sims/sim_data_1.csv') %>% as.matrix()

set.seed(850340)
index <- sample(1:nrow(dat), 750, replace = FALSE)
train <- dat[index,]
test <- dat[-index,]

apply(train, 2, sd)

# Model
start_time <- Sys.time()
out_44 <- GMMM(x = train, k = 4, sigma_2 = 0.5, m = 0, lambda_2 = 1)
end_time <- Sys.time()
end_time - start_time

out_4 <- GMMM(x = train, k = 4, sigma_2 = 1, m = 0, lambda_2 = 1)
#save(out_4, file = "./out_4.RData")
#load("./out_4.RData")

# Data Cleaning
out_mu <- out_4$TheStateRecord[grep("mu", names(out_4$TheStateRecord))]
out_mu <- array(as.numeric(unlist(out_mu)), dim=c(4, ncol(train), out_4$StepCount))
dim(out_mu)

# Remove Burn in
out_mu_burn <- out_mu[,,101:1001]

out_theta <- out_4$TheStateRecord[grep("theta", names(out_4$TheStateRecord))]
out_theta <- array(as.numeric(unlist(out_theta)), dim=c(nrow(train), 4, out_4$StepCount))
dim(out_theta)

# Remove Burn in
out_theta_burn <- out_theta[,,101:1001]

# Don't care about z
#out_z <- out_4$TheStateRecord[grep("z", names(out_4$TheStateRecord))]
#out_z <- array(as.numeric(unlist(out_z)), dim=c(nrow(train), ncol(train), out_4$StepCount))
#dim(out_z)

# Remove Burn in
#out_z_burn <- out_z[,,101:1001]

out_lp <- out_4$TheStateRecord[grep("dev", names(out_4$TheStateRecord))] %>% as_vector() %>% 
  cbind(., 1:out_4$StepCount) %>% 
  as_tibble() %>% rename(deviance = 1, iteration = 2)

# Results
out_lp %>% 
  ggplot(aes(x = iteration, y = deviance)) +
  geom_line() +
  theme_minimal() +
  labs(x = "Iterations", y = "Deviance") +
  ylim(-5000000, -3750000)

# pdf("gmmm_logjoint.pdf")
# out_lp %>% 
#   ggplot(aes(x = iteration, y = deviance)) +
#   geom_line() +
#   theme_minimal() +
#   labs(x = "Iterations", y = "Deviance") +
#   ylim(-4500000, -3850000)
# dev.off()

# Burn-in dropped!
mu_mc <- apply(out_mu_burn, c(1,2), mean)
theta_mc <- apply(out_theta_burn, c(1,2), mean)   
# Does average sum to one?
rowSums(theta_mc) # Yes!

# Don't care about z
#z_mc <- apply(out_z_burn, c(1,2), mean)

# Burn-in dropped!
mu_mc_sd <- apply(out_mu_burn, c(1,2), function(x) sd(as.vector(x)))
theta_mc_sd <- apply(out_theta_burn, c(1,2), function(x) sd(as.vector(x)))
# Don't care about z
#z_mc_sd <- apply(out_z_burn, c(1,2), function(x) sd(as.vector(x)))

# Mu
mu_mean <- t(mu_mc) %>% as_tibble(.) %>% 
  rename(ave1 = 1, ave2 = 2, ave3 = 3, ave4 = 4) %>% 
  cbind(., chemicals = colnames(dat)) %>% 
  pivot_longer(cols = starts_with("ave"), 
               names_to = c("pattern"),
               values_to = "mean") %>% 
  mutate(pattern = substring(pattern, 4))


mu_sd <- t(mu_mc_sd) %>% as_tibble(.) %>% 
  rename(sd1 = 'V1', sd2 = 'V2', sd3 = 'V3', sd4 = 'V4') %>% 
  cbind(., chemicals = colnames(dat)) %>% 
  pivot_longer(cols = starts_with("sd"), 
               names_to = c("pattern"),
               values_to = "sd") %>% 
  mutate(pattern = substring(pattern, 3))
  
mu_full <- full_join(mu_mean, mu_sd, by = c("chemicals", "pattern")) %>% 
  mutate(class = case_when(grepl("pcb", chemicals) == TRUE ~ "PCBs",
                           grepl("BDE", chemicals) == TRUE ~ "PBDEs",
                           grepl("^M", chemicals) == TRUE ~ "Phthalates",
                           grepl("pcb", chemicals) == TRUE ~ "PCBs",
                           grepl("pcb", chemicals) == TRUE ~ "PCBs",
                           chemicals %in% c("tcs", "b_pb", "bp_3", "bpa", 
                                            "dcp_24", "dcp_25", "m_pb", "p_pb") ~ "Phenols"))

mu_full %>% 
  ggplot(aes(x = chemicals, y = mean, fill = class)) +
  geom_col() +
  facet_wrap(.~ pattern) +
  theme_classic() +
  labs(x = "Patterns", y = "Average chemical concentration", fill = "Chemical class") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# pdf("mu_plot.pdf", width = 15)
# mu_full %>% 
#   ggplot(aes(x = chemicals, y = mean, fill = class)) +
#   geom_col() +
#   facet_wrap(.~ pattern) +
#   theme_classic() +
#   labs(x = "Patterns", y = "Average chemical concentration", fill = "Chemical class") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = "bottom")
# dev.off()
  
# Theta
theta_mean <- theta_mc %>% as_tibble(.) %>% 
  rename(ave1 = 1, ave2 = 2, ave3 = 3, ave4 = 4) %>% 
  cbind(., person = row_number(train)) %>% 
  pivot_longer(cols = starts_with("ave"), 
               names_to = c("pattern"),
               values_to = "mean") %>% 
  mutate(pattern = substring(pattern, 4))

theta_sd <- theta_mc_sd %>% as_tibble(.) %>% 
  rename(sd1 = 'V1', sd2 = 'V2', sd3 = 'V3', sd4 = 'V4') %>% 
  cbind(., person = row_number(train)) %>% 
  pivot_longer(cols = starts_with("sd"), 
               names_to = c("pattern"),
               values_to = "sd") %>% 
  mutate(pattern = substring(pattern, 3))

theta_full <- full_join(theta_mean, theta_sd, by = c("person", "pattern"))

# Randomly chosen
theta_full %>% 
  filter(person %in% c(2667, 3546, 33629, 29136)) %>% 
  ggplot(aes(x = pattern, y = mean, fill = pattern)) +
  geom_col() +
  facet_wrap(.~ as.character(person)) +
  theme_classic() +
  labs(x = "Patterns", y = "Proportion of exposure contributed") +
  theme(legend.position = "none")

# pdf("theta_plot.pdf")
# theta_full %>% 
#   filter(person %in% c(2667, 3546, 33629, 29136)) %>% 
#   ggplot(aes(x = pattern, y = mean, fill = pattern)) +
#   geom_col() +
#   facet_wrap(.~ as.character(person)) +
#   theme_classic() +
#   labs(x = "Patterns", y = "Proportion of exposure contributed") +
#   theme(legend.position = "none")
# dev.off()
  
# Test set
dim(mu_mc)
dim(theta_mc)  
dim(test)  


test_out <- test %*% t(mu_mc) %>% 
  as_tibble(.)

test_out <- test_out / rowSums(test_out)
  





