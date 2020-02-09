library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# The Stan model as a string.
model_string <- "
data {
// Number of data points
int n1;
int n2;
// Number of successes
int y1[n1];
int y2[n2];
}

parameters {
real<lower=0, upper=1> theta1;
real<lower=0, upper=1> theta2;
}

model {  
theta1 ~ beta(1, 1);
theta2 ~ beta(1, 1);
y1 ~ bernoulli(theta1);
y2 ~ bernoulli(theta2); 
}

generated quantities {
}
"

y1 <- c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0)
y2 <- c(0, 0, 1, 1, 1, 0, 1, 1, 1, 0)
data_list <- list(y1 = y1, y2 = y2, n1 = length(y1), n2 = length(y2))

# Compiling and producing posterior samples from the model.
stan_samples <- stan(model_code = model_string, data = data_list)

# Plotting and summarizing the posterior distribution
stan_samples
plot(stan_samples)
