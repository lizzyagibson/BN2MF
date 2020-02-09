data {
  int N;
  vector[N] x;
  vector[N] y;
}

parameters {
  real beta;
  real beta_0;
  real sigma;
}

model {
  // prior
  beta ~ normal(1, 1);
  sigma ~ gamma(1, 1);
  
  // likelihood
  y ~ normal(beta_0 + x * beta, sigma);
  
}

generated quantities {
  vector[N] y_pred;
  for (i in 1:N)
    y_pred[i] = normal_rng(beta_0 + x[i] * beta, sigma);
  
}
