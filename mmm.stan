// Mixed membership model in STAN

data {

  // data
  int  <lower=1> K;         // number of mixture components, topics
  int  <lower=1> N;         // number of individuals
  int  <lower=1> P;         // num chemicals
  real <lower=0> X[N, P];   // chemical concentrations

  // prior
  vector<lower=0>[K] alpha; // topic prior, column vector
  vector<lower=0>[P] beta;  // chemical prior, column vector
}

parameters {

  real <lower=0> phi[P,K];    // chemical dist for topic k, PxK
  real <lower=0> theta[K,N];  // topic dist for individual m :: mixing proportions, KxN
  real <lower=0> sigma[P];    // variance of original data
}

transformed parameters {
  real <lower=0> mu[N,P];

for (n in 1:N) {
  for (p in 1:P) {
    mu[n,p] = (phi[p,n]*theta[p,n])';
}}}

model {

for (n in 1:N)
  theta[N] ~ lognormal(alpha, beta);  // prior on theta

for (k in 1:K)
  phi[k] ~ lognormal(alpha, beta);    // prior on phi, non-negative and continuous dist

for (p in 1:P)
  sigma[p] ~ gamma(0,0); // prior on variance

  for (n in 1:N) {
    for (p in 1:P) {
    real out[K];
      for (k in 1:K) {
      
      out[k] = log(normal_lpdf(X[n,p] | mu[n,p], sigma[p]));
      // for h in 
      // recall bayes rule p(param | data ) = p (data | param) * p(param) / prob(data)
      // we want to esitmate param by sample optimizing over a f(param), 
      // prob(data) := prob(data | prior), constant so we drop it
      // prob(param ) := prop(param | prior), keep it
      // prob(data | param) := prob(data | param=p , Prob(param=p | prior )) 
      // prob(param | data,prior )   ===BAYES RULE == prob(data | param, prior) * prob(param | prior ) / prob(data | prior)

      target += log_sum_exp(out); // likelihood is p(chemical | theta, phi)
}
}
}
}
// log-sum-of-exponents function is used to stabilize the numerical arithmetic.

generated quantities {
}
