// Mixed membership model in STAN

data {
  // data
  int  <lower=1> K;         // number of mixture components, topics
  int  <lower=1> N;         // number of individuals
  int  <lower=1> P;         // num chemicals
  real <lower=0> X[N, P];   // chemical concentrations

  // prior
  real<lower=0> mu1; // topic prior, column vector
  real<lower=0> sigma1;
  real<lower=0> mu2; // chemical prior, column vector
  real<lower=0> sigma2;  // 
}

parameters {
  positive_ordered [K] scores[N];    // topic dist for individual m :: mixing proportions, NxK
  vector <lower=0> [K] topics[P];             // chemical dist for topic k, PxK
}

model {
  for (n in 1:N) {
    scores[n] ~ lognormal(mu1, sigma1);    // prior on scores, non-negative and continuous dist
  }

  for (p in 1:P) {
    topics[p] ~ lognormal(mu2, sigma2);  // prior on topics
  }

  for (n in 1:N) {
    for (p in 1:P) 
      real out[1];
      out = log(scores[n]) + log(topics[p]);
      target += log_sum_exp(out); // likelihood is p(chemical | topics, scores)
}
}
// log-sum-of-exponents function is used to stabilize the numerical arithmetic.
// for (m in 1:M) {
//   for (v in 1:V) {
//     real gamma[K];
//       for (k in 1:K) {
//         gamma[k] = log(theta[m,k]) + log(phi[k,v]); // prob(word | theta, phi) = theta * phi
//       }
//     target += log_sum_exp(gamma); // likelihood;
//   }
// }

generated quantities {
}
