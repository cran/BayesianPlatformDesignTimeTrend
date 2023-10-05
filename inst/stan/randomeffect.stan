data {
  int<lower=1> N;  // number of patients
  int<lower=1> K;  // number of treatment arms
  int<lower=1> groupmax;  // number of time points
  array[N] int<lower=1, upper=K> X;  // predictor variable
  array[N] int<lower=0, upper=1> Y;  // response variable
  array[N] int<lower=1> group;
}

parameters {
  real beta_0;  // intercept
  vector[K] beta;  // coefficients for predictor variable
  vector[groupmax] alpha;  // random effects for each time point
  real<lower=0> gamma;  // precision parameter for alpha
}

model {
  // priors
  beta_0 ~ normal(0, 1.8);
  beta ~ normal(0, 1.8);
  alpha[1] ~ normal(0,0.01);
  alpha[2] ~ normal(0, sqrt(1/gamma));
  alpha[3:groupmax] ~ normal(2*alpha[2:(groupmax-1)] - alpha[1:(groupmax-2)], sqrt(1/gamma));
  gamma ~ gamma(0.1, 0.01);

  // likelihood
    // for(n in 1:N) {
    //   y[n] ~ bernoulli(inv_logit(beta_0 + alpha[groupmax+1-group[n]] +  x[n]*beta));
    // }

  for (n in 1:N) {
    vector[K] pi;
    for (k in 1:K) {
      pi[k] = inv_logit(beta_0 + beta[k] * (X[n] == k) + alpha[groupmax-group[n]+1]);
    }
    Y[n] ~ bernoulli(pi[X[n]]);
  }
}
