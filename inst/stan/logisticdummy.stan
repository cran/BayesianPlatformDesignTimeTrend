functions {
}
data {
  real beta0_prior_mu;
  real beta1_prior_mu;
  real <lower=0> beta0_prior_sigma;
  real <lower=0> beta1_prior_sigma;
  int <lower=0> beta0_nu;
  int <lower=0> beta1_nu;
  
  int<lower=1> N;  // total number of observations
  int y[N];  // response variable
  int<lower=1> K;  // number of population-level effects
  int <lower=0, upper=K> z[N]; // arm applied on trial n
   matrix[N, K] x;  // population-level design matrix
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] xc;  // centered version of X without an intercept
  vector[Kc] means_x;  // column means of X before centering
  for (i in 2:K) {
    means_x[i - 1] = mean(x[, i]);
    xc[, i - 1] = x[, i] - means_x[i - 1];
  }
}
parameters {
  vector[K-1] beta1;  // population-level effects
  real beta0;  // temporary intercept for centered predictors
}
transformed parameters {
}
model {
  // likelihood including constants
  // priors including constants
  target += student_t_lpdf(beta1 | beta1_nu, beta1_prior_mu, beta1_prior_sigma);
  target += student_t_lpdf(beta0 | beta0_nu, beta0_prior_mu, beta0_prior_sigma);
  // target += normal_lpdf(beta0 | beta0_prior_mu, beta0_prior_sigma);
  // target += normal_lpdf(beta1 | beta1_prior_mu, beta1_prior_sigma);
  target += bernoulli_logit_glm_lpmf(y | xc, beta0, beta1);
    // y ~ bernoulli_logit(xdummy0*beta0+xdummy1*beta1);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = beta0 - dot_product(means_x, beta1);
  // // additionally sample draws from priors
  // real prior_b = normal_rng(0,10);
  // real prior_Intercept = normal_rng(0,10);
}
