data {
  int<lower=1> K; // number of  arms
  int<lower=0> N; // number of trials
  int<lower=1, upper=K> z[N]; // arm on trial n
  int<lower=0, upper=1> y[N]; // outcome on trial n
  real <lower=0> pistar;
  int <lower=0> pess;
}
parameters {
  vector<lower=0, upper=1>[K] theta; // arm return prob
}
model {
  theta ~ beta(pess*pistar,pess*(1-pistar)); //Prior of theta
  y ~ bernoulli(theta[z]); //Assuming i.i.d. for each arm.
  //Each time we apply arm k, the outcome follows the same distribution and
  //does not depend on last observation
}
generated quantities {
  simplex[K] times_to_be_best; // phi_{n} = times_to_be_best[n] = sum(phi_{n,k})=1
  //for k in 1...K which is the length of this simplex.
  {
    real best_prob = max(theta);
    for (k in 1 : K) { // Due to 1000 burn-in samples, and 3000 iterations,
    // there will be 2000 times comparison for each round n.
      times_to_be_best[k] = theta[k] >= best_prob;
    }
    times_to_be_best /= sum(times_to_be_best);
  }
}
