
data {
  int<lower=0> n;
  vector[n] steps;
  vector[n-1] turns;
}


parameters {
  real<lower= -pi(), upper = pi()> mu;
  real<lower=0> kappa;
  real<lower=0> scale;
  real<lower=0> shape;
}


model {
  mu ~ normal(0,1);
  kappa ~ normal(0,10);
  scale ~ normal(0,1);
  shape ~ normal(1,1);
  steps ~ weibull(shape, scale);
  turns ~ von_mises(mu, kappa);
}

