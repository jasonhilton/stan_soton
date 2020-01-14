data {
  int <lower=0> J;
  real y[J];
  real <lower=0> sigma[J];
}

parameters {
  real mu;
  real<lower=0> tau;
  vector[J] theta_tilde;
}


model {
  vector[J] theta;
  theta = mu + theta_tilde * tau;
  mu ~ normal(0, 5);
  tau ~ cauchy(0, 5);
  theta_tilde ~ normal(0,1);
  y ~ normal(theta, sigma);
}
