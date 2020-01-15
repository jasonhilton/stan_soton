data{
  int N;
  int k;
  matrix[N, k] X;
  vector[N] y;
}

// The parameters accepted by the model. Our model
// accepts  parameters 'alpha', 'beta' and 'sigma'.
parameters{
  real alpha;
  vector[k] beta;
  real<lower=0> sigma;
}

// The model to be estimated. 
model{
    y ~ normal(alpha + X * beta, sigma);
}

