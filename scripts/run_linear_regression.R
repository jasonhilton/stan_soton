library(rstan)



x1 <- seq(-5,5, 0.1)
x2 <- rnorm(length(x1))
X <- cbind(x1,x2)
Beta <- c(1, -3) 
alpha <- 10
sigma_e  <- 0.5 
y <- rnorm(length(x1), 
           alpha +  X %*% Beta, 
           sigma_e)

stan_input_data <- list(
  N=dim(X)[1],
  k=dim(X)[2],
  X=X,
  y=y
)

## running stan command here
stan_out <- stan(file="stan/linear_regression.stan",
                 data=stan_input_data,
                 iter=2000,
                 chains=2,
                 cores=2,
                 thin=1)
## traceplots
stan_trace(stan_out)
## ac plots
stan_ac(stan_out)
# rhat plots
stan_rhat(stan_out)

## posterior predictive checking using generating quantities
library(bayesplot)
y_rep <- as.matrix(stan_out, "y_rep")
y_real <- stan_input_data$y
ppc_dens_overlay(y_real, y_rep) + theme_bw(base_size = 16)
