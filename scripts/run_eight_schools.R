library(rstan)
library(bayesplot)


## Taken from https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

# Example based on https://mc-stan.org/users/documentation/case-studies/divergences_and_bias.html



# Example based on https://mc-stan.org/users/documentation/case-studies/divergences_and_bias.html


stan_out <- stan(file="stan/eight_schools.stan",
                 data= schools_dat,
                 iter=2000,
                 cores=2,
                 chains=2)



np <- nuts_params(stan_out)
lp <- log_posterior(stan_out)
bayesplot::mcmc_nuts_divergence(np, lp)

all_pars <- as.array(stan_out)

bayesplot::mcmc_scatter(all_pars, pars = c("theta[1]", "tau"), np=np)

bayesplot::mcmc_scatter(all_pars, pars = c("theta[1]", "tau"), 
                        transformations=list(tau=log), np=np)


