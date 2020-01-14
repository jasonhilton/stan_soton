library(rstan)
library(bayesplot)


## Taken from https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

# Example based on https://mc-stan.org/users/documentation/case-studies/divergences_and_bias.html


# see  ?bayesplot::mcmc_nuts_divergence