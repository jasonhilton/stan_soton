library(rstan)

args = commandArgs(trailingOnly=TRUE)

model_name <- args[1]

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

## traceplots

## ac plots

# rhat plots

## generating quantities

