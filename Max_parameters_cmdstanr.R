setwd("C:/Users/SmithAC/Documents/GitHub/bbsStanBayes")

library(cmdstanr)

N = 250000

x = rnorm(N)

y = x+rnorm(N,0,0.3)

stan_data <- list(N = N,
                  y = y,
                  x = x)

mod <- "models/simple_regression.stan"
model <- cmdstan_model(mod)


stanfit <- model$sample(
  data=stan_data,
  refresh=200,
  chains=4, 
  iter_sampling=1000,
  iter_warmup=1000,
  parallel_chains = 4,
  output_dir = "output",
  output_basename = "simple_regression_fit")




