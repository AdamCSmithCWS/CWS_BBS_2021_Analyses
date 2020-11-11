## building a Stan version of the bbsBayes models

library(bbsBayes)
library(tidyverse)
library(rstan)
library(shinystan)



# load and stratify CASW data ---------------------------------------------
species = "Pacific Wren"
strat = "bbs_usgs"
model = "gamye"

strat_data = stratify(by = strat)
jags_data = prepare_jags_data(strat_data = strat_data,
                             species_to_run = species,
                             model = model,
                             n_knots = 7,
                             min_year = 1999)





stan_data = jags_data[c("ncounts",
                         "nstrata",
                         "nobservers",
                         "count",
                         "strat",
                         "obser",
                         "year",
                         "firstyr",
                         "fixedyear",
                        "nonzeroweight")]
stan_data[["nyears"]] <- max(jags_data$year)
stan_data[["max_nobservers"]] <- max(jags_data$nobservers)
stan_data[["nknots_year"]] <- jags_data$nknots
stan_data[["year_basispred"]] <- jags_data$X.basis


mod.file = "models/gamye.stan"

parms = c("sdnoise",
          "sdyear",
          "sdobs",
          "beta_p",
          "sdbeta",
          "strata_p",
          "sdstrata",
          "BETA",
          "STRATA",
          # "n",
          # "nsmooth",
          "eta")

## compile model
slope_model = stan_model(file=mod.file)

## run sampler on model, data
stime = system.time(slope_stanfit <-
                      sampling(slope_model,
                               data=stan_data,
                               verbose=TRUE, refresh=100,
                               chains=3, iter=2000,
                               warmup=1000,
                               cores = 3,
                               pars = parms,
                               control = list(adapt_delta = 0.8,
                                              max_treedepth = 15)))


launch_shinystan(slope_stanfit) 



