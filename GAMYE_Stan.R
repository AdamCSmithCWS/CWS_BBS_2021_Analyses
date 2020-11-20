## building a Stan version of the bbsBayes models

library(bbsBayes)
library(tidyverse)
library(rstan)
library(shinystan)



# load and stratify CASW data ---------------------------------------------
species = "Mourning Dove"
strat = "bbs_cws"
model = "gamye"

strat_data = stratify(by = strat)
jags_data = prepare_jags_data(strat_data = strat_data,
                             species_to_run = species,
                             model = model,
                             #min_year = 1999,
                             n_knots = 13)





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
          "sdyear_gam",
          "sdobs",
          "beta_p",
          "sdbeta",
          "strata_p",
          "sdstrata",
          "BETA",
          "STRATA",
           "n",
           "nsmooth",
          "eta")

## compile model
slope_model = stan_model(file=mod.file)

## run sampler on model, data
stime = system.time(slope_stanfit <-
                      sampling(slope_model,
                               data=stan_data,
                               verbose=TRUE, refresh=100,
                               chains=3, iter=500,
                               warmup=400,
                               cores = 3,
                               pars = parms,
                               control = list(adapt_delta = 0.8,
                                              max_treedepth = 13)))


paste(stime[[3]]/3600,"hours")
launch_shinystan(slope_stanfit) 



