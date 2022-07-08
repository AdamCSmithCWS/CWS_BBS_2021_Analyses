library(bbsBayes)
library(tidyverse)
library(cmdstanr)

bbs_data <- stratify(by = "bbs_usgs")


setwd("C:/GitHub/bbsStanBayes")


source("Functions/prepare-data-alt.R")


species <- "Barn Swallow"


species_f <- gsub(species,pattern = " ",replacement = "_")


sp_data <- prepare_data(bbs_data,
                        species_to_run = species,
                        model = "gamye",
                        min_max_route_years = 10,
                        basis = "mgcv")


stan_data <- sp_data

stan_data[["stratify_by"]] <- NULL
stan_data[["alt_data"]] <- NULL
stan_data[["model"]] <- NULL



mod.file = "models/gamye_bbs_CV.stan"

## compile model
model <- cmdstan_model(mod.file)

out_base <- paste(species_f,sp_data$model,"BBS",sep = "_")


init_def <- function(){ list(noise_raw = rnorm(ncounts*stan_data$use_pois,0,0.1),
                             strata_raw = rnorm(stan_data$nstrata,0,0.1),
                             STRATA = 0,
                             nu = 10,
                             sdstrata = runif(1,0.01,0.1),
                             eta = 0,
                             yeareffect_raw = matrix(rnorm(stan_data$nstrata*nyears,0,0.1),nrow = stan_data$nstrata,ncol = stan_data$nyears),
                             obs_raw = rnorm(stan_data$nobservers,0,0.1),
                             ste_raw = rnorm(stan_data$nsites,0,0.1),
                             sdnoise = runif(1,0.3,1.3),
                             sdobs = runif(1,0.01,0.1),
                             sdste = runif(1,0.01,0.2),
                             sdbeta = runif(stan_data$nstrata,0.01,0.1),
                             sdBETA = runif(1,0.01,0.1),
                             sdyear = runif(stan_data$nstrata,0.01,0.1),
                             BETA_raw = rnorm(stan_data$nknots_year,0,0.1),
                             beta_raw = matrix(rnorm(stan_data$nknots_year*stan_data$nstrata,0,0.01),nrow = stan_data$nstrata,ncol = stan_data$nknots_year))}




stanfit <- model$sample(
  data=stan_data,
  refresh=100,
  chains=3, iter_sampling=1000,
  iter_warmup=1000,
  parallel_chains = 3,
  #pars = parms,
  adapt_delta = 0.95,
  max_treedepth = 14,
  seed = 123,
  init = init_def,
  output_dir = output_dir,
  output_basename = out_base)


#stanfit1 <- as_cmdstan_fit(files = paste0(output_dir,csv_files))

loo_out <- stanfit$loo()


save(list = c("stanfit","stan_data","csv_files",
              "out_base","loo_out"),
     file = paste0(output_dir,"/",out_base,"_fit.RData"))
}




