library(bbsBayes)
library(tidyverse)
library(cmdstanr)

bbs_data <- stratify(by = "bbs_usgs")


setwd("C:/GitHub/bbsStanBayes")


source("Functions/prepare-data-alt.R")
source("Functions/neighbours_define_alt.R")


species <- "Pacific Wren"


species_f <- gsub(species,pattern = " ",replacement = "_")


sp_data <- prepare_data(bbs_data,
                        species_to_run = species,
                        model = "slope",
                        min_max_route_years = 10)


stan_data <- sp_data





# Spatial neighbourhoods --------------------------------------------------

base_strata_map <- bbsBayes::load_map(stratify_by = stan_data[["stratify_by"]])

alt_df <- stan_data[["alt_data"]][[1]]
strata_df <- alt_df %>% 
  select(strat,strat_name) %>% 
  distinct() %>% 
  arrange(strat)

realized_strata_map <- base_strata_map %>% 
  inner_join(.,strata_df,by = c("ST_12" = "strat_name"))




  
  
neighbours <- neighbours_define(real_strata_map = realized_strata_map, #sf map of strata
                                strat_link_fill = 10000, #distance to fill if strata are not connected
                                buffer = TRUE,
                                convex_hull = FALSE,
                                plot_neighbours = TRUE,
                                species = "",
                                plot_dir = "maps/",
                                plot_file = "_strata_map",
                                save_plot_data = TRUE,
                                voronoi = FALSE,
                                nn_fill = FALSE,
                                add_map = NULL,
                                strat_indicator = "strat",
                                island_link_dist_factor = 1.2 #consider nearest strata neighbours if distances are within this factor of each other, when linking otherwise isolated islands of strata
                                )



stan_data[["N_edges"]] = neighbours$N_edges
stan_data[["node1"]] <- neighbours$node1
stan_data[["node2"]] <- neighbours$node2


stan_data[["stratify_by"]] <- NULL
stan_data[["model"]] <- NULL




mod.file = "models/slope_spatial_bbs_CV.stan"

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
                             sdBETA = runif(1,0.01,0.1),
                             sdyear = runif(stan_data$nstrata,0.01,0.1),
                             BETA = rnorm(1,0,0.1),
                             beta_raw = rnorm(stan_data$nstrata,0,0.1))}




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




