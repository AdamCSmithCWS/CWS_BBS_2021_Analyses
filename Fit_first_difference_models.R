library(bbsBayes)
library(tidyverse)
library(cmdstanr)

bbs_data <- stratify(by = "bbs_cws")


#setwd("C:/GitHub/bbsStanBayes")

setwd("C:/Users/SmithAC/Documents/GitHub/bbsStanBayes")


model_sel <- "firstdiff"
Non_hierarchical <- FALSE


species <- "Golden-winged Warbler"
species_f <- gsub(species,pattern = " ",replacement = "_") # species name without spaces


# # Initial fit using JAGS for comparison of some key parameters (see parameters_to_save argument in run_model()) ----------------------------------------------------------
jags_data <- prepare_data(strat_data = bbs_data,
                          species_to_run = species,
                          model = model_sel,
                          min_max_route_years = 2,
                          heavy_tailed = TRUE)

# model_to_file(model = model_sel,heavy_tailed = TRUE,filename = "temp_firstdif.R")
# model_to_file(model = model_sel,heavy_tailed = FALSE,filename = "temp_firstdif2.R")

jagsfit <- bbsBayes::run_model(jags_data = jags_data,
                               parameters_to_save = c("n","nslope",
                                                      "yeareffect","STRATA",
                                                      "sdobs","sdbeta","eta",
                                                      "sdnoise","nu",
                                                      "sdyear"),
                               parallel = TRUE,
                               modules = "glm",
                               model_file_path = "temp_firstdif.R")


inds <- generate_indices(jagsfit, jags_data,
                         stratify_by = "bbs_cws",
                         regions = c("stratum","continental","prov_state"))


plots <- plot_indices(inds,add_observed_means = TRUE,add_number_routes = TRUE)

pdf("temp.pdf")
for(i in 1:length(plots)){
  print(plots[[i]])
}

dev.off()

# 
# save(list = c("jagsfit","jags_data","species"),
#      file = paste("output/saved_bbsBayes_fit_data_",species_f,"_",model,".RData",sep = "_"))
# 



# Stan models -------------------------------------------------------------

## spatial versions of both teh slope and gamye exist for the Stan models and can be fit with this script 
fit_spatial <- TRUE # TRUE = spatial sharing of information and FALSE = non-spatial sharing

## the bbsBayes prepare_data function doesn't create all of the objects required for the Stan versions of the models
## this source() call over-writes the bbsBayes function prepare_data()
source("Functions/prepare-data-Stan.R")
if(fit_spatial){
source("Functions/neighbours_define.R") # function to generate spatial neighbourhoods to add to the spatial applications of the models
}



sp_data <- prepare_data(bbs_data,
                        species_to_run = species,
                        model = model_sel,
                        min_max_route_years = 2)


stan_data <- sp_data

# Spatial neighbourhoods --------------------------------------------------
if(fit_spatial){
  
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
                                species = species,
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
}#end of if fit_spatial



# extra list elements not required by Stan
# temporarily save them as objects then return them to the stan_data list after model fitting (below)
tmp_stratify_by <- stan_data[["stratify_by"]]  
tmp_model <- stan_data[["model"]]
tmp_alt_data <- stan_data[["alt_data"]]


stan_data[["stratify_by"]] <- NULL 
stan_data[["model"]] <- NULL
stan_data[["alt_data"]] <- NULL

stan_data[["zero_betas"]] <- rep(0,stan_data$nstrata)
stan_data[["Iy1"]] <- c((stan_data$fixedyear-1):1)
stan_data[["nIy1"]] <- length(stan_data[["Iy1"]])

stan_data[["Iy2"]] <- c((stan_data$fixedyear+1):stan_data$nyears)
stan_data[["nIy2"]] <- length(stan_data[["Iy2"]])


if(fit_spatial){
  
mod.file = "models/first_difference_spatial_bbs_CV.stan"
out_base <- paste(species_f,sp_data$model,"Spatial","BBS",sep = "_") # text string to identify the saved output from the Stan process unique to species and model, but probably something the user wants to control

}else{
  if(Non_hierarchical){
  mod.file = "models/first_difference_NonHier_bbs_CV.stan"
  out_base <- paste(species_f,sp_data$model,"NonHier_BBS",sep = "_")
  }else{
  mod.file = "models/first_difference_bbs_CV.stan"
  out_base <- paste(species_f,sp_data$model,"BBS",sep = "_")
  }
}

## compiles Stan model (this is only necessary if the model has been changed since it was last run on this machine)
model <- cmdstan_model(mod.file, stanc_options = list("Oexperimental"))

output_dir <- "output/" # Stan writes output to files as it samples. This is great because it's really stable, but the user needs to think about where to store that output

### this init_def is something that the JAGS versions don't need. It's a function definition, so perhaps something we'll have to build
### into the fit_model function
### the slightly annoying thing is that it's model-specific, so we'll need to have a few versions of it
init_def <- function(){ list(noise_raw = rnorm(stan_data$ncounts*stan_data$use_pois,0,0.1),
                             strata_raw = rnorm(stan_data$nstrata,0,0.1),
                             STRATA = 0,
                             nu = 10,
                             sdstrata = runif(1,0.01,0.1),
                             eta = 0,
                             obs_raw = rnorm(stan_data$nobservers,0,0.1),
                             ste_raw = rnorm(stan_data$nsites,0,0.1),
                             sdnoise = runif(1,0.3,1.3),
                             sdobs = runif(1,0.01,0.1),
                             sdste = runif(1,0.01,0.2),
                             sdbeta = runif(1,0.01,0.1),
                             sdBETA = runif(1,0.01,0.1),
                             BETA_raw = rnorm((stan_data$nyears-1),0,0.1),
                             beta_raw = matrix(rnorm(stan_data$nstrata*(stan_data$nyears-1),0,0.1),nrow = stan_data$nstrata,ncol = stan_data$nyears-1))}

if(Non_hierarchical){
  init_def <- function(){ list(noise_raw = rnorm(stan_data$ncounts*stan_data$use_pois,0,0.1),
                               strata_raw = rnorm(stan_data$nstrata,0,0.1),
                               STRATA = 0,
                               nu = 10,
                               sdstrata = runif(1,0.01,0.1),
                               eta = 0,
                               obs_raw = rnorm(stan_data$nobservers,0,0.1),
                               ste_raw = rnorm(stan_data$nsites,0,0.1),
                               sdnoise = runif(1,0.3,1.3),
                               sdobs = runif(1,0.01,0.1),
                               sdste = runif(1,0.01,0.2),
                               sdbeta = runif(stan_data$nstrata,0.01,0.1),
                               beta_raw = matrix(rnorm(stan_data$nstrata*(stan_data$nyears-1),0,0.1),nrow = stan_data$nstrata,ncol = stan_data$nyears-1))}
  
}


stanfit <- model$sample(
  data=stan_data,
  refresh=100,
  chains=3, iter_sampling=1000,
  iter_warmup=1000,
  parallel_chains = 3,
  #pars = parms,
  adapt_delta = 0.8,
  max_treedepth = 12,
  seed = 123,
  init = init_def,
  output_dir = output_dir,
  output_basename = out_base)

# shinystan::launch_shinystan(shinystan::as.shinystan(stanfit))


# loo_out <- stanfit$loo()


fit_summary <- stanfit$summary()

stan_data[["stratify_by"]] <- tmp_stratify_by 
stan_data[["model"]] <- tmp_model
stan_data[["alt_data"]] <- tmp_alt_data
stan_data[["strat_name"]] <- tmp_alt_data$strat_name

save(list = c("stanfit","stan_data",
              "out_base",
              "fit_summary"),
     file = paste0(output_dir,"/",out_base,"_Stan_fit.RData"))





