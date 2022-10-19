library(bbsBayes)
library(tidyverse)
library(cmdstanr)

stratified_data <- stratify(by = "bbs_cws")


allspecies.eng = stratified_data$species_strat$english
allspecies.fre = stratified_data$species_strat$french
allspecies.num = stratified_data$species_strat$sp.bbs

allspecies.file = str_replace_all(str_replace_all(allspecies.eng,"[:punct:]",replacement = ""),
                                  "\\s",replacement = "_")

###################################################
# Analysis by Species X Model Combination
###################################################


nspecies = length(allspecies.eng)

name_simpl_function <- function(y){
  new <- gsub("[[:punct:]]", x = y, replacement = "")
  new <- gsub("\\s", x = new, replacement = "_")
  return(new)
}

nrecs_sp <- stratified_data$bird_strat %>% 
  group_by(AOU) %>% 
  summarise(num_counts = n(),
            num_routes = length(unique(Route))) %>% 
  left_join(.,stratified_data$species_strat,by = c("AOU" = "sp.bbs")) %>% 
  filter(num_counts > 200,
         !grepl("unid",english)) %>% 
  mutate(species_file = name_simpl_function(english),
         grouping = rep_len(1:9,length.out = length(unique(species_file)))) %>% 
  arrange(num_counts)


# split species groups that can't be separated in the early years ---------
splitters = c("Clark's Grebe","Western Grebe","Alder Flycatcher","Willow Flycatcher")
split_miny = c(1990,1990,1978,1978)
names(split_miny) <- splitters


save(list = c("nrecs_sp","splitters","stratified_data","split_miny"),
     file = "species_lists.RData")



# Run in separate R sessions ----------------------------------------------

GG <- 1

library(bbsBayes)
library(tidyverse)
library(cmdstanr)

consider_spatial <- TRUE


#setwd("C:/GitHub/bbsStanBayes")
setwd("C:/Users/SmithAC/Documents/GitHub/bbsStanBayes")

load("species_lists.RData")
model = "gamye"


species_to_run <- nrecs_sp %>% 
  filter(grouping == GG)

source("Functions/prepare-data-Stan.R")
if(consider_spatial){
  source("Functions/neighbours_define.R") # function to generate spatial neighbourhoods to add to the spatial applications of the models
}
output_dir <- "F:/bbsStanBayes/output/" # Stan writes output to files as it samples. This is great because it's really stable, but the user needs to think about where to store that output


for(jj in (1:nrow(species_to_run))){

species <- as.character(species_to_run[jj,"english"])
species_f <- as.character(species_to_run[jj,"species_file"])

## replaces the bbsBayes prepare_dta function because it includes additional infor required for Stan models


## Temporarily prepare the data for the species to see how many strata are likely
start_year <- NULL


if(consider_spatial){
  sp_data <- prepare_data(strat_data = stratified_data,
                          species_to_run = species,
                          model = model,
                          min_year = start_year,
                          min_max_route_years = 2,
                          basis = "mgcv")
if(sp_data$nstrata > 5){
  fit_spatial <- TRUE # TRUE = spatial sharing of information and FALSE = non-spatial sharing
  
}else{
  fit_spatial <- FALSE # TRUE = spatial sharing of information and FALSE = non-spatial sharing
}
}else{
  fit_spatial <- FALSE
}

if(fit_spatial){
  out_base <- paste(species_f,model,"Spatial","BBS",sep = "_") # text string to identify the saved output from the Stan process unique to species and model, but probably something the user wants to control
  
}else{
  out_base <- paste(species_f,model,"BBS",sep = "_")
  
} 

if(!file.exists(paste0(output_dir,"/",out_base,"_Stan_fit.RData")) |
   species %in% splitters){
  
 if(species %in% splitters){start_year <- split_miny[species]}
  
  sp_data <- prepare_data(strat_data = stratified_data,
                          species_to_run = species,
                          model = model,
                          min_year = start_year,
                          min_max_route_years = 2,
                          basis = "mgcv")
  
  
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
tmp_strat_name <- stan_data[["strat_name"]]

stan_data[["stratify_by"]] <- NULL 
stan_data[["model"]] <- NULL
stan_data[["alt_data"]] <- NULL
stan_data[["strat_name"]] <- NULL



if(fit_spatial){
  
mod.file = paste0("models/",model,"_spatial_bbs_CV.stan")
out_base <- paste(species_f,model,"Spatial","BBS",sep = "_") # text string to identify the saved output from the Stan process unique to species and model, but probably something the user wants to control

init_def <- function(){ list(noise_raw = rnorm(stan_data$ncounts*stan_data$use_pois,0,0.1),
                             strata_raw = rnorm(stan_data$nstrata,0,0.1),
                             STRATA = 0,
                             nu = 10,
                             sdstrata = runif(1,0.01,0.1),
                             eta = 0,
                             yeareffect_raw = matrix(rnorm(stan_data$nstrata*stan_data$nyears,0,0.1),nrow = stan_data$nstrata,ncol = stan_data$nyears),
                             obs_raw = rnorm(stan_data$nobservers,0,0.1),
                             ste_raw = rnorm(stan_data$nsites,0,0.1),
                             sdnoise = runif(1,0.03,1.3),
                             sdobs = runif(1,0.01,0.1),
                             sdste = runif(1,0.01,0.2),
                             #sdbeta = runif(stan_data$nstrata,0.01,0.1),
                             sdbeta = runif(1,0.01,0.1),
                             sdBETA = runif(1,0.01,0.1),
                             sdyear = runif(stan_data$nstrata,0.01,0.1),
                             BETA_raw = rnorm(stan_data$nknots_year,0,0.1),
                             beta_raw = matrix(rnorm(stan_data$nknots_year*stan_data$nstrata,0,0.01),nrow = stan_data$nstrata,ncol = stan_data$nknots_year))}



}else{
  mod.file = paste0("models/",model,"_bbs_CV.stan")
  out_base <- paste(species_f,model,"BBS",sep = "_")
  
  init_def <- function(){ list(noise_raw = rnorm(stan_data$ncounts*stan_data$use_pois,0,0.1),
                               strata_raw = rnorm(stan_data$nstrata,0,0.1),
                               STRATA = 0,
                               nu = 10,
                               sdstrata = runif(1,0.01,0.1),
                               eta = 0,
                               yeareffect_raw = matrix(rnorm(stan_data$nstrata*stan_data$nyears,0,0.1),nrow = stan_data$nstrata,ncol = stan_data$nyears),
                               obs_raw = rnorm(stan_data$nobservers,0,0.1),
                               ste_raw = rnorm(stan_data$nsites,0,0.1),
                               sdnoise = runif(1,0.03,1.3),
                               sdobs = runif(1,0.01,0.1),
                               sdste = runif(1,0.01,0.2),
                               sdbeta = runif(stan_data$nstrata,0.01,0.1),
                               #sdbeta = runif(1,0.01,0.1),
                               sdBETA = runif(1,0.01,0.1),
                               sdyear = runif(stan_data$nstrata,0.01,0.1),
                               BETA_raw = rnorm(stan_data$nknots_year,0,0.1),
                               beta_raw = matrix(rnorm(stan_data$nknots_year*stan_data$nstrata,0,0.01),nrow = stan_data$nstrata,ncol = stan_data$nknots_year))}
  
  
  
}

## compiles Stan model (this is only necessary if the model has been changed since it was last run on this machine)
stan_model <- cmdstan_model(mod.file, stanc_options = list("Oexperimental"))


### this init_def is something that the JAGS versions don't need. It's a function definition, so perhaps something we'll have to build
### into the fit_model function
### the slightly annoying thing is that it's model-specific, so we'll need to have a few versions of it






stanfit <- stan_model$sample(
  data=stan_data,
  refresh=100,
  chains=3, 
  iter_warmup=1200,
  iter_sampling=1000,
  parallel_chains = 3,
  #pars = parms,
  adapt_delta = 0.95,
  max_treedepth = 13,
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
stan_data[["strat_name"]] <- tmp_strat_name

save(list = c("stanfit","stan_data",
              "out_base",
              "fit_summary"),
     file = paste0(output_dir,"/",out_base,"_Stan_fit.RData"))


print(paste("Completed",out_base))
}else{
  print(paste("Skipped",out_base,"already run"))
  
}
}



