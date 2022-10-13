## Fitting the BYM model to 1995 - 2021 BBS data on Azure



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


save(list = c("nrecs_sp","splitters","stratified_data"),
     file = "species_lists.RData")



# Prepare Data ------------------------------------------------------------


consider_spatial <- FALSE



#setwd("C:/GitHub/bbsStanBayes")

model = "gamye"

# 
# species_to_run <- nrecs_sp %>% 
#   filter(grouping == GG)

species_to_run <- nrecs_sp %>% 
  filter(english %in% c("American Robin","Barn Swallow"))

source("Functions/prepare-data-Stan.R")
if(consider_spatial){
  source("Functions/neighbours_define.R") # function to generate spatial neighbourhoods to add to the spatial applications of the models
}
 output_dir <- "D:/bbsStanBayes/output/" # Stan writes output to files as it samples. This is great because it's really stable, but the user needs to think about where to store that output
# 

 
 # parallel function -------------------------------------------------------
 
 run_model_azure <- function(species){
   
   species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)
   
   
   
   
   mod_file <- cmdstanr::write_stan_file(code = mod_string)
   
   
   
   
   
   
   ## compiles Stan model (this is only necessary if the model has been changed since it was last run on this machine)
   stan_model <- cmdstan_model(mod_file, stanc_options = list("Oexperimental"))
   
   
   ### this init_def is something that the JAGS versions don't need. It's a function definition, so perhaps something we'll have to build
   ### into the fit_model function
   ### the slightly annoying thing is that it's model-specific, so we'll need to have a few versions of it
   
   
   
   
   
   
   stanfit <- stan_model$sample(
     data=stan_data,
     refresh=100,
     chains=3, iter_sampling=1000,
     iter_warmup=1200,
     parallel_chains = 3,
     #pars = parms,
     adapt_delta = 0.8,
     max_treedepth = 13,
     seed = 123,
     init = init_def)
   
   # shinystan::launch_shinystan(shinystan::as.shinystan(stanfit))
   
   
   # loo_out <- stanfit$loo()
   
   
   tmp <- stanfit$draws()
   tmp <- stanfit$sampler_diagnostics()
   
   ret <- list(stanfit = stanfit)
   
   # save(list = c("stanfit","summ"),
   #      file = paste0(output_dir,"/",out_base,"_stanfit.RData"))
   return(ret)
 }
 
 
 
 # Load the doAzureParallel library 
 library(doAzureParallel) 
 
 setCredentials("credentials.json")

 # # generate your cluster in the cloud; this takes a few minutes
 cluster <- makeCluster("cluster.json") 
 # #this .json file currently includes instructions to use a rocker Bayesian image
 # 
 # 
 # # Register your parallel backend 
  registerDoAzureParallel(cluster) 
 
 # Check that the nodes are running 
 #getDoParWorkers() 
 
 # clusters <- vector(mode = "list",length = 2)
 async_jobs <- vector(mode = "list",length = 2) 
 
 for(jj in nrow(species_to_run):1){
   
   species <- as.character(species_to_run[jj,"english"])
   species_f <- as.character(species_to_run[jj,"species_file"])
   species_job <- gsub(pattern = " ",species,replacement = "")
   ## replaces the bbsBayes prepare_dta function because it includes additional infor required for Stan models
   
   
   
   if(consider_spatial){
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
   

     sp_data <- prepare_data(strat_data = stratified_data,
                             species_to_run = species,
                             model = model,
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
     
     
     stan_data[["stratify_by"]] <- NULL 
     stan_data[["model"]] <- NULL
     stan_data[["alt_data"]] <- NULL
     
     
    saveRDS(stan_data,file = paste0("Data/Stan_data_",out_base,".rds")) 
    
    stan_data2 <- stan_data
    stan_data2[["stratify_by"]] <- tmp_stratify_by 
    stan_data2[["model"]] <- tmp_model
    stan_data2[["alt_data"]] <- tmp_alt_data
    stan_data2[["strat_name"]] <- tmp_alt_data$strat_name
    
    saveRDS(stan_data2,file = paste0("Data/Stan_data2_",out_base,".rds")) 
    
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
    
    
    mod_string <- read_file(mod.file)
    


# AZURE Run ---------------------------------------------------------------





opt <- list(wait = FALSE)#, #you can set this to false if you want to free-up your local machine, but you then have to frequently check back with the cluster to download the results
            #enableCloudCombine = TRUE, #this combines the results into a single R-object list before downloading
            #job = ) 


t1 = Sys.time()

async_jobs[[jj]] <- foreach(i = species,
                     .packages = c("cmdstanr"),
                     .errorhandling = "pass", #this option passes any errors back to the console
                     .options.azure = opt) %dopar% {
                       
                       # This code is executed, in parallel, across your cluster.
                       run_model_azure(species = i)
                       
                     }



}



results[[jj]] <- try(getJobResult(async_jobs[[jj]]),silent = TRUE)

while(any(class(results) == "try-error")){
  results <- try(getJobResult(async_job),silent = TRUE)
  Sys.sleep(time = 3600*2)  
  
}


output_dir <- "output"

for(species in species_list){
  i <- which(species_list == species)
  stanfit <- results[[i]][["stanfit"]]
  
  species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)
  
  spp <- "_BYM_"
  
  out_base <- paste0(species_f,spp,firstYear,"_",lastYear)
  
  # stanfit$save_object(dir = output_dir, basename = out_base, timestamp = FALSE, random = FALSE)
  stanfit$save_object(file = paste0(output_dir,"/",out_base,"_stanfit.RDS"))
  # 
  # save(list = c("summ"),
  #      file = paste0(output_dir,"/",out_base,"_summ.RData"))
  
}
#results

t2 = Sys.time()
t2-t1
### when finished, be sure to run this stopCluster function, so you avoid paying for idle virtual machines
stopCluster(cluster)






# Post cluster summaries --------------------------------------------------
library(tidyverse)
library(cmdstanr)

species_list <- c("Chestnut-collared Longspur",
                  "Thick-billed Longspur")
firstYear = 1995
lastYear = 2021

source("Functions/posterior_summary_functions.R")

output_dir <- "output"

for(species in species_list){
  
  species_f <- gsub(gsub(species,pattern = " ",replacement = "_",fixed = T),pattern = "'",replacement = "",fixed = T)
  
  spp <- "_BYM_"
  
  out_base <- paste0(species_f,spp,firstYear,"_",lastYear)
  
  # stanfit$save_object(dir = output_dir, basename = out_base, timestamp = FALSE, random = FALSE)
  stanfit <- readRDS(paste0(output_dir,"/",out_base,"_stanfit.RDS"))
  
  print(stanfit$time())
  
  # slopes <- posterior_samples(fit = stanfit,
  #                             parm = "beta",
  #                             dims = "ste") 
  # intercepts <- posterior_samples(fit = stanfit,
  #                             parm = "alpha",
  #                             dims = "ste")
  # 
  # slope_sum <- posterior_sums(slopes,dims = "ste")
  
}












