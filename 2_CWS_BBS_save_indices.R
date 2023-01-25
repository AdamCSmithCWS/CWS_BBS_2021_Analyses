## plotting

### NOTE: this script imports some functions that overwrite some of the standard bbsBayes functions to allow the rest of the package to work with Stan output
## prepare_data(), plot_indices(), generate_indices(), generate_trends(), and extract_index_data() (the last is called within generate_indices())


library(bbsBayes)
library(tidyverse)
library(cmdstanr)
library(patchwork)
library(doParallel)
library(foreach)

# Basic requirements: to generate and save
#   both official languages for outputs?
#   indices - both kinds - output rds of species-named lists
#   trends for short and long-term - output rds of species-named lists
#   trajectory plots - diagnostic output pdf - save higher level plots in list of ggplots for combined printing
#   trajectory plots - save all plots in list of ggplots for combined printing
#   trend maps - save all plots in list of ggplots for printing
#   rolling-trends - output rds of species-named lists

#
#   abundance maps? are these useful anymore in a world of eBird?
    
  
# Set a quality threshold for regions with extreme upper limits to the abundance distribution (BARS, CLSW, etc.)
# consider excluding strata with no data for first x-years
# use the MDPI for cis and ci width


# Generate trend uncertainty maps that show the variation in trend estimates across ~50 samples from the posterior
# or at least a map that plots the CIs of the trends

source("functions/posterior_summary_functions.R")
source("functions/generate-indices-alt.R")
source("functions/plot-indices-alt.R")
source("functions/extract-index-data-alt.R")
source("functions/generate-trends-alt.R")


setwd("C:/Users/SmithAC/Documents/GitHub/bbsStanBayes")

strat_sel <- "bbs_cws"


model_sel = "gamye"
Non_hierarchical <- FALSE

fit_spatial <- FALSE # TRUE = spatial sharing of information and FALSE = non-spatial sharing
source("Functions/prepare-data-Stan.R")
if(fit_spatial){
  source("Functions/neighbours_define.R") # function to generate spatial neighbourhoods to add to the spatial applications of the models
}
#output_dir <- "output/" # Stan writes output to files as it samples. This is great because it's really stable, but the user needs to think about where to store that output
output_dir <- "F:/bbsStanBayes/output" # Stan writes output to files as it samples. This is great because it's really stable, but the user needs to think about where to store that output



YYYY = 2021 ## BBS data version year
short_time = 10 #length of time for short-term trend


# load previous coverage data -----------------------------------------------------------

lastyear = read.csv("C:/Users/SmithAC/Documents/GitHub/BBS_Summaries/2017estimates/All BBS trends 2017 w reliab.csv",stringsAsFactors = F)
covs = lastyear[,c("sp","species","geo.area","trendtype","trendtime","startyear","reliab.cov")]

covs <- covs[which((covs$trendtime == "full")),]

oldregs = read.csv("C:/Users/SmithAC/Documents/GitHub/BBS_Summaries/old region names.csv",stringsAsFactors = F)

covs = merge(covs,oldregs,by = "geo.area")



# regions to estimate -----------------------------------------------------

regs_to_estimate <- c("continental","national","prov_state","bcr","stratum","bcr_by_country")



load("species_lists.RData") # loads objects created at the beginning of the script Fit_gamye_models_cws.R

species_to_run <- nrecs_sp 

# sp_tmp <- c("American Robin",
#             "Barn Swallow",
#             "House Sparrow",
#             "American Crow",
#             "Carolina Wren",
#             "Greater Roadrunner")
# 
# species_to_run <- filter(species_to_run,
#                          english %in% sp_tmp)
# 
# Species loop to generate the indices in parallel ------------------------------------------------------------
n_cores <- 6#length(provs)
cluster <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cluster)


fullrun <- foreach(jj = rev(1:nrow(species_to_run)),
                   .packages = c("posterior","tidyverse","bbsBayes"),
                   .inorder = FALSE,
                   .errorhandling = "pass") %dopar%
  {
    

#for(jj in rev(1:nrow(species_to_run))[c(1:20)]){
  
  
  source("functions/posterior_summary_functions.R")
  source("functions/generate-indices-alt.R")
  source("functions/plot-indices-alt.R")
  source("functions/extract-index-data-alt.R")
  source("functions/generate-trends-alt.R")
  
  
  source("Functions/prepare-data-Stan.R")
  if(fit_spatial){
    source("Functions/neighbours_define.R") # function to generate spatial neighbourhoods to add to the spatial applications of the models
  }
 
  species <- as.character(species_to_run[jj,"english"])
  species_f <- as.character(species_to_run[jj,"species_file"])
  espece <- as.character(species_to_run[jj,"french"])
  
  aou <- as.character(species_to_run[jj,"AOU"])
  

if(fit_spatial){
  
  out_base <- paste(species_f,model_sel,"Spatial","BBS",sep = "_") # text string to identify the saved output from the Stan process unique to species and model, but probably something the user wants to control
  
}else{
  out_base <- paste(species_f,model_sel,"BBS",sep = "_")
  if(Non_hierarchical){
    out_base <- paste(species_f,model_sel,"NonHier_BBS",sep = "_")
  }
}

  

# load model fit ----------------------------------------------------------

  
if(!file.exists(paste0(output_dir,"/",out_base,"_Stan_fit.RData"))){next}

load(paste0(output_dir,"/",out_base,"_Stan_fit.RData"))

if(model_sel == "gamye"){alt_n <- "nsmooth"}
if(model_sel == "slope"){alt_n <- "nslope"}

if(model_sel %in% c("firstdiff","gam")){alt_n <- NA}
if(is.null(stan_data$strat_name)){
  stan_data$strat_name <- stan_data$alt_data$strat_name
}


# estimate indices --------------------------------------------------------


ind <- generate_indices(jags_mod = stanfit,
                        jags_data = stan_data,
                        backend = "Stan",
                        stratify_by = strat_sel,
                        alternate_n = "n",
                        max_backcast = 15,
                        n_obs_backcast = 2,
                        #drop_exclude = TRUE,
                        regions = regs_to_estimate)



if(!is.na(alt_n)){
  
inds <- generate_indices(jags_mod = stanfit,
                         jags_data = stan_data,
                         backend = "Stan",
                         stratify_by = strat_sel,
                         alternate_n = alt_n,
                         max_backcast = 15,
                         n_obs_backcast = 2,
                         #drop_exclude = TRUE,
                         regions = regs_to_estimate)
}else{
  inds <- ind
}


saveRDS(ind,file = paste0(output_dir,"/temp_rds_storage/",species_f,"_indices.RDS"))
saveRDS(inds,file = paste0(output_dir,"/temp_rds_storage/",species_f,"_indices_smooth.RDS"))


} #temp end loop


stopCluster(cl = cluster)





