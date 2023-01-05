## testing for completeness




library(bbsBayes)
library(tidyverse)
library(cmdstanr)

consider_spatial <- FALSE


#setwd("C:/GitHub/bbsStanBayes")
setwd("C:/Users/SmithAC/Documents/GitHub/bbsStanBayes")

load("species_lists.RData")
model = "gamye"


source("Functions/prepare-data-Stan.R")
if(consider_spatial){
  source("Functions/neighbours_define.R") # function to generate spatial neighbourhoods to add to the spatial applications of the models
}
output_dir <- "F:/bbsStanBayes/output/" # Stan writes output to files as it samples. This is great because it's really stable, but the user needs to think about where to store that output

species_remaining <- NULL

for(jj in rev(1:nrow(nrecs_sp))){
  
  species <- as.character(nrecs_sp[jj,"english"])
  species_f <- as.character(nrecs_sp[jj,"species_file"])
  
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
  
  if(!file.exists(paste0(output_dir,"/",out_base,"_Stan_fit.RData"))){
    print(species)
    species_remaining <- c(species_remaining,species)
  } 
  
  
}

n_rem <- nrecs_sp %>% 
  filter(english %in% species_remaining)
  
  
     