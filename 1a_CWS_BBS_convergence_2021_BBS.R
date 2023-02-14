## plotting

### NOTE: this script depends on a development version of bbsBayes at the following repo: https://github.com/AdamCSmithCWS/bbsBayes/tree/testing_Stan 
## the specific changes in that branch are in the two functions
## generate_indices() and extract_index_data() (the second is called within generate_indices())

setwd("C:/Users/SmithAC/Documents/GitHub/bbsStanBayes")

strat_sel <- "bbs_cws"


library(tidyverse)
library(cmdstanr)


load("species_lists.RData") # loads objects created at the beginning of the script Fit_gamye_models_cws.R

model_sel = "gamye"
Non_hierarchical <- FALSE

species_to_run <- nrecs_sp 

fit_spatial <- FALSE # TRUE = spatial sharing of information and FALSE = non-spatial sharing
source("Functions/prepare-data-Stan.R")
if(fit_spatial){
  source("Functions/neighbours_define.R") # function to generate spatial neighbourhoods to add to the spatial applications of the models
}
#output_dir <- "output/" # Stan writes output to files as it samples. This is great because it's really stable, but the user needs to think about where to store that output
output_dir <- "F:/bbsStanBayes/output" # Stan writes output to files as it samples. This is great because it's really stable, but the user needs to think about where to store that output

pfail <- function(x,criteria,less_than = TRUE){
  if(less_than){
    y <- length(which(x < criteria))/length(x)
  }else{
    y <- length(which(x > criteria))/length(x)
  }
  return(round(y,3))
 
}
fail_params <- NULL

fail_summary <- NULL

for(jj in rev(1:nrow(species_to_run))){
  
  species <- as.character(species_to_run[jj,"english"])
  species_f <- as.character(species_to_run[jj,"species_file"])
  



 
if(fit_spatial){
  
  out_base <- paste(species_f,model_sel,"Spatial","BBS",sep = "_") # text string to identify the saved output from the Stan process unique to species and model, but probably something the user wants to control
  
}else{
  out_base <- paste(species_f,model_sel,"BBS",sep = "_")
  if(Non_hierarchical){
    out_base <- paste(species_f,model_sel,"NonHier_BBS",sep = "_")
  }
}

if(!file.exists(paste0(output_dir,"/",out_base,"_Stan_fit.RData"))){next}

load(paste0(output_dir,"/",out_base,"_Stan_fit.RData"))

failt <- fit_summary %>% 
  filter((ess_bulk < 200 | rhat > 1.05)) %>% 
  mutate(species = species)

fail_sum <- fit_summary %>% 
  mutate(parameter_group = gsub(pattern = "[[:punct:]]|[[:digit:]]",
                                x = variable,
                                replacement = "")) %>% 
  group_by(parameter_group) %>% 
  summarise(pfail_ess_bulk = pfail(ess_bulk,criteria = 200),
            pfail_rhat = pfail(rhat,criteria = 1.05,less_than = FALSE)) %>% 
  mutate(species = species)

fail_summary <- bind_rows(fail_summary,fail_sum)
fail_params <- bind_rows(fail_params,failt)

write.csv(fail_params,"output/temp_failed_parameters.csv",
          row.names = FALSE)
write.csv(fail_summary,"output/temp_fail_summary.csv",
          row.names = FALSE)
print(paste(round(jj/length(species_to_run$english),3),species,"max_rhat =",max(failt$rhat,na.rm = TRUE)))
}


max_fail <- fail_summary %>% 
  group_by(species) %>% 
  summarise(max_rhat = max(pfail_rhat),
            max_ess = max(pfail_ess_bulk)) %>% 
  arrange(-max_rhat)
re_run_species <- max_fail %>% 
  filter(max_rhat > 0)
saveRDS(re_run_species,"temp_species_to_rerun.rds")





# extracting SD parameter estimates ---------------------------------------



sd_year<- NULL
sd_other <- NULL


for(jj in rev(1:nrow(species_to_run))){
  
  species <- as.character(species_to_run[jj,"english"])
  species_f <- as.character(species_to_run[jj,"species_file"])

  if(fit_spatial){
    
    out_base <- paste(species_f,model_sel,"Spatial","BBS",sep = "_") # text string to identify the saved output from the Stan process unique to species and model, but probably something the user wants to control
    
  }else{
    out_base <- paste(species_f,model_sel,"BBS",sep = "_")
    if(Non_hierarchical){
      out_base <- paste(species_f,model_sel,"NonHier_BBS",sep = "_")
    }
  }
  
  if(!file.exists(paste0(output_dir,"/",out_base,"_Stan_fit.RData"))){next}
  
  load(paste0(output_dir,"/",out_base,"_Stan_fit.RData"))
  
  sd_tmp <- fit_summary %>% 
    filter(grepl("sd",variable)) %>% 
    mutate(species = species)
  
  sd_year_tmp <- sd_tmp %>% 
    filter(grepl("year",variable))
  
  sd_other_tmp <- sd_tmp %>% 
    filter(!grepl("year",variable))
  
  sd_year <- bind_rows(sd_year,sd_year_tmp)
  sd_other <- bind_rows(sd_other,sd_other_tmp)
  

  write.csv(sd_year,"output/temp_sd_year.csv",
            row.names = FALSE)
  write.csv(sd_other,"output/temp_sd_other.csv",
            row.names = FALSE)
  
}
