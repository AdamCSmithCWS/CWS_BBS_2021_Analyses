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

setwd("C:/Users/SmithAC/Documents/GitHub/bbsStanBayes")

source("functions/posterior_summary_functions.R")
source("functions/generate-indices-alt.R")
source("functions/plot-indices-alt.R")
source("functions/extract-index-data-alt.R")
source("functions/generate-trends-alt.R")
source("functions/web_trends.R")
source("functions/reliability.R")


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

covs <- covs %>% 
  mutate(trendtype = stringr::str_to_sentence(trendtype))


# regions to estimate -----------------------------------------------------

regs_to_estimate <- c("continental","national","prov_state","bcr","stratum")



load("species_lists.RData") # loads objects created at the beginning of the script Fit_gamye_models_cws.R

species_to_run <- nrecs_sp 

CV_threshold <- function(m,ci,thresh = 100){
  y <- ifelse(ci/m > thresh,TRUE,FALSE)
  return(y)
}




# reliability category definitions ----------------------------------------

prec_cuts = c(abs(2*((0.7^(1/20))-1)),
              abs(2*((0.5^(1/20))-1)))*100 
names(prec_cuts) <- c("High","Medium")

cov_cuts = c(0.5,0.25)
names(cov_cuts) <- c("High","Medium")

pool_cuts = c(0.33,0.1)
names(pool_cuts) <- c("High","Medium")

backcast_cuts = c(0.90,0.75)
names(backcast_cuts) <- c("High","Medium")
















# Species loop to save indices and generate trends ------------------------------------------------------------
n_cores <- 20#length(provs)
cluster <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cluster)


fullrun <- foreach(jj = (1:nrow(species_to_run)),
                   .packages = c("posterior","tidyverse","bbsBayes",
                                 "patchwork"),
                   .inorder = FALSE,
                   .errorhandling = "pass") %dopar%
  {
    
    
   # for(jj in rev(1:nrow(species_to_run))[21:31]){
    
    
    source("functions/posterior_summary_functions.R")
    source("functions/generate-indices-alt.R")
    source("functions/plot-indices-alt.R")
    source("functions/extract-index-data-alt.R")
    source("functions/generate-trends-alt.R")
    source("functions/generate-map-alt.R")
    source("functions/web_trends.R")
    source("functions/reliability.R")

    

    
    species <- as.character(species_to_run[jj,"english"])
    species_f <- as.character(species_to_run[jj,"species_file"])
    espece <- as.character(species_to_run[jj,"french"])
    species_f_bil <- gsub(paste(species,espece),pattern = "[[:space:]]|[[:punct:]]",
                                          replacement = "_")
    aou <- as.character(species_to_run[jj,"AOU"])
    
    
    if(fit_spatial){
      
      out_base <- paste(species_f,model_sel,"Spatial","BBS",sep = "_") # text string to identify the saved output from the Stan process unique to species and model, but probably something the user wants to control
      
    }else{
      out_base <- paste(species_f,model_sel,"BBS",sep = "_")
      if(Non_hierarchical){
        out_base <- paste(species_f,model_sel,"NonHier_BBS",sep = "_")
      }
    }
    
    
  
    
    if(!file.exists(paste0(output_dir,"/",out_base,"_Stan_fit.RData"))){next}
    

    if(model_sel == "gamye"){alt_n <- "nsmooth"}
    if(model_sel == "slope"){alt_n <- "nslope"}
    

    
    # load indices --------------------------------------------------------
    
    
    ind <- read_rds(paste0(output_dir,"/temp_rds_storage/",species_f,"_indices.RDS"))
    inds <- read_rds(paste0(output_dir,"/temp_rds_storage/",species_f,"_indices_smooth.RDS"))
    
    
   

        fy <- max(1970,min(inds$data_summary$Year)) 
        if(species %in% c("Alder Flycatcher","Willow Flycatcher")){
          fy <- 1978 #5 years after the split 
        }
        if(species %in% c("Clark's Grebe","Western Grebe","Eurasian Collared-Dove")){
          fy <- 1990 #5 years after the split and first year EUCD observed on > 3 BBS routes
        }
        if(species == "Cave Swallow"){
          fy = 1985
        }
        
        
        
        
        trends_70 <- generate_trends(inds,
                                     Min_year = fy,
                                     quantiles = c(0.025, 0.05, 0.10, 0.25, 0.75, 0.9, 0.95, 0.975),
                                     prob_decrease = c(0,25,30,50),
                                     prob_increase = c(0,33,100))
        
        trends_short <- generate_trends(inds,
                                        Min_year = (YYYY-short_time),
                                        quantiles = c(0.025, 0.05, 0.10, 0.25, 0.75, 0.9, 0.95, 0.975),
                                        prob_decrease = c(0,25,30,50),
                                        prob_increase = c(0,33,100))
        
        
        tl <- trends_70 %>% 
          mutate(species = species,
                 espece = espece,
                 bbs_num = as.integer(aou),
                 Trend_Time = "Long-term",
                 For_web = for_web_func(Strata_included))
        
        
        
        ts <- trends_short %>% 
          mutate(species = species,
                 espece = espece,
                 bbs_num = as.integer(aou),
                 Trend_Time = "Short-term",
                 For_web = for_web_func(Strata_included))
        
        
       
        
        
        # insert last year's coverage estimates -----------------------------------
        
        covsp = covs %>% 
          filter(sp == aou) %>% 
          select(sp,new.area,trendtype,reliab.cov)
        

        
        # insert reliability categories---------------------------------------------------
        
        tall <- bind_rows(tl,ts) %>% 
          left_join(.,covsp,by = c("bbs_num" = "sp",
                                   "Region_alt" = "new.area",
                                   "Trend_Time" = "trendtype")) %>% 
          mutate(precision = reliab_func_prec(Width_of_95_percent_Credible_Interval),
                 coverage = reliab_func_cov(reliab.cov),
                 backcast_reliab = reliab_func_backcast(backcast_flag),
                 reliability = reliability_func(precision,coverage,backcast_reliab)) %>% 
          mutate(across(where(is.numeric),~signif(.,3)))
        
        
  
        write.csv(tall,file = paste0("trends/Trends_by_species/",species_f_bil,"_trends.csv"),
                  row.names = FALSE)
        
        
        
        map <- generate_map(trends_70,select = TRUE,stratify_by = strat_sel,
                            species = paste(species,espece))

        
        
        map_25 <- generate_map(trends_70,select = TRUE,stratify_by = strat_sel,
                            species = paste("25 Percentile trends",species,espece),
                            heat_column = "Trend_Q0.25")
        map_75 <- generate_map(trends_70,select = TRUE,stratify_by = strat_sel,
                               species = paste("75 Percentile trends",species,espece),
                               heat_column = "Trend_Q0.75")
        
        # map_10 <- generate_map(trends_70,select = TRUE,stratify_by = strat_sel,
        #                        species = paste("10 Percentile trends",species,espece),
        #                        heat_column = "Trend_Q0.1")
        # map_90 <- generate_map(trends_70,select = TRUE,stratify_by = strat_sel,
        #                        species = paste("90 Percentile trends",species,espece),
        #                        heat_column = "Trend_Q0.9")
        
        mapshort <- generate_map(trends_short,select = TRUE,stratify_by = strat_sel,
                                 species = paste(species,espece))

        mapshort_25 <- generate_map(trends_short,select = TRUE,stratify_by = strat_sel,
                               species = paste("25 Percentile trends",species,espece),
                               heat_column = "Trend_Q0.25")
        mapshort_75 <- generate_map(trends_short,select = TRUE,stratify_by = strat_sel,
                               species = paste("75 Percentile trends",species,espece),
                               heat_column = "Trend_Q0.75")
        
        # mapshort_10 <- generate_map(trends_short,select = TRUE,stratify_by = strat_sel,
        #                             species = paste("10 Percentile trends",species,espece),
        #                             heat_column = "Trend_Q0.1")
        # mapshort_90 <- generate_map(trends_short,select = TRUE,stratify_by = strat_sel,
        #                             species = paste("90 Percentile trends",species,espece),
        #                             heat_column = "Trend_Q0.9")
        
        pdf(file = paste0("trends/Trend_maps/",species_f_bil,"_trend_maps.pdf"),
            width = 11,
            height = 8.5)
        print(map)
        print(map_25 / map_75)
        print(mapshort)
        print(mapshort_25 / mapshort_75)
        dev.off()
    
        traj_out <- vector("list",2)
        names(traj_out) <- c("Long-term","Short-term")
          traj_out[[1]] <- map
          traj_out[[2]] <- mapshort
          
  
        saveRDS(traj_out,file = paste0(output_dir,"/temp_rds_storage/",species_f,"_trend_maps.RDS"))
        
        

# export indices ----------------------------------------------------------

        ind_out <- ind$data_summary %>% 
          filter(Year >= fy) %>% 
          mutate(Trend_Time = "Long-term",
                 For_web = for_web_func(Strata_included),
                 species = species,
                 espece = espece,
                 bbs_num = as.integer(aou)) %>% 
          mutate(across(where(is.numeric),~signif(.,3)))
        
        write.csv(ind_out,file = paste0("indices/full/",species_f_bil,"_annual_indices.csv"),
                  row.names = FALSE)
        
        inds_out <- inds$data_summary %>% 
          filter(Year >= fy) %>% 
          mutate(Trend_Time = "Long-term",
                 For_web = for_web_func(Strata_included),
                 species = species,
                 espece = espece,
                 bbs_num = as.integer(aou),
                 index_type = "Smooth") %>% 
          mutate(across(where(is.numeric),~signif(.,3)))
        
        write.csv(inds_out,file = paste0("indices/smooth/",species_f_bil,"_smoothed_annual_indices.csv"),
                  row.names = FALSE)
        
  } #temp end loop


stopCluster(cl = cluster)













# All species high-level plots --------------------------------------------

species_to_run <- species_to_run %>% 
  arrange(-AOU)

pdf(file = paste0("Figures/BBS_High_level_summary_",YYYY,".pdf"),
    height = 9,
    width = 17)

for(jj in (1:nrow(species_to_run))){
  
  
  species <- as.character(species_to_run[jj,"english"])
  species_f <- as.character(species_to_run[jj,"species_file"])
  espece <- as.character(species_to_run[jj,"french"])
  species_f_bil <- gsub(paste(species,espece),pattern = "[[:space:]]|[[:punct:]]",
                        replacement = "_")
  aou <- as.character(species_to_run[jj,"AOU"])
  

tmaps <- readRDS(paste0(output_dir,"/temp_rds_storage/",species_f,"_trend_maps.RDS"))

trajs <- readRDS(paste0(output_dir,"/temp_rds_storage/",species_f,"_highlevel_trajs.RDS"))

design <- "
124
135
"

layt <- trajs[[1]] + trajs[[2]] + trajs[[3]] +
  tmaps[[1]] + tmaps[[2]] +
  plot_layout(design = design,
              guides = "collect") 

print(layt)

print(species)

}

dev.off()

