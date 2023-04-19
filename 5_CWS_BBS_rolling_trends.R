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


# regions to estimate -----------------------------------------------------

regs_to_estimate <- c("continental","national","prov_state","bcr")



load("species_lists.RData") # loads objects created at the beginning of the script Fit_gamye_models_cws.R

species_to_run <- nrecs_sp 


# 
# 
# sp_tmp <- c( "Brown-headed Cowbird",
#              "Clark's Nutcracker",
#              "Western Wood-Pewee",
#              "Eastern Towhee",
#              "Common Grackle",
#              "Song Sparrow",
#              "Townsend's Solitaire",
#              "Black-throated Green Warbler",
#              "Common Yellowthroat",
#              "American Robin",
#              "Hermit Thrush",
#              "American Wigeon",
#              "Bufflehead",
#              "Broad-tailed Hummingbird",
#              "Common Nighthawk",
#              "Downy Woodpecker",
#              "Acorn Woodpecker",
#              "Spotted Sandpiper",
#              "Wood Stork",
#              "Horned Grebe",
#              "Double-crested Cormorant",
#              "House Wren",
#              "Rock Wren",
#              "Townsend's Warbler",
#              "Savannah Sparrow",
#              "American Robin",
#              "Barn Swallow",
#              "House Sparrow",
#              "American Crow",
#              "Carolina Wren",
#              "Greater Roadrunner")
# 
# species_to_run <- filter(species_to_run,
#                          english %in% sp_tmp)
# 


# Species loop to save indices and generate trends ------------------------------------------------------------
n_cores <- 15#length(provs)
cluster <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cluster)


fullrun <- foreach(jj = rev(1:nrow(species_to_run)),
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
    source("functions/web_trends.R")
    source("functions/reliability.R")
    source("functions/generate-map-alt.R")
    

    
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
        
        
        
        # rolling trends ----------------------------------------------------------
        
        
        starts <- c(seq(fy,(YYYY-short_time),by = 1))
        
        roll_trends_out <- NULL
        pdf(paste0("trends/rolling_trend_maps/",species_f_bil,"_rolling_trend_map.pdf"))
        for(dd in starts){
          trends_10temp <- generate_trends(inds,
                                           Min_year = dd,
                                           Max_year = dd+10,
                                           prob_decrease = c(0,25,30,50),
                                           prob_increase = c(0,33,100))
          roll_trends_out <- bind_rows(roll_trends_out,trends_10temp)
          
          
          map_tmp <- generate_map(trends_10temp,select = TRUE,
                                  stratify_by = strat_sel,
                                  species = paste(species,espece),
                                  exclude_backcast = TRUE,
                                  alpha_backcast = FALSE)
         
          print(map_tmp)
        }
        dev.off()
        
        roll_trends_out <- roll_trends_out %>% 
          mutate(species = species,
                 espece = espece,
                 bbs_num = aou,
                 Trend_Time = "rolling-short")
        readr::write_excel_csv(roll_trends_out,file = paste0("trends/rolling_trend_estimates/",species_f_bil,"_rolling_trends.csv"))
        
        roll_trends_out <- roll_trends_out %>% 
          filter(Region_type %in% c("continental","national","prov_state"))
        
        thresh30 = (0.7^(1/short_time)-1)*100
        thresh50 = (0.5^(1/short_time)-1)*100
        
        
        pdf(paste0("trends/rolling_trend_graphs/",species_f_bil,"_rolling_trends.pdf"),
            width = 11,
            height = 8.5)
        
        for(rr in unique(roll_trends_out$Region_alt)){
          rttmp <- roll_trends_out %>% 
            filter(Region_alt == rr)
          rollTrend <- rttmp %>% 
            filter(End_year == YYYY) %>% 
            select(Trend,prob_decrease_30_percent,prob_decrease_50_percent)
          pth_30_labs <- paste0(round(rollTrend[,"prob_decrease_30_percent"],2)*100,
                             "% probability of 30% decrease")
          pth_50_labs <- paste0(round(rollTrend[,"prob_decrease_50_percent"],2)*100,
                               "% probability of 50% decrease")
          rtp <- ggplot(data = rttmp,
                      aes(y = Trend,x = End_year))+
            geom_errorbar(aes(ymin = Trend_Q0.025,
                              ymax = Trend_Q0.975),
                          width = 0,alpha = 0.2)+
            geom_errorbar(aes(ymin = Trend_Q0.25,
                              ymax = Trend_Q0.75),
                          width = 0,alpha = 0.6)+
            geom_point(aes(alpha = backcast_flag))+
            scale_alpha_continuous(range = c(0.1,1))+
            guides(alpha = "none")+
            geom_hline(yintercept = thresh30,colour = "darkorange")+
            geom_hline(yintercept = thresh50,colour = "darkred")+
            geom_hline(yintercept = 0)+
            labs(title = paste(species,espece,"rolling",short_time,"year trends",rr),
                 subtitle = paste("Based on trend in",YYYY,":",pth_30_labs,"and",pth_50_labs))+
            xlab(paste("Ending year of",short_time,"trend"))+
            ylab(paste(short_time,"year trends"))+
            theme_bw()
            
        print(rtp)
  
        
        } #temp end loop

        dev.off()
        
        
  }

stopCluster(cl = cluster)




