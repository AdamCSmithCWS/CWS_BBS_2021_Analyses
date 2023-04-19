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

regs_to_estimate <- c("continental","national","prov_state","bcr","stratum")



load("species_lists.RData") # loads objects created at the beginning of the script Fit_gamye_models_cws.R

species_to_run <- nrecs_sp 

CV_threshold <- function(m,ci,thresh = 100){
  y <- ifelse(ci/m > thresh,TRUE,FALSE)
  return(y)
}

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
# # 
# species_to_run <- filter(species_to_run,
#                          english %in% sp_tmp)
# 
# 










# Species loop to plot diagnostic trajectories ------------------------------------------------------------
n_cores <- 15#
cluster <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cluster)


fullrun <- foreach(jj = rev(1:nrow(species_to_run)),
                   .packages = c("posterior","tidyverse","bbsBayes",
                                 "patchwork"),
                   .inorder = FALSE,
                   .errorhandling = "pass") %dopar%
  {
    

#for(jj in rev(1:nrow(species_to_run))){
  
  
  source("functions/posterior_summary_functions.R")
  source("functions/generate-indices-alt.R")
  source("functions/plot-indices-alt.R")
  source("functions/extract-index-data-alt.R")
  source("functions/generate-trends-alt.R")
  
  
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

  

# load model fit ----------------------------------------------------------

  
if(!file.exists(paste0(output_dir,"/",out_base,"_Stan_fit.RData"))){next}

# load(paste0(output_dir,"/",out_base,"_Stan_fit.RData"))
# 
if(model_sel == "gamye"){alt_n <- "nsmooth"}
if(model_sel == "slope"){alt_n <- "nslope"}

# if(model_sel %in% c("firstdiff","gam")){alt_n <- NA}
# if(is.null(stan_data$strat_name)){
#   stan_data$strat_name <- stan_data$alt_data$strat_name
# }


# load indices --------------------------------------------------------


ind <- read_rds(paste0(output_dir,"/temp_rds_storage/",species_f,"_indices.RDS"))
inds <- read_rds(paste0(output_dir,"/temp_rds_storage/",species_f,"_indices_smooth.RDS"))


# species diagnostic trajectory plots -------------------------------------


trajs <- plot_indices(inds,
                      species = paste(species,espece),
                      ci_width = 0.9,
                      add_observed_means = TRUE,
                      add_number_routes = TRUE,
                      add_adjusted_means = TRUE,
                      add_extras = TRUE,
                      title_size = 12,
                      axis_title_size = 10,
                      axis_text_size = 10)

trajshort <- plot_indices(inds,
                          species = paste(species,espece),
                          add_observed_means = TRUE,
                          add_number_routes = TRUE,
                          min_year = ((YYYY-short_time)-5), # add five years to the short-term trajectories
                          ci_width = 0.9,
                          add_adjusted_means = TRUE,
                          add_extras = TRUE,
                          title_size = 12,
                          axis_title_size = 10,
                          axis_text_size = 10)



pdf(file = paste0("Figures/diagnostic_trajectories/",species_f_bil,"_diagnostic_trajectories.pdf"),width = 11,height = 8.5)

for(i in names(trajs)){
  t1 <- trajs[[i]] 
  t2 <- trajshort[[i]]
  
  # st <- str_trunc(t1$labels$title, width = 8,
  #                 side = "left",
  #                 ellipsis = "")
  if(!is.na(alt_n)){
    n1 <- ind$data_summary %>% 
      mutate(Reg_traj = gsub(Region_alt,pattern = "[[:punct:]]|[[:space:]]",replacement = "")) %>% 
      filter(Reg_traj == gsub(i,pattern = "[[:punct:]]|[[:space:]]",replacement = "")) #,
    #Region_type == "stratum"
    t1 <- t1 +
      geom_ribbon(data = n1, aes(x = Year,y = Index,ymin = Index_q_0.05,ymax = Index_q_0.95),
                  fill = grey(0.5),alpha = 0.2)+
      geom_line(data = n1, aes(x = Year,y = Index),
                colour = grey(0.5),)
    
    n2 <- ind$data_summary %>% 
      mutate(Reg_traj = gsub(Region_alt,pattern = "[[:punct:]]|[[:space:]]",replacement = "")) %>% 
      filter(Reg_traj == gsub(i,pattern = "[[:punct:]]|[[:space:]]",replacement = ""),
             #Region_type == "stratum",
             Year >= ((YYYY-short_time)-5))
    t2 <- t2 +
      geom_ribbon(data = n2, aes(x = Year,y = Index,ymin = Index_q_0.05,ymax = Index_q_0.95),
                  fill = grey(0.5),alpha = 0.2)+
      geom_line(data = n2, aes(x = Year,y = Index),
                colour = grey(0.5))
  }
  print(t1 + t2)
  

}

dev.off()  # close diagnostic trajectory plotting


# 
# saveRDS(ind,file = paste0(output_dir,"/temp_rds_storage/",species_f,"_indices.RDS"))
# saveRDS(inds,file = paste0(output_dir,"/temp_rds_storage/",species_f,"_indices_smooth.RDS"))


} #temp end loop

# 
# stopCluster(cl = cluster)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Species loop to plot simple trajectories ------------------------------------------------------------
n_cores <- 15#length(provs)
cluster <- makeCluster(n_cores, type = "PSOCK")
registerDoParallel(cluster)


fullrun <- foreach(jj = rev(1:nrow(species_to_run)),
                   .packages = c("posterior","tidyverse","bbsBayes",
                                 "patchwork"),
                   .inorder = FALSE,
                   .errorhandling = "pass") %dopar%
  {
    
    
    #for(jj in rev(1:nrow(species_to_run))){
    
    
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
    
    
    
    # load model fit ----------------------------------------------------------
    
    
    if(!file.exists(paste0(output_dir,"/",out_base,"_Stan_fit.RData"))){next}
    
    #load(paste0(output_dir,"/",out_base,"_Stan_fit.RData"))
    
    if(model_sel == "gamye"){alt_n <- "nsmooth"}
    if(model_sel == "slope"){alt_n <- "nslope"}
    
    # if(model_sel %in% c("firstdiff","gam")){alt_n <- NA}
    # if(is.null(stan_data$strat_name)){
    #   stan_data$strat_name <- stan_data$alt_data$strat_name
    # }
    # 
    
    # load indices --------------------------------------------------------
    
    
    ind <- read_rds(paste0(output_dir,"/temp_rds_storage/",species_f,"_indices.RDS"))
    inds <- read_rds(paste0(output_dir,"/temp_rds_storage/",species_f,"_indices_smooth.RDS"))
    
    
   

        fy <- 1970 
        if(species %in% c("Alder Flycatcher","Willow Flycatcher")){
          fy <- 1978 #5 years after the split 
        }
        if(species %in% c("Clark's Grebe","Western Grebe","Eurasian Collared-Dove")){
          fy <- 1990 #5 years after the split and first year EUCD observed on > 3 BBS routes
        }
        if(species == "Cave Swallow"){
          fy = 1985
        }
        
     
      
    
    trajs <- plot_indices(ind,
                          species = paste(species,espece, sep = " / "),
                          ci_width = 0.9,
                          min_year = fy, 
                          title_size = 10,
                          axis_title_size = 8,
                          axis_text_size = 8)
    
    trajshort <- plot_indices(ind,
                              species = "",
                              min_year = ((YYYY-short_time)-5), # add five years to the short-term trajectories
                              ci_width = 0.9,
                              title_size = 10,
                              axis_title_size = 8,
                              axis_text_size = 8)
    
    
    traj_out <- vector("list",3)
    names(traj_out) <- c("Continental","Canada","United_States_of_America")
    
    pdf(file = paste0("Figures/trajectories/",species_f_bil,"_trajectories.pdf"),width = 11,height = 8.5)
    for(i in names(trajs)){
      t1 <- trajs[[i]] 
      t2 <- trajshort[[i]]
      
      # st <- str_trunc(t1$labels$title, width = 8,
      #                 side = "left",
      #                 ellipsis = "")
      if(!is.na(alt_n)){
        n1 <- inds$data_summary %>% 
          mutate(Reg_traj = gsub(Region_alt,pattern = "[[:punct:]]|[[:space:]]",replacement = "")) %>% 
          filter(Reg_traj == gsub(i,pattern = "[[:punct:]]|[[:space:]]",replacement = ""),
                 Year >= fy) #,
        #Region_type == "stratum"
        t1 <- t1 +
          geom_ribbon(data = n1, aes(x = Year,y = Index,ymin = Index_q_0.05,ymax = Index_q_0.95),
                      fill = grey(0.5),alpha = 0.2)+
          geom_line(data = n1, aes(x = Year,y = Index),
                    colour = grey(0.5))
        
        n2 <- inds$data_summary %>% 
          mutate(Reg_traj = gsub(Region_alt,pattern = "[[:punct:]]|[[:space:]]",replacement = "")) %>% 
          filter(Reg_traj == gsub(i,pattern = "[[:punct:]]|[[:space:]]",replacement = ""),
                 #Region_type == "stratum",
                 Year >= ((YYYY-short_time)-5))
        t2 <- t2 +
          geom_ribbon(data = n2, aes(x = Year,y = Index,ymin = Index_q_0.05,ymax = Index_q_0.95),
                      fill = grey(0.5),alpha = 0.2)+
          geom_line(data = n2, aes(x = Year,y = Index),
                    colour = grey(0.5))
      }
      print(t1 + t2)
      
      
      if(i %in% c("Continental","Canada","United_States_of_America")){
        traj_out[[i]] <- t1
      }
    }
    
    dev.off()  # 
    
    
    saveRDS(traj_out,file = paste0(output_dir,"/temp_rds_storage/",species_f,"_highlevel_trajs.RDS"))
    
    
    # 
    # saveRDS(ind,file = paste0(output_dir,"/temp_rds_storage/",species_f,"_indices.RDS"))
    # saveRDS(inds,file = paste0(output_dir,"/temp_rds_storage/",species_f,"_indices_smooth.RDS"))
    
    
  } #temp end loop


stopCluster(cl = cluster)




