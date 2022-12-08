## plotting

### NOTE: this script depends on a development version of bbsBayes at the following repo: https://github.com/AdamCSmithCWS/bbsBayes/tree/testing_Stan 
## the specific changes in that branch are in the two functions
## generate_indices() and extract_index_data() (the second is called within generate_indices())
source("functions/posterior_summary_functions.R")
source("functions/generate-indices-alt.R")


setwd("C:/Users/SmithAC/Documents/GitHub/bbsStanBayes")

strat_sel <- "bbs_cws"


library(bbsBayes)
library(tidyverse)
library(cmdstanr)
library(patchwork)


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


# pdf(paste0("Figures/temp_trajectories2.pdf"),
#     width = 11,
#     height = 8.5)

jj = 131 #connecticut warbler
regs_to_estimate <- c("stratum","continental","national","prov_state","bcr")

for(jj in 1:nrow(species_to_run)){
  
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

if(model_sel == "gamye"){alt_n <- "nsmooth"}
if(model_sel == "slope"){alt_n <- "nslope"}

if(model_sel %in% c("firstdiff","gam")){alt_n <- NA}
if(is.null(stan_data$strat_name)){
  stan_data$strat_name <- stan_data$alt_data$strat_name
}
ind <- generate_indices(jags_mod = stanfit,
                        jags_data = stan_data,
                        backend = "Stan",
                        stratify_by = strat_sel,
                        alternate_n = "n",
                        regions = regs_to_estimate)

# extract_observer_route_effects <- function(jags_mod = stan_fit,
#                                            jags_dat = stan_data){
#   
#   
#   
# }


if(!is.na(alt_n)){
  
inds <- generate_indices(jags_mod = stanfit,
                         jags_data = stan_data,
                         backend = "Stan",
                         stratify_by = strat_sel,
                         alternate_n = alt_n,
                         regions = regs_to_estimate)
}else{
  inds <- ind
}

## add option to choose which observed means get added and whether to add the mean effects
## add option to choose which observed means get added and whether to add the mean effects
## add option to choose which observed means get added and whether to add the mean effects
## add option to choose which observed means get added and whether to add the mean effects
## add option to choose which observed means get added and whether to add the mean effects
## add option to choose which observed means get added and whether to add the mean effects
## add option to choose which observed means get added and whether to add the mean effects
## add option to choose which observed means get added and whether to add the mean effects




trajs <- plot_indices(inds,
                      species = species,
                      add_observed_means = TRUE,
                      add_number_routes = TRUE)

trajshort <- plot_indices(inds,
                      species = species,
                      add_observed_means = TRUE,
                      add_number_routes = TRUE,
                      min_year = 2004)

pdf(file = paste0("Figures/",out_base,".pdf"),width = 11,height = 8.5)
for(i in names(trajs)){
  t1 <- trajs[[i]] 
  t2 <- trajshort[[i]]
  
  # st <- str_trunc(t1$labels$title, width = 8,
  #                 side = "left",
  #                 ellipsis = "")
  if(!is.na(alt_n)){
  n1 <- ind$data_summary %>% 
    mutate(Reg_traj = gsub(Region_alt,pattern = "[[:punct:]]",replacement = "")) %>% 
    filter(Reg_traj == gsub(i,pattern = "[[:punct:]]",replacement = "")) #,
  #Region_type == "stratum"
  t1 <- t1 +
    geom_ribbon(data = n1, aes(x = Year,y = Index,ymin = Index_q_0.025,ymax = Index_q_0.975),
                fill = grey(0.5),alpha = 0.2)+
    geom_line(data = n1, aes(x = Year,y = Index),
              colour = grey(0.5))
  
  n2 <- ind$data_summary %>% 
    filter(Region_alt == gsub(i,pattern = "_",replacement = "-"),
           #Region_type == "stratum",
           Year >= 2004)
  t2 <- t2 +
    geom_ribbon(data = n2, aes(x = Year,y = Index,ymin = Index_q_0.025,ymax = Index_q_0.975),
                fill = grey(0.5),alpha = 0.2)+
    geom_line(data = n2, aes(x = Year,y = Index),
              colour = grey(0.5))
  }
  print(t1 + t2)
  
}



trends <- generate_trends(inds)
trends_short <- generate_trends(inds,Min_year = 2011)
map <- generate_map(trends,select = TRUE,stratify_by = strat_sel,species = species)
mapshort <- generate_map(trends_short,select = TRUE,stratify_by = strat_sel,species = species)
print(map + mapshort)

starts <- c(1970,seq(1976,2011,by = 5))
maps <- vector("list",length = length(starts))
names(maps) <- paste(starts)

for(dd in starts){
trends_10temp <- generate_trends(inds,Min_year = dd,Max_year = dd+10)
maps[[paste(dd)]] <- generate_map(trends_10temp,select = TRUE,stratify_by = strat_sel,species = species)
print(maps[[paste(dd)]])
}

dev.off()


}





## doesn't yet work
#gfacet <- geofacet_plot(inds,select = TRUE,trends = trends,add_observed_means = TRUE,stratify_by = strat_sel)


