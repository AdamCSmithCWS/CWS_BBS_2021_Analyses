## plotting

library(bbsBayes)
library(tidyverse)
library(cmdstanr)
library(patchwork)


load("species_lists.RData")
model = "gamye"


species_to_run <- nrecs_sp 

fit_spatial <- TRUE # TRUE = spatial sharing of information and FALSE = non-spatial sharing
source("Functions/prepare-data-alt.R")
if(fit_spatial){
  source("Functions/neighbours_define_alt.R") # function to generate spatial neighbourhoods to add to the spatial applications of the models
}
output_dir <- "output/" # Stan writes output to files as it samples. This is great because it's really stable, but the user needs to think about where to store that output


pdf(paste0("Figures/temp_trajectories.pdf"),
    width = 11,
    height = 8.5)

for(jj in 1:nrow(species_to_run)){
  
  species <- as.character(species_to_run[jj,"english"])
  species_f <- as.character(species_to_run[jj,"species_file"])
  

output_dir <- "output"




if(fit_spatial){
  
  out_base <- paste(species_f,model,"Spatial","BBS",sep = "_") # text string to identify the saved output from the Stan process unique to species and model, but probably something the user wants to control
  
}else{
  out_base <- paste(species_f,model,"BBS",sep = "_")
  
}

if(!file.exists(paste0(output_dir,"/",out_base,"_Stan_fit.RData"))){next}

load(paste0(output_dir,"/",out_base,"_Stan_fit.RData"))
explore_observers <- FALSE

if(explore_observers){
sdobs = fit_summary %>% filter(variable == "sdobs")

obs_effs <- fit_summary %>% 
  filter(grepl("obs_raw",variable)) %>% 
  mutate(effect_estimate = mean*sdobs$mean,
         effect_mult = exp(effect_estimate),
         observer = row_number()) %>% 
  select(observer,effect_estimate,effect_mult)

obs_sum <- stan_data$alt_data[[1]] %>% #original data frame
  group_by(ObsN,observer) %>% 
  summarise(mean_count = mean(count),
            log_mean_count = log(mean_count),
            n_surveys = n(),
            strats = median(strat),
            n_strat = length(unique(strat))) %>% 
  left_join(.,obs_effs,by = "observer")

obs_sum_p <- ggplot(data = obs_sum,aes(x = effect_mult,y = strats,colour = n_surveys))+
  geom_point(alpha = 0.4)+
  labs(title = species)+
  scale_colour_viridis_c()

print(obs_sum_p)


weird_obs <- obs_sum %>% filter(effect_mult > 20)
obs_weird <- stan_data$alt_data[[1]] %>% 
  filter(ObsN %in% weird_obs$ObsN)

}
inds <- generate_indices(jags_mod = stanfit,
                         jags_data = stan_data,
                         backend = "Stan",
                         stratify_by = "bbs_usgs",
                         alternate_n = "nsmooth")


trajs <- plot_indices(inds,
                      species = species,
                      add_observed_means = TRUE,
                      add_number_routes = TRUE)

trajshort <- plot_indices(inds,
                      species = species,
                      add_observed_means = TRUE,
                      add_number_routes = TRUE,
                      min_year = 2004)

#pdf(file = "temp.pdf",width = 11,height = 8.5)
for(i in 1:length(trajs)){
  print(trajs[[i]] + trajshort[[i]])
}


trends <- generate_trends(inds)
trends_short <- generate_trends(inds,Min_year = 2004)

map <- generate_map(trends,select = TRUE,stratify_by = "bbs_usgs",species = species)
mapshort <- generate_map(trends_short,select = TRUE,stratify_by = "bbs_usgs",species = species)
print(map + mapshort)


}

dev.off()





## doesn't yet work
#gfacet <- geofacet_plot(inds,select = TRUE,trends = trends,add_observed_means = TRUE,stratify_by = "bbs_usgs")


