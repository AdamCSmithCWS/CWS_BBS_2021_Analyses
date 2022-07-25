## plotting

library(bbsBayes)
library(tidyverse)
library(cmdstanr)

species <- "Golden-winged Warbler"


species_f <- gsub(species,pattern = " ",replacement = "_")
model <- "gamye"
out_base <- paste(species_f,model,"BBS",sep = "_")
output_dir <- "output"

load(paste0(output_dir,"/",out_base,"_fit.RData"))


inds <- generate_indices(jags_mod = stanfit,
                         jags_data = stan_data,
                         backend = "Stan",
                         stratify_by = "bbs_usgs",
                         alternate_n = "nsmooth")


trajs <- plot_indices(inds,add_observed_means = TRUE,add_number_routes = TRUE)
pdf(paste0("Figures/",species_f,"_temp_trajectories.pdf"),
    width = 8,
    height = 8)
for(i in 1:length(trajs)){
  print(trajs[[i]])
}
dev.off()


trends <- generate_trends(inds)
trends_short <- generate_trends(inds,Min_year = 2004)

map <- generate_map(trends,select = TRUE,stratify_by = "bbs_usgs",species = species)
print(map)
mapshort <- generate_map(trends_short,select = TRUE,stratify_by = "bbs_usgs",species = species)
print(mapshort)


## doesn't yet work
gfacet <- geofacet_plot(inds,select = TRUE,trends = trends,add_observed_means = TRUE,stratify_by = "bbs_usgs")


