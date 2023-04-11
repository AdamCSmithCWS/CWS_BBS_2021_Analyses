### bbsBayes2 post fit plotting with bbsBayes2

library(bbsBayes2)
library(tidyverse)
library(patchwork)

load("species_lists.RData")
model_sel = "gamye"
strat = "bbs_cws"
output_dir <- "F:/bbsStanBayes/output/" # Stan writes output to files as it samples. This is great because it's really stable, but the user needs to think about where to store that output
output_dir <- "output/"
species_to_run <- nrecs_sp

for(jj in 1:nrow(species_to_run)){

  # prep data with bbsBayes2 ------------------------------------------------
species <- as.character(species_to_run[jj,"english"])

species_f <- as.character(species_to_run[jj,"species_file"])
start_year <- NULL

if(species %in% splitters){start_year <- split_miny[species]}

s <- stratify(by = strat, species = species)
p <- prepare_data(s,
                  min_year = start_year,
                  min_max_route_years = 2)
md <- prepare_model(p,
              model = model_sel)

m <- md

# load the fitted model from 2021 CWS -------------------------------------

out_base <- paste(species_f,model_sel,"BBS",sep = "_")

load(paste0(output_dir,"/",out_base,"_Stan_fit.RData"))

# combine the fitted model with prepped data

#tmp <- stanfit$summary()

m[["model_fit"]] <- stanfit

saveRDS(m,file = paste0("output/saved_bbsBayes2_",species_f,".rds"))  

}



