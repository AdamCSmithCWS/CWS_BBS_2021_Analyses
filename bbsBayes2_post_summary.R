### bbsBayes2 post fit plotting with bbsBayes2

library(bbsBayes2)
setwd("C:/GitHub/bbsStanBayes")
load("species_lists.RData")





model_sel = "gamye"
strat = "bbs_cws"
output_dir <- "F:/bbsStanBayes/output/" # Stan writes output to files as it samples. This is great because it's really stable, but the user needs to think about where to store that output
output_dir <- "output/"
species_to_run <- nrecs_sp
sps <- c("Common Yellowthroat",
         "Song Sparrow",
         "Downy Woodpecker",
         "Common Nighthawk",
         "Savannah Sparrow",
         "Black-throated Green Warbler",
         "American Robin")
wtemp <- which(species_to_run$english %in% sps)

for(jj in wtemp[6:7]){ #:nrow(species_to_run)){
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







# confirm that bbsBayes2 functions work -----------------------------------



species <- "Song Sparrow"
jj <- which(nrecs_sp$english == species)
species_f <- as.character(nrecs_sp[jj,"species_file"])

m <- readRDS(paste0("output/saved_bbsBayes2_",species_f,".rds"))


ind <- generate_indices(m,
                        regions = "continent")
trajs <- plot_indices(ind)
print(trajs[["continent"]])

inds <- generate_indices(m,
                         alternate_n = "nsmooth")

starts <- c(1966,
            1978,
            1990,
            2010,
            2021)

maps <- vector(mode = "list",
               length = length(starts)-1)
for(y in 1:(length(starts)-1)){
  fy <- starts[y]
  ly <- starts[y+1]

tr <- generate_trends(inds,
                      min_year = fy,
                      max_year = ly)
maps[[y]] <- plot_map(tr)
}

