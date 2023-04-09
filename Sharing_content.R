setwd("C:/GitHub/bbsStanBayes")


library(bbsBayes2)
library(tidyverse)
library(patchwork)
load("species_lists.RData")





species <- "Wood Thrush"
jj <- which(nrecs_sp$english == species)
species_f <- as.character(nrecs_sp[jj,"species_file"])

m <- readRDS(paste0("output/saved_bbsBayes2_",species_f,".rds"))


ind <- generate_indices(m,
                        regions = "continent")
trajs <- plot_indices(ind,
                      title = FALSE)
print(trajs[["continent"]])

inds <- generate_indices(m,
                         alternate_n = "nsmooth")

starts <- c(1966,
            1977,
            1990,
            2011,
            2021) #identify some change points


# CAlculate trends between change points
# map trends and pair with continental trajectory
# mark start and end points on continental trajectory
# create shareable sized graphics (no names)

#create combined graphic with name? or two versions of graphics (with and without names)

maps <- vector(mode = "list",
               length = length(starts)-1)
for(y in 1:(length(starts)-1)){
  fy <- starts[y]
  ly <- starts[y+1]
  
  tr <- generate_trends(inds,
                        min_year = fy,
                        max_year = ly)
  map <- plot_map(tr,
                  title = FALSE)
  indt <- ind$indices %>% 
    filter(year %in% c(fy,ly)) %>% 
    mutate(yup = index_q_0.025,
           ydown = 0)
  traj <- trajs[[1]] +
    geom_errorbar(data = indt,
              aes(x = year, ymin = ydown,ymax = yup),
              colour = grey(0.5),
              linewidth = 0.5,
              width = 0)
    
  mapt <- traj + map + plot_layout(nrow = 2,heights = c(2,4))
  mag = 3
  png(filename = paste0("Sharing/",species_f,"_",fy,"-",ly,".png"),
      res = 72 *mag,
      height = 1000*mag,
      width = 600*mag)
  print(mapt)
  dev.off()
  maps[[y]] <- mapt
  
}

