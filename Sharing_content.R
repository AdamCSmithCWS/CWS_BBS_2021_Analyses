setwd("C:/GitHub/bbsStanBayes")


library(bbsBayes2)
library(tidyverse)
library(patchwork)
library(sf)
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

base <- bbsBayes2::load_map("bbs_cws") %>% 
  filter(strata_name %in% inds$indices$region) %>% 
  st_bbox()

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
                  title = FALSE)+
    coord_sf(xlim = c(base["xmin"],base["xmax"]),
                    ylim = c(base["ymin"],base["ymax"]))+
    labs(title = paste("Trends",fy,"-",ly),
         subtitle = "North American Breeding Bird Survey")+
    theme(title = element_text(size = 20))
  
  indt <- ind$indices %>% 
    filter(year %in% c(fy,ly)) %>% 
    mutate(yup = index_q_0.025,
           ydown = 0)
  traj <- trajs[[1]] +
    geom_errorbar(data = indt,
              aes(x = year, ymin = ydown,ymax = yup),
              colour = grey(0.5),
              linewidth = 0.5,
              width = 0)+
    labs(title = paste(fy,"-",ly))
    
  mapt <- traj + map + plot_layout(nrow = 2,heights = c(2,4))
  mag = 3
  png(filename = paste0("Sharing/",species_f,"_",fy,"-",ly,".png"),
      res = 72 *mag,
      height = 900*mag,
      width = 600*mag)
  print(mapt)
  dev.off()
  maps[[y]] <- map
  
}

fy <- starts[1]
ly <- 2021

tr <- generate_trends(inds,
                      min_year = fy,
                      max_year = ly)
map <- plot_map(tr,
                title = FALSE)+
  coord_sf(xlim = c(base["xmin"],base["xmax"]),
           ylim = c(base["ymin"],base["ymax"]))+
  labs(title = paste(species,"Trends",fy,"-",ly),
       subtitle = "North American Breeding Bird Survey")+
  theme(title = element_text(size = 20))

indt <- ind$indices %>% 
  filter(year %in% c(fy,ly)) %>% 
  mutate(yup = index_q_0.025,
         ydown = 0)
traj <- trajs[[1]] +
  geom_errorbar(data = indt,
                aes(x = year, ymin = ydown,ymax = yup),
                colour = grey(0.5),
                linewidth = 0.5,
                width = 0)+
  labs(title = paste(species,fy,"-",ly))

mapt <- traj + map + plot_layout(nrow = 2,heights = c(2,4))
mag = 3
png(filename = paste0("Sharing/",species_f,"_",fy,"-",ly,".png"),
    res = 72 *mag,
    height = 900*mag,
    width = 600*mag)
print(mapt)
dev.off()



indt <- ind$indices %>% 
  filter(year %in% starts) %>% 
  mutate(yup = index_q_0.025,
         ydown = 0)
traj <- trajs[[1]] +
  geom_errorbar(data = indt,
                aes(x = year, ymin = ydown,ymax = yup),
                colour = grey(0.5),
                linewidth = 0.5,
                width = 0)+
  labs(title = paste(species,"BBS",fy,"-",ly))

map1 <- maps[[1]] + theme(legend.position = "none") +
  labs(title = paste(starts[[1]],starts[[2]],sep = "-"),
       subtitle = "")
map2 <- maps[[2]] + theme(legend.position = "none")+
  labs(title = paste(starts[[2]],starts[[3]],sep = "-"),
       subtitle = "")
map3 <- maps[[3]] + theme(legend.position = "none")+
  labs(title = paste(starts[[3]],starts[[4]],sep = "-"),
       subtitle = "")
map4 <- maps[[4]] + theme(legend.position = "none")+
  labs(title = paste(starts[[4]],starts[[5]],sep = "-"),
       subtitle = "")

mapt <- traj + map1 + map2 + map3 + map4 + 
plot_layout(design = c("
                                   1111
                                   2345
                                 "))
mag = 3
png(filename = paste0("Sharing/",species_f,"_",fy,"-",ly,".png"),
    res = 72 *mag,
    height = 500*mag,
    width = 600*mag)
print(mapt)
dev.off()

