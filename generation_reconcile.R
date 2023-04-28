
library(bbsBayes)
library(tidyverse)

strat <- stratify("bbs_usgs")
sp_full <- strat$species_strat %>% 
  mutate(sci_name = paste(genus,species)) # the full BBS species list from the 2021 data release



gens = read.csv("alt_data/cobi13486-sup-0004-tables4.csv") # supplement from Bird et al.

gens <- gens %>% select(Scientific_name,
                        GenLength)

spslist <- left_join(sp_full,gens,by = c("sci_name" = "Scientific_name"))



# reconciling scientific names --------------------------------------------

sps_nomatch <- spslist[which(is.na(spslist$GenLength)),"english"]

fullgensnames = read.csv("alt_data/cobi13486-sup-0001-tables1.csv") # an additional supplement from Bird et al. with alternate common names
sps_altmatch <- sps_nomatch[sps_nomatch %in% unique(fullgensnames$Common_name)]
for(s1 in sps_altmatch){
  scin1 = fullgensnames[which(fullgensnames$Common_name == s1),"Scientific_name"]
  spslist[which(spslist$english == s1),"GenLength"] <- gens[which(gens$Scientific_name == scin1),"GenLength"]
}



# Manual adjustments for unmatched species --------------------------------



spslist[which(spslist$english == "American Three-toed Woodpecker"),"GenLength"] <- gens[which(gens$Scientific_name == "Picoides tridactylus"),"GenLength"]
spslist[which(spslist$english == "Black-necked Stilt"),"GenLength"] <- gens[which(gens$Scientific_name == "Himantopus himantopus"),"GenLength"]
spslist[which(spslist$english == "Green Heron"),"GenLength"] <- gens[which(gens$Scientific_name == "Butorides striata"),"GenLength"]
spslist[which(spslist$english == "Spotted Dove"),"GenLength"] <- gens[which(gens$Scientific_name == "Spilopelia chinensis"),"GenLength"]
spslist[which(spslist$english == "Island Scrub-Jay"),"GenLength"] <- gens[which(gens$Scientific_name == "Aphelocoma californica"),"GenLength"]
spslist[which(spslist$english == "Woodhouse's Scrub-Jay"),"GenLength"] <- gens[which(gens$Scientific_name == "Aphelocoma californica"),"GenLength"]
spslist[which(spslist$english == "Black Oystercatcher"),"GenLength"] <- gens[which(gens$Scientific_name == "Haematopus ater"),"GenLength"]
spslist[which(spslist$english == "Hoary Redpoll"),"GenLength"] <- gens[which(gens$Scientific_name == "Acanthis flammea"),"GenLength"]
spslist[which(spslist$english == "Redpoll (Common/Hoary)"),"GenLength"] <- gens[which(gens$Scientific_name == "Acanthis flammea"),"GenLength"]
spslist[which(spslist$english == "Mexican Duck"),"GenLength"] <- gens[which(gens$Scientific_name == "Anas platyrhynchos"),"GenLength"]

spslist[which(spslist$english == "Blue Grouse (Dusky/Sooty)"),"GenLength"] <- mean(gens[grepl(gens$Scientific_name, pattern = "Dendragapus"),"GenLength"])
spslist[which(spslist$english == "Western Grebe (Clark's/Western)"),"GenLength"] <- mean(gens[grepl(gens$Scientific_name, pattern = "Aechmophorus"),"GenLength"])
spslist[which(spslist$english == "Neotropic Cormorant"),"GenLength"] <- gens[which(gens$Scientific_name == "Nannopterum brasilianus"),"GenLength"]
spslist[which(spslist$english == "Sapsuckers (Yellow-bellied/Red-naped/Red-breasted/Williamson's)"),"GenLength"] <- mean(gens[grepl(gens$Scientific_name, pattern = "Sphyrapicus"),"GenLength"])
spslist[which(spslist$english == "Northern Oriole (Bullock's/Baltimore)"),"GenLength"] <- mean(gens[which(gens$Scientific_name == "Icterus bullockiorum" |
                                                                                                           gens$Scientific_name == "Icterus galbula"),"GenLength"])
spslist[which(spslist$english == "(Blue Goose) Snow Goose"),"GenLength"] <- gens[which(gens$Scientific_name == "Anser caerulescens"),"GenLength"]
spslist[which(spslist$english == "(Black Brant) Brant"),"GenLength"] <- gens[which(gens$Scientific_name == "Branta bernicla"),"GenLength"]
spslist[which(spslist$english == "(Great White Heron) Great Blue Heron"),"GenLength"] <- gens[which(gens$Scientific_name == "Ardea herodias"),"GenLength"]
spslist[which(spslist$english == "(Harlan's Hawk) Red-tailed Hawk"),"GenLength"] <- gens[which(gens$Scientific_name == "Buteo jamaicensis"),"GenLength"]
spslist[which(spslist$english == "(Yellow-shafted Flicker) Northern Flicker"),"GenLength"] <- gens[which(gens$Scientific_name == "Colaptes auratus"),"GenLength"]
spslist[which(spslist$english == "(Red-shafted Flicker) Northern Flicker"),"GenLength"] <- gens[which(gens$Scientific_name == "Colaptes auratus"),"GenLength"]
spslist[which(spslist$english == "Traill's Flycatcher (Alder/Willow)"),"GenLength"] <-  mean(gens[which(gens$Scientific_name == "Empidonax alnorum" |
                                                                                                          gens$Scientific_name == "Empidonax traillii"),"GenLength"])
spslist[which(spslist$english == "Western Flycatcher (Cordilleran/Pacific-slope)"),"GenLength"] <-  mean(gens[which(gens$Scientific_name == "Empidonax difficilis" |
                                                                                                          gens$Scientific_name == "Empidonax occidentalis"),"GenLength"])


spslist[which(spslist$english == "Solitary Vireo (Blue-headed/Cassin's)"),"GenLength"] <- mean(gens[which(gens$Scientific_name == "Vireo cassinii" |
                                                                                                            gens$Scientific_name == "Vireo solitarius"),"GenLength"])

spslist[grepl(spslist$english,pattern = "Dark-eyed Junco"),"GenLength"] <- gens[which(gens$Scientific_name == "Junco hyemalis"),"GenLength"]
spslist[grepl(spslist$english,pattern = "Yellow-rumped Warbler"),"GenLength"] <- gens[which(gens$Scientific_name == "Setophaga coronata"),"GenLength"]
spslist[grepl(spslist$english,pattern = "American Crow") &
          !grepl(spslist$english,pattern = "unid"),"GenLength"] <- gens[which(gens$Scientific_name == "Corvus brachyrhynchos"),"GenLength"]



write_csv(spslist,
          "alt_data/full_bbs_species_list_w3_generations.csv")


