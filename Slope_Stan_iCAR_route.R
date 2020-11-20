## building a Stan version of the bbsBayes models

library(bbsBayes)
library(tidyverse)
library(rstan)
library(shinystan)
library(sf)
library(spdep)
# library(ggforce)
# library(tidybayes)
source("functions/mungeCARdata4stan.R")


# load and stratify CASW data ---------------------------------------------
species = "Pacific Wren"
#species = "Barred Owl"
strat = "bbs_usgs"
model = "slope"

strat_data = stratify(by = strat)
jags_data = prepare_jags_data(strat_data = strat_data,
                             species_to_run = species,
                             model = model,
                             n_knots = 10,
                             min_year = 1999)


stan_data = jags_data[c("ncounts",
                         "nstrata",
                         "nobservers",
                         "count",
                         "strat",
                         "obser",
                         "year",
                         "firstyr",
                         "fixedyear",
                        "nonzeroweight")]
stan_data[["nyears"]] <- max(jags_data$year)
stan_data[["max_nobservers"]] <- max(jags_data$nobservers)





# spatial neighbourhood define --------------------------------------------


#load strata map
laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system



locat = system.file("maps",
                    package="bbsBayes")
map.file = "BBS_USGS_strata"

strata_map = sf::read_sf(dsn = locat,
                         layer = map.file)

strata_map = st_transform(strata_map,crs = laea)

real_strata_map = filter(strata_map,ST_12 %in% unique(jags_data$strat_name))


strata_list <- data.frame(ST_12 = unique(jags_data$strat_name),
                          strat = unique(jags_data$strat))

real_strata_map <- inner_join(real_strata_map,strata_list, by = "ST_12")


# generate neighbourhoods -------------------------------------------------

coords = st_coordinates(st_centroid(real_strata_map))


# nb_dbpoly <- spdep::poly2nb(real_strata_map,row.names = real_strata_map$strat)
# 
# plot(real_strata_map,max.plot = 1,reset = F)
# plot(nb_dbpoly,coords,add = T,col = "red")



### currently using 2 nearest neighbours to define the spatial relationships
## many regions may  have > 2 neighbours because of the symmetry condition
nb_db <- spdep::knn2nb(spdep::knearneigh(coords,k = 2),row.names = real_strata_map$strat,sym = TRUE)

plot(real_strata_map,max.plot = 1,reset = F)
plot(nb_db,coords,add = T,col = "red")

nb_info = spdep::nb2WB(nb_db)

### re-arrange GEOBUGS formated nb_info into appropriate format for Stan model
car_stan_dat <- mungeCARdata4stan(adjBUGS = nb_info$adj,
                                  numBUGS = nb_info$num)



stan_data[["N_edges"]] = car_stan_dat$N_edges
stan_data[["node1"]] = car_stan_dat$node1
stan_data[["node2"]] = car_stan_dat$node2


mod.file = "models/slope_iCAR.stan"

parms = c("sdnoise",
          "sdyear",
          "sdobs",
          "beta_p",
          "sdbeta",
          "strata_p",
          "sdstrata",
          "BETA",
          "STRATA",
           "n",
          # "nsmooth",
          "beta",
          "eta")

## compile model
slope_model = stan_model(file=mod.file)

## run sampler on model, data
stime = system.time(slope_stanfit <-
                      sampling(slope_model,
                               data=stan_data,
                               verbose=TRUE, refresh=100,
                               chains=3, iter=500,
                               warmup=400,
                               cores = 3,
                               pars = parms,
                               control = list(adapt_delta = 0.8,
                                              max_treedepth = 17)))


paste(stime[[3]]/3600,"hours")

launch_shinystan(slope_stanfit) 



