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
source("functions/prepare-jags-data-alt.R") ## small alteration of the bbsBayes function
## changes captured in a commit on Nov 20, 2020


# load and stratify CASW data ---------------------------------------------
species = "Pacific Wren"
#species = "Barred Owl"
strat = "bbs_cws"
model = "slope"

strat_data = stratify(by = strat)



jags_data = prepare_jags_data(strat_data = strat_data,
                             species_to_run = species,
                             model = model,
                             #n_knots = 10,
                             min_year = 1999)


stan_data = jags_data[c("ncounts",
                         #"nstrata",
                         #"nobservers",
                         "count",
                         #"strat",
                         #"obser",
                         "year",
                         "firstyr",
                         "fixedyear")]
stan_data[["nyears"]] <- max(jags_data$year)
stan_data[["observer"]] <- as.integer(factor((jags_data$ObsN)))
stan_data[["nobservers"]] <- max(stan_data$observer)


# spatial neighbourhood define --------------------------------------------
laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system

locat = system.file("maps",
                    package = "bbsBayes")
map.file = "BBS_CWS_strata"

strata_map = read_sf(dsn = locat,
                     layer = map.file)
strata_map = st_transform(strata_map,crs = laea)

real_strata_map = filter(strata_map,ST_12 %in% unique(jags_data$strat_name))

strata_list <- data.frame(ST_12 = unique(jags_data$strat_name),
                          strat = unique(jags_data$strat))


real_strata_map <- inner_join(real_strata_map,strata_list, by = "ST_12")
strata_bounds <- st_union(real_strata_map)
strata_bounds_buf = st_buffer(strata_bounds,dist = 100000)


jags_data[["routeF"]] <- as.integer(factor((jags_data$route)))

route_map = unique(data.frame(route = jags_data$route,
                              routeF = jags_data$routeF,
                              strat = jags_data$strat_name,
                              Latitude = jags_data$Latitude,
                              Longitude = jags_data$Longitude))


# reconcile duplicate spatial locations -----------------------------------
dups = which(duplicated(route_map[,c("Latitude","Longitude")]))
if(length(dups) > 0){
  route_map[dups,"Latitude"] <- route_map[dups,"Latitude"]+0.01 #=0.01 decimal degrees ~ 1km
  route_map[dups,"Longitude"] <- route_map[dups,"Longitude"]+0.01 #=0.01 decimal degrees ~ 1km
  
}
route_map = st_as_sf(route_map,coords = c("Longitude","Latitude"))
st_crs(route_map) <- 4269 #NAD83 commonly used by US federal agencies
#load strata map


route_map = st_transform(route_map,crs = laea)

# generate neighbourhoods -------------------------------------------------

# coords = st_coordinates(route_map)


# nb_dbpoly <- spdep::poly2nb(real_strata_map,row.names = real_strata_map$strat)
# 
# plot(real_strata_map,max.plot = 1,reset = F)
# plot(nb_dbpoly,coords,add = T,col = "red")


# Voronoi polygons from route locations -----------------------------------
box <- st_as_sfc(st_bbox(route_map))

v <- st_cast(st_voronoi(st_union(route_map), envelope = box))

vint = st_sf(st_cast(st_intersection(v,strata_bounds_buf),"POLYGON"))
vintj = st_join(vint,route_map,join = st_contains)
vintj = arrange(vintj,routeF)

nb_db = poly2nb(vintj,row.names = vintj$route,queen = FALSE)


### currently using 2 nearest neighbours to define the spatial relationships
## many regions may  have > 2 neighbours because of the symmetry condition
# nb_db <- spdep::knn2nb(spdep::knearneigh(coords,k = 4),row.names = route_map$route,sym = TRUE)
cc = st_coordinates(st_centroid(vintj))
plot(nb_db,cc,col = "red")

ggp = ggplot(data = route_map)+
  geom_sf(data = vintj,alpha = 0.3)+
  geom_sf(aes(col = strat))+
  geom_sf_text(aes(label = routeF),size = 3,alpha = 0.3)
print(ggp)

# wca = which(grepl(route_map$strat,pattern = "-CA-",fixed = T))
# wak = which(grepl(route_map$strat,pattern = "-AK-",fixed = T))
# 
# nb2[[wak]]


nb_info = spdep::nb2WB(nb_db)

### re-arrange GEOBUGS formated nb_info into appropriate format for Stan model
car_stan_dat <- mungeCARdata4stan(adjBUGS = nb_info$adj,
                                  numBUGS = nb_info$num)



stan_data[["N_edges"]] = car_stan_dat$N_edges
stan_data[["node1"]] = car_stan_dat$node1
stan_data[["node2"]] = car_stan_dat$node2
stan_data[["route"]] = jags_data$routeF
stan_data[["nroutes"]] = max(jags_data$routeF)

mod.file = "models/slope_iCAR_route.stan"

parms = c("sdnoise",
          "sdyear",
          "sdobs",
          "sdbeta",
          "alpha",
          "sdstrata",
          "BETA",
          "ALPHA",
          "beta",
          "eta",
          "log_lik")

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



