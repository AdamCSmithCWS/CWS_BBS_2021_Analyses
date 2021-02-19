## building a Stan version of the bbsBayes models

library(bbsBayes)
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE,javascript = FALSE)
library(shinystan)
library(sf)
library(spdep)


source("functions/mungeCARdata4stan.R") ## function to modify the BUGS formatted spatial neighbourhood data to the required format for the Stan iCAR model
source("functions/prepare-jags-data-alt.R") ## function to modify the BUGS formatted spatial neighbourhood data to the required format for the Stan iCAR model


# spatial data load -------------------------------------------------------

laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system

locat = system.file("maps",
                    package = "bbsBayes")
map.file = "BBS_CWS_strata" # the most standard BBS stratification, used here just to provide some spatial limits

# reading in the strata spatial file (GIS shapefile)
strata_map = read_sf(dsn = locat,
                     layer = map.file)
strata_map = st_transform(strata_map,crs = laea) #reprojecting the geographic coordinate file to an equal area projection




# load and stratify CASW data ---------------------------------------------
species = "Wood Thrush"
strat = "bbs_cws"
model = "gamye"

strat_data = stratify(by = strat)
jags_data = prepare_jags_data(strat_data = strat_data,
                             species_to_run = species,
                             model = model,
                             #min_year = 1999,
                             n_knots = 13)


str_link <- unique(data.frame(strat = jags_data$strat,
                              strat_name = jags_data$strat_name))
str_link <- str_link %>% arrange(strat)

real_strata_map <- inner_join(strata_map,str_link,by = c("ST_12" = "strat_name")) %>% 
  arrange(strat)
# generate neighbourhoods -------------------------------------------------
strat_link_fill = 100000 #distance between strata to fill if isolated strata remain

reg_bounds <- st_union(real_strata_map)
reg_bounds_buf = st_buffer(reg_bounds,dist = strat_link_fill)

centres = suppressWarnings(st_centroid(real_strata_map))
coords = st_coordinates(centres)


box <- st_as_sfc(st_bbox(centres))

v <- st_cast(st_voronoi(st_union(centres), envelope = box))

vint = st_sf(st_cast(st_intersection(v,reg_bounds_buf),"POLYGON"))
#vint <- bind_cols(centres,v)
vintj = st_join(vint,centres,join = st_contains)
vintj = arrange(vintj,strat)

nb_db = poly2nb(vintj,row.names = vintj$strat,queen = FALSE)


### currently using 2 nearest neighbours to define the spatial relationships
## many regions may  have > 2 neighbours because of the symmetry condition
# nb_db <- spdep::knn2nb(spdep::knearneigh(coords,k = 4),row.names = route_map$route,sym = TRUE)
cc = suppressWarnings(st_coordinates(st_centroid(vintj)))
#

ggp = ggplot(data = real_strata_map)+
  geom_sf(data = vintj,alpha = 0.3)+ 
  geom_sf(aes(col = strat))+
  geom_sf_text(aes(label = strat),size = 3,alpha = 0.3)+
  theme(legend.position = "none")
pdf(file = paste0("Figures/",species,"strata_connections ",strat_link_fill/1000,".pdf"))
plot(nb_db,cc,col = "pink")
text(labels = rownames(cc),cc ,pos = 2)
print(ggp)
dev.off()


# wca = which(grepl(route_map$strat,pattern = "-CA-",fixed = T))
# wak = which(grepl(route_map$strat,pattern = "-AK-",fixed = T))
# 
# nb2[[wak]]


nb_info = spdep::nb2WB(nb_db)










### re-arrange GEOBUGS formated nb_info into appropriate format for Stan model
car_stan_dat <- mungeCARdata4stan(adjBUGS = nb_info$adj,
                                  numBUGS = nb_info$num)




stan_data = jags_data[c("ncounts",
                         "nstrata",
                         "count",
                         "strat",
                         "year",
                         "firstyr",
                        "nonzeroweight")]
stan_data[["observer"]] <- as.integer(factor(jags_data$ObsN))
stan_data[["nobservers"]] <- max(stan_data[["observer"]])

stan_data[["route"]] <- as.integer(factor(jags_data$route))
stan_data[["nroutes"]] <- max(stan_data[["route"]])
strat_route <- unique(data.frame(strat = stan_data$strat,
                                 route = stan_data$route))
strat_route <- strat_route %>% arrange(strat)

nroutes_strata <- table(strat_route$strat)
maxnroutes_strata <- max(nroutes_strata)
stan_data[["nroutes_strata"]] <- nroutes_strata
stan_data[["maxnroutes_strata"]] <- maxnroutes_strata

rte_mat <- matrix(data = 0,
                  nrow = stan_data$nstrata,
                  ncol = maxnroutes_strata)
for(i in 1:stan_data$nstrata){
  rte_mat[i,1:nroutes_strata[i]] <- strat_route[which(strat_route$strat == i),"route"]
}

stan_data[["rte_mat"]] <- rte_mat

stan_data[["nyears"]] <- max(jags_data$year)
stan_data[["nknots_year"]] <- jags_data$nknots
stan_data[["year_basispred"]] <- jags_data$X.basis

stan_data[["N_edges"]] = car_stan_dat$N_edges
stan_data[["node1"]] = car_stan_dat$node1
stan_data[["node2"]] = car_stan_dat$node2



mod.file = "models/gamye_iCAR.stan"

parms = c("sdnoise",
          "sdyear",
          "sdBETA",
          "sdobs",
          "sdrte",
          "beta",
          "sdbeta",
          "strata",
          "sdstrata",
          "eta",
          "nu",
          "BETA",
          "STRATA",
           "n",
           "nsmooth",
          "eta")

## compile model
slope_model = stan_model(file=mod.file)

## run sampler on model, data
slope_stanfit <- sampling(slope_model,
                               data=stan_data,
                               verbose=TRUE, refresh=100,
                               chains=3, iter=500,
                               warmup=400,
                               cores = 3,
                               pars = parms,
                               control = list(adapt_delta = 0.8,
                                              max_treedepth = 14))




save(list = c("slope_stanfit","stan_data","jags_data","model"),
     file = paste0("output/",species,"_gamye_iCAR.RData"))

launch_shinystan(slope_stanfit) 
