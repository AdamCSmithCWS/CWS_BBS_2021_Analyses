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
model = "gamye"

strat_data = stratify(by = strat)
jags_data = prepare_jags_data(strat_data = strat_data,
                             species_to_run = species,
                             model = model,
                             n_knots = 7,
                             min_year = 1999)


stan_data = jags_data[c("ncounts",
                         "nstrata",
                         "nobservers",
                         "count",
                         "strat",
                         #"obser",
                        #"fixedyear",
                        #"nonzeroweight",
                         "year",
                         "firstyr"
                        )]
stan_data[["nyears"]] <- max(jags_data$year)
stan_data[["sum_observers"]] <- sum(jags_data$nobservers)
stan_data[["obser"]] <- as.integer(factor(paste(jags_data$strat,jags_data$obser,sep = "_")))
stan_data[["nknots_year"]] <- jags_data$nknots
stan_data[["year_basispred"]] <- jags_data$X.basis


### needed if derived annual index values are included
# stan_data[["max_nobservers"]] <- max(jags_data$nobservers)
# 
# obsstr = unique(data.frame(obser_o = jags_data$obser,
#                            strat = jags_data$strat,
#                            obser = as.integer(factor(paste(jags_data$strat,jags_data$obser,sep = "_")))))
# obsstr = obsstr[order(obsstr$strat,obsstr$obser_o),]
# 
# obs_mat = matrix(1,nrow = jags_data$nstrata,ncol = stan_data$max_nobservers)
# 
# for(s in 1:jags_data$nstrata){
#   obs_mat[s,1:stan_data$nobservers[s]] <- obsstr[which(obsstr$strat == s),"obser"]
# }
#stan_data[["obs_mat"]] = obs_mat





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


mod.file = "models/gamye_iCAR.stan"

parms = c("sdnoise",
          "sdyear",
          "sdobs",
          "beta",
          "sdbeta",
          "strata",
          "sdstrata",
          "BETA",
          "STRATA",
          "sdBETA_gam",
          # "n",
          # "nsmooth",
          "eta",
          "log_lik")

## compile model
slope_model = stan_model(file=mod.file)

## run sampler on model, data
stime = system.time(slope_stanfit <-
                      sampling(slope_model,
                               data=stan_data,
                               verbose=TRUE, refresh=100,
                               chains=3, iter=600,
                               warmup=500,
                               cores = 3,
                               pars = parms,
                               control = list(adapt_delta = 0.8,
                                              max_treedepth = 15)))

save(list = c("slope_stanfit","species","mod.file","model","strat"),
     file = "output/Pacific Wren_iCAR_GAMYE_t4tail.RData")


get_elapsed_time(slope_stanfit)/3600 ## in hours
# warmup   sample
# chain:1 9.186139 1.975275
# chain:2 8.263056 1.897203
# chain:3 8.954222 2.016814

library(loo)
library(tidyverse)

log_lik_1 <- extract_log_lik(slope_stanfit, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik_1), cores = 10)
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 10)
print(loo_1)

# 
# Computed from 300 by 4171 log-likelihood matrix
# 
# Estimate    SE
# elpd_loo  -7936.8  86.6
# p_loo      1079.0  24.5
# looic     15873.6 173.3
# ------
#   Monte Carlo SE of elpd_loo is NA.
# 
# Pareto k diagnostic values:
#   Count Pct.    Min. n_eff
# (-Inf, 0.5]   (good)     2835  68.0%   15        
# (0.5, 0.7]   (ok)        817  19.6%   12        
# (0.7, 1]   (bad)       456  10.9%   6         
# (1, Inf)   (very bad)   63   1.5%   1         
# See help('pareto-k-diagnostic') for details.


doy = ((jags_data$month-4)*30+jags_data$day)
plot(loo_1$pointwise[,"influence_pareto_k"],log(stan_data$count+1))
plot(loo_1$pointwise[,"influence_pareto_k"],doy)
plot(doy,log(stan_data$count+1))



loo2 = data.frame(loo_1$pointwise)

loo2$flag = cut(loo2$influence_pareto_k,breaks = c(0,0.5,0.7,1,Inf))
dts = data.frame(count = stan_data$count,
                 obser = stan_data$obser,
                 strat = stan_data$strat,
                 year = stan_data$year)
loo2 = cbind(loo2,dts)

plot(log(loo2$count+1),loo2$influence_pareto_k)

obserk = loo2 %>% group_by(obser) %>% 
  summarise(n = n(),
            mean_k = mean(influence_pareto_k),
            max_k = max(influence_pareto_k),
            sd_k = sd(influence_pareto_k),
            mean_looic = mean(looic),
            mean_ploo = mean(p_loo))
plot(obserk$n,obserk$max_k)
plot(obserk$n,obserk$mean_k)
plot(obserk$n,obserk$sd_k)
plot(obserk$n,obserk$mean_looic)
plot(obserk$n,obserk$mean_ploo)


yeark = loo2 %>% group_by(year) %>% 
  summarise(n = n(),
            mean_k = mean(influence_pareto_k),
            q90 = quantile(influence_pareto_k,0.9),
            max_k = max(influence_pareto_k),
            sd_k = sd(influence_pareto_k),
            strat = mean(strat),
            sd = sd(strat))
plot(yeark$year,yeark$max_k)
plot(yeark$year,yeark$mean_k)
plot(yeark$year,yeark$sd_k)
plot(yeark$year,yeark$q90)

stratk = loo2 %>% group_by(strat) %>% 
  summarise(n = n(),
            mean_k = mean(influence_pareto_k),
            q90_k = quantile(influence_pareto_k,0.9),
            max_k = max(influence_pareto_k),
            sd_k = sd(influence_pareto_k),
            strat = mean(strat),
            sd = sd(strat))
plot(stratk$strat,stratk$max_k)
plot(stratk$n,stratk$mean_k)

plot(stratk$strat,stratk$mean_k)
plot(stratk$strat,stratk$sd_k)
plot(stratk$strat,stratk$q90_k)



launch_shinystan(slope_stanfit) 



