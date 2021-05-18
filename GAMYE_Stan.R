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
source("functions/neighbours_define.R")

# spatial data load -------------------------------------------------------

laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system

locat = system.file("maps",
                    package = "bbsBayes")
map.file = "BBS_CWS_strata" # the most standard BBS stratification, used here just to provide some spatial limits

# reading in the strata spatial file (GIS shapefile)
strata_map = read_sf(dsn = locat,
                     layer = map.file)
strata_map = st_transform(strata_map,crs = laea) #reprojecting the geographic coordinate file to an equal area projection


strat = "bbs_cws"
model = "gamye"

strat_data = stratify(by = strat)

# load and stratify species data ---------------------------------------------
species = "Bobolink"

start_year = 1970

jags_data = prepare_jags_data(strat_data = strat_data,
                             species_to_run = species,
                             model = model,
                             min_year = start_year,
                             n_knots = 10)


str_link <- unique(data.frame(strat = jags_data$strat,
                              strat_name = jags_data$strat_name))
str_link <- str_link %>% arrange(strat)

realized_strata_map <- inner_join(strata_map,str_link,by = c("ST_12" = "strat_name")) %>% 
  arrange(strat)
# generate neighbourhoods -------------------------------------------------
buff_dist <- 10000




### re-arrange GEOBUGS formated nb_info into appropriate format for Stan model
car_stan_dat <- neighbours_define(real_strata_map = realized_strata_map,
                                   strat_link_fill = buff_dist,
                                   species = species,
                                  voronoi = FALSE)




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

nroutes_strata <- as.integer(table(strat_route$strat))
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




ncounts = stan_data$ncounts
nroutes = stan_data$nroutes
nstrata = stan_data$nstrata
nyears = stan_data$nyears
nknots_year = stan_data$nknots_year
nobservers = stan_data$nobservers



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
          #"nu",
          "BETA",
          "STRATA",
           "n",
           "nsmooth",
          "eta",
          "log_lik")


# # cmdStanR ----------------------------------------------------------------
mod.file = "models/gamye_iCAR.stan"

library(cmdstanr)
## compile model
slope_model <- cmdstan_model(mod.file)

init_def <- function(){ list(noise_raw = rnorm(ncounts,0,0.1),
                            strata_raw = rnorm(nstrata,0,0.1),
                            STRATA = 0,
                            eta = 0,
                            yeareffect_raw = matrix(rnorm(nstrata*nyears,0,0.1),nrow = nstrata,ncol = nyears),
                            obs_raw = rnorm(nobservers,0,0.1),
                            rte_raw = rnorm(nroutes,0,0.1),
                            sdnoise = 0.2,
                            sdobs = 0.1,
                            sdrte = 0.2,
                            sdbeta = runif(nknots_year,0.01,0.1),
                            sdBETA = 0.1,
                            sdyear = runif(nstrata,0.01,0.1),
                            nu = 10,
                            BETA_raw = rnorm(nknots_year,0,0.1),
                            beta_raw = matrix(rnorm(nknots_year*nstrata,0,0.01),nrow = nstrata,ncol = nknots_year))}



## run sampler on model, data
# data_file <- "tmpdata/tmp_data.json"
# write_stan_json(stan_data,file = data_file)
slope_stanfit <- slope_model$sample(
                          data=stan_data,
                          refresh=25,
                          chains=3, iter_sampling=400,
                          iter_warmup=500,
                          parallel_chains = 3,
                          #pars = parms,
                          adapt_delta = 0.8,
                          max_treedepth = 14,
                          seed = 123,
                          init = init_def)



species_file = gsub(pattern = "([[:punct:]]|[[:blank:]])","",species)

save(list = c("slope_stanfit","stan_data","jags_data","model"),
     file = paste0("output/cmdStan_",species_file,"ten_yr__gamye_iCAR.RData"))

slope_stanfit$cmdstan_diagnose()


# export to csv and read in as rstan --------------------------------------

slope_stanfit$save_output_files(dir = "output",
                                basename = paste0(species_file,"_cmdStan_out_",nyears))
csv_files <- dir("output/",pattern = paste0(species,"_cmdStan_out_",nyears),full.names = TRUE)

#sl_rstan <- As.mcmc.list(read_stan_csv(csv_files))


# mod_sum = slope_stanfit$summary(variables = "n")



loo_sum = slope_stanfit$loo()
ploo = as.data.frame(loo_sum$pointwise)
ploo$k_cat <- cut(ploo$influence_pareto_k,breaks = c(-Inf,0.5,0.7,1,Inf),levels = c("good","ok","bad","very_bad"))
df = data.frame(count = stan_data$count,
                route = stan_data$route,
                observer = stan_data$observer,
                strat = stan_data$strat,
                year = stan_data$year)
ploo = bind_cols(ploo,df)
nyear_sqrt <- ceiling(sqrt(length(unique(ploo$year))))
k_count = ggplot(data = ploo,aes(x = count,y = influence_pareto_k,colour = k_cat))+
  geom_point()+
  facet_wrap(~year,nrow = nyear_sqrt,scales = "fixed")+
geom_abline(slope = 0,intercept = 0.7,colour = "red")
print(k_count)

looic_count = ggplot(data = ploo,aes(x = count,y = looic,colour = k_cat))+
  geom_point()+
  facet_wrap(~year,nrow = nyear_sqrt,scales = "fixed")
print(looic_count)


k_looic = ggplot(data = ploo,aes(x = looic,y = influence_pareto_k,colour = k_cat))+
  geom_point()+
  facet_wrap(~year,nrow = nyear_sqrt,scales = "fixed")+
  geom_abline(slope = 0,intercept = 0.7,colour = "red")
print(k_looic)

looic_strat = ggplot(data = ploo,aes(x = strat,y = looic,colour = k_cat))+
  geom_point()+
  facet_wrap(~year,nrow = nyear_sqrt,scales = "fixed")
print(looic_strat)

strat_loo_sum <- ploo %>% group_by(strat) %>% 
  summarise(m_k = mean(influence_pareto_k),
            m_loo = mean(looic),
            m_ploo = mean(p_loo),
            md_k = median(influence_pareto_k),
            md_loo = median(looic),
            md_ploo = median(p_loo))

map_loo <- inner_join(realized_strata_map,strat_loo_sum,by = "strat")
plot_map_loo = ggplot(data = map_loo,aes(fill = md_loo))+
  geom_sf()
print(plot_map_loo)

plot_map_k = ggplot(data = map_loo,aes(fill = md_k))+
  geom_sf()
print(plot_map_k)

plot_map_ploo = ggplot(data = map_loo,aes(fill = md_ploo))+
  geom_sf()
print(plot_map_ploo)



source("functions/cmdstanr_gather.R")
# index and trend plotting ------------------------------------------------

nsmooth_samples <- gather_samples(fit = slope_stanfit,
                                  parm = "nsmooth",
                                  dims = c("s","y")) %>% 
  left_join(.,str_link,by = c("s" = "strat")) %>% 
  mutate(year = y+(start_year-1))


n_samples <- gather_samples(fit = slope_stanfit,
                                  parm = "n",
                                  dims = c("s","y")) %>% 
  left_join(.,str_link,by = c("s" = "strat")) %>% 
  mutate(year = y+(start_year-1))



# STratum level indices ---------------------------------------------------


ind_sm <- nsmooth_samples %>% group_by(s,strat_name,year) %>% 
  summarise(mean = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(type = "smooth")

ind_full <- n_samples %>% group_by(s,strat_name,year) %>% 
  summarise(mean = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975))%>% 
  mutate(type = "full")


inds <- bind_rows(ind_sm,ind_full)

nstrata <- max(inds$s)
nf <- ceiling(sqrt(nstrata))

ind_fac <- ggplot(data = inds,aes(x = year,y = mean))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = type),alpha = 0.2)+
  geom_line(aes(colour = type))+
  scale_y_continuous(limits = c(0,NA))+
  facet_wrap(~strat_name,nrow = nf,scales = "free")
  
print(ind_fac)


# continental indices -----------------------------------------------------

a_weights <- as.data.frame(realized_strata_map) %>% 
  rename(s = strat) %>% 
  mutate(area = AREA/sum(AREA)) %>% 
  select(s,area)


I_sm <- nsmooth_samples %>% left_join(.,a_weights,by = "s") %>% 
  mutate(.value = .value*area) %>% 
  group_by(.draw,year) %>% 
  summarise(.value = sum(.value)) %>% 
  group_by(year) %>% 
  summarise(mean = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(type = "smooth")

I_full <- n_samples %>% left_join(.,a_weights,by = "s") %>% 
  mutate(.value = .value*area) %>% 
  group_by(.draw,year) %>% 
  summarise(.value = sum(.value)) %>% 
  group_by(year) %>% 
  summarise(mean = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>%
  mutate(type = "full")


Is <- bind_rows(I_sm,I_full)


I_plot <- ggplot(data = Is,aes(x = year,y = mean))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = type),alpha = 0.2)+
  geom_line(aes(colour = type))+
  scale_y_continuous(limits = c(0,NA))

print(I_plot)





# Trends ------------------------------------------------------------------




tr_func <- function(samples = nsmooth_samples,
                    scale = "strat_name",
                    start_year = 2000,
                    end_year = 2019){
  
 
  nyear = end_year - start_year
  
  #samples$reg <- as.character(samples[,scale])
  p_py <- function(i1,i2,ny = nyear){
    ppy <- 100*((i2/i1)^(1/(ny))-1)
    return(ppy)
  }
  
  p_ch <- function(i1,i2){
    pch <- 100*((i2/i1)-1)
    return(pch)
  }
  
  
  tr_sm <- samples %>% filter(year %in% c(start_year,end_year)) %>% 
    pivot_wider(.,id_cols = any_of(c(scale,"strat_name",
                                     ".variable",".chain",".iteration",".draw")),
                names_from = year,
                values_from = .value,
                names_prefix = "Y") %>% 
    rename_with(.,~gsub(pattern = paste0("Y",start_year),"Y1",.x,fixed = TRUE))%>% 
    rename_with(.,~gsub(pattern = paste0("Y",end_year),"Y2",.x,fixed = TRUE)) %>% 
    rename_with(.,~gsub(pattern = scale,"reg",.x,fixed = TRUE)) %>% 
    mutate(ch = p_ch(Y1,Y2),
           t = p_py(Y1,Y2)) %>% 
    group_by(reg) %>% 
    summarise(mean_trend = mean(t),
              lci_trend = quantile(t,0.025),
              uci_trend = quantile(t,0.975),
              mean_ch = mean(ch),
              lci_ch = quantile(ch,0.025),
              uci_ch = quantile(ch,0.975)) %>% 
    mutate(start_year = start_year,
           end_year = end_year) %>% 
    rename_with(.,~gsub(pattern = "reg",scale,.x,fixed = TRUE))
  
  return(tr_sm)
  
  
}







trend_90 = tr_func(nsmooth_samples,start_year = 1990)
trend_09 = tr_func(nsmooth_samples,start_year = 2009)


trend_map <- function(base_map = realized_strata_map,
                      trends = trend_90){
  breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
  labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
  labls = paste0(labls, " %")
  map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                   "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
  names(map_palette) <- labls
  
  tmap <- left_join(base_map,trends,by = c("ST_12" = "strat_name")) %>% 
    mutate(Tplot = cut(mean_trend,breaks = c(-Inf, breaks, Inf),labels = labls))
  
  fyear = unique(trends$start_year)
  lyear = unique(trends$end_year)
  tplot <- ggplot(data = tmap)+
    geom_sf(aes(fill = Tplot)) +
    scale_colour_manual(values = map_palette, 
                        aesthetics = c("fill"),
                        guide = ggplot2::guide_legend(reverse=TRUE),
                        name = paste0("Trend\n",fyear,"-",lyear))
  
  
  return(tplot)
  
}

m1 = trend_map(trends = trend_09)
print(m1)


m2 = trend_map(trends = trend_90)
print(m2)

# RStan -------------------------------------------------------------------

mod.file = "models/gamye_iCAR.stan"

## compile model
slope_model = stan_model(file=mod.file)

## run sampler on model, data
slope_stanfit <- sampling(slope_model,
                               data=stan_data,
                               verbose=TRUE, refresh=25,
                               chains=3, 
                          iter=800,
                               warmup=500,
                          cores = 3,
                               pars = parms,
                               control = list(adapt_delta = 0.8,
                                              max_treedepth = 14),
                          init = init_def)




save(list = c("slope_stanfit","stan_data","jags_data","model"),
     file = paste0("output/",species,"full__gamye_iCAR.RData"))

loo_sum_rstan = loo(slope_stanfit)


loo_sum_rstan

launch_shinystan(slope_stanfit) 




n_draws = gather_draws(slope_stanfit,n[s,y])
nsmooth_draws = gather_draws(slope_stanfit,nsmooth[s,y])














