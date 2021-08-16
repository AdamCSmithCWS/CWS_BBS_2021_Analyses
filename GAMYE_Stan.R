## building a Stan version of the bbsBayes models

library(bbsBayes)
library(tidyverse)

# library(rstan)
# rstan_options(auto_write = TRUE,javascript = FALSE)
#library(shinystan)
library(sf)
library(spdep)
library(cmdstanr)



source("functions/mungeCARdata4stan.R") ## function to modify the BUGS formatted spatial neighbourhood data to the required format for the Stan iCAR model
source("functions/prepare-jags-data-alt.R") ## 
source("functions/neighbours_define.R")

# spatial data load -------------------------------------------------------

laea = st_crs("+proj=laea +lat_0=40 +lon_0=-95") # Lambert equal area coord reference system

locat = system.file("maps",
                    package = "bbsBayes")
map.file = "BBS_USGS_strata" # the most standard BBS stratification, used here just to provide some spatial limits

# reading in the strata spatial file (GIS shapefile)
strata_map = read_sf(dsn = locat,
                     layer = map.file)
strata_map = st_transform(strata_map,crs = laea) #reprojecting the geographic coordinate file to an equal area projection

if(any(grepl(names(strata_map),pattern = "AREA_1"))){
  strata_map <- strata_map %>% 
    rename(AREA = AREA_1)
}




strat = "bbs_usgs"
model = "gamye"

strat_data = stratify(by = strat)

# load and stratify species data ---------------------------------------------

allspecies.eng = strat_data$species_strat$english

cdn_trends <- read.csv("c:/BBS_summaries/2019 BBS trends for website.csv")

cdn_t_sp = unique(cdn_trends$species)

cdn_t_sp <- cdn_t_sp[-grep(pattern = "unid.",cdn_t_sp)]

output_dir <- "G:/bbsStanBayes/output"



for(species in rev(cdn_t_sp)){
  
  species_file = gsub(pattern = "([[:punct:]]|[[:blank:]])","",species)
  
  spchk <- dir(output_dir,pattern = species_file,all.files = FALSE)
  
  if(length(spchk) > 0){next}
  
  
shiny_explore = FALSE # set to TRUE to automatically launch shinystan on model finish

#first_year = 1970

jags_data = prepare_jags_data(strat_data = strat_data,
                             species_to_run = species,
                             model = model,
                             min_max_route_years = 2,
                             #min_year = first_year,
                             n_knots = 14)

first_year = min(jags_data$r_year)

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
                          refresh=200,
                          chains=3, iter_sampling=1000,
                          iter_warmup=1000,
                          parallel_chains = 3,
                          #pars = parms,
                          adapt_delta = 0.8,
                          max_treedepth = 14,
                          seed = 123,
                          init = init_def)




out_base <- paste0(species_file,"_",nyears)

# export to csv and read in as rstan --------------------------------------
slope_stanfit$save_output_files(dir = output_dir,
                                basename = out_base,
                                random = FALSE,
                                timestamp = FALSE)

csv_files <- dir(output_dir,pattern = out_base,full.names = TRUE)



save(list = c("slope_stanfit","stan_data","jags_data","csv_files","realized_strata_map"),
     file = paste0(output_dir,"/",out_base,"_gamye_iCAR.RData"))




# export to csv and read in as rstan --------------------------------------






}# end of species loop


#sl_rstan <- As.mcmc.list(read_stan_csv(csv_files))
if(shiny_explore){
  
  
  csv_files <- dir(output_dir,pattern = out_base,full.names = TRUE)
  
  
  
sl_rstan <- rstan::read_stan_csv(csv_files)
launch_shinystan(as.shinystan(sl_rstan))

loo_stan = loo(sl_rstan)
}
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

pdf(paste0("figures/",species_file,"loo_plots.pdf"),
    width = 11,
    height = 8.5)


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
            m_looic = mean(looic),
            m_ploo = mean(p_loo),
            md_k = median(influence_pareto_k),
            md_looic = median(looic),
            md_ploo = median(p_loo))

map_loo <- inner_join(realized_strata_map,strat_loo_sum,by = "strat")
plot_map_loo = ggplot(data = map_loo,aes(fill = md_looic))+
  geom_sf()
print(plot_map_loo)

plot_map_k = ggplot(data = map_loo,aes(fill = md_k))+
  geom_sf()
print(plot_map_k)

plot_map_ploo = ggplot(data = map_loo,aes(fill = md_ploo))+
  geom_sf()
print(plot_map_ploo)

dev.off()

source("functions/posterior_summary_functions.R")
source("functions/trend_function.R")

# index and trend plotting ------------------------------------------------

nsmooth_samples <- posterior_samples(fit = sl_rstan,
                                  parm = "nsmooth",
                                  dims = c("s","y")) %>% 
  left_join(.,str_link,by = c("s" = "strat")) %>% 
  mutate(year = y+(first_year-1))


n_samples <- posterior_samples(fit = sl_rstan,
                                  parm = "n",
                                  dims = c("s","y")) %>% 
  left_join(.,str_link,by = c("s" = "strat")) %>% 
  mutate(year = y+(first_year-1))



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

raw <- data.frame(count = stan_data$count,
                  strat = stan_data$strat,
                  year = stan_data$year+(first_year-1)) %>% 
  left_join(.,str_link,by = "strat")

raw_means <- raw %>% group_by(year,strat,strat_name) %>% 
  summarise(mean_count = mean(count),
            uqrt_count = quantile(count,0.75),
            n_surveys = n())



ind_fac <- ggplot(data = inds,aes(x = year,y = mean))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = type),alpha = 0.2)+
  geom_line(aes(colour = type))+
  geom_point(data = raw_means,aes(x = year,y = mean_count,fill = n_surveys),alpha = 0.2,inherit.aes = FALSE)+
  scale_y_continuous(limits = c(0,NA))+
  facet_wrap(~strat_name,nrow = nf,scales = "free")
  
#print(ind_fac)


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
  labs(title = paste(species,"survey wide trajectory"))+
  scale_y_continuous(limits = c(0,NA))

pdf(file = paste0("figures/",species_file,"trajectories.pdf"))

print(I_plot)


for(st in str_link$strat_name){
  indst = filter(inds,strat_name == st)
  raw_mt = filter(raw_means,strat_name == st)
  nsur = filter(raw,strat_name == st)
  ind_fac <- ggplot(data = indst,aes(x = year,y = mean))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = type),alpha = 0.2)+
  geom_line(aes(colour = type))+
  geom_point(data = raw_mt,
             aes(x = year,y = mean_count),
             alpha = 0.1,colour = "blue",
             size = 1,inherit.aes = FALSE)+
    geom_dotplot(data = nsur,aes(x = year),
                 inherit.aes = FALSE,binwidth = 1,
                 colour = grey(0.5),
                 fill = grey(0.5),
                 alpha = 0.1,
                 method = "histodot",dotsize = 0.5)+
  scale_y_continuous(limits = c(0,NA))+
  labs(title = st)
print(ind_fac)
}
dev.off()







# Trends ------------------------------------------------------------------









trend_90 = tr_func(nsmooth_samples,start_year = 1990)
trend_09 = tr_func(nsmooth_samples,start_year = 2009)
trend_70 = tr_func(nsmooth_samples,start_year = 1970)


pdf(paste0("Figures/",species_file,"trend_maps_spatial_GAMYE.pdf"),
    width = 11,
    height = 8.5)

for(yy in c(1970,1990,2009)){

  trends_t = tr_func(nsmooth_samples,start_year = yy)
  
  m1 = trend_map(trends = trends_t)
print(m1)


}
for(yy in c(seq(1970,2000,by = 10),2009)){
  trends_t = tr_func(nsmooth_samples,start_year = yy,
                     end_year = yy+10)
  
  m1 = trend_map(trends = trends_t)
  print(m1)
  
  
}

dev.off()

# RStan -------------------------------------------------------------------
# 
# mod.file = "models/gamye_iCAR.stan"
# 
# ## compile model
# slope_model = stan_model(file=mod.file)
# 
# ## run sampler on model, data
# slope_stanfit <- sampling(slope_model,
#                                data=stan_data,
#                                verbose=TRUE, refresh=25,
#                                chains=3, 
#                           iter=800,
#                                warmup=500,
#                           cores = 3,
#                                pars = parms,
#                                control = list(adapt_delta = 0.8,
#                                               max_treedepth = 14),
#                           init = init_def)
# 
# 
# 
# 
# save(list = c("slope_stanfit","stan_data","jags_data","model"),
#      file = paste0("output/",species,"full__gamye_iCAR.RData"))
# 
# loo_sum_rstan = loo(slope_stanfit)
# 
# 
# loo_sum_rstan
# 
# launch_shinystan(slope_stanfit) 
# 
# 
# 
# 
# n_draws = gather_draws(slope_stanfit,n[s,y])
# nsmooth_draws = gather_draws(slope_stanfit,nsmooth[s,y])
# 













