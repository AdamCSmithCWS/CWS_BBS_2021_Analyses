## building a Stan version of the bbsBayes models

library(bbsBayes)
library(tidyverse)
library(shinystan)
library(sf)
library(spdep)
library(cmdstanr)
library(posterior)

source("functions/mungeCARdata4stan.R") ## function to modify the BUGS formatted spatial neighbourhood data to the required format for the Stan iCAR model
source("functions/prepare-jags-data-alt.R") ## function to modify the BUGS formatted spatial neighbourhood data to the required format for the Stan iCAR model
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


strat = "bbs_usgs"
model = "gamye"

strat_data = stratify(by = strat)

# load and stratify species data ---------------------------------------------
species = "Northwestern Crow"

shiny_explore = FALSE # set to TRUE to automatically launch shinystan on model finish

first_year = 1990

jags_data = prepare_jags_data(strat_data = strat_data,
                             species_to_run = species,
                             model = model,
                             min_max_route_years = 3,
                             #n_knots = 4,
                             min_year = first_year)


str_link <- unique(data.frame(strat = jags_data$strat,
                              strat_name = jags_data$strat_name))
str_link <- str_link %>% arrange(strat)

realized_strata_map <- inner_join(strata_map,str_link,by = c("ST_12" = "strat_name")) %>% 
  arrange(strat)


stan_data = jags_data[c("ncounts",
                         "nstrata",
                         "count",
                         "strat",
                         "year",
                         "firstyr",
                        "nonzeroweight")]
 stan_data[["observer"]] <- as.integer(factor(paste0(jags_data[["ObsN"]],"-",jags_data[["route"]])))
  stan_data[["nobservers"]] <- max(stan_data[["observer"]])

#stan_data[["route"]] <- as.integer(factor(jags_data$route))
#stan_data[["nroutes"]] <- max(stan_data[["route"]])
strat_routeobs <- unique(data.frame(strat = stan_data$strat,
                                 obs = stan_data$observer))
strat_routeobs <- strat_routeobs %>% arrange(strat)

nobservers_strata <- as.integer(table(strat_routeobs$strat))
maxnobservers_strata <- max(nobservers_strata)
stan_data[["nobservers_strata"]] <- nobservers_strata
stan_data[["maxnobservers_strata"]] <- maxnobservers_strata

rteobs_mat <- matrix(data = 0,
                  nrow = stan_data$nstrata,
                  ncol = maxnobservers_strata)
for(i in 1:stan_data$nstrata){
  rteobs_mat[i,1:nobservers_strata[i]] <- strat_routeobs[which(strat_routeobs$strat == i),"obs"]
}

stan_data[["rteobs_mat"]] <- rteobs_mat

stan_data[["nyears"]] <- max(jags_data$year)
stan_data[["nknots_year"]] <- jags_data$nknots
stan_data[["year_basispred"]] <- jags_data$X.basis



ncounts = stan_data$ncounts
#nroutes = stan_data$nroutes
nstrata = stan_data$nstrata
nyears = stan_data$nyears
nknots_year = stan_data$nknots_year
nobservers = stan_data$nobservers





# # cmdStanR ----------------------------------------------------------------
mod.file = "models/gamye_jagsmatch_inform_prior_year.stan"

## compile model
model <- cmdstan_model(mod.file)

init_def <- function(){ list(noise_raw = rnorm(ncounts,0,0.1),
                            strata_raw = rnorm(nstrata,0,0.1),
                            STRATA = 0,
                            eta = 0,
                            yeareffect_raw = matrix(rnorm(nstrata*nyears,0,0.1),nrow = nstrata,ncol = nyears),
                            obs_raw = rnorm(nobservers,0,0.1),
                            #rte_raw = rnorm(nroutes,0,0.1),
                            sdnoise = 0.2,
                            sdobs = 0.1,
                            #sdrte = 0.2,
                            sdbeta = runif(nknots_year,0.01,0.1),
                            sdBETA = 0.1,
                            # logsdyear = runif(nstrata,-0.01,0.1),
                            sdyear = runif(nstrata,0,0.2),
                            #nu = 10,
                            BETA_raw = rnorm(nknots_year,0,0.1),
                            beta_raw = matrix(rnorm(nknots_year*nstrata,0,0.01),nrow = nstrata,ncol = nknots_year))}


species_file = gsub(pattern = "([[:punct:]]|[[:blank:]])","",species)
out_base <- paste0(species_file,"_",nyears)
## run sampler on model, data
# data_file <- "tmpdata/tmp_data.json"
# write_stan_json(stan_data,file = data_file)
model_stanfit <- model$sample(
                          data=stan_data,
                          refresh=100,
                          chains=3, iter_sampling=1000,
                          iter_warmup=1000,
                          parallel_chains = 3,
                          #pars = parms,
                          adapt_delta = 0.8,
                          max_treedepth = 14,
                          seed = 123,
                          init = init_def,
                          output_dir = "output",
                          output_basename = out_base)



summr = as_draws_df(model_stanfit$draws(format = "draws_list" ,
                                      variables = c("BETA",
                                                    "beta",
                                                    "sdbeta",
                                                    "sdBETA"))) %>% 
  summary()

save(list = c("stan_data","jags_data","model_stanfit"),
     file = paste0("output/cmdStan_",species_file,"_",nyears,"_gamye_inform_prior.RData"))



# export to csv and read in as rstan --------------------------------------

csv_files <- dir("output",pattern = out_base,full.names = TRUE)
csv_files <- csv_files[grepl(csv_files,pattern = ".csv",fixed = TRUE)]

#sl_rstan <- As.mcmc.list(read_stan_csv(csv_files))
if(shiny_explore){
sl_rstan <- rstan::read_stan_csv(csv_files)
launch_shinystan(as.shinystan(sl_rstan))
}
# mod_sum = slope_stanfit$summary(variables = "n")


# 
# loo_sum = slope_stanfit$loo()
# ploo = as.data.frame(loo_sum$pointwise)
# ploo$k_cat <- cut(ploo$influence_pareto_k,breaks = c(-Inf,0.5,0.7,1,Inf),levels = c("good","ok","bad","very_bad"))
# df = data.frame(count = stan_data$count,
#                 route = stan_data$route,
#                 observer = stan_data$observer,
#                 strat = stan_data$strat,
#                 year = stan_data$year)
# ploo = bind_cols(ploo,df)
# nyear_sqrt <- ceiling(sqrt(length(unique(ploo$year))))
# 
# pdf(paste0("figures/",species_file,"loo_plots_standard.pdf"),
#     width = 11,
#     height = 8.5)
# 
# 
# k_count = ggplot(data = ploo,aes(x = count,y = influence_pareto_k,colour = k_cat))+
#   geom_point()+
#   facet_wrap(~year,nrow = nyear_sqrt,scales = "fixed")+
# geom_abline(slope = 0,intercept = 0.7,colour = "red")
# print(k_count)
# 
# looic_count = ggplot(data = ploo,aes(x = count,y = looic,colour = k_cat))+
#   geom_point()+
#   facet_wrap(~year,nrow = nyear_sqrt,scales = "fixed")
# print(looic_count)
# 
# 
# k_looic = ggplot(data = ploo,aes(x = looic,y = influence_pareto_k,colour = k_cat))+
#   geom_point()+
#   facet_wrap(~year,nrow = nyear_sqrt,scales = "fixed")+
#   geom_abline(slope = 0,intercept = 0.7,colour = "red")
# print(k_looic)
# 
# looic_strat = ggplot(data = ploo,aes(x = strat,y = looic,colour = k_cat))+
#   geom_point()+
#   facet_wrap(~year,nrow = nyear_sqrt,scales = "fixed")
# print(looic_strat)
# 
# strat_loo_sum <- ploo %>% group_by(strat) %>% 
#   summarise(m_k = mean(influence_pareto_k),
#             m_looic = mean(looic),
#             m_ploo = mean(p_loo),
#             md_k = median(influence_pareto_k),
#             md_looic = median(looic),
#             md_ploo = median(p_loo))
# 
# map_loo <- inner_join(realized_strata_map,strat_loo_sum,by = "strat")
# plot_map_loo = ggplot(data = map_loo,aes(fill = md_looic))+
#   geom_sf()
# print(plot_map_loo)
# 
# plot_map_k = ggplot(data = map_loo,aes(fill = md_k))+
#   geom_sf()
# print(plot_map_k)
# 
# plot_map_ploo = ggplot(data = map_loo,aes(fill = md_ploo))+
#   geom_sf()
# print(plot_map_ploo)
# 
# dev.off()

source("functions/posterior_summary_functions.R")
source("functions/trend_function.R")

# index and trend plotting ------------------------------------------------

nsmooth_samples <- posterior_samples(fit = model_stanfit,
                                  parm = "nsmooth",
                                  dims = c("s","y")) %>% 
  left_join(.,str_link,by = c("s" = "strat")) %>% 
  mutate(year = y+(first_year-1))


n_samples <- posterior_samples(fit = model_stanfit,
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


# 
# ind_fac <- ggplot(data = inds,aes(x = year,y = mean))+
#   geom_ribbon(aes(ymin = lci,ymax = uci,fill = type),alpha = 0.2)+
#   geom_line(aes(colour = type))+
#   geom_point(data = raw_means,aes(x = year,y = mean_count,fill = n_surveys),alpha = 0.2,inherit.aes = FALSE)+
#   scale_y_continuous(limits = c(0,NA))+
#   facet_wrap(~strat_name,nrow = nf,scales = "free")
#   
#print(ind_fac)


# continental indices -----------------------------------------------------

a_weights <- as.data.frame(realized_strata_map) %>% 
  rename(s = strat) %>% 
  mutate(area = AREA_1/sum(AREA_1)) %>% 
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

pdf(file = paste0("figures/",species_file,"trajectories_standard.pdf"))

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





# bbsBayes JAGS verions of same model -------------------------------------

param = c("n","nsmooth","sdbeta",
          "beta","BETA","sdBETA")


gamye_sdpriors = run_model(jags_data = jags_data,
                           parameters_to_save = c(param),
                           n_burnin = 20000,
                           parallel = TRUE,
                           n_thin = 40,
                           n_iter = 40000,
                           model_file_path = "models/gamye_sdpriors.R",
                           modules = NULL)
my_sso <- shinystan::launch_shinystan(
  shinystan::as.shinystan(gamye_sdpriors$samples, 
                          model_name = "gamye_sdpriors"))

inds = generate_indices(gamye_sdpriors,jags_data = jags_data)
indsmooth = generate_indices(gamye_sdpriors,
                             jags_data = jags_data,
                             alternate_n = "nsmooth")
plots = plot_indices(inds)
library(tidyverse)
pdf("gamye_sdpriors_JAGS_indices.pdf")
for(j in 1:length(plots)){
  tmpn <- gsub(names(plots)[j],pattern = "_", replacement = "-")
  tmp = indsmooth$data_summary %>% filter(Region == tmpn)
  tmpp <- plots[[j]]+
    geom_ribbon(data = tmp,aes(x = Year,ymax = Index_q_0.95,ymin = Index_q_0.025),
                alpha = 0.2,fill = "red")+
    geom_line(data= tmp,aes(x = Year,y = Index),colour = "red")
  print(tmpp)
}
dev.off()



# Trends ------------------------------------------------------------------










trend_90 = tr_func(nsmooth_samples,start_year = 1990)
trend_09 = tr_func(nsmooth_samples,start_year = 2009)
trend_70 = tr_func(nsmooth_samples,start_year = 1970)


pdf(paste0("Figures/",species_file,"trend_maps_standard_GAMYE.pdf"),
    width = 11,
    height = 8.5)

for(yy in c(1970,1990,2009)){

  trends_t = tr_func(nsmooth_samples,start_year = yy)
  
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













