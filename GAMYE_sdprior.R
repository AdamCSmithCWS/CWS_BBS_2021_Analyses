
library(bbsBayes)
species = "Fox Sparrow"
model = "gamye"
strat = "bbs_usgs"
library(shinystan)

strat_data = stratify(by = strat)
param = c("n","nsmooth","sdbeta",
          "beta","BETA","sdBETA")




# 
# jags_data = prepare_data(strat_data = strat_data,species_to_run = species,model = model,n_knots = 13)
# 
# mod = run_model(jags_data = jags_data,parameters_to_save = c("n","n3"),
#                 parallel = TRUE)
# fin_val <- get_final_values(mod)
# 
# 
# mod2 = run_model(jags_data = jags_data,parameters_to_save = c("n","n3"),
#                 parallel = TRUE,
#                 inits = fin_val,n_burnin = 0,
#                 n_thin = 30,
#                 n_iter = 20000)
# 
# 
# 
# 


# alternate basis fxn -----------------------------------------------------

source("functions/prepare-jags-data-alt.R")



jags_data_alt <- prepare_jags_data_alt(strat_data = strat_data,
                                       species_to_run = species,
                                       model = model,
                                       n_knots = 13,
                                       heavy_tailed = TRUE,
                                       min_max_route_years = 3)




gamye_sdpriors = run_model(jags_data = jags_data_alt,
                            parameters_to_save = c(param),
                            n_burnin = 20000,
                            parallel = TRUE,
                            n_thin = 40,
                            n_iter = 40000,
                            model_file_path = "gamye_sdpriors.R",
                            modules = NULL)
my_sso <- shinystan::launch_shinystan(
  shinystan::as.shinystan(gamye_sdpriors$samples, 
                          model_name = "gamye_sdpriors"))

inds = generate_indices(mod_alt_constr7,jags_data = jags_data_alt)
indsmooth = generate_indices(mod_alt_constr7,
                             jags_data = jags_data_alt,
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
