## plotting

### NOTE: this script depends on a development version of bbsBayes at the following repo: https://github.com/AdamCSmithCWS/bbsBayes/tree/testing_Stan 
## the specific changes in that branch are in the two functions
## generate_indices() and extract_index_data() (the second is called within generate_indices())


strat_sel <- "bbs_cws"


library(bbsBayes)
library(tidyverse)
library(cmdstanr)
library(patchwork)
source("Functions/animated_maps_function.R")
source("Functions/generate_map_data.R")
## overwriting two of the bbsBayes functions with versions that accept Stan output
source("Functions/extract-index-data-alt.R")
source("Functions/generate-indices-alt.R")


load("species_lists.RData") # loads objects created at the beginning of the script Fit_gamye_models_cws.R

models_sel = c("gamye","gamye_Spatial",
               "firstdiff_Spatial","firstdiff_NonHier",
               "firstdiff")

models_sel = c("gamye_Spatial",
               "firstdiff_Spatial")

species <- "Golden-winged Warbler"


output_dir <- "output" # Stan writes output to files as it samples. This is great because it's really stable, but the user needs to think about where to store that output
#output_dir <- "F:/bbsStanBayes/output" # Stan writes output to files as it samples. This is great because it's really stable, but the user needs to think about where to store that output


regs_to_estimate <- c("stratum","continental","national","prov_state","bcr")


  jj <- which(nrecs_sp[,"english"] == species)
  species_f <- as.character(nrecs_sp[jj,"species_file"])
  #species_f <- gsub(pattern = "[[:punct:]]|[[:blank:]]",replacement = "",species)
  trends_roll_out <- NULL
  
  trends_ann_out <- NULL
  
  ind_all_out <- NULL
  inds_all_out <- NULL
  trends_short_long <- NULL
  
for(model_sel in models_sel){
 
  out_base <- paste(species_f,model_sel,"BBS",sep = "_") # text string to identify the saved output from the Stan process unique to species and model, but probably something the user wants to control
  

if(!file.exists(paste0(output_dir,"/",out_base,"_Stan_fit.RData"))){next}

load(paste0(output_dir,"/",out_base,"_Stan_fit.RData"))

alt_n <- ifelse(grepl(model_sel, pattern = "gamye"),"nsmooth",NA)


# if(grepl(pattern = "firstdiff",model_sel)){alt_n <- NA}

ind <- generate_indices(jags_mod = stanfit,
                        jags_data = stan_data,
                        backend = "Stan",
                        stratify_by = strat_sel,
                        alternate_n = "n",
                        regions = regs_to_estimate)


if(!is.na(alt_n)){
  
inds <- generate_indices(jags_mod = stanfit,
                         jags_data = stan_data,
                         backend = "Stan",
                         stratify_by = strat_sel,
                         alternate_n = alt_n,
                         regions = regs_to_estimate)
inds_out_tmp <- inds$data_summary %>% 
  mutate(model = model_sel,
         species = species)

inds_all_out <- bind_rows(inds_all_out,inds_out_tmp)

}else{
  inds <- ind
}

ind_out_tmp <- ind$data_summary %>% 
  mutate(model = model_sel,
         species = species)

ind_all_out <- bind_rows(ind_all_out,ind_out_tmp)

trajs <- plot_indices(inds,
                      species = species,
                      add_observed_means = TRUE,
                      add_number_routes = TRUE)

trajshort <- plot_indices(inds,
                      species = species,
                      add_observed_means = TRUE,
                      add_number_routes = TRUE,
                      min_year = 2004)

pdf(file = paste0("Figures/",out_base,".pdf"),width = 11,height = 8.5)
for(i in names(trajs)){
  t1 <- trajs[[i]] 
  t2 <- trajshort[[i]]
  
  # st <- str_trunc(t1$labels$title, width = 8,
  #                 side = "left",
  #                 ellipsis = "")
  if(!is.na(alt_n)){
  n1 <- ind$data_summary %>% 
    mutate(Reg_traj = gsub(Region_alt,pattern = "[[:punct:]]",replacement = "")) %>% 
    filter(Reg_traj == gsub(i,pattern = "[[:punct:]]",replacement = "")) #,
  #Region_type == "stratum"
  t1 <- t1 +
    geom_ribbon(data = n1, aes(x = Year,y = Index,ymin = Index_q_0.025,ymax = Index_q_0.975),
                fill = grey(0.5),alpha = 0.2)+
    geom_line(data = n1, aes(x = Year,y = Index),
              colour = grey(0.5))
  
  n2 <- ind$data_summary %>% 
    filter(Region_alt == gsub(i,pattern = "_",replacement = "-"),
           #Region_type == "stratum",
           Year >= 2004)
  t2 <- t2 +
    geom_ribbon(data = n2, aes(x = Year,y = Index,ymin = Index_q_0.025,ymax = Index_q_0.975),
                fill = grey(0.5),alpha = 0.2)+
    geom_line(data = n2, aes(x = Year,y = Index),
              colour = grey(0.5))
  }
  print(t1 + t2)
  
}



trends <- generate_trends(inds)
trends_short <- generate_trends(inds,Min_year = 2011)
map <- generate_map(trends,select = TRUE,stratify_by = strat_sel,species = species)
mapshort <- generate_map(trends_short,select = TRUE,stratify_by = strat_sel,species = species)
print(map + mapshort)

trends <- trends %>% 
  mutate(model = model_sel,
         species = species,
         trend_type = "long-term")

trends_short <- trends_short %>% 
  mutate(model = model_sel,
         species = species,
         trend_type = "short-term")

trends_short_long <- bind_rows(trends_short_long,trends)
trends_short_long <- bind_rows(trends_short_long,trends_short)




starts <- c(1970,seq(1976,2011,by = 5))
maps <- vector("list",length = length(starts))
names(maps) <- paste(starts)

for(dd in starts){
  trends_10temp <- generate_trends(inds,Min_year = dd,Max_year = dd+10)
  maps[[paste(dd)]] <- generate_map(trends_10temp,select = TRUE,stratify_by = strat_sel,species = species)
  print(maps[[paste(dd)]])
}



# Rolling 10-year trends --------------------------------------------------

starts <- c(1970:2011)
if(!is.na(alt_n)){
  
  inds2 <- generate_indices(jags_mod = stanfit,
                           jags_data = stan_data,
                           backend = "Stan",
                           stratify_by = strat_sel,
                           alternate_n = alt_n,
                           regions = c("continental","national","prov_state"))
  
}else{
  inds2 <- generate_indices(jags_mod = stanfit,
                            jags_data = stan_data,
                            backend = "Stan",
                            stratify_by = strat_sel,
                            regions = c("continental","national","prov_state"))
}

trends_roll <- NULL
for(dd in starts){
  trends_10temp <- generate_trends(inds2,Min_year = dd,Max_year = dd+10)
  trends_roll <- bind_rows(trends_roll,trends_10temp)

}
trends_roll <- trends_roll %>% 
  mutate(model = model_sel,
         species = species)

trends_roll_out <- bind_rows(trends_roll_out,trends_roll)

t_roll <- trends_roll %>% 
  filter(Region_type %in% c("continental","national"))

t_roll_plot <- ggplot(data = t_roll,
                      aes(x = End_year,
                          y = Trend))+
  geom_errorbar(aes(ymin = Trend_Q0.025,
                    ymax = Trend_Q0.975),
                width = 0,alpha = 0.2)+
  geom_errorbar(aes(ymin = Trend_Q0.25,
                    ymax = Trend_Q0.75),
                width = 0,alpha = 0.2,size = 1.5)+
  geom_point()+
  labs(title = paste(species,model_sel,"Rolling ten-year trends through time"))+
  ylab("Trend %/year")+
  xlab("End of ten-year trend")+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = 100*(0.7^(1/10)-1),colour = "darkorange")+
  geom_hline(yintercept = 100*(0.5^(1/10)-1),colour = "darkred")+
  theme_bw()+
  facet_wrap(facets = vars(Region_alt),nrow = 1)

print(t_roll_plot)

dev.off()
# Animated annual trend maps

# 
# animated_trend_map(fit = stanfit,
#                                rawdata = stan_data,
#                                stratification = strat_sel,
#                                alt_n = alt_n,
#                                firstYear = NULL,
#                                lastYear = NULL,
#                    trend_length = 10,
#                                res_mag = 3,
#                                dir_out = "Figures/",
#                    file_name_prefix = paste0(model_sel,"_"),
#                                species = species)
# 

}

  trends_roll_out 
  
  trends_ann_out 
  
  ind_all_out
  inds_all_out
  
  

# overplot the trends by state country and survey-wide --------------

  provs <- c("AB","SK","BC","YT","NT","MB","NU",
             "ON","QC","NB","NS","PE","NL")
  # prov_arrange <- load_map(stratify_by = "state") %>% 
  #   as.data.frame() %>% 
  #   mutate(country = ifelse(ST_12 %in% provs,"Canada","United States")) %>% 
  #   rename(Region = ST_12)
  # 
  trend_longstate <- trends_short_long %>% 
    filter(Region_type == "prov_state",
           model == "gamye_Spatial",
           trend_type == "long-term") %>% 
    mutate(country = ifelse(Region %in% provs,"Canada","United States")) %>% 
    group_by(country) %>% 
    mutate(Region_order = reorder(Region,Trend,mean)) %>% 
    select(Region,Region_order,country)
  
  region_order <- c("Continental","CA","US",levels(trend_longstate$Region_order))
  
  allregs_t <- trends_short_long %>% 
    filter(Region_type %in% c("prov_state","continental","national"),
           model != "firstdiff_NonHier") %>% 
    mutate(trend_time = paste(Start_year,End_year,sep = "-"),
           Region_order = factor(Region,ordered = TRUE,levels = region_order)) %>% 
    arrange(Region_order) 
  
  contsp <- ggplot(data = allregs_t,aes(x = Region_order,y = Trend,colour = model))+
    geom_errorbar(aes(ymin = Trend_Q0.025,ymax = Trend_Q0.975),
                  width = 0,alpha = 0.5,
                  position = position_dodge(width = 0.5))+
    geom_errorbar(aes(ymin = Trend_Q0.25,ymax = Trend_Q0.75),
                  width = 0,alpha = 0.75,size = 0.8,
                  position = position_dodge(width = 0.5))+
    geom_point(position = position_dodge(width = 0.5),size = 0.85)+
    xlab("")+
    geom_hline(yintercept = 0)+
    geom_hline(yintercept = 100*(0.7^(1/10)-1),colour = "darkorange")+
    geom_hline(yintercept = 100*(0.5^(1/10)-1),colour = "darkred")+
    theme_bw()+
    coord_flip()+
    scale_colour_viridis_d(begin = 0.3,end = 0.9)+
    guides(colour = guide_legend(title = "Alternative Trend Models",reverse = TRUE))+
    facet_wrap(facet = vars(trend_time),ncol = 2,
               scales = "free_x")
 
  pdf(paste0("Figures/",species_f,"long_short_term_trends.pdf"),
      width = 8,
      height = 4.5) 
  print(contsp)
dev.off()  
  


# overplot the trajectories by same groupings -----------------------------------

allregs_i <- ind_all_out %>% 
  filter(Region_type %in% c("prov_state","continental","national"),
         model != "firstdiff_NonHier") %>% 
  mutate(Region_order = factor(Region,ordered = TRUE,levels = region_order),
         obs_mean = ifelse(nrts == 0,NA,obs_mean)) %>% 
  arrange(Region_order) 

provs_regs_i <- allregs_i %>% 
  filter(Region_type == "prov_state")

national_regs_i <- allregs_i %>% 
  filter(Region_type != "prov_state")

conttraj <- ggplot(data = national_regs_i,aes(x = Year,y = Index,group = model))+
  # geom_ribbon(aes(ymin = Index_q_0.025,ymax = Index_q_0.975,
  #                 fill = model),
  #               alpha = 0.2)+
  geom_ribbon(aes(ymin = Index_q_0.25,ymax = Index_q_0.75,
                    fill = model),
                alpha = 0.1)+
  geom_line(aes(colour = model))+
  xlab("")+
  theme_bw()+
  theme(legend.position = "bottom")+
  scale_y_continuous(limits = c(0,NA))+
  scale_colour_viridis_d(begin = 0.3,end = 0.9,aesthetics = c("fill","colour"),
                         guide = guide_legend(title = "Alternative Trend Models",reverse = TRUE))+
  facet_wrap(facet = vars(Region_order),ncol = 3,
             scales = "free_y")

pdf(paste0("Figures/",species_f,"_continental_national_trajectories.pdf"),
    width = 8,
    height = 4.5) 
print(conttraj)
dev.off()  

obs_means <- provs_regs_i %>% 
  filter(model == "gamye")


provtraj <- ggplot(data = provs_regs_i,aes(x = Year,y = Index,group = model))+
  # geom_ribbon(aes(ymin = Index_q_0.025,ymax = Index_q_0.975,
  #                 fill = model),
  #               alpha = 0.2)+
  geom_ribbon(aes(ymin = Index_q_0.25,ymax = Index_q_0.75,
                  fill = model),
              alpha = 0.1)+
  geom_line(aes(colour = model))+
  geom_point(data = obs_means,
             aes(x = Year,y = obs_mean,alpha = nrts),inherit.aes = FALSE)+
  xlab("")+
  theme_bw()+
  theme(legend.position = "bottom")+
  scale_y_continuous(limits = c(0,NA))+
  scale_colour_viridis_d(begin = 0.3,end = 0.9,aesthetics = c("fill","colour"),
                         guide = guide_legend(title = "Alternative Trend Models",reverse = TRUE))+
  facet_wrap(facet = vars(Region_order),ncol = 4,
             scales = "free_y")

pdf(paste0("Figures/",species_f,"_state_prov_trajectories.pdf"),
    width = 11,
    height = 8) 
print(provtraj)
dev.off()  



## panel rolling trends

t_roll <- trends_roll_out %>% 
  filter(Region_type %in% c("continental","national"),
         model != "firstdiff_NonHier") %>% 
  mutate(Region_order = factor(Region,ordered = TRUE,levels = c("Continental","CA","US"))) %>% 
  arrange(Region_order) 



t_roll_plot <- ggplot(data = t_roll,
                      aes(x = End_year,
                          y = Trend,
                          colour = model))+
  geom_errorbar(aes(ymin = Trend_Q0.025,
                    ymax = Trend_Q0.975),
                width = 0,alpha = 0.2,
                position = position_dodge(width = 0.2))+
  geom_errorbar(aes(ymin = Trend_Q0.25,
                    ymax = Trend_Q0.75),
                width = 0,alpha = 0.3,
                position = position_dodge(width = 0.2))+
  geom_point(position = position_dodge(width = 0.2))+
  labs(title = paste(species,model_sel,"Rolling ten-year trends through time"))+
  ylab("Trend %/year")+
  xlab("End of ten-year trend")+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = 100*(0.7^(1/10)-1),colour = "darkorange")+
  geom_hline(yintercept = 100*(0.5^(1/10)-1),colour = "darkred")+
  theme_bw()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(begin = 0.3,end = 0.9,aesthetics = c("fill","colour"),
                         guide = guide_legend(title = "Alternative Trend Models",reverse = TRUE))+
  facet_wrap(facets = vars(Region_order),nrow = 1)

#print(t_roll_plot)

pdf(paste0("Figures/",species_f,"_continental_national_rolling_trends.pdf"),
    width = 11,
    height = 4.5) 
print(t_roll_plot)
dev.off()  




## panel rolling trends by prov state

t_roll <- trends_roll_out %>% 
  filter(Region_type %in% c("prov_state"),
         model != "firstdiff_NonHier") %>% 
  mutate(Region_order = factor(Region,ordered = TRUE,levels = region_order)) %>% 
  arrange(Region_order) 



t_roll_plot <- ggplot(data = t_roll,
                      aes(x = End_year,
                          y = Trend,
                          colour = model))+
  geom_errorbar(aes(ymin = Trend_Q0.025,
                    ymax = Trend_Q0.975),
                width = 0,alpha = 0.2,
                position = position_dodge(width = 0.2))+
  geom_errorbar(aes(ymin = Trend_Q0.25,
                    ymax = Trend_Q0.75),
                width = 0,alpha = 0.3,
                position = position_dodge(width = 0.2))+
  geom_point(position = position_dodge(width = 0.2))+
  labs(title = paste(species,model_sel,"Rolling ten-year trends through time"))+
  ylab("Trend %/year")+
  xlab("End of ten-year trend")+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = 100*(0.7^(1/10)-1),colour = "darkorange")+
  geom_hline(yintercept = 100*(0.5^(1/10)-1),colour = "darkred")+
  theme_bw()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(begin = 0.3,end = 0.9,aesthetics = c("fill","colour"),
                         guide = guide_legend(title = "Alternative Trend Models",reverse = TRUE))+
  facet_wrap(facets = vars(Region_order),nrow = 3,
             scales = "free_y")

#print(t_roll_plot)

pdf(paste0("Figures/",species_f,"_state_prov_rolling_trends.pdf"),
    width = 11,
    height = 8.5) 
print(t_roll_plot)
dev.off()  



