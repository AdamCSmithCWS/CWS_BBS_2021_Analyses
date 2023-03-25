
## compile trend and index files for 
### 1 - website
### 2 - State of Canada's Birds
YYYY <- 2021

library(bbsBayes)
library(tidyverse)


source("functions/mapping.R")
source("functions/loess_func.R")


# Compile all trends and indices ------------------------------------------------------

load("species_lists.RData") # loads objects created at the beginning of the script Fit_gamye_models_cws.R

species_to_run <- nrecs_sp 

trends <- NULL  
indices <- NULL 
smooth_indices <- NULL
  
for(jj in 1:nrow(species_to_run)){
  species <- as.character(species_to_run[jj,"english"])
  species_f <- as.character(species_to_run[jj,"species_file"])
  espece <- as.character(species_to_run[jj,"french"])
  species_f_bil <- gsub(paste(species,espece),pattern = "[[:space:]]|[[:punct:]]",
                        replacement = "_")
  aou <- as.character(species_to_run[jj,"AOU"])
  
  
  
smooth_indices_1 <- read.csv(paste0("indices/smooth/",species_f_bil,"_smoothed_annual_indices.csv"))
indices_1 <- read.csv(paste0("indices/full/",species_f_bil,"_annual_indices.csv"))
trends_1 <- read.csv(paste0("trends/Trends_by_species/",species_f_bil,"_trends.csv"))

trends <- bind_rows(trends,trends_1)  
indices <- bind_rows(indices,indices_1)  
smooth_indices <- bind_rows(smooth_indices,smooth_indices_1)  

  
  
print(round(jj/nrow(species_to_run),2))
  
}# end species loop

numeric_cols_trends <- names(trends)[which(as.character(lapply(trends,typeof)) %in% c("double"))]

trends_round <- trends %>% 
  mutate(across(all_of(numeric_cols_trends),~signif(.,3)))

saveRDS(trends,"output/alltrends.rds")
write.csv(trends_round, paste0("website/All_",YYYY,"_BBS_trends.csv"),row.names = F)



numeric_cols_indices <- names(indices)[which(as.character(lapply(indices,typeof)) %in% c("double"))]

indices_round <- indices %>% 
  mutate(across(all_of(numeric_cols_indices),~signif(.,3)))

saveRDS(indices,"output/allindices.rds")
write.csv(indices_round, paste0("website/All_",YYYY,"_BBS_indices.csv"),row.names = F)

numeric_cols_smooth_indices <- names(smooth_indices)[which(as.character(lapply(smooth_indices,typeof)) %in% c("double"))]

smooth_indices_round <- smooth_indices %>% 
  mutate(across(all_of(numeric_cols_indices),~signif(.,3)))
saveRDS(smooth_indices,"output/allsmooth_indices.rds")
write.csv(smooth_indices_round, paste0("website/All_",YYYY,"_BBS_smooth_indices.csv"),row.names = F)



# Website trends ----------------------------------------------------------


web <- trends_round %>% 
  filter(For_web == TRUE,
         Region != "Continental",
         Region_type != "bcr") %>% 
  mutate(prob_decrease_0_25_percent = prob_decrease_0_percent-prob_decrease_25_percent,
         prob_decrease_25_50_percent = prob_decrease_0_percent - (prob_decrease_0_25_percent + prob_decrease_50_percent),
         prob_increase_0_33_percent = prob_increase_0_percent-prob_increase_33_percent,
         prob_increase_33_100_percent = prob_increase_0_percent - (prob_increase_0_33_percent + prob_increase_100_percent),
         mapfile = paste(bbs_num,Region_alt,Trend_Time,"map.png",sep = "_"))




# generate maps for CWS website -------------------------------------------

#loading some base maps for the website plots
canmap <- rgdal::readOGR(dsn = system.file("maps",
                                           package = "bbsBayes"),
                         layer = "BBS_CWS_strata",
                         verbose = FALSE)

canmap@data$country <- substr(canmap@data$ST_12,1,2)

canmap <- canmap[canmap$country == "CA",]

basmap = rgdal::readOGR(dsn = system.file("maps",
                                          package = "bbsBayes"),
                        layer = "BBS_USGS_strata",
                        verbose = FALSE)
basmap@data$country <- substr(basmap@data$ST_12,1,2)

basmap <- basmap[basmap$country == "CA",]

## looping through species to create the maps in folder webmaps
for(sp in unique(web$bbs_num)){
  dft <- web %>% 
    filter(bbs_num == sp,
           Trend_Time == "Long-term")
  generate_web_maps(dft,
                    canmap = canmap)
  dft <- web %>% 
    filter(bbs_num == sp,
           Trend_Time == "Short-term")
  generate_web_maps(dft,
                    canmap = canmap)
  
}
## maps have Canadian regions only and show the regions included in each trend





clout = c("bbs_num",
          "species",
          "espece",
          "Region_alt",
          "Trend_Time",
          "Start_year",
          "End_year",
          "Trend",
          "Trend_Q0.025",
          "Trend_Q0.975",
          "reliability",
          "Width_of_95_percent_Credible_Interval",
          "reliab.cov",
          "backcast_flag",
          "prob_decrease_0_percent",
          "prob_increase_0_percent",
          "prob_decrease_50_percent",
          "prob_decrease_25_50_percent",
          "prob_decrease_0_25_percent",
          "prob_increase_0_33_percent",
          "prob_increase_33_100_percent",
          "prob_increase_100_percent",
          "Percent_Change",
          "Percent_Change_Q0.025",
          "Percent_Change_Q0.975",
          "Number_of_Routes",
          "Strata_included",
          "Strata_excluded",
          "mapfile")

clnms = c("sp","species","espece","geo.area","trendtype",
          "startyear","endyear","trend",
          "llimit","ulimit","reliab.over",
          "reliab.prec","reliab.cov","reliab.pool",
          "p.decrease","p.increase","p.d50","pd50.25",
          "pd25.0","pi0.33","pi33.100","pi100",
          "percent.change","percent.change.llimit",
          "percent.change.ulimit","nroutesduringtrend",
          "strata.inc","st.excl.long","mapfile")

web = web[,clout]
names(web) = clnms

web_species <- read.csv("alt_data/BBS_AvianCore.csv")
names_match <- web %>% 
  select(species,espece,bbs_num) %>% 
  distinct()

miss_bbs_num <- names_match %>% 
  select(bbs_num,species) %>% 
  left_join(.,
            web_species,
            by = c("bbs_num" = "bbsNumber"),
            multiple = "all") %>% 
  arrange(bbs_num) %>% 
  filter(is.na(commonNameE))


miss_english_names <- names_match %>% 
  select(bbs_num,species,espece) %>% 
  left_join(.,
            web_species,
            by = c("species" = "commonNameE"),
            multiple = "all") %>% 
  arrange(bbs_num) %>% 
  filter(is.na(bbsNumber))



write.csv(web, paste0("website/",YYYY," BBS trends for website.csv"),row.names = F)

# Indices for website -----------------------------------------------------



webi_short <- indices_round %>% 
  filter(For_web == TRUE,
         Region != "Continental",
         Region_type != "bcr",
         Year > (YYYY-11)) %>% 
  mutate(Trend_Time = "Short-term")

webi <- indices_round %>% 
  filter(For_web == TRUE,
               Region != "Continental",
               Region_type != "bcr") %>% 
  bind_rows(.,webi_short)


clouti =  c("bbs_num",
            "species",
            "espece",
            "Region_alt",
            "Trend_Time",
            "Year",
            "Index",
            "Index_q_0.05",
            "Index_q_0.95")
clnmsi = c("sp","species","espece","geo.area","trendtype",
           "year","an.index",
           "llimit","ulimit")

webi = webi[,clouti]
names(webi) <- clnmsi


write.csv(webi, paste0("website/",YYYY," BBS indices for website.csv"),row.names = F)

webi_short = webi_short[,clouti]
names(webi_short) <- clnmsi

write.csv(webi_short, paste0("website/",YYYY," Short-term only BBS indices for website.csv"),row.names = F)


# SOCB Trends -------------------------------------------------------------

## export all trends that include some area within Canada, as well as the US trends
## see note from M-A, Teams chat Jan 26, 2023
# ok, here we go: we want the ability to see trends for Canada (of course), but also continental, US, BCR, prov/terr, and the intersections of BCR and prov/terr, because that will feed into the goals (ugh)

# trends <- readRDS("output/alltrends.rds")
# indices <- readRDS("output/allindices.rds")
#smooth_indices <- readRDS("output/allsmooth_indices.rds")

trends_out <- trends_round %>% 
   filter((For_web == TRUE | Region %in% c("Continental","US"))) #%>% 
  # mutate(prob_decrease_0_25_percent = prob_decrease_0_percent-prob_decrease_25_percent,
  #        prob_decrease_25_50_percent = prob_decrease_0_percent - (prob_decrease_0_25_percent + prob_decrease_50_percent),
  #        prob_increase_0_33_percent = prob_increase_0_percent-prob_increase_33_percent,
  #        prob_increase_33_100_percent = prob_increase_0_percent - (prob_increase_0_33_percent + prob_increase_100_percent))
socb_headings <- c("results_code",
                   "version",
                   "area_code",
                   "species_code",
                   "species_id",
                   "season",
                   "period",
                   "years",
                   "year_start",
                   "year_end",
                   "trnd",
                   "index_type",
                   "lower_ci",
                   "upper_ci",
                   "stderr",
                   "model_type",
                   "model_fit",
                   "percent_change",
                   "percent_change_low",
                   "percent_change_high",
                   "prob_decrease_0",
                   "prob_decrease_25",
                   "prob_decrease_30",
                   "prob_decrease_50",
                   "prob_increase_0",
                   "prob_increase_33",
                   "prob_increase_100",
                   "reliability",
                   "precision_num",
                   "precision_cat",
                   "coverage_num",
                   "coverage_cat",
                   "goal",
                   "goal_lower",
                   "sample_size",
                   "sample_total",
                   "subtitle",
                   "prob_LD",
                   "prob_MD",
                   "prob_LC",
                   "prob_MI",
                   "prob_LI")

trend_headings_match <- c("",
                          "",
                          "Region_alt",
                          "",
                          "espece",
                          "species",
                          "Trend_Time",
                          "",
                          "Start_year",
                          "End_year",
                          "Trend",
                          "",
                          "Trend_Q0.025",
                          "Trend_Q0.975",
                          "",
                          "",
                          "",
                          "Percent_Change",
                          "Percent_Change_Q0.025",
                          "Percent_Change_Q0.975",
                          "prob_decrease_0_percent",
                          "prob_decrease_25_percent",
                          "prob_decrease_30_percent",
                          "prob_decrease_50_percent",
                          "prob_increase_0_percent",
                          "prob_increase_33_percent",
                          "prob_increase_100_percent",
                          "reliability",
                          "Width_of_95_percent_Credible_Interval",
                          "precision",
                          "reliab.cov",
                          "coverage",
                          "",
                          "",
                          "Mean_Number_of_Routes",
                          "Number_of_Routes",
                          "",
                          "",
                          "",
                          "",
                          "",
                          "")

socb_headings_extract <- data.frame(socb = socb_headings,
                                    trend = trend_headings_match)

trends_select <- trend_headings_match[-which(trend_headings_match == "")]

socb_headings_select <- socb_headings_extract %>% 
  filter(trend != "")

trends_socb <- trends_out %>% 
  select(any_of(trends_select))

if(any(names(trends_socb) != socb_headings_select[,2])) stop("Stop columns don't match")

write.csv(trends_socb,paste0("website/BBS_",YYYY,"_trends_for_socb.csv"))

socb_headings_extract <- socb_headings_extract %>% 
  rename(NatureCountsTrendsSample = socb,
         BBS_trend_headers = trend)

write.csv(socb_headings_extract,
          "website/linking_columns_naturecounts_bbs.csv")




# SOCB indices ------------------------------------------------------------

indices_socb <- indices_round %>% 
  filter((For_web == TRUE | Region %in% c("Continental","US"))) %>% 
  group_by(species,Region,Trend_Time) %>% 
  mutate(LOESS_index = loess_func(Index,Year)) %>% 
  select(Year,Region,Index,Trend_Time,
         Index_q_0.05,Index_q_0.95,
         species,espece,bbs_num,
         LOESS_index) %>% 
  rename(upper_ci = Index_q_0.95,
         lower_ci = Index_q_0.05) 

write.csv(indices_socb,
          file = paste0("website/BBS_",YYYY,"_annual_indices_for_socb.csv"))

