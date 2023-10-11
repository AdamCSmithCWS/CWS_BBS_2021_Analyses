
## compile trend and index files for 
### 1 - website
### 2 - State of Canada's Birds
YYYY <- 2021

webmaps <- FALSE # set to true if needing to create all map images for ECCC website

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
  mutate(across(all_of(numeric_cols_trends),~signif(.,3))) %>% 
  relocate(species,espece,bbs_num,Trend_Time)

saveRDS(trends_round,"output/alltrends.rds")
readr::write_excel_csv(trends_round, paste0("website/All_",YYYY,"_BBS_trends.csv"))

trends_high_level <- trends_round %>% 
  filter(Region_type %in% c("national","continental"))
readr::write_excel_csv(trends_high_level, 
          paste0("website/","continent_national_trends_cws_BBS_",YYYY,".csv"))

numeric_cols_indices <- names(indices)[which(as.character(lapply(indices,typeof)) %in% c("double"))]

indices_round <- indices %>% 
  mutate(across(all_of(numeric_cols_indices),~signif(.,3))) %>% 
  relocate(species,espece,bbs_num,Trend_Time)

saveRDS(indices_round,"output/allindices.rds")
readr::write_excel_csv(indices_round, paste0("website/All_",YYYY,"_BBS_indices.csv"))

indices_high_level <- indices_round %>% 
  filter(Region_type %in% c("national","continental"))
readr::write_excel_csv(indices_high_level, 
                       paste0("website/","continent_national_indicess_cws_BBS_",YYYY,".csv"))



numeric_cols_smooth_indices <- names(smooth_indices)[which(as.character(lapply(smooth_indices,typeof)) %in% c("double"))]

smooth_indices_round <- smooth_indices %>% 
  mutate(across(all_of(numeric_cols_indices),~signif(.,3))) %>% 
  relocate(species,espece,bbs_num,Trend_Time)
saveRDS(smooth_indices_round,"output/allsmooth_indices.rds")
readr::write_excel_csv(smooth_indices_round, 
                       paste0("website/All_",YYYY,"_BBS_smooth_indices.csv"))



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

if(nrow(miss_bbs_num) > 0){
  warning("At least one bbs number is missing from Avian Core")
  
  print(paste("Avian core is missing",
              paste(miss_bbs_num$bbs_num,
                    collapse = ", ")))
  web <- web %>% 
    filter(bbs_num %in% web_species$bbsNumber)
}


miss_english_names <- names_match %>% 
  select(bbs_num,species,espece) %>% 
  left_join(.,
            web_species,
            by = c("species" = "commonNameE"),
            multiple = "all") %>% 
  arrange(bbs_num) %>% 
  filter(is.na(bbs_num))


if(nrow(miss_english_names) > 0){
  warning("At least one species name is missing from Avian Core")
  web <- web %>% 
    filter(species %in% web_species$commonNameE)

}


if(webmaps){

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

map_dup_test <- any(duplicated(web$mapfile))

test_map <- any(!file.exists(paste0("website/webmaps/",web$mapfile)))

if(test_map){
  stop("At least one map is missing")
}

maps <- list.files("website/webmaps/",
                   pattern = ".png")
test_map_extra <- any((web$mapfile %in% maps) == FALSE)

if(test_map_extra){
  warning("There are extra maps in the webmaps folder")
}

}#end webmaps

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


readr::write_excel_csv(web, paste0("website/",YYYY," BBS trends for website.csv"))











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

webi <- webi %>% 
  filter(bbs_num %in% web_species$bbsNumber)

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



tshort = web %>% 
  filter(trendtype == "Short-term")
tlong= web %>% 
  filter(trendtype == "Long-term")
index_short <- webi %>% 
  filter(year == YYYY-10,
         trendtype == "Short-term")
index_long <- webi %>% 
  filter(year == YYYY-10,
         trendtype == "Long-term")
  
if(nrow(index_long) != nrow(tlong) |
   nrow(index_short) != nrow(tshort)){
warning("The number of indices and trends don't match \n explore index_trend_test")  

index_trend_test <- webi %>% 
  filter(year == YYYY-10) %>% 
  select(species,trendtype,geo.area,an.index) %>% 
  full_join(.,web,
            by = c("species","trendtype","geo.area"))
}

readr::write_excel_csv(webi, paste0("website/",YYYY," BBS indices for website.csv"))

# webi_short = webi_short[,clouti]
# names(webi_short) <- clnmsi
# 
# readr::write_excel_csv(webi_short, paste0("website/",YYYY," Short-term only BBS indices for website.csv"))













# SOCB Trends -------------------------------------------------------------

## export all trends that include some area within Canada, as well as the US trends
## see note from M-A, Teams chat Jan 26, 2023
# ok, here we go: we want the ability to see trends for Canada (of course), but also continental, US, BCR, prov/terr, and the intersections of BCR and prov/terr, because that will feed into the goals (ugh)

 trends_round <- readRDS("output/alltrends.rds")
trends_3g <- readRDS("output/all_3_generation_trends.rds") %>% 
  mutate(Trend_Time = "Three-generations")
trends_round <- bind_rows(trends_round,
                          trends_3g)
# indices_round <- readRDS("output/allindices.rds")
# smooth_indices_round <- readRDS("output/allsmooth_indices.rds")


# tlong <- trends_round %>% 
#   filter(Trend_Time == "Long-term") %>% 
#   group_by(species) %>% 
#   summarise(earliest = min(Start_year),
#             latest = max(Start_year))

trends_out <- trends_round %>% 
   filter((For_web == TRUE | Region %in% c("Continental","US"))) %>% 
  mutate(prob_LD = prob_decrease_50_percent,
         prob_MD = prob_decrease_25_percent - prob_decrease_50_percent,
         prob_LC = (prob_decrease_0_percent-prob_decrease_25_percent)+(prob_increase_0_percent-prob_increase_33_percent) ,
         prob_MI = prob_increase_33_percent - prob_increase_100_percent,
         prob_LI = prob_increase_100_percent)

test_probs <- trends_out %>% 
  mutate(prob_test = prob_LD+prob_MD+prob_LC+prob_MI+prob_LI)
  
if(any(round(test_probs$prob_test,2) != 1)){stop("probabilites of change categories don't sum properly")}


trends_out <- trends_out %>%
  mutate(years = paste(Start_year,End_year,sep = "-"),
         results_code = "BBS",
         season = "breeding",
         version = YYYY,
         species_id = "",
         area_code = ifelse(Region_type == "prov_state",Region,Region_alt),
         area_code = gsub(area_code,pattern = "United States of America",
                          replacement = "USA"),
         model_type = "GAMYE") %>% 
  rename(species_code = bbs_num,
         species_name = species,
         period = Trend_Time,
         year_start = Start_year,
         year_end = End_year,
         trnd = Trend,
         lower_ci = Trend_Q0.025,
         upper_ci = Trend_Q0.975,
         percent_change = Percent_Change,
         percent_change_low = Percent_Change_Q0.025,
         percent_change_high = Percent_Change_Q0.975,
         prob_decrease_0 = prob_decrease_0_percent,
         prob_decrease_25 = prob_decrease_25_percent,
         prob_decrease_30 = prob_decrease_30_percent,
         prob_decrease_50 = prob_decrease_50_percent,
         prob_increase_0 = prob_increase_0_percent,
         prob_increase_33 = prob_increase_33_percent,
         prob_increase_100 = prob_increase_100_percent,
         confidence = reliability,
         precision_num = Width_of_95_percent_Credible_Interval,
         precision_cat = precision,
         coverage_num = reliab.cov,
         coverage_cat = coverage,
         sample_size = Mean_Number_of_Routes,
         sample_total = Number_of_Routes,
         prob_LD = prob_LD,
         prob_MD = prob_MD,
         prob_LC = prob_LC,
         prob_MI = prob_MI,
         prob_LI = prob_LI)




trends_socb <- trends_out %>% 
  relocate(results_code,
           season,
           version,
           model_type,
           area_code,
         species_code,
         species_id,
         species_name,
         period,
         years,
         year_start,
         year_end,
         trnd,
         lower_ci,
         upper_ci,
         percent_change,
         percent_change_low,
         percent_change_high,
         prob_decrease_0,
         prob_decrease_25,
         prob_decrease_30,
         prob_decrease_50,
         prob_increase_0,
         prob_increase_33,
         prob_increase_100,
         confidence,
         precision_num,
         precision_cat,
         coverage_num,
         coverage_cat,
         sample_size,
         sample_total,
         prob_LD,
         prob_MD,
         prob_LC,
         prob_MI,
         prob_LI)

readr::write_excel_csv(trends_socb,
          paste0("website/BBS_",YYYY,"_trends_for_socb.csv"))

# readr::write_excel_csv(trends_socb[1:300,],
#           paste0("website/sample_BBS_",YYYY,"_trends_for_socb.csv"))
# 




# 
# 
# socb_headings <- c("results_code",
#                    "version",
#                    "area_code",
#                    "species_code",
#                    "species_id",
#                    "season",
#                    "period",
#                    "years",
#                    "year_start",
#                    "year_end",
#                    "trnd",
#                    "index_type",
#                    "lower_ci",
#                    "upper_ci",
#                    "stderr",
#                    "model_type",
#                    "model_fit",
#                    "percent_change",
#                    "percent_change_low",
#                    "percent_change_high",
#                    "prob_decrease_0",
#                    "prob_decrease_25",
#                    "prob_decrease_30",
#                    "prob_decrease_50",
#                    "prob_increase_0",
#                    "prob_increase_33",
#                    "prob_increase_100",
#                    "reliability",
#                    "precision_num",
#                    "precision_cat",
#                    "coverage_num",
#                    "coverage_cat",
#                    "goal",
#                    "goal_lower",
#                    "sample_size",
#                    "sample_total",
#                    "subtitle",
#                    "prob_LD",
#                    "prob_MD",
#                    "prob_LC",
#                    "prob_MI",
#                    "prob_LI")
# 
# trend_headings_match <- c("",
#                           "",
#                           "Region_alt",
#                           "bbs_num",
#                           "species",
#                           "Trend_Time",
#                           "",
#                           "Start_year",
#                           "End_year",
#                           "Trend",
#                           "",
#                           "Trend_Q0.025",
#                           "Trend_Q0.975",
#                           "",
#                           "",
#                           "",
#                           "Percent_Change",
#                           "Percent_Change_Q0.025",
#                           "Percent_Change_Q0.975",
#                           "prob_decrease_0_percent",
#                           "prob_decrease_25_percent",
#                           "prob_decrease_30_percent",
#                           "prob_decrease_50_percent",
#                           "prob_increase_0_percent",
#                           "prob_increase_33_percent",
#                           "prob_increase_100_percent",
#                           "reliability",
#                           "Width_of_95_percent_Credible_Interval",
#                           "precision",
#                           "reliab.cov",
#                           "coverage",
#                           "",
#                           "",
#                           "Mean_Number_of_Routes",
#                           "Number_of_Routes",
#                           "",
#                           "prob_LD",
#                           "prob_MD",
#                           "prob_LC",
#                           "prob_MI",
#                           "prob_LI")
# 
# socb_headings_extract <- data.frame(socb = socb_headings,
#                                     trend = trend_headings_match)
# 
# trends_select <- trend_headings_match[-which(trend_headings_match == "")]
# 
# socb_headings_select <- socb_headings_extract %>% 
#   filter(trend != "")
# 
# trends_socb <- trends_out %>% 
#   select(any_of(trends_select))
# 
# if(any(names(trends_socb) != socb_headings_select[,2])) stop("Stop columns don't match")
# 
# names(trends_socb) <- socb_headings_select[,1]
# readr::write_excel_csv(trends_socb,paste0("website/BBS_",YYYY,"_trends_for_socb.csv"))
# 
# socb_headings_extract <- socb_headings_extract %>% 
#   rename(NatureCountsTrendsSample = socb,
#          BBS_trend_headers = trend)
# 
# readr::write_excel_csv(socb_headings_extract,
#           "website/linking_columns_naturecounts_bbs.csv")
# 



# SOCB indices ------------------------------------------------------------

indices_round <- readRDS("output/allindices.rds")
smooth_indices_round <- readRDS("output/allsmooth_indices.rds")
smooth_join <- smooth_indices_round %>% 
  select(species,Region,Region_type,Trend_Time,
         Year,Index) %>% 
  rename(smooth_index = Index)

indices_socb <- indices_round %>% 
  filter((For_web == TRUE | Region %in% c("Continental","US"))) %>% 
  inner_join(.,smooth_join,
             by = c("species",
                    "Region",
                    "Region_type",
                    "Trend_Time",
                    "Year")) %>% 
  group_by(species,Region,Region_type,Trend_Time) %>% 
  mutate(LOESS_index = loess_func(Index,Year),
         area_code = ifelse(Region_type == "prov_state",Region,Region_alt),
         area_code = gsub(area_code,pattern = "United States of America",
                          replacement = "USA"),
         results_code = "BBS",
         season = "breeding",
         version = YYYY,
         model_type = "GAMYE") %>% 
  ungroup() %>% 
  rename(species_code = bbs_num,
    species_id = species,
    index = Index,
    year = Year,
         period = Trend_Time,
         upper_ci = Index_q_0.95,
         lower_ci = Index_q_0.05) %>% 
  select(-c(Index_q_0.025,
            Index_q_0.975,
            Region_type,
            Region)) %>% 
  relocate(results_code,
           season,
           version,
           model_type,
           area_code,
           year,
           period, 
           species_code,
           species_id,
           index,
           upper_ci,
           lower_ci,
           LOESS_index,
           smooth_index)

readr::write_excel_csv(indices_socb,
          file = paste0("website/BBS_",YYYY,"_annual_indices_for_socb.csv"))

# readr::write_excel_csv(indices_socb[sample(1:nrow(indices_socb),100,FALSE),],"sample_indices_output_bbs.csv")

