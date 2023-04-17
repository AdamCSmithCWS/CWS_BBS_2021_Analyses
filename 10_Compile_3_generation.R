
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

  
for(jj in 1:nrow(species_to_run)){
  species <- as.character(species_to_run[jj,"english"])
  species_f <- as.character(species_to_run[jj,"species_file"])
  espece <- as.character(species_to_run[jj,"french"])
  species_f_bil <- gsub(paste(species,espece),pattern = "[[:space:]]|[[:punct:]]",
                        replacement = "_")
  aou <- as.character(species_to_run[jj,"AOU"])
  
  
  
trends_1 <- read.csv(paste0("trends/Three_generation/",species_f_bil,"_trends.csv"))

trends <- bind_rows(trends,trends_1)  

  
  
print(round(jj/nrow(species_to_run),2))
  
}# end species loop

numeric_cols_trends <- names(trends)[which(as.character(lapply(trends,typeof)) %in% c("double"))]

trends_round <- trends %>% 
  mutate(across(all_of(numeric_cols_trends),~signif(.,3))) %>% 
  relocate(species,espece,bbs_num,Trend_Time)

saveRDS(trends_round,"output/all_trends3_generation_short.rds")
write.csv(trends_round, paste0("website/All_",YYYY,"_BBS_trends_3_generation_short.csv"),row.names = F)

trends_high_level <- trends_round %>% 
  filter(Region_type %in% c("national","continental"))
write.csv(trends_high_level, paste0("website/","continent_national_3_generation_short_trends_cws_BBS_",YYYY,".csv"),row.names = F)



trends_round_short <- trends_round %>% 
  filter(Trend_Time == "Short-term")

saveRDS(trends_round_short,"output/all_3_generation_trends.rds")
write.csv(trends_round_short, paste0("website/All_",YYYY,"_BBS_short-term_3_generation_trends.csv"),row.names = F)


# 
# 
# 
# 
# # SOCB Trends -------------------------------------------------------------
# 
# ## export all trends that include some area within Canada, as well as the US trends
# ## see note from M-A, Teams chat Jan 26, 2023
# # ok, here we go: we want the ability to see trends for Canada (of course), but also continental, US, BCR, prov/terr, and the intersections of BCR and prov/terr, because that will feed into the goals (ugh)
# 
#  trends_round <- readRDS("output/alltrends.rds")
# # indices_round <- readRDS("output/allindices.rds")
# # smooth_indices_round <- readRDS("output/allsmooth_indices.rds")
# 
# 
# # tlong <- trends_round %>% 
# #   filter(Trend_Time == "Long-term") %>% 
# #   group_by(species) %>% 
# #   summarise(earliest = min(Start_year),
# #             latest = max(Start_year))
# 
# trends_out <- trends_round %>% 
#    filter((For_web == TRUE | Region %in% c("Continental","US"))) %>% 
#   mutate(prob_LD = prob_decrease_50_percent,
#          prob_MD = prob_decrease_25_percent - prob_decrease_50_percent,
#          prob_LC = (prob_decrease_0_percent-prob_decrease_25_percent)+(prob_increase_0_percent-prob_increase_33_percent) ,
#          prob_MI = prob_increase_33_percent - prob_increase_100_percent,
#          prob_LI = prob_increase_100_percent)
# 
# test_probs <- trends_out %>% 
#   mutate(prob_test = prob_LD+prob_MD+prob_LC+prob_MI+prob_LI)
#   
# if(any(round(test_probs$prob_test,2) != 1)){stop("probabilites of change categories don't sum properly")}
# 
# 
# trends_out <- trends_out %>%
#   mutate(years = paste(Start_year,End_year,sep = "-"),
#          results_code = "BBS",
#          season = "breeding",
#          version = YYYY,
#          species_id = "",
#          area_code = ifelse(Region_type == "prov_state",Region,Region_alt),
#          area_code = gsub(area_code,pattern = "United States of America",
#                           replacement = "USA"),
#          model_type = "GAMYE") %>% 
#   rename(species_code = bbs_num,
#          species_name = species,
#          period = Trend_Time,
#          year_start = Start_year,
#          year_end = End_year,
#          trnd = Trend,
#          lower_ci = Trend_Q0.025,
#          upper_ci = Trend_Q0.975,
#          percent_change = Percent_Change,
#          percent_change_low = Percent_Change_Q0.025,
#          percent_change_high = Percent_Change_Q0.975,
#          prob_decrease_0 = prob_decrease_0_percent,
#          prob_decrease_25 = prob_decrease_25_percent,
#          prob_decrease_30 = prob_decrease_30_percent,
#          prob_decrease_50 = prob_decrease_50_percent,
#          prob_increase_0 = prob_increase_0_percent,
#          prob_increase_33 = prob_increase_33_percent,
#          prob_increase_100 = prob_increase_100_percent,
#          confidence = reliability,
#          precision_num = Width_of_95_percent_Credible_Interval,
#          precision_cat = precision,
#          coverage_num = reliab.cov,
#          coverage_cat = coverage,
#          sample_size = Mean_Number_of_Routes,
#          sample_total = Number_of_Routes,
#          prob_LD = prob_LD,
#          prob_MD = prob_MD,
#          prob_LC = prob_LC,
#          prob_MI = prob_MI,
#          prob_LI = prob_LI)
# 
# 
# 
# 
# trends_socb <- trends_out %>% 
#   relocate(results_code,
#            season,
#            version,
#            model_type,
#            area_code,
#          species_code,
#          species_id,
#          species_name,
#          period,
#          years,
#          year_start,
#          year_end,
#          trnd,
#          lower_ci,
#          upper_ci,
#          percent_change,
#          percent_change_low,
#          percent_change_high,
#          prob_decrease_0,
#          prob_decrease_25,
#          prob_decrease_30,
#          prob_decrease_50,
#          prob_increase_0,
#          prob_increase_33,
#          prob_increase_100,
#          confidence,
#          precision_num,
#          precision_cat,
#          coverage_num,
#          coverage_cat,
#          sample_size,
#          sample_total,
#          prob_LD,
#          prob_MD,
#          prob_LC,
#          prob_MI,
#          prob_LI)
# 
# write.csv(trends_socb,
#           paste0("website/BBS_",YYYY,"_trends_for_socb.csv"),
#           row.names = FALSE)
# 
# # write.csv(trends_socb[1:300,],
# #           paste0("website/sample_BBS_",YYYY,"_trends_for_socb.csv"),
# #           row.names = FALSE)
# # 
# 
# 
# 
# 
# # 
# # 
# # socb_headings <- c("results_code",
# #                    "version",
# #                    "area_code",
# #                    "species_code",
# #                    "species_id",
# #                    "season",
# #                    "period",
# #                    "years",
# #                    "year_start",
# #                    "year_end",
# #                    "trnd",
# #                    "index_type",
# #                    "lower_ci",
# #                    "upper_ci",
# #                    "stderr",
# #                    "model_type",
# #                    "model_fit",
# #                    "percent_change",
# #                    "percent_change_low",
# #                    "percent_change_high",
# #                    "prob_decrease_0",
# #                    "prob_decrease_25",
# #                    "prob_decrease_30",
# #                    "prob_decrease_50",
# #                    "prob_increase_0",
# #                    "prob_increase_33",
# #                    "prob_increase_100",
# #                    "reliability",
# #                    "precision_num",
# #                    "precision_cat",
# #                    "coverage_num",
# #                    "coverage_cat",
# #                    "goal",
# #                    "goal_lower",
# #                    "sample_size",
# #                    "sample_total",
# #                    "subtitle",
# #                    "prob_LD",
# #                    "prob_MD",
# #                    "prob_LC",
# #                    "prob_MI",
# #                    "prob_LI")
# # 
# # trend_headings_match <- c("",
# #                           "",
# #                           "Region_alt",
# #                           "bbs_num",
# #                           "species",
# #                           "Trend_Time",
# #                           "",
# #                           "Start_year",
# #                           "End_year",
# #                           "Trend",
# #                           "",
# #                           "Trend_Q0.025",
# #                           "Trend_Q0.975",
# #                           "",
# #                           "",
# #                           "",
# #                           "Percent_Change",
# #                           "Percent_Change_Q0.025",
# #                           "Percent_Change_Q0.975",
# #                           "prob_decrease_0_percent",
# #                           "prob_decrease_25_percent",
# #                           "prob_decrease_30_percent",
# #                           "prob_decrease_50_percent",
# #                           "prob_increase_0_percent",
# #                           "prob_increase_33_percent",
# #                           "prob_increase_100_percent",
# #                           "reliability",
# #                           "Width_of_95_percent_Credible_Interval",
# #                           "precision",
# #                           "reliab.cov",
# #                           "coverage",
# #                           "",
# #                           "",
# #                           "Mean_Number_of_Routes",
# #                           "Number_of_Routes",
# #                           "",
# #                           "prob_LD",
# #                           "prob_MD",
# #                           "prob_LC",
# #                           "prob_MI",
# #                           "prob_LI")
# # 
# # socb_headings_extract <- data.frame(socb = socb_headings,
# #                                     trend = trend_headings_match)
# # 
# # trends_select <- trend_headings_match[-which(trend_headings_match == "")]
# # 
# # socb_headings_select <- socb_headings_extract %>% 
# #   filter(trend != "")
# # 
# # trends_socb <- trends_out %>% 
# #   select(any_of(trends_select))
# # 
# # if(any(names(trends_socb) != socb_headings_select[,2])) stop("Stop columns don't match")
# # 
# # names(trends_socb) <- socb_headings_select[,1]
# # write.csv(trends_socb,paste0("website/BBS_",YYYY,"_trends_for_socb.csv"))
# # 
# # socb_headings_extract <- socb_headings_extract %>% 
# #   rename(NatureCountsTrendsSample = socb,
# #          BBS_trend_headers = trend)
# # 
# # write.csv(socb_headings_extract,
# #           "website/linking_columns_naturecounts_bbs.csv")
# # 
# 
# 
# 
# # SOCB indices ------------------------------------------------------------
# 
# indices_round <- readRDS("output/allindices.rds")
# 
# indices_socb <- indices_round %>% 
#   filter((For_web == TRUE | Region %in% c("Continental","US"))) %>% 
#   group_by(species,Region,Region_type,Trend_Time) %>% 
#   mutate(LOESS_index = loess_func(Index,Year),
#          area_code = ifelse(Region_type == "prov_state",Region,Region_alt),
#          area_code = gsub(area_code,pattern = "United States of America",
#                           replacement = "USA")) %>% 
#   ungroup() %>% 
#   rename(species_code = bbs_num,
#     species_id = species,
#     index = Index,
#     year = Year,
#          period = Trend_Time,
#          upper_ci = Index_q_0.95,
#          lower_ci = Index_q_0.05) %>% 
#   select(-c(Index_q_0.025,
#             Index_q_0.975,
#             Region_type,
#             Region)) %>% 
#   relocate(area_code,
#            year,
#            period, 
#            species_code,
#            species_id,
#            index,
#            upper_ci,
#            lower_ci,
#            LOESS_index)
# 
# write.csv(indices_socb,
#           file = paste0("website/BBS_",YYYY,"_annual_indices_for_socb.csv"),
#           row.names = FALSE)
# 
# # write.csv(indices_socb[sample(1:nrow(indices_socb),100,FALSE),],"sample_indices_output_bbs.csv",row.names = FALSE)

