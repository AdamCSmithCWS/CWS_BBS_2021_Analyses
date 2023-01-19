


# Cross-validation folds definition ---------------------------------------

cv_folds <- function(orig_df = NULL, # original full observational dataframe of BBS data
                     K = 10,  # number of folds
                     fold_groups = "obs_n", #primary grouping factor alternative is Route_Factored
                     na_singles = TRUE,
                     strata_group = "strata_name",
                     first_year = "first_year"){  
         
  require(tidyverse)             
if(is.null(orig_df)){
  stop("orig_df must be supplied - dataframe with all observations and with column matching fold_groups")
  return(NULL)
}
  

# fold_groups define ------------------------------------------------------
# fold_groups is the critical grouping factor for the leave future out cross-validation
# each fold-k identifies a test-set of observations that include all future observations from
# a randomly selected set of groups in fold_group (e.g., all future observations for a set of observers)
# balanced across all routes and strata
  if(!fold_groups %in% names(orig_df)){
    stop(paste("column", fold_groups,
               "not found in dataframe supplied to orig_df"))
  }
# ID training and testing -------------------------------------------------
# each training-fold (inverse of those in fold-k) must include all values of fold_groups
# Challenging for fold_groups with only 1-year of data - they can't be used in any of the 
## testing sets
  
  # split observers in to K groups, based on fold_groups
  out_df <- orig_df %>% 
    rename_with(.,~ gsub(pattern = fold_groups,replacement = "grouping",.x)) %>% 
    rename_with(.,~ gsub(pattern = strata_group,replacement = "Stratum_Factored",.x)) %>% 
    rename_with(.,~ gsub(pattern = first_year,replacement = "first_year",.x)) 
  
  n_groups <- length(unique(out_df$grouping))
  
  #makes a dataframe that assigns each group to one of the K folds
  #after sorting by Stratum, so that each stratum is represented in each fold
grps <- out_df %>% 
  distinct(grouping,Stratum_Factored) %>% 
  arrange(Stratum_Factored) %>% 
  distinct(grouping) %>% 
  mutate(fold = rep(1:K,length.out = n_groups))
  

# optional exclude all groups with no replication from testing sets -------
if(na_singles){
  n_events_by_group <- out_df %>% 
    group_by(grouping,first_year) %>% 
    summarise(n_events = n(),.groups = "drop") %>% 
    group_by(grouping) %>% 
    summarise(n_cats = n(),
              only_first_year = ifelse(n_cats == 1,TRUE,FALSE),
              n_events_cv = sum(n_events))
  
  out_df <- out_df %>% 
    left_join(.,n_events_by_group,by = "grouping") %>% 
    left_join(.,grps,by = "grouping") %>% 
    group_by(grouping) %>% 
    mutate(fold2 = fold-first_year, #bumps first-year observations to the previous fold
           fold = ifelse(fold2 < 1,K,fold2), #wraps around the list of K values for fold == 1 and first_year == 1
           fold = ifelse(only_first_year,NA,fold)) %>% #removes groups with only first-year observer route combinations
    rename_with(.,~ gsub(replacement = fold_groups,pattern = "grouping",.x)) %>% 
    rename_with(.,~ gsub(replacement = strata_group,pattern = "Stratum_Factored",.x)) %>% 
    rename_with(.,~ gsub(replacement = first_year,pattern = "first_year",.x)) %>% 
  select(-fold2)
  
}else{
  
out_df <- out_df %>% 
  left_join(.,grps,by = "grouping") %>% 
  group_by(grouping) %>% 
  mutate(fold2 = fold+first_year, #bumps first-year observations to the next fold
         fold = ifelse(fold2 > K,1,fold2)) %>% 
  rename_with(.,~ gsub(replacement = fold_groups,pattern = "grouping",.x)) %>% 
  rename_with(.,~ gsub(replacement = strata_group,pattern = "Stratum_Factored",.x)) %>% 
  rename_with(.,~ gsub(replacement = first_year,pattern = "first_year",.x)) %>% 
  select(-fold2)
}

return(out_df)

}
