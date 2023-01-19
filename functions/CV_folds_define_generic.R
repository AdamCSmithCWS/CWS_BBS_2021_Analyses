


# Cross-validation folds definition ---------------------------------------

cv_folds <- function(orig_df = NULL, # original full observational dataframe of BBS data
                     K = 10,  # number of folds
                     fold_groups = "Observer_Factored",
                     na_singles = TRUE){  
                      
if(is.null(orig_df)){
  stop("orig_df must be supplied")
  return(NULL)
}
# ID training and testing -------------------------------------------------
# each training-fold (inverse of those in fold-k) must include all values of fold_groups
# Challenging for fold_groups with only 1-year of data - they can't be used in any of the 
## testing sets
  
  # split observers in to K groups, based on fold_groups
  out_df <- orig_df %>% 
    rename_with(.,~ gsub(pattern = fold_groups,replacement = "grouping",.x)) 
    
  n_groups <- length(unique(out_df$grouping))
  
  #makes a dataframe that assigns each group to one of the K folds
  #after sorting by Stratum, so that each stratum is represented in each fold
grps <- out_df %>% 
  distinct(grouping,Stratum_Factored) %>% 
  arrange(Stratum_Factored) %>% 
  distinct(grouping) %>% 
  mutate(fold = rep(1:K,length.out = n_groups))
  
if(na_singles){
  n_events_by_group <- out_df %>% 
    group_by(grouping,First_Year) %>% 
    summarise(n_events = n(),.groups = "drop") %>% 
    group_by(grouping) %>% 
    summarise(n_cats = n(),
              only_first_year = ifelse(n_cats == 1,TRUE,FALSE),
              n_events_cv = sum(n_events))
  
  out_df <- out_df %>% 
    left_join(.,n_events_by_group,by = "grouping") %>% 
    left_join(.,grps,by = "grouping") %>% 
    group_by(grouping) %>% 
    mutate(fold2 = fold-First_Year, #bumps first-year observations to the previous fold
           fold = ifelse(fold2 < 1,K,fold2), #wraps around the list of K values for fold == 1 and first_year == 1
           fold = ifelse(only_first_year,NA,fold)) %>% 
    rename_with(.,~ gsub(replacement = fold_groups,pattern = "grouping",.x)) %>% 
    select(-fold2)
  
}else{
  
out_df <- out_df %>% 
  left_join(.,grps,by = "grouping") %>% 
  group_by(grouping) %>% 
  mutate(fold2 = fold+First_Year, #bumps first-year observations to the next fold
         fold = ifelse(fold2 > K,1,fold2)) %>% 
  rename_with(.,~ gsub(replacement = fold_groups,pattern = "grouping",.x)) %>% 
  select(-fold2)
}

return(out_df)

}
