


# Cross-validation folds definition ---------------------------------------

cv_folds <- function(orig_df = NULL, # original full observational dataframe of BBS data
                     K = 10,  # number of folds
                     fold_groups = "Observer_Factored"){  # set to FALSE to have folds move forwards in time
  
  if(is.null(orig_df)){
    stop("orig_df must be supplied")
    return(NULL)
  }
  # ID training and testing -------------------------------------------------
  # each fold must include all strata and routes
  # easy to do if first-year observations are always kept in the training fold
  # possibly problematic because no way to compare the predictions for first years
  # for each stratum, split observers in to K groups
  
  orig_df <- orig_df %>% 
    rename_with(.,~ gsub(pattern = fold_groups,replacement = "grouping",.x)) 
  
  n_groups <- length(unique(orig_df$grouping))
  
  #makes a dataframe that assigns each group to one of the K folds
  #after sorting by Stratum, so that each stratum is represented in each fold
  grps <- orig_df %>% 
    distinct(grouping,Stratum_Factored) %>% 
    arrange(Stratum_Factored) %>% 
    distinct(grouping) %>% 
    mutate(fold = rep(1:K,length.out = n_groups))
  
  orig_df <- orig_df %>% 
    left_join(.,grps,by = "grouping") %>% 
    group_by(grouping) %>% 
    mutate(fold2 = fold+First_Year,
           fold = ifelse(fold2 > K,1,fold2)) %>% 
    rename_with(.,~ gsub(replacement = fold_groups,pattern = "grouping",.x))
  
  
  
  out_df <- orig_df 
  
  return(out_df)
  
}