
dim_ext <- function(dim = 1,
                     var = "",
                     cl = "Parameter",
                     dat = NULL){
  ##3 function to extract the indicator values from cmdstanr output
  require(stringr)
  
  pat = paste0("(?<=",var,"\\[")
  
  if(dim > 1){
    for(j in 1:(dim-1)){
      
      pat2 = paste0(pat,")[:digit:]+")
      cl2 = str_extract(unlist(dat[,cl]),pattern = pat2)
      
      d = max(nchar(cl2))
      
      pat = paste0(pat,"[:digit:]{1,",d,"}[:punct:]")
    }
  }
  
  
  pat = paste0(pat,")[:digit:]+")
  dds = as.integer(str_extract(unlist(dat[,cl]),pattern = pat))
  return(dds)
  
}


### function to generate the same tidy output as gather-draws in tidybayes package
## dims should be a character vector defining the dimensions of the parameter
## e.g., parm = "nsmooth", dims = c("stratum","year"),
gather_samples <- function(fit = cmdstanfit,
                        parm = "nsmooth",
                        dims = NULL){
  require(posterior)
  samples <- as_draws_df(fit$draws(variables = c(parm)))
  if(length(dims) > 0){
    parm_ex <- paste0(parm,"\\[")
  }else{
    parm_ex <- parm
  }
  
  plong <- suppressWarnings(samples %>% pivot_longer(
    cols = matches(parm_ex,ignore.case = FALSE),
    names_to = c(".variable"),
    values_to = ".value",
    values_drop_na = TRUE
  )) 
  
  for(dn in 1:length(dims)){
    dd = dims[dn]
    plong[,dd] = dim_ext(dim = dn,
                          var = parm,
                          cl = ".variable",
                          dat = plong)
    
  }
  
  plong <- plong %>% mutate(.variable = parm)
  return(plong)
  
}

