for_web_func <- function(x){
  y = rep(FALSE,length(x))
  for(j in 1:length(x)){
    sts = unlist(strsplit(x[j],split = " ; "))
    if(any(stringr::str_starts(sts,pattern = "CA-"))){
      y[j] <- TRUE
    }
  }
  return(y)
}




reliability_func <- function(x,y,z){
  o = rep("High",length(x))
  for(j in 1:length(x)){
    sts = c(x[j],y[j],z[j])
    if(any(grepl(sts,pattern = "edium"))){
      o[j] <- "Medium"
    }
    if(any(grepl(sts,pattern = "ow"))){
      o[j] <- "Low"
    }
  }
  
  o <- factor(o,levels = c("High","Medium","Low"),ordered = TRUE)
  return(o)
}



