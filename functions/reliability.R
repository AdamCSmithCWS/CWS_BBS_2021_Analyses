reliab_func_prec <- function(x){
  y = rep("Low",length(x))
  y[which(x <= prec_cuts["Medium"])] <- "Medium"
  y[which(x < prec_cuts["High"])] <- "High"
  y = factor(y,levels = c("Low","Medium","High"),ordered = T)
  
  return(y)
}

reliab_func_cov <- function(x){
  y = rep(NA,length(x))
  y[which(x < cov_cuts["Medium"])] <- "Low"
  y[which(x >= cov_cuts["Medium"])] <- "Medium"
  y[which(x > cov_cuts["High"])] <- "High"
  y = factor(y,levels = c("Low","Medium","High"),ordered = T)
  return(y)
}

reliab_func_pool <- function(x){
  y = rep(NA,length(x))
  y[which(x < pool_cuts["Medium"])] <- "Low"
  y[which(x >= pool_cuts["Medium"])] <- "Medium"
  y[which(x > pool_cuts["High"])] <- "High"
  y = factor(y,levels = c("Low","Medium","High"),ordered = T)
  
  return(y)
}

reliab_func_backcast <- function(x){
  y = rep(NA,length(x))
  y[which(x < backcast_cuts["Medium"])] <- "Low"
  y[which(x >= backcast_cuts["Medium"])] <- "Medium"
  y[which(x > backcast_cuts["High"])] <- "High"
  y = factor(y,levels = c("Low","Medium","High"),ordered = T)

  return(y)
}

# reliab_func_backcast <- function(x){
#   y = rep("Low",length(x))
#   y[which(x <= backcast_cuts["Medium"])] <- "Medium"
#   y[which(x < backcast_cuts["High"])] <- "High"
#   y = factor(y,levels = c("Low","Medium","High"),ordered = T)
#   
#   return(y)
# }

