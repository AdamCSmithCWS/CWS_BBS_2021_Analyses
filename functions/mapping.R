
generate_web_maps <- function(df,
                              canmap = canmap){
  
  dfstrata = df[which(df$Region_type == "stratum"),]
  for(j in 1:nrow(df)){
    if(df[j,"For_web"]){
      
      mapname <- df[j,"mapfile"]
      
      sts = unlist(strsplit(df[j,"Strata_included"],split = " ; "))
     
      mapo = map_f(st = dfstrata,
                   stinc = sts,
                   map = canmap)
      
      png(filename = paste0("website/WebMaps/",mapname),
          bg = "white",width = 480, height = 320)
      print(mapo)
      dev.off()
      
      
      
    }
    
  }
  
  
  
  return(invisible(NULL))
  
}

map_f <- function(st, #dataframe of trends
                  map = canmap,
                  slope = F,
                  stinc = sts){
  
  
  stplot = st[which(st$Region %in% c(stinc)),]

  breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
  labls = c(paste0("< ",breaks[1]),
            paste0(breaks[-c(length(breaks))],
                   ":",
                   breaks[-c(1)]),
            paste0("> ",breaks[length(breaks)]))
  labls = paste0(labls, " %")
  #labls[length(labls)+1] <- "X"
  map@data$row_num <- 1:nrow(map@data)
  map@data <- merge(map@data, stplot, by.x = "ST_12", by.y = "Region", all.x = T,sort = F)
  map@data <- map@data[order(map@data$row_num), ]
  if(slope){
    map@data$Trend <- as.numeric(as.character(map@data$Slope_Trend))
  }else{
    map@data$Trend <- as.numeric(as.character(map@data$Trend))
  }
  
  map@data$Trend <- cut(map@data$Trend, breaks = c(-Inf, breaks,Inf),
                        labels = labls,
                        ordered_result = T)
  
  map@data <- subset(map@data, select = c(Trend))
  
  map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                   "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
  names(map_palette) <- labls
  
  return(
    sp::spplot(map, 
               col.regions = map_palette,
               edge.col = grey(0.5),
               sp.layout = list(basmap, edge.col = grey(0.5),
                                fill="transparent", first=FALSE),
               par.settings = list(axis.line = list(col = 'transparent')))
  )
}
