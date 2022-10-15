generate_map_data <- function(trend = NULL,
         select = FALSE,
         stratify_by = NULL,
         slope = FALSE,
         species = "",
         col_viridis = FALSE)
{
  # Silly things to remove "visible binding" note in check
  Trend <- NULL
  rm(Trend)
  Tplot <- NULL
  rm(Tplot)
  
  if(select){
    trend = trend[which(trend$Region_type == "stratum"),]
  }
  if (is.null(stratify_by))
  {
    stop("No stratification specified."); return(NULL)
  }
  
  if(isFALSE(is.element(stratify_by, c("state", "bcr", "latlong", "bbs_cws", "bbs_usgs"))))
  {
    stop("Invalid stratification specified, choose one of state, bcr, latlong, bbs_cws, or bbs_usgs"); return(NULL)
  }
  
  if (is.null(trend))
  {
    stop("Argument trend is empty, or stratification of model does not match stratify_by argument"); return(NULL)
  }
  
  
  fyr = min(trend$Start_year)
  lyr = min(trend$End_year)
  
  # map <- sf::read_sf(dsn = system.file("maps",
  #                                      package = "bbsBayes"),
  #                    layer = maps[[stratify_by]],
  #                    quiet = TRUE)
  map <- load_map(stratify_by = stratify_by)
  
  
  breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
  labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
  labls = paste0(labls, " %")
  
  if(slope){
    trend$Tplot <- as.numeric(as.character(trend$Slope_Trend))
  }else{
    trend$Tplot <- as.numeric(as.character(trend$Trend))
  }
  trend$Tplot <- cut(trend$Tplot,breaks = c(-Inf, breaks, Inf),labels = labls)
  # tlev <- levels(trend$Tplot)
  # trend$Tplot <- factor(trend$Tplot,levels = c(tlev))
  trend$ST_12 = trend$Region
  map = dplyr::left_join(x = map,y = trend,by = "ST_12")
  
  if(species != ""){
    ptit = paste(species,"trends",fyr,"-",lyr)
  }else{
    ptit = ""
  }
  
  if (col_viridis)
  {
    map_palette <- c("#fde725", "#dce319", "#b8de29", "#95d840", "#73d055", "#55c667",
                     "#238a8d", "#2d708e", "#39568c", "#453781", "#481567")
  }else
  {
    map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                     "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
  }
  
  names(map_palette) <- labls
  
  # mp.plot = ggplot2::ggplot()+
  #   ggplot2::geom_sf(data = map,ggplot2::aes(fill = Tplot),colour = grey(0.4),size = 0.1)+
  #   ggplot2::theme_minimal()+
  #   ggplot2::ylab("")+
  #   ggplot2::xlab("")+
  #   ggplot2::labs(title = ptit)+
   mp_theme <- ggplot2::theme(legend.position = "right", line = ggplot2::element_line(size = 0.4),
                   rect = ggplot2::element_rect(size = 0.1),
                   axis.text = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank())
   mp_colour_scale <- ggplot2::scale_colour_manual(values = map_palette, aesthetics = c("fill"),
                                 guide = ggplot2::guide_legend(reverse=TRUE),
                                 name = paste0("Trend\n",fyr,"-",lyr))

  ret_list <- list(map_palette = map_palette,
                   map = map,
                   title = ptit,
                    mp_theme = mp_theme,
                    mp_colour_scale = mp_colour_scale
  )
  return(ret_list)
}