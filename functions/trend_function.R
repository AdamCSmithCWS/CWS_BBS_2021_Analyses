tr_func <- function(samples = nsmooth_samples,
                    scale = "strat_name",
                    start_year = 2000,
                    end_year = 2019){
  
  
  nyear = end_year - start_year
  
  #samples$reg <- as.character(samples[,scale])
  p_py <- function(i1,i2,ny = nyear){
    ppy <- 100*((i2/i1)^(1/(ny))-1)
    return(ppy)
  }
  
  p_ch <- function(i1,i2){
    pch <- 100*((i2/i1)-1)
    return(pch)
  }
  
  
  tr_sm <- samples %>% filter(year %in% c(start_year,end_year)) %>% 
    pivot_wider(.,id_cols = any_of(c(scale,"strat_name",
                                     ".variable",".chain",".iteration",".draw")),
                names_from = year,
                values_from = .value,
                names_prefix = "Y") %>% 
    rename_with(.,~gsub(pattern = paste0("Y",start_year),"Y1",.x,fixed = TRUE))%>% 
    rename_with(.,~gsub(pattern = paste0("Y",end_year),"Y2",.x,fixed = TRUE)) %>% 
    rename_with(.,~gsub(pattern = scale,"reg",.x,fixed = TRUE)) %>% 
    mutate(ch = p_ch(Y1,Y2),
           t = p_py(Y1,Y2)) %>% 
    group_by(reg) %>% 
    summarise(mean_trend = mean(t),
              lci_trend = quantile(t,0.025),
              uci_trend = quantile(t,0.975),
              mean_ch = mean(ch),
              lci_ch = quantile(ch,0.025),
              uci_ch = quantile(ch,0.975)) %>% 
    mutate(start_year = start_year,
           end_year = end_year) %>% 
    rename_with(.,~gsub(pattern = "reg",scale,.x,fixed = TRUE))
  
  return(tr_sm)
  
  
}




trend_map <- function(base_map = realized_strata_map,
                      trends = trend_90){
  breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
  labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
  labls = paste0(labls, " %")
  map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
                   "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
  names(map_palette) <- labls
  
  tmap <- left_join(base_map,trends,by = c("ST_12" = "strat_name")) %>% 
    mutate(Tplot = cut(mean_trend,breaks = c(-Inf, breaks, Inf),labels = labls))
  
  fyear = unique(trends$start_year)
  lyear = unique(trends$end_year)
  tplot <- ggplot(data = tmap)+
    geom_sf(aes(fill = Tplot)) +
    labs(title = paste(species,"trends spatial GAMYE"))+
    scale_colour_manual(values = map_palette, 
                        aesthetics = c("fill"),
                        guide = ggplot2::guide_legend(reverse=TRUE),
                        name = paste0("Trend\n",fyear,"-",lyear))
  
  
  return(tplot)
  
}



