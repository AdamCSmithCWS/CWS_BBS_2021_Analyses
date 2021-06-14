tr_func <- function(samples = nsmooth_samples,
                    scale = "strat_name",
                    start_year = 2000,
                    end_year = 2019,
                    slope = FALSE,
                    nsig = 3){
  
  
  nyear = end_year - start_year
  
  #samples$reg <- as.character(samples[,scale])
  #end-point %/year
  p_py <- function(i,ny = nyear){
    ppy <- 100*((i[2]/i[1])^(1/(ny))-1)
    return(ppy)
  }
  
  # end-point %-change function
  p_ch <- function(i){
    pch <- 100*((i[2]/i[1])-1)
    return(pch)
  }
  
  # slope %-change function
  p_ch_sl <- function(tt,ny = nyear){
    (((((tt)/100)+1)^ny)-1)*100
  }
  
  # slope trend %/year
  p_py_sl = function(i){
    il = log(i)
    n = length(il)
    sy = sum(il)
    sx = sum(1:n)
    ssx = sum((1:n)^2)
    sxy = sum(il*(1:n))
    b = (n*sxy - sx*sy)/(n*ssx - sx^2)
    py_sl = 100*(exp(b)-1)
    return(py_sl)
  }

  
  
  
  tr_sm <- samples %>% filter(year %in% c(start_year,end_year))%>% 
    arrange(year) %>% 
    rename_with(.,~gsub(pattern = scale,"reg",.x,fixed = TRUE))  %>%
    group_by(.draw,.variable,reg) %>% 
    summarise(ch = p_ch(.value),
              t = p_py(.value),ny = nyear)%>% 
    group_by(reg) %>% 
    summarise(mean_trend = signif(mean(t),nsig),
              median_trend = signif(median(t),nsig),
              lci_trend = signif(quantile(t,0.025),nsig),
              uci_trend = signif(quantile(t,0.975),nsig),
              mean_ch = signif(mean(ch),nsig),
              median_ch = signif(median(ch),nsig),
              lci_ch = signif(quantile(ch,0.025),nsig),
              uci_ch = signif(quantile(ch,0.975),nsig)) %>% 
    mutate(start_year = start_year,
           end_year = end_year) %>% 
    rename_with(.,~gsub(pattern = "reg",scale,.x,fixed = TRUE)) 
  
  
  if(slope){

  tr_slope <- samples %>% filter(year %in% c(start_year:end_year))%>% 
    arrange(year) %>% 
    rename_with(.,~gsub(pattern = scale,"reg",.x,fixed = TRUE))  %>%
    group_by(.draw,.variable,reg) %>% 
    summarise(sl_t = p_py_sl(.value))%>% 
    group_by(reg) %>% 
  summarise(mean_slope_trend = signif(mean(sl_t),nsig),
            median_slope_trend = signif(median(sl_t),nsig),
            lci_slope_trend = signif(quantile(sl_t,0.025),nsig),
            uci_slope_trend = signif(quantile(sl_t,0.975),nsig),
            mean_ch_slope_trend = signif(p_ch_sl(mean_slope_trend,nyear),nsig),
            median_ch_slope_trend = signif(p_ch_sl(median_slope_trend,nyear),nsig),
            lci_ch_slope_trend = signif(p_ch_sl(lci_slope_trend,nyear),nsig),
            uci_ch_slope_trend = signif(p_ch_sl(uci_slope_trend,nyear),nsig))  
  

  tr_sm <- tr_sm %>% rename_with(.,~gsub(pattern = scale,"reg",.x,fixed = TRUE)) %>% 
    left_join(.,tr_slope,by = "reg") %>% 
    rename_with(.,~gsub(pattern = "reg",scale,.x,fixed = TRUE))
  }
  
  tr_sm <- tr_sm %>% 
  return(tr_sm)
  
  
}



# 
# trend_map <- function(base_map = realized_strata_map,
#                       trends = trend_90){
#   breaks <- c(-7, -4, -2, -1, -0.5, 0.5, 1, 2, 4, 7)
#   labls = c(paste0("< ",breaks[1]),paste0(breaks[-c(length(breaks))],":", breaks[-c(1)]),paste0("> ",breaks[length(breaks)]))
#   labls = paste0(labls, " %")
#   map_palette <- c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf",
#                    "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695")
#   names(map_palette) <- labls
#   
#   tmap <- left_join(base_map,trends,by = c("ST_12" = "strat_name")) %>% 
#     mutate(Tplot = cut(mean_trend,breaks = c(-Inf, breaks, Inf),labels = labls))
#   
#   fyear = unique(trends$start_year)
#   lyear = unique(trends$end_year)
#   tplot <- ggplot(data = tmap)+
#     geom_sf(aes(fill = Tplot)) +
#     labs(title = paste(species,"trends spatial GAMYE"))+
#     scale_colour_manual(values = map_palette, 
#                         aesthetics = c("fill"),
#                         guide = ggplot2::guide_legend(reverse=TRUE),
#                         name = paste0("Trend\n",fyear,"-",lyear))
#   
#   
#   return(tplot)
#   
# }



