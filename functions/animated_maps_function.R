animated_trend_map <- function(fit = stanfit,
                               rawdata = stan_data,
                               stratification = strat_sel,
                               alt_n = NA,
                               firstYear = NULL,
                               lastYear = NULL,
                               trend_length = 10,
                               res_mag = 3,
                               dir_out = "",
                               file_name_prefix = "_",
                               species = species){
  
require(gifski)
  require(gganimate)
require(magick)
require(patchwork)

  species_f <- gsub("[[:punct:]]","",gsub(" ","",species))

if(!is.na(alt_n)){
  
  inds2 <- generate_indices(jags_mod = fit,
                            jags_data = rawdata,
                            backend = "Stan",
                            stratify_by = stratification,
                            alternate_n = alt_n,
                            max_backcast = 1)
  
  inds1 <- generate_indices(jags_mod = fit,
                            jags_data = rawdata,
                            backend = "Stan",
                            stratify_by = stratification,
                            max_backcast = 1)
  
}else{
  inds2 <- generate_indices(jags_mod = fit,
                            jags_data = rawdata,
                            backend = "Stan",
                            stratify_by = stratification,
                            max_backcast = 1)
  
  inds1 <- inds2
}

  if(is.null(firstYear)){ firstYear <- min(inds2$data_summary$Year)}
  if(is.null(lastYear)){ lastYear <- max(inds2$data_summary$Year)}
  
  starts <- c(firstYear:(lastYear-trend_length))
  
  
  ann_trends_maps <- NULL
  
for(j in 1:length(starts)){
  dd <- starts[j]
  trends_anntemp <- generate_trends(inds2,Min_year = dd,Max_year = dd+trend_length)
  ann_trends_map_data <- generate_map_data(trends_anntemp,select = TRUE,stratify_by = strat_sel,
                                        species = paste0(species," ",model_sel))
  mptmp <- ann_trends_map_data$map %>% 
    mutate(Start_year = dd,
           End_year = dd+trend_length,
           backcast_flag = ifelse(backcast_flag > 0.5,backcast_flag,0))
  ann_trends_maps <- bind_rows(ann_trends_maps,mptmp)
  
}
  map_palette <- ann_trends_map_data$map_palette
 
   
  mp.plot <- ggplot2::ggplot()+
    ggplot2::geom_sf(data = ann_trends_maps,ggplot2::aes(fill = Tplot,alpha = backcast_flag),colour = grey(0.4),size = 0.1)+
    ggplot2::theme_minimal()+
    ggplot2::ylab("")+
    ggplot2::xlab("")+
    ggplot2::labs(subtitle = species)+
    ggplot2::theme(legend.position = "right", line = ggplot2::element_line(size = 0.4),
                                                    rect = ggplot2::element_rect(size = 0.1),
                                                    axis.text = ggplot2::element_blank(),
                                                    axis.line = ggplot2::element_blank())+
    ggplot2::scale_colour_manual(values = map_palette, aesthetics = c("fill"),
                                                  guide = ggplot2::guide_legend(reverse=TRUE),
                                 name = paste0(trend_length,"-year Trend"),na.value = "white")+
    ggplot2::scale_alpha(range = c(0,1))+
    ggplot2::guides(alpha = "none")
  
 
  mp.plot2 <- mp.plot+
    transition_states(End_year,
                      transition_length = 1,
                      state_length = 1)+ 
    ease_aes()+ 
    ggtitle('Trends ending in {closest_state}')
  
  #mp.plot2

  animate(mp.plot2,
          duration = 60,
          res = 72*res_mag,
          width = 480*res_mag,
          height = 380*res_mag)
  
  anim_save(filename = paste0(dir_out,"animated_map2_",file_name_prefix,species_f,".gif"))
  
  


}
