animated_trend_map <- function(fit = stanfit,
                               rawdata = stan_data,
                               stratification = strat_sel,
                               alt_n = NA,
                               firstYear = NULL,
                               lastYear = NULL,
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
  
  starts <- c(firstYear:(lastYear-1))
  
  pl_Indstt <- plot_indices(inds1)
  pl_Inds <- pl_Indstt$Continental
  imgs <- vector(mode = "character",length = length(starts))
  if(!is.na(alt_n)){
    tmpi <- inds2$data_summary %>% 
      filter(Region == "Continental")
    
    pl_Inds <- pl_Inds+
      geom_line(data = tmpi,
                aes(x = Year,y = Index))+
      geom_ribbon(data = tmpi,
                  aes(x = Year,y = Index,ymin = Index_q_0.025,ymax = Index_q_0.975),
                  alpha = 0.2)+
      ylab("Relative Abundance")
  }
  pl_Inds <- pl_Inds+
    labs(title = paste("Survey-wide trajectory",species))+
    xlab("")
  
  #print(pl_Inds)
    
  
  
  
  # ggplot(data = Indices_all,aes(x = true_year,y = median))+
  #   geom_ribbon(aes(ymin = lci,ymax = uci,fill = version),alpha = 0.1)+
  #   geom_line(aes(colour = version))+
  #   geom_point(data = isel,size = 2,colour = "black")+
  #   theme_bw()+
  #   xlab("")+
  #   ylab("Mean annual prediction")+
  #   theme(legend.position = "none")

  # dir_out_temp <- paste0(dir_out,"temp/")
  # if(!dir.exists(dir_out_temp)){
  #   dir.create(dir_out_temp)
  # }
  dir_out_temp <- tempdir()
  ann_trends_maps <- NULL
  
for(j in 1:length(starts)){
  dd <- starts[j]
  trends_anntemp <- generate_trends(inds2,Min_year = dd,Max_year = dd+1)
  ann_trends_map_data <- generate_map_data(trends_anntemp,select = TRUE,stratify_by = strat_sel,
                                        species = paste0(species," ",model_sel))
  mptmp <- ann_trends_map_data$map %>% 
    mutate(Start_year = dd,
           End_year = dd+1)
  ann_trends_maps <- bind_rows(ann_trends_maps,mptmp)
  
}
  map_palette <- ann_trends_map_data$map_palette
  # map_palette <- c("#FFFFFF",map_palette)
  # names(map_palette)[1] <- NA
 
   
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
                                 name = paste0("Annual Trend"),na.value = "white")+
    ggplot2::guides(alpha = "none")
  
 
  mp.plot2 <- mp.plot+
    transition_states(End_year,
                      transition_length = 3,
                      state_length = 1)+ 
    ease_aes()+ 
    ggtitle('One-year trend {closest_state}')
  
  #mp.plot2

  animate(mp.plot2,
          duration = 60,
          res = 72*res_mag,
          width = 480*res_mag,
          height = 380*res_mag)
  
  anim_save(filename = paste0(dir_out,"animated_map2_",file_name_prefix,species_f,".gif"))
  
  
#   
#   # angif <- ggplot(data = ttmd)+
# #   geom_sf(data = prov_state,alpha = 0,
# #           colour = grey(0.8),inherit.aes = FALSE)+
# #   geom_sf(aes(fill = trend_plot))+
# #   coord_sf(xlim = xb,
# #            ylim = yb)+
# #   theme_bw()+
# #   xlab("")+
# #   ylab("")+
# #   labs(title = paste("Annual Trends in space",yy,"-",yy2))+
# #   scale_colour_manual(values = map_palette_s, 
# #                       aesthetics = c("fill"),
# #                       guide = guide_legend(reverse=TRUE),
# #                       name = paste0("Trend"))
# 
# isel <- inds2$data_summary %>% 
#   filter(Year == dd+1,
#          Region == "Continental")
# 
# 
# pl_Inds2 <- pl_Inds +
#   geom_point(data = isel,aes(x = Year, y = Index),
#              size = 3)
#   
# #print(pl_Inds2)
# 
# 
# pout <- pl_Inds2 + ann_trends_maps +
#   plot_layout(ncol = 1,
#               nrow = 2,
#               heights = c(2,4))
# 
# 
# imgs[j] <- paste0(dir_out_temp,
#                   species_f,j,
#                   ".png")
# 
# png(filename = paste0(dir_out_temp,
#                       species_f,j,
#                       ".png"),
#     res = 72*res_mag,
#     width = 480*res_mag,
#     height = 550*res_mag)
# 
# print(pout)
# dev.off()
# 
# 
# }
# 
# 
# ## list file names and read in
# img_list <- lapply(imgs, image_read)
# 
# ## join the images together
# img_joined <- image_join(img_list)
# 
# ## animate at 2 frames per second
# img_animated <- image_animate(img_joined, fps = 2)
# 
# ## view animated image
# #img_animated
# 
# ## save to disk
# image_write(image = img_animated,
#             path = paste0(dir_out,"animated_map_",file_name_prefix,species_f,".gif"))


}
