plot_indices <- function(indices_list = NULL,
         ci_width = 0.95,
         min_year = NULL,
         max_year = NULL,
         species = "",
         title_size = 20,
         axis_title_size = 18,
         axis_text_size = 16,
         line_width = 1,
         add_observed_means = FALSE,
         add_number_routes = FALSE,
         add_adjusted_means = FALSE,
         add_extras = FALSE){

  Year <- NULL; rm(Year)
  Index <- NULL; rm(Index)
  Stratum <- NULL; rm(Stratum)
  obs_mean <- NULL; rm(obs_mean)
  lci <- NULL; rm(lci)
  uci <- NULL; rm(uci)
  
  indices = indices_list$data_summary
  
  lq = (1-ci_width)/2
  uq = ci_width+lq
  lqc = paste0("Index_q_",lq)
  uqc = paste0("Index_q_",uq)
  
  
  indices$lci = indices[,lqc]
  indices$uci = indices[,uqc]
  
  cl = "#39568c"
  
  plot_list <- list()
  
  if (!is.null(min_year))
  {
    indices <- indices[which(indices$Year >= min_year), ]
  }
  
  if(!is.null(max_year))
  {
    indices <- indices[which(indices$Year <= max_year), ]
  }
  
  mny = min(indices$Year)
  mxy = max(indices$Year)
  yys = pretty(seq(mny, mxy))
  yys = c(yys[-length(yys)],mxy)
  
  
  
  plot_index <- 1
  for (i in unique(indices$Region_alt))
  {
    to_plot <- indices[which(indices$Region_alt == i), ]
    
    if(add_number_routes){
      
      if(max(to_plot$nrts) > 200){
        ncby_y = ceiling(to_plot$nrts/50)
        annot = c("each dot ~ 50 routes")
      }else{
        if(max(to_plot$nrts) > 100){
          ncby_y = ceiling(to_plot$nrts/10)
          annot = c("each dot ~ 10 routes")
        }else{
          ncby_y = to_plot$nrts
          annot = c("each dot = 1 route")
        }
      }
      
      names(ncby_y) <- to_plot$Year
      dattc = data.frame(Year = rep(as.integer(names(ncby_y)),times = ncby_y))
    }
    
    
    if(add_observed_means){
      
      annotobs = to_plot[4,c("obs_mean","Year")]
      
      p <- ggplot2::ggplot() +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(),
                       axis.line = element_line(colour = "black"),
                       plot.title = ggplot2::element_text(size = title_size),
                       axis.title = ggplot2::element_text(size = axis_title_size),
                       axis.text = ggplot2::element_text(size = axis_text_size)) +
        ggplot2::labs(title = paste(species, " ", i, sep = ""),
                      x = "Year",
                      y = "Annual index of abundance (mean count)") +
        ggplot2::geom_point(data = to_plot,ggplot2::aes(x = Year,y = obs_mean),colour = grDevices::grey(0.6),
                            alpha = 0.4,
                            shape = 21)+
        ggplot2::geom_line(data = to_plot, ggplot2::aes(x = Year, y = Index), colour = cl, size = line_width) +
        ggplot2::geom_ribbon(data = to_plot, ggplot2::aes(x = Year, ymin = lci, ymax = uci),fill = cl,alpha = 0.3)+
        ggplot2::scale_x_continuous(breaks = yys)+
        ggplot2::scale_y_continuous(limits = c(0,NA))+
        ggplot2::annotate(geom = "text",x = annotobs$Year,y = annotobs$obs_mean,label = "",colour = grDevices::grey(0.6))
      
      if(add_number_routes){
        
        p <- p + ggplot2::geom_dotplot(data = dattc,mapping = ggplot2::aes(x = Year),drop = TRUE,binaxis = "x", stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = FALSE,fill = grDevices::grey(0.6),colour = grDevices::grey(0.6),alpha = 0.2,dotsize = 0.3)+
          ggplot2::annotate(geom = "text",x = min(dattc$Year)+5,y = 0,label = annot,alpha = 0.4,colour = grDevices::grey(0.6))
        
      }
      
      if(add_adjusted_means){
       
        p <- p + 
          ggplot2::geom_point(data = to_plot,
                            ggplot2::aes(x = Year,y = obs_mean_scaled),
                            colour = grDevices::grey(0.6))

      }
      
      if(add_extras){
        upy <- max(to_plot$uci,na.rm = TRUE)/5
        to_plot <- to_plot %>% 
          mutate(mean_obs = (mean_obs*upy)+upy,
                 mean_site = (mean_site*upy)+upy)
        p <- p + 
          ggplot2::geom_point(data = to_plot,
                              ggplot2::aes(x = Year,y = mean_obs),alpha = 0.2,
                              colour = "darkred")+
          ggplot2::geom_point(data = to_plot,
                              ggplot2::aes(x = Year,y = mean_site),alpha = 0.3,
                              colour = "darkorange")+
          ggplot2::geom_hline(yintercept = upy,alpha = 0.1,colour = grey(0.5))+
          labs(subtitle = "Red = route-effects, Orange = observer effects")
          
          
      }
      
    }else{
      
      p <- ggplot2::ggplot() +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.background = ggplot2::element_blank(),
                       axis.line = element_line(colour = "black"),
                       plot.title = ggplot2::element_text(size = title_size),
                       axis.title = ggplot2::element_text(size = axis_title_size),
                       axis.text = ggplot2::element_text(size = axis_text_size)) +
        ggplot2::labs(title = paste(species, " ", i, sep = ""),
                      x = "Year",
                      y = "Annual index of abundance (mean count)") +
        ggplot2::geom_line(data = to_plot, ggplot2::aes(x = Year, y = Index), colour = cl, size = line_width) +
        ggplot2::geom_ribbon(data = to_plot, ggplot2::aes(x = Year, ymin = lci, ymax = uci),fill = cl, alpha = 0.3)+
        ggplot2::scale_x_continuous(breaks = yys)+
        ggplot2::scale_y_continuous(limits = c(0,NA))
      if(add_number_routes){
        
        p <- p + ggplot2::geom_dotplot(data = dattc,mapping = ggplot2::aes(x = Year),drop = TRUE,binaxis = "x", stackdir = "up",method = "histodot",binwidth = 1,width = 0.2,inherit.aes = FALSE,fill = grDevices::grey(0.6),colour = grDevices::grey(0.6),alpha = 0.2,dotsize = 0.3)+
          ggplot2::annotate(geom = "text",x = min(dattc$Year)+5,y = 0,label = annot,alpha = 0.4,colour = grDevices::grey(0.6))
      }
      
    }
    plot_list[[stringr::str_replace_all(paste(i),
                                        "[[:punct:]\\s]+",
                                        "_")]] <- p
    plot_index <- plot_index + 1
  }
  
  return(plot_list)
}