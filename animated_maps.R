
library(gifski)
library(magick)
library(patchwork)

angif <- ggplot(data = ttmd)+
  geom_sf(data = prov_state,alpha = 0,
          colour = grey(0.8),inherit.aes = FALSE)+
  geom_sf(aes(fill = trend_plot))+
  coord_sf(xlim = xb,
           ylim = yb)+
  theme_bw()+
  xlab("")+
  ylab("")+
  labs(title = paste("Annual Trends in space",yy,"-",yy2))+
  scale_colour_manual(values = map_palette_s, 
                      aesthetics = c("fill"),
                      guide = guide_legend(reverse=TRUE),
                      name = paste0("Trend"))

isel <- Indices_all %>% 
  filter(true_year == yy2,
         version == "smooth")

pl_Inds <- ggplot(data = Indices_all,aes(x = true_year,y = median))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = version),alpha = 0.1)+
  geom_line(aes(colour = version))+
  geom_point(data = isel,size = 2,colour = "black")+
  theme_bw()+
  xlab("")+
  ylab("Mean annual prediction")+
  theme(legend.position = "none")+
  labs(title = paste("Survey-wide trajectory",species,dd))

pout <-  pl_Inds + angif +
  plot_layout(ncol = 1,
              nrow = 2,
              heights = c(1,4))
mg = 4
imgs[j] <- paste0(dir_out,
                  species_f,j,
                  ".png")

png(filename = paste0(dir_out,
                      species_f,j,
                      ".png"),
    res = 72*mg,
    width = 480*mg,
    height = 550*mg)

print(pout)
dev.off()


}



## list file names and read in
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 2)

## view animated image
#img_animated

## save to disk
image_write(image = img_animated,
            path = paste0("Figures/animated_map_index_",species_f,".gif"))

