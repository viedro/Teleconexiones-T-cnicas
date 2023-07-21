setwd("C:/Users/ASUS/Desktop/TF_Tecnicas/Tf_tecnicas/teleconexiones")

library(rasterVis)
library(raster)
library(rgl)
library(rgdal)
library(elevatr)

library(sf)



shp <- shapefile("data/raw/Sudamérica.shp")
shp <- shp[shp$PAÍS == "Perú", ]

elevation <- get_elev_raster(shp, z = 6)

plot(elevation)


crs(elevation) <- crs(shp)
#############
recor <- raster::crop(elevation,shp[shp$PAÍS == "Perú", ]) %>%
  raster::mask(shp[shp$PAÍS == "Perú", ])


plot(recor)

myColorkey <- list(tri.lower = F, tri.upper = TRUE, ## where the colors change
                   labels=list(labels=c(seq(0,6000,500)), ##what to print
                               at=seq(0,6000,500)))



library(ggplot2)
library(metR)
peru2 <- as.data.frame(recor,xy=T)


peru <- st_read("data/raw/Sudamérica.shp")
peru <- peru[peru$PAÍS == "Perú",]
colnames(peru2)[3] <- "COR"


ggplot(data=peru)+
  #geom_tile(data=grid,aes(x,y,fill=corr))+
  geom_contour_fill(data=peru2,aes(x,y,z=COR),breaks=seq(0,6000,500))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),legend.key.height = unit(4.5,"cm"))+
  scale_fill_gradientn(colours=(paletteer_c("grDevices::terrain.colors",length(seq(0,6000,500))-1)),breaks=seq(0,6000,500),limits=c(0,6000))+
  geom_sf(fill="transparent",col="black",lwd=0.6)+
  ggtitle(paste0("Mapa Topográfico del Perú "))+
  xlab("Longitud")+
  ylab("Latitud")+
  guides(fill=guide_colorsteps(title="m"))+
  coord_sf(xlim=c(-81,-69),ylim=c(-17.55,-0.8))

ggsave(paste0("figures/","map.png"), width = 8.5, height = 10,units = 'in',dpi=500)





