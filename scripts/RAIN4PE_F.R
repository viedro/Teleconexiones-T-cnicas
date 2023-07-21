setwd("C:/Users/ASUS/Desktop/TF_Tecnicas/Tf_tecnicas/teleconexiones")

library(raster )

#Donde se tiene los datos RAIN4PE
nc <- brick("data/raw/RAIN4PE_daily_0.1d_1981_2015_v1.0.nc")

dates <- getZ(nc) 
range(dates)
m_y <-  format(as.Date(dates), "%Y-%m")
m_y

s <- stackApply(nc, m_y, fun=sum)
s

library(rasterVis)

levelplot(s[[1]],margin=F)

library(SPEI)
###Calculando el índice SPI usando la función calc
SPI <- calc(s, fun= function(x, scale = 12, distribution = "Gamma", na.rm=T,...) as.numeric((spi(x,scale= scale,na.rm=na.rm,...))$fitted))

shp <- shapefile("data/raw/Sudamérica.shp")

###Filtrando Peru
shp <- shp[shp$PAÍS == "Perú" ,]
shp
###Estableciendo sistema de refererencia al rasterbrick de los indices
crs(SPI) <- crs(shp)
SPI


levelplot(SPI[[12]],margin=F)



library(tidyverse)
###Recortando con el shapefile de Perú
recor <- raster::crop(SPI,shp) %>%
  raster::mask(shp)
levelplot(recor[[12]],margin=F)


SPI_arr <- raster::as.array(recor)

str(SPI_arr)



dates <- as.Date(paste0(unique(m_y),"-01"))
range(dates)


###Extrayendo las coordenadas de latitud y longitud del raster anterior
lon <- raster::xFromCol(recor, 1:ncol(recor))
lat <- raster::yFromRow(recor, 1:nrow(recor))

###Vemos las longitudes de cada eje
c(length(lon),length(lat))


SPI_arr <- aperm(SPI_arr, c(2,1,3))
###Corregido
str(SPI_arr)

lat <- rev(lat)  

#Inversión también de latitudes en el arreglo
SPI_arr <- SPI_arr[,ncol(SPI_arr):1 ,]

head(lat)


library(rsoi)



SPI_arr <- SPI_arr[,,seq(12,12*35,12)]

####################################Correlacion




functions <- list("download_oni","download_soi","download_aao","download_ao","download_pdo","download_dmi","download_mei","download_npgo")
name_index <- toupper(substr(functions,10,nchar(functions)))

graph <- list()
for (k in 1:length(functions)){
  result <- do.call(functions[[k]], list(use_cache = FALSE))
  result <- result[,c("Date","Year",name_index[k])]
  result <- result[result$Date>=as.Date("1981-01-01") & result$Date<=as.Date("2015-12-31"),]
  
  index <- result[,c("Year",name_index[k])]
  colnames(index) <- c("Year","name_index")
  index <- index  %>% group_by(Year) %>%
    summarise(name_index = mean(name_index))
  index <- data.frame(index)
  index <- index[,"name_index"]
  #####COrrelación
  
  c.matrix <- matrix(NA,length(lon),length(lat))
  t.matrix <- matrix(NA,length(lon),length(lat))
  
  
  for (i in 1:length(lon)) {
    for (j in 1:length(lat)) {
      c.matrix[i,j] <- cor(SPI_arr[i,j,],index,method="spearman", use = "pairwise.complete.obs")
      if (all(is.na(SPI_arr[i,j,]))) {
        t.matrix[i, j] <- NA
      }else{
        t.matrix[i, j] <-cor.test(SPI_arr[i,j,],index,method="spearman")$p.value
      }
    }
  }
  
  
  grid <- expand.grid(x=lon, y=lat)
  grid$corr <- as.vector(c.matrix)
  grid$pval <- as.vector(t.matrix)
  
  sig <- subset(grid[, c(1, 2, 4)], pval < 0.05)
  
  
  graph[[k]] <- ggplot(data=peru)+
    #geom_tile(data=grid,aes(x,y,fill=corr))+
    geom_contour_fill(data=grid,aes(x,y,z=corr),breaks=seq(-1,1,0.1))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5,face = "bold"),legend.key.height = unit(4.5,"cm"))+
    scale_fill_gradientn(colours=rev(paletteer_c("ggthemes::Red-Blue-White Diverging",length(seq(-1,1,0.1))-1)),breaks=seq(-1,1,0.1),limits=c(-1,1))+
    geom_sf(fill="transparent",col="black",lwd=0.6)+
    geom_point(data=sig,aes(x,y),size=0.00001)+
    ggtitle(paste0("Correlación de Spearman SPI-12 RAIN4PE e Índice ",name_index[k] ," , pval<0.05"))+
    xlab("Longitud")+
    ylab("Latitud")+
    guides(fill=guide_colorsteps(title="Rho"))+
    coord_sf(xlim=c(-81,-69),ylim=c(-17.55,-0.8))
  
  ggsave(paste0("figures/",name_index[k],"_RAIN4PE.png"), width = 8.5, height = 10,units = 'in',dpi=500)
  
}

###############Tendencia

library(trend)
c.matrix <- matrix(NA,length(lon),length(lat))
t.matrix <- matrix(NA,length(lon),length(lat))


for (i in 1:length(lon)) {
  for (j in 1:length(lat)) {
    if (all(is.na(SPI_arr[i,j,]))){
      c.matrix[i, j] <- NA
      t.matrix[i,j] <- NA
  }else{
    #c.matrix[i,j] <- MannKendall(SPI_arr[i,j,])$tau
    c.matrix[i, j] <- mk.test(SPI_arr[i,j,])$estimates[3]
    t.matrix[i, j] <- mk.test(SPI_arr[i,j,])$p.value
  }

    
  }}




grid <- expand.grid(x=lon, y=lat)
grid$corr <- as.vector(c.matrix)
grid$pval <- as.vector(t.matrix)


#sig <- subset(grid[, c(1, 2, 4)], pval < 0.05) 
#sig <- SpatialPointsDataFrame(coords = sig[, c(1, 2)], data = sig)
#sig



library(ggplot2)
library(plotly)
library(sf)
library(metR)
library(paletteer)

sig <- subset(grid[, c(1, 2, 4)], pval < 0.05) 


ggplot(data=peru)+
  #geom_tile(data=grid,aes(x,y,fill=corr))+
  geom_contour_fill(data=grid,aes(x,y,z=corr),breaks=seq(-1,1,0.1))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),legend.key.height = unit(4.5,"cm"))+
  scale_fill_gradientn(colours=rev(paletteer_c("ggthemes::Classic Red-White-Black",length(seq(-1,1,0.1))-1)),breaks=seq(-1,1,0.1),limits=c(-1,1))+
  geom_sf(fill="transparent",col="black",lwd=0.6)+
  geom_point(data=sig,aes(x,y),size=0.00001)+
  ggtitle("Tendencia de Mann Kendall SPI-12 RAIN4PE , pval<0.05")+
  xlab("Longitud")+
  ylab("Latitud")+
  guides(fill=guide_colorsteps(title="Rho"))+
  coord_sf(xlim=c(-81,-69),ylim=c(-17.55,-0.8))

ggsave(paste0("figures/","Kendall_RAIN4PE.png"), width = 8.5, height = 10,units = 'in',dpi=500)

###Para correlacion spi pisco y rain4pe

SPI2 <- SPI_arr[2:127,,]

