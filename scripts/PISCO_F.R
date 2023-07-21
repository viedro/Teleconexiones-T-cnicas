setwd("C:/Users/ASUS/Desktop/TF_Tecnicas/Tf_tecnicas/teleconexiones")

library(raster)
#####Donde pp.nc es la version estable de PISCOP mensual para todo el Peru (1981-2016)
nc <-brick("data/raw/pp.nc ")

###Paquete para el calculo del SPI
library(SPEI)
###Calculando el índice SPI usando la función calc
SPI <- calc(nc, fun= function(x, scale = 12, distribution = "Gamma", na.rm=T,...) as.numeric((spi(x,scale= scale,na.rm=na.rm,...))$fitted))

###Abriendo Shapefile de todo sudamerica
shp <- shapefile("data/raw/Sudamérica.shp")

###Filtrando Peru
shp <- shp[shp$PAÍS == "Perú" ,]
shp
###Estableciendo sistema de refererencia al rasterbrick de los indices
crs(SPI) <- crs(shp)
SPI

###Paquete para el grafciado de datos raster o netcdf
library(rasterVis)

###Visualizando los maximos o minimos del indice desde 1981 al 2016
max(maxValue(SPI),na.rm=T)
min(minValue(SPI),na.rm=T)

library(tidyverse)
###Recortando con el shapefile de Perú
recor <- raster::crop(SPI,shp) %>%
  raster::mask(shp)

#Convirtiendo a array
SPI_arr <- raster::as.array(recor)

str(SPI_arr)

dates <- getZ(nc)
dates <-as.Date("1960-01-01") %m+% months(dates-0.5)
head(dates)

###Extrayendo las coordenadas de latitud y longitud del raster anterior
lon <- raster::xFromCol(recor, 1:ncol(recor))
lat <- raster::yFromRow(recor, 1:nrow(recor))

###Vemos las longitudes de cada eje
c(length(lon),length(lat))

#Corrigiendo el orden del array
SPI_arr <- aperm(SPI_arr, c(2,1,3))
###Corregido
str(SPI_arr)


###La latitud no va de menor a mayor
lat <- rev(lat)  

#Inversión también de latitudes en el arreglo
SPI_arr <- SPI_arr[,ncol(SPI_arr):1 ,]

head(lat)

SPI_arr <- SPI_arr[,,seq(12,12*36,12)]

library(rsoi)
library(ggplot2)
library(plotly)
library(sf)
library(metR)
library(paletteer)

####Automatizando graficas de correlacion
functions <- list("download_oni","download_soi","download_aao","download_ao","download_pdo","download_dmi","download_mei","download_npgo")
name_index <- toupper(substr(functions,10,nchar(functions)))

graph <- list()
for (k in 1:length(functions)){
  result <- do.call(functions[[k]], list(use_cache = FALSE))
  result <- result[,c("Date","Year",name_index[k])]
  result <- result[result$Date>=as.Date("1981-01-01") & result$Date<=as.Date("2016-12-31"),]
  
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
    ggtitle(paste0("Correlación de Spearman SPI-12 PISCO e Índice ",name_index[k] ," , pval<0.05"))+
    xlab("Longitud")+
    ylab("Latitud")+
    guides(fill=guide_colorsteps(title="Rho"))+
    coord_sf(xlim=c(-81,-69),ylim=c(-17.55,-0.8))
  
  ggsave(paste0("figures/",name_index[k],"_PISCO.png"), width = 8.5, height = 10,units = 'in',dpi=500)
  
}


###############Tendencia de Mann Kendall

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



peru <- st_read("data/raw/Sudamérica.shp")
peru <- peru[peru$PAÍS == "Perú",]
sig <- subset(grid[, c(1, 2, 4)], pval < 0.05) 

ggplot(data=peru)+
  #geom_tile(data=grid,aes(x,y,fill=corr))+
  geom_contour_fill(data=grid,aes(x,y,z=corr),breaks=seq(-1,1,0.1))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),legend.key.height = unit(4.5,"cm"))+
  scale_fill_gradientn(colours=rev(paletteer_c("ggthemes::Classic Red-White-Black",length(seq(-1,1,0.1))-1)),breaks=seq(-1,1,0.1),limits=c(-1,1))+
  geom_sf(fill="transparent",col="black",lwd=0.6)+
  geom_point(data=sig,aes(x,y),size=0.00001)+
  ggtitle("Tendencia de Mann Kendall SPI-12 PISCO , pval<0.05")+
  xlab("Longitud")+
  ylab("Latitud")+
  guides(fill=guide_colorsteps(title="Rho"))+
  coord_sf(xlim=c(-81,-69),ylim=c(-17.55,-0.8))

ggsave(paste0("figures/","Kendall_PISCO.png"), width = 8.5, height = 10,units = 'in',dpi=500)


####Correlación SPI PISCO y RAIN4PE
###Del final del script de RAIN4PE
SPI2


c.matrix <- matrix(NA,length(lon),length(lat))
t.matrix <- matrix(NA,length(lon),length(lat))

for (i in 1:length(lon)) {
  for (j in 1:length(lat)) {
    c.matrix[i,j] <- cor(SPI_arr[i,j,1:35],SPI2[i,j,] ,method="spearman", use = "pairwise.complete.obs")
    if (all(is.na(SPI_arr[i,j,1:35]))) {
      t.matrix[i, j] <- NA
    }else{
      t.matrix[i, j] <-cor.test(SPI_arr[i,j,1:35],SPI2[i,j,] ,method="spearman")$p.value
    }
  }
}



grid <- expand.grid(x=lon, y=lat)
grid$corr <- as.vector(c.matrix)
grid$pval <- as.vector(t.matrix)


library(ggplot2)
library(plotly)
library(sf)
library(metR)
library(paletteer)

peru <- st_read("data/raw/Sudamérica.shp")
peru <- peru[peru$PAÍS == "Perú",]
sig <- subset(grid[, c(1, 2, 4)], pval < 0.05) 



ggplot(data=peru)+
  #geom_tile(data=grid,aes(x,y,fill=corr))+
  geom_contour_fill(data=grid,aes(x,y,z=corr),breaks=seq(-1,1,0.1))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),legend.key.height = unit(4.5,"cm"))+
  scale_fill_gradientn(colours=rev(paletteer_c("ggthemes::Red-Blue-White Diverging",length(seq(-1,1,0.1))-1)),breaks=seq(-1,1,0.1),limits=c(-1,1))+
  geom_sf(fill="transparent",col="black",lwd=0.6)+
  geom_point(data=sig,aes(x,y),size=0.00001)+
  ggtitle("Correlación de Spearman SPI-12 PISCO y RAIN4PE , pval<0.05")+
  xlab("Longitud")+
  ylab("Latitud")+
  guides(fill=guide_colorsteps(title="Rho"))+
  coord_sf(xlim=c(-81,-69),ylim=c(-17.55,-0.8))

ggsave(paste0("figures/","COR_spi_PISCO_RAIN4PE.png"), width = 8.5, height = 10,units = 'in',dpi=500)





