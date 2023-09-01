library(hdf5r)
library(data.table)
library(raster)
library(ggplot2)
library(rgdal)
source("R/GEDI_FUN.R")

extent.poly<-rgdal::readOGR("Data/extent/Zambezi_poly.shp")

#get a list of GEDI L2B files
L2B.ls<-list.files("Data/GEDIL2B", recursive = TRUE, full.names = TRUE)
print(L2B.ls[1:4])

#Pick one!
gedi.file<-L2B.ls[1]

#read in the GEDI L2B file
L2B.h5<-readLevel2B(gedi.file)

#Let's look at the GEDI data structure
list.datasets(L2B.h5[["BEAM0000"]])

#now we can read all GEDI files into a list and combine

#we use a lapply function to do this as it can easily be made parallel
L2B.ls.all<-lapply(L2B.ls, function(x){
  L2B.m<-NULL
  L2B.h5<-readLevel2B(x)
  try(L2B.m<-getLevel2BVPM(L2B.h5), silent = TRUE)
  return(L2B.m)
})

L2B.all<-do.call(rbind,L2B.ls.all)

#Make a spatial object
L2B.sp<-SpatialPointsDataFrame(coords=cbind(L2B.all$longitude_lastbin,L2B.all$latitude_lastbin),
                               data=L2B.all,
                               proj4string = crs(raster()))
ext.ex<-extract(extent.poly, L2B.sp)

L2B.sp@data$ext.flag<-ext.ex[,2]
L2B.sp<-L2B.sp[!is.na(L2B.sp$ext.flag),]
L2B.sp$shot_number<-as.numeric(L2B.sp$shot_number)

writeOGR(L2B.sp, "output/GEDI_L2B.shp", layer='shot_number', driver="ESRI Shapefile")


L2B.q1.sp<-L2B.sp[L2B.sp$l2b_quality_flag==1,]

writeOGR(L2B.q1.sp, "output/GEDI_L2B_high_quality.shp", layer='shot_number', driver="ESRI Shapefile")



L2B.pavd.layers.ls.all<-lapply(L2B.ls, function(x){
  L2B.pavd<-NULL
  try(L2B.h5<-readLevel2B(x), silent = TRUE)
  
  try(L2B.pavd<-getLevel2BPAVDProfile(L2B.h5), silent = TRUE)
  return(L2B.pavd)
})

L2B.pavd.layers.all<-do.call(rbind,L2B.pavd.layers.ls.all)

L2B.pavd.layers.sp<-SpatialPointsDataFrame(coords=cbind(L2B.pavd.layers.all$lon_lowestmode,
                                            L2B.pavd.layers.all$lat_lowestmode),
                               data=L2B.pavd.layers.all,
                               proj4string = crs(raster()))
ext.ex<-extract(extent.poly, L2B.pavd.layers.sp)

L2B.pavd.layers.sp@data$ext.flag<-ext.ex[,2]
L2B.pavd.layers.sp<-L2B.pavd.layers.sp[!is.na(L2B.pavd.layers.sp$ext.flag),]
L2B.pavd.layers.sp$shot_number<-as.numeric(L2B.pavd.layers.sp$shot_number)

writeOGR(L2B.pavd.layers.sp, "output/GEDI_L2B_layers.shp", layer='shot_number', driver="ESRI Shapefile")


L2B.pavd.layers.sp<-L2B.pavd.layers.sp[L2B.pavd.layers.sp$l2b_quality_flag==1,]

r<-raster::raster(ext = extent(L2B.pavd.layers.sp), resolution = 0.01)

pavd.5m.r<-rasterize(L2B.pavd.layers.sp, r, fun=function(x, na.rm=TRUE) mean(x), field='pavd_z0_5m')
pavd.10m.r<-rasterize(L2B.pavd.layers.sp, r, fun=function(x, na.rm=TRUE) mean(x), field='pavd_z5_10m')
pavd.15m.r<-rasterize(L2B.pavd.layers.sp, r, fun=function(x, na.rm=TRUE) mean(x), field='pavd_z10_15m')
pavd.20m.r<-rasterize(L2B.pavd.layers.sp, r, fun=function(x, na.rm=TRUE) mean(x), field='pavd_z15_20m')
pavd.25m.r<-rasterize(L2B.pavd.layers.sp, r, fun=function(x, na.rm=TRUE) mean(x), field='pavd_z20_25m')
pavd.30m.r<-rasterize(L2B.pavd.layers.sp, r, fun=function(x, na.rm=TRUE) mean(x), field='pavd_z25_30m')

pavd.s<-stack(pavd.5m.r,
              pavd.10m.r,
              pavd.15m.r,
              pavd.20m.r,
              pavd.25m.r,
              pavd.30m.r)
names(pavd.s)<-c("pavd_5m",
                 "pavd_10m",
                 "pavd_15m",
                 "pavd_20m",
                 "pavd_25m",
                 "pavd_30m")
# animate(pavd.s)
plot(pavd.s)




#we use a lapply function to do this as it can easily be made parallel
L2B.pavd.ls.all<-lapply(L2B.ls, function(x){
  L2B.m<-NULL
  try(L2B.h5<-readLevel2B(x), silent = TRUE)
  try(L2B.m<-getLevel2BVPM(L2B.h5), silent = TRUE)
  
  if(!is.null(L2B.m)){
    
    L2B.pavd<-getLevel2BPAVDProfile(L2B.h5)
    L2B.pavd.long<-getLevel2BPAVDlong(L2B.pavd, L2B.m)
    
    L2B.pavd.long<-L2B.pavd.long[L2B.pavd.long$pavd_z!=-9999,]
    return(L2B.pavd.long)
  }
})

L2B.pavd.all<-do.call(rbind,L2B.pavd.ls.all)

L2B.pavd.all$z.bin<-ceiling(L2B.pavd.all$rh100/100/5)*5

L2B.pavd.long.ag<-aggregate(pavd_z~Z+z.bin, FUN="mean",L2B.pavd.all)
L2B.pavd.long.ag$sd<-aggregate(pavd_z~Z+z.bin, FUN="sd",L2B.pavd.all)[,3]
L2B.pavd.long.ag$sd[is.na(L2B.pavd.long.ag$sd)]<-0

ggplot(L2B.pavd.long.ag[L2B.pavd.long.ag$z.bin<=40&
                          L2B.pavd.long.ag$pavd_z>=0  ,], aes(Z, pavd_z))+
  geom_path(data=L2B.pavd.all[L2B.pavd.all$z.bin<=40&
                                L2B.pavd.all$pavd_z>=0,],aes(group=shot_number),linewidth=0.1)+
  geom_ribbon(aes(x=Z, ymin=pavd_z-sd, ymax=pavd_z+sd),alpha=0.4, fill="forestgreen")+
  geom_path(size=1, color="forestgreen")+
  facet_wrap(~reorder(paste0(z.bin-5," - ",z.bin," m"),z.bin), nrow=2)+
  coord_flip()+
  theme_bw()+
  theme(panel.grid =element_blank())

