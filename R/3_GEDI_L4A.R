library(hdf5r)
library(data.table)
library(raster)
library(ggplot2)
library(rgdal)
source("R/GEDI_FUN.R")

L4A.ls<-list.files("Data/GEDIL4A", pattern="h5",
                   recursive = TRUE, full.names = TRUE)
print(L4A.ls[1:4])

#Pick one!
gedi.file<-L4A.ls[1]

#read in the GEDI L2B file
L4A.h5<-readLevel4A(gedi.file)

#Let's look at a subset of the GEDI data structure
list.datasets(L4A.h5[["BEAM0000"]])

L4A.m<-getLevel4A(level4A = L4A.h5)


L4A.ls.all<-lapply(L4A.ls, function(x){
  L4A.m<-NULL
  L4A.h5<-readLevel4A(x)
  try(L4A.m<-getLevel4A(L4A.h5), silent = TRUE)
  return(L4A.m)
})

L4A.all<-do.call(rbind,L4A.ls.all)

#Make a spatial object
L4A.sp<-SpatialPointsDataFrame(coords=cbind(L4A.all$lon_lowestmode,L4A.all$lat_lowestmode),
                               data=L4A.all,
                               proj4string = crs(raster()))

extent.poly<-rgdal::readOGR("Data/extent/Zambezi_poly.shp")
ext.ex<-extract(extent.poly, L4A.sp)
# clip.col<-ncol(L4A.sp)

L4A.sp@data$ext.flag<-ext.ex[,2]
L4A.sp<-L4A.sp[!is.na(L4A.sp$ext.flag),]
# L4A.sp<-L4A.sp[,1:clip.col]
L4A.sp$shot_number<-as.numeric(L4A.sp$shot_number)

writeOGR(L4A.sp, "output/GEDI_L4A.shp", layer='shot_number', driver="ESRI Shapefile")


L4A.q1.sp<-L4A.sp[L4A.sp$degrade_flag==0&
                    L4A.sp$l2_quality_flag==1&
                    L4A.sp$l4_quality_flag==1,]
hist(L4A.q1.sp$agbd)
hist(L4A.q1.sp$sensitivity)
hist(L4A.q1.sp$agbd_pi_lower)

 
# min(L4A.q1.sp@data$agbd)
# max(L4A.q1.sp@data$agbd)
# mean(L4A.q1.sp@data$agbd)
# quantile(L4A.q1.sp@data$agbd, 0.99)

L4A.q1.sp<-L4A.q1.sp[(L4A.q1.sp$agbd>=0 )&
                    (L4A.q1.sp$agbd<=quantile(L4A.q1.sp@data$agbd, 0.99)),]


writeOGR(L4A.q1.sp, "output/GEDI_L4A_high_quality.shp", layer='shot_number', driver="ESRI Shapefile")

r<-raster::raster(ext = extent(L4A.q1.sp), resolution = 0.01)

agb.r<-rasterize(L4A.q1.sp, r, fun=function(x, na.rm=TRUE) mean(x), field='agbd')
agb.sd.r<-rasterize(L4A.q1.sp, r, fun=function(x, na.rm=TRUE) sd(x), field='agbd')

agb.s<-stack(agb.r,agb.sd.r)
names(agb.s)<-c("agb","agb_sd")
plot(agb.s)
plot(agb.r,agb.sd.r)

gedi.agb.global<-raster("Data/GEDIL4B/GEDI04_B_MW019MW138_02_002_05_R01000M_MU_clipped.tif")

gedi.agb.global.t<-projectRaster(gedi.agb.global, r)

gedi.agb.global<-crop(gedi.agb.global.t, r)

agb.r<-rasterize(L4A.q1.sp, gedi.agb.global, fun=function(x, na.rm=TRUE) mean(x), field='agbd')
gedi.agb.global[is.na(agb.r)]<-NA
agb.r[is.na(gedi.agb.global)]<-NA

plot(stack(agb.r,gedi.agb.global))

plot(agb.r,gedi.agb.global)
abline(0,1,col="blue")

plot(stack(agb.r, gedi.agb.global))

GLAD.CHM<-raster("Data/Aux/GLAD_GEDI.tif")

glad.chm.ex<-extract(GLAD.CHM, L4A.q1.sp, buffer = 25/2/111139, fun='mean',
                na.rm=T) # 25/2 m radius of GEDI footprint / ~111,139 m per degree of lat/lon 


L4A.q1.sp$GLAD.chm<-glad.chm.ex

ggplot(L4A.q1.sp@data[L4A.q1.sp$agbd<1000,], aes(GLAD.chm, agbd))+
  # geom_errorbar(aes(x=GLAD.chm, ymin=agbd-(agbd_pi_upper-agbd),ymax=agbd_pi_upper), alpha=0.2)+
  geom_point(size=0.5,alpha=0.2)+
  facet_wrap(~selected_algorithm)+
  geom_smooth()

ggplot(L4A.q1.sp@data[L4A.q1.sp$agbd<1000&
                        L4A.q1.sp$selected_algorithm==2,
                      # &
                      # L4A.q1.sp$agbd>80&
                      # L4A.q1.sp$GLAD.chm>5
], aes(GLAD.chm, agbd))+
  # geom_errorbar(aes(x=GLAD.chm, ymin=agbd-(agbd_pi_upper-agbd),ymax=agbd_pi_upper), alpha=0.2)+
  geom_point(size=0.5,alpha=0.2)+
  # facet_wrap(~selected_algorithm)+
  geom_smooth()


srtm<-raster("Data/Aux/Mangrove_hmax95_Mozambique.tiff")

srtm.chm.ex<-extract(srtm, L4A.q1.sp, buffer = 25/2/111139, fun='mean',
                     na.rm=T) # 25/2 m radius of GEDI footprint / ~111,139 m per degree of lat/lon 


L4A.q1.sp$srtm.chm<-srtm.chm.ex

ggplot(L4A.q1.sp@data, aes(L4A.q1.sp@data$srtm.chm, L4A.q1.sp@data$agbd))+
  # geom_errorbar(aes(x=GLAD.chm, ymin=agbd-(agbd_pi_upper-agbd),ymax=agbd_pi_upper), alpha=0.2)+
  geom_point(size=0.5,alpha=0.2)+
  facet_wrap(~selected_algorithm)+
  geom_smooth()

ggplot(L4A.q1.sp@data[L4A.q1.sp$agbd<1000&
                        L4A.q1.sp$selected_algorithm==1,], aes(GLAD.chm, agbd))+
  # geom_errorbar(aes(x=GLAD.chm, ymin=agbd-(agbd_pi_upper-agbd),ymax=agbd_pi_upper), alpha=0.2)+
  geom_point(size=0.5,alpha=0.2)+
  # facet_wrap(~selected_algorithm)+
  geom_smooth()


AGB.ALS<-raster("Data/ALS/agb_map.tif")
AGB.ALS[AGB.ALS>800]<-NA
AGB.dd<-raster::projectRaster(AGB.ALS, crs=crs(raster()), res=1/111139*7, method="ngb")
plot(AGB.dd)

agb.ex<-extract(AGB.dd, L4A.q1.sp, buffer = 25/2/111139, fun='mean',
                na.rm=T) # 25/2 m radius of GEDI footprint / ~111,139 m per degree of lat/lon 
L4A.q1.sp$agb.als<-as.data.frame(agb.ex)[,1]
L4A.q1.sp<-L4A.q1.sp[!is.na(L4A.q1.sp$agb.als),]

ggplot(L4A.q1.sp@data[L4A.q1.sp$agbd<1000,], aes(agb.als, agbd))+
  geom_errorbar(aes(x=agb.als, ymin=agbd-(agbd_pi_upper-agbd),ymax=agbd_pi_upper), alpha=0.2)+
  geom_point(size=0.5,alpha=0.2)+
  facet_wrap(~selected_algorithm)+
  geom_smooth(method="lm", formula=y~x-1)+
  geom_abline(slope=1,intercept=0)


