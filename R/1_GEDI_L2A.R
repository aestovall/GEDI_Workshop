library(hdf5r)
library(data.table)
library(raster)
library(ggplot2)
source("R/GEDI_FUN.R")

#get a list of GEDI L2A files
L2A.ls<-list.files("Data/GEDIL2A", recursive = TRUE, full.names = TRUE)

#Pick one!
gedi.file<-L2A.ls[1]

#read in the GEDI L2A file
L2A.h5<-readLevel2A(gedi.file)

#Let's look at the GEDI data structure
list.datasets(L2A.h5[["BEAM0000"]])

#now we can read all GEDI files into a list and combine

#we use a lapply function to do this as it can easily be made parallel
L2A.ls.all<-lapply(L2A.ls, function(x){
  L2A.m<-NULL
  L2A.h5<-readLevel2A(x)
  try(L2A.m<-getLevel2A(L2A.h5), silent = TRUE)
  return(L2A.m)
})

#lets actually do that! Make it parallel:
library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)

L2A.ls.all<-
  parLapply(cl, L2A.ls, function(x){
  
  library(hdf5r)
  library(data.table)
  library(raster)
  
  source("R/GEDI_FUN.R")
  
  L2A.m<-NULL
  try(L2A.h5<-readLevel2A(x), silent = TRUE)
  try(L2A.m<-getLevel2A(L2A.h5), silent = TRUE)
  return(L2A.m)
})
stopCluster(cl)

#combine the list into a single long data frame with the rbind command
L2A.all<-do.call(rbind,L2A.ls.all)

#Make a spatial object from the DF for subsetting and extraction. Don't forget the CRS!
L2A.sp<-SpatialPointsDataFrame(coords=cbind(L2A.all$lon_lowestmode,L2A.all$lat_lowestmode),
                               data=L2A.all,
                               proj4string = crs(raster()))
#we bring in our extent polygon and use it for subsetting.
extent.poly<-rgdal::readOGR("Data/extent/Zambezi_poly.shp")

#Why do we use extract? The clip function will not accept geometry, 
# so the square extent is the only thing it will clip to. If you want to use 
# complex geometry, use this method.
ext.ex<-extract(extent.poly, L2A.sp)

#add extracted values to the GEDI data
L2A.sp@data$ext.flag<-ext.ex[,2]

#now subset by the identified extent flag
L2A.sp<-L2A.sp[!is.na(L2A.sp$ext.flag),]

#we also need to convert the shot number column to an integer for writing purposes
L2A.sp$shot_number<-as.numeric(L2A.sp$shot_number)

#write the file!
writeOGR(L2A.sp, "GEDI_L2A.shp", layer='shot_number', driver="ESRI Shapefile")

#We can explore the GEDI data in QGIS to determine how to QAQC


###### QAQC Subsetting ############

#now let's subset based on our observations
L2A.q1.sp<-L2A.sp[L2A.sp$degrade_flag==0&
                    L2A.sp$quality_flag==1,]

#write the file!
writeOGR(L2A.q1.sp, "GEDI_L2A_high_quality.shp", layer='shot_number', driver="ESRI Shapefile")

#now we can take a look at some of our GEDI values
gedi.sub<-L2A.q1.sp@data[L2A.q1.sp@data$delta_time <= (104334330.891843348741531) &
                           L2A.q1.sp@data$delta_time >= (104334330.891843348741531-100) &
                           L2A.q1.sp@data$beam == "BEAM0101",]

#where are we? Let's look at an interactive map
library(leaflet)
pal <- colorNumeric(
  palette = "Blues",
  domain = gedi.sub$rh100)

m1<-leaflet() %>%
  addCircleMarkers(gedi.sub$lon_lowestmode,
                   gedi.sub$lat_lowestmode,
                   radius = 1,
                   opacity = 1,
                   color = pal(gedi.sub$rh100))  %>%
  addScaleBar(options = list(imperial = FALSE)) %>%
  addProviderTiles(providers$Esri.WorldImagery)  %>%
  addLegend(colors = "blue", labels= c("GEDI Shots"),
            title ="GEDI Level2A") 
m1

library(ggplot2)
ggplot(gedi.sub,
       aes(x=delta_time,y=elev_lowestmode))+
  geom_ribbon(aes(ymin=elev_lowestmode,ymax=rh100+elev_lowestmode),
              fill="forestgreen", alpha=0.5)+
  geom_path(aes(y=rh100+elev_lowestmode), color="forestgreen")+
  geom_path(color="black", linewidth=0.4)+
  theme_bw()+theme(panel.grid = element_blank())


#Let's compare to another 7 m resolution canopy height product from ALS
CHM.dd<-raster("Data/ALS/CHM_dd.tif")
chm.ex<-raster::extract(CHM.dd, 
                        L2A.q1.sp,
                        buffer = 25/2/111139, fun='mean')
#add extracted values to the GEDI data
L2A.q1.sp$CHM<-chm.ex

gedi.sub<-L2A.q1.sp@data[L2A.q1.sp@data$delta_time <= (104334330.891843348741531) &
                           L2A.q1.sp@data$delta_time >= (104334330.891843348741531-10) &
                           L2A.q1.sp@data$beam == "BEAM0101",]

library(ggplot2)
ggplot(gedi.sub,
       aes(x=delta_time,y=0))+
  geom_ribbon(aes(ymin=0,ymax=rh98),
              fill="forestgreen", alpha=0.5)+
  geom_path(aes(y=rh98), color="forestgreen")+
  geom_path(color="black", linewidth=0.4)+
  geom_path(aes(y=CHM),color="purple")
  theme_bw()+theme(panel.grid = element_blank())


L2A.all.sub<-L2A.q1.sp[!is.na(L2A.q1.sp$CHM)&
                         L2A.q1.sp$quality_flag==1&
                         L2A.q1.sp$sensitivity>0.95,]

library(ggplot2)
ggplot(L2A.all.sub@data, aes(rh98, CHM))+
  geom_point(aes(color=sensitivity))+
  facet_wrap(~selected_algorithm)+
  geom_abline(slope=1,intercept=0)


L2A.all.sub<-L2A.all[!is.na(L2A.all$CHM)&
                       L2A.all$quality_flag==1
                     # L2A.all$degrade_flag==0
                     ,]

plot(L2A.all.sub$rh98,L2A.all.sub$CHM)
abline(0,1, col="blue")




GLAD.CHM<-raster("Data/Aux/GLAD/GLAD_GEDI.tif")
GLAD.CHM.ex<-raster::extract(GLAD.CHM, 
                             L2A.q1.sp,
                             buffer = 25/2/111139, fun='mean')

L2A.q1.sp$GLAD.CHM<-GLAD.CHM.ex

lm(rh100~GLAD.CHM-1, L2A.q1.sp@data[L2A.q1.sp$GLAD.CHM<30&
                                      L2A.q1.sp$GLAD.CHM>4&
                                      L2A.q1.sp$rh100>=5&
                                      L2A.q1.sp$rh100<=30&
                                      L2A.q1.sp$selected_algorithm==2,])

GLAD.CHM[GLAD.CHM>40]<-NA
plot(GLAD.CHM*1.411)

chm.df<-as.data.frame(rasterToPoints(CHM.dd))
chm.ag.mean<-raster::rasterize(chm.df[,1:2], GLAD.CHM, field=chm.df$chm_map, fun=mean)
chm.ag.sd<-raster::rasterize(chm.df[,1:2], GLAD.CHM, field=chm.df$chm_map, fun=sd)
plot(stack(chm.ag.mean,chm.ag.sd))

chm_comp<-stack(GLAD.CHM,chm.ag.mean)
chm_comp.df<-as.data.frame(rasterToPoints(chm_comp))
colnames(chm_comp.df)[4]<-"chm.ag"

ggplot(chm_comp.df[!is.na(chm_comp.df$chm.ag&
                            chm_comp.df$chm.ag<35),], 
       aes(chm.ag,GLAD_GEDI*1.411))+
  geom_jitter(size=0.2,height=0.7, alpha=0.05)+
  # geom_hex(bins=100)+
  geom_smooth()+
  geom_abline(slope=1,intercept=0)

ggplot(L2A.q1.sp@data[L2A.q1.sp$GLAD.CHM<30&
                        L2A.q1.sp$GLAD.CHM>4&
                        L2A.q1.sp$rh100>=5&
                        L2A.q1.sp$rh100<=30 ,], aes(GLAD.CHM,rh100))+
  geom_jitter(aes(color=sensitivity), size=0.5,width=0.5)+
  # geom_hex(bins=100)+
  facet_wrap(~selected_algorithm)+
  geom_smooth(method="lm", formula=y~x-1)+
  geom_abline(slope=1,intercept=0)
