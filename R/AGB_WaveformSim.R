library(devtools)
# devtools::install_git("https://github.com/carlos-alberto-silva/rGEDI", dependencies = TRUE)
# install.packages('rGEDI')

library(rGEDI)
library(lidR)
library(plot3D)
library(rgdal)
proj.utm<-crs(rgdal::readOGR(list.files("Data/ALS/index",pattern="shp", full.names = TRUE)))

outdir<-"Data/ALS/plots"
las.files<-list.files(outdir, full.names = TRUE)

lapply(las.files, function(x){
  # Specifying the path to ALS data
  lasfile_plot <- x
  # Reading ALS file
  las_plot<-readLAS(lasfile_plot)
  
  # Extracting plot center geolocations
  xcenter_plot = (las_plot@header$`Max X`+las_plot@header$`Min X`)/2
  ycenter_plot = (las_plot@header$`Max Y`+las_plot@header$`Min Y`)/2
  
  # Simulating GEDI full-waveform
  wf_plot<-gediWFSimulator(input=lasfile_plot,
                           output=gsub(".las",".h5",lasfile_plot),
                           coords = c(xcenter_plot, ycenter_plot),
                           waveID=lasfile_plot)
  
  wf_plot_metrics<-gediWFMetrics(input=wf_plot,
                                 outRoot=file.path(getwd(), "plot"),
                                 rhRes = 1)
  
  metrics<-wf_plot_metrics
  rownames(metrics)<-c("plot")
  head(metrics[,1:8])
  
  if(x==las.files[1]) write.table(metrics,"GEDI_L2A_sim_metrics.txt", append = TRUE) else{
    write.table(metrics,"GEDI_L2A_sim_metrics.txt", append = TRUE, col.names = FALSE)
  }
  
  wf_plot_metrics_noise<-gediWFMetrics(input=wf_plot,
                                       outRoot=file.path(getwd(), "plot"),
                                       linkNoise= c(3.0103,0.95),
                                       maxDN= 4096,
                                       sWidth= 0.5,
                                       varScale= 3,
                                       rhRes = 1)
  close(wf_plot)
  metrics_noise<-rbind(wf_plot_metrics_noise)
  rownames(metrics_noise)<-c("plot")
  head(metrics_noise[,1:8])
  
  if(x==las.files[1]) write.table(metrics_noise,"GEDI_L2A_sim_noise_metrics.txt") else{
    write.table(metrics_noise,"GEDI_L2A_sim_noise_metrics.txt", append=TRUE, col.names = FALSE)
  }
  
})

library(ggplot2)
library(data.table)
library(reshape2)
library(reshape)
library(raster)
library(randomForest)
df<-fread("GEDI_L2A_sim_metrics.txt")
df$filename

#getting plot from file name
df$plot<-gsub(".h5","",
              gsub("/Users/atticusstovall/Documents/R/GEDI_tutorial/Data/ALS/plots/","",df$filename))

#Read plot data. This must be modified depending on input data format.
plots<-read.csv("Data/field_data/Zambezi_2013_Subplot_Data.csv")

#Make a spatial object
plots.sp<-SpatialPointsDataFrame(coords=cbind(plots$Longitude,plots$Latitude)[!is.na(plots$Longitude),1:2],
                                 data=plots[!is.na(plots$Longitude),1:32],
                                 proj4string = crs(raster()))

#Reproject plots to lat/lon
plots.sp.t<-sp::spTransform(plots.sp, proj.utm)

df$AGB<-plots.sp.t$Total.AGB[as.numeric(df$plot)]

df<-df[df$AGB!=9999,]

plot(df$`rhGauss 100`, df$AGB)

colnames(df)
metrics_all.m<-df[,c("AGB","filename",paste0(("rhGauss "), seq(0,100, by=1)),
                     "waveEnergy","FHD","cover")]

install.packages('leaps')
library(leaps)
library(leaps)

AGBreg<-regsubsets(AGB~., nbest=1, nvmax=5, 
                   data=metrics_all.m[,-c('filename')],
                   really.big = T)
summary(AGBreg)
plot(AGBreg)

mlr.m<-lm(AGB~`rhGauss 100`+`rhGauss 0`+`rhGauss 98`+`rhGauss 2`+waveEnergy, metrics_all.m)
summary(mlr.m)

#check RMSE?
sqrt(mean((predict(mlr.m)-metrics_all.m$AGB)^2))/mean(metrics_all.m$AGB)
sqrt(mean((predict(mlr.m)-metrics_all.m$AGB)^2))

#how does the model fit?
plot(predict(mlr.m), metrics_all.m$AGB)
abline(0,1)

plot(df$`true top`,df$`rhGauss 100`)
abline(0,1)

sqrt(mean((df$`true top`-df$`rhGauss 100`)^2))/mean(df$`true top`)
sqrt(mean((df$`true top`-df$`rhGauss 100`)^2))

# rf.m<-randomForest(AGB~., do.trace=TRUE,data=metrics_all.m,
#                    ntree=1000)
# plot(rf.m)
# names(importance(rf.m)[rev(order(importance(rf.m))),])[1:10]
# randomForest::varImpPlot(rf.m)
# 
# rf.m<-randomForest(AGB~zq90+zmax+zq55, do.trace=TRUE,data=metrics_all.m)

library(caret)
sqrt(mean((predict(rf.m)-metrics_all.m$AGB)^2))/mean(metrics_all.m$AGB)
sqrt(mean((predict(rf.m)-metrics_all.m$AGB)^2))

plot(predict(rf.m), metrics_all.m$AGB)
abline(0,1)



