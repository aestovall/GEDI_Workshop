
L4A.q1.sp<-rgdal::readOGR("output/GEDI_L4A_high_quality.shp")

tdx<-raster("Data/Aux/TDX/Zambezi_TDX_CHM.tiff")
tdx.ex<-extract(tdx, L4A.q1.sp, buffer = 25/2/111139, fun='mean',
                na.rm=T) # 25/2 m radius of GEDI footprint / ~111,139 m per degree of lat/lon 
L4A.q1.sp$tdx<-tdx.ex

ggplot(L4A.q1.sp@data, aes(tdx, agbd))+
  geom_errorbar(aes(x=tdx, ymin=agbd-(agbd_pi_upper-agbd),ymax=agbd_pi_upper), 
                alpha=0.05, col="orange")+
  geom_point(size=0.5,alpha=0.2)+
  facet_wrap(~selected_algorithm)+
  geom_smooth()+
  coord_cartesian(ylim=c(0,400))+
  theme_bw()+
  theme(panel.grid = element_blank())

#let's create our model subset so we can make an AGB model
L4A.q1.sub<-L4A.q1.sp@data[L4A.q1.sp@data$agbd<=250&
                             L4A.q1.sp@data$agbd>0&
                             !is.na(L4A.q1.sp@data$agbd)&
                             !is.na(L4A.q1.sp@data$tdx)&
                             L4A.q1.sp@data$selected_algorithm==2,]


#fit a non-linear model
nls.m<-nls(agbd~a*tdx^b+c, start=list(a=0.1,b=1, c=50), L4A.q1.sub)
summary(nls.m)

#does a height metric look non-linear?
plot(L4A.q1.sub$tdx, L4A.q1.sub$agbd)

#look at the model
lines(0:35,predict(nls.m, newdata=data.frame(tdx=0:35)), col="red")

#check error
sqrt(mean((predict(nls.m)-L4A.q1.sub$agbd)^2))/mean(L4A.q1.sub$agbd)
sqrt(mean((predict(nls.m)-L4A.q1.sub$agbd)^2))


tdx.agb<-coef(nls.m)[1]*tdx^coef(nls.m)[2]+coef(nls.m)[3]
#import our ALS AGB raster
als.agb<-raster("Data/ALS/agb_map.tif")

#reproject the ALS raster to the same resolution as the TDX data
als.agb.t<-projectRaster(als.agb, tdx.agb, method='ngb')

#remove NA values
tdx.agb[is.na(als.agb.t)]<-NA

hist(als.agb.t-tdx.agb)

#look at the difference!
plot(als.agb.t-tdx.agb)

writeRaster(tdx.agb-als.agb.t, "output/agb_als_tdx_diff.tif")

#look at the difference!
plot(als.agb.t-tdx.agb)

agb.s.df<-as.data.frame(rasterToPoints(stack(als.agb.t,tdx.agb)))
colnames(agb.s.df)[3:4]<-c("als","tdx")
agb.s.df<-na.exclude(agb.s.df)

m<-lm(als~tdx, agb.s.df)
summary(m)

tdx.agb.unbiased<-2.695e+01+tdx.agb*2.021

writeRaster(tdx.agb.unbiased-als.agb.t, "output/agb_als_tdx_unbiased_diff.tif")

library(viridis)
ggplot(agb.s.df[2.695e+01+agb.s.df$tdx*2.021>=90, ], aes(als,tdx))+
  geom_hex(bins=250, color=NA)+
  coord_fixed(xlim=c(0,300),ylim=c(0,300))+
  geom_abline(intercept = 0, slope=1, color="blue")+
  scale_fill_viridis(option='magma')+
  scale_color_viridis(option='magma')+
  theme(panel.background = element_rect(fill="black"),
        panel.grid = element_blank())

ggplot(agb.s.df[2.695e+01+agb.s.df$tdx*2.021>=90, ], 
       aes(als,(2.695e+01+tdx*2.021)))+
  geom_hex(bins=250, color=NA)+
  coord_fixed(xlim=c(0,300),ylim=c(0,300))+
  geom_abline(intercept = 0, slope=1, color="blue")+
  scale_fill_viridis(option='magma')+
  scale_color_viridis(option='magma')+
  theme(panel.background = element_rect(fill="black"),
        panel.grid = element_blank())
