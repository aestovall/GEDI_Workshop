---
editor_options: 
  markdown: 
    wrap: 72
---

# Become a GEDI Master

## Installation of rGEDI

``` r
#The CRAN version:
install.packages("rGEDI")

#The development version:
library(devtools)
devtools::install_git("https://github.com/carlos-alberto-silva/rGEDI", dependencies = TRUE)

# loading rGEDI package
library(rGEDI)
```

## Find GEDI data within your study area (GEDI finder tool)

``` r
# Study area boundary box coordinates
ul_lat<- -44.0654
lr_lat<- -44.17246
ul_lon<- -13.76913
lr_lon<- -13.67646
```

Better yet, I like to import a raster or shapefile of the area of
interest to grab the bounding coordinates like this:

``` r
#pick your country of interest
countries<-rnaturalearth::ne_countries()

#read in a shapefile of interest
country<-readOGR("D:/US_conterminous_extent.shp")

#extent to get GEDI data
ext_bb<-extent(country)

# Study area boundary box coordinates (less manual. Woo!)
ul_lat<- ext_bb[c(4)]
lr_lat<- ext_bb[c(3)]
ul_lon<- ext_bb[c(1)]
lr_lon<- ext_bb[c(2)]
```

Note the date range of the data you desire. This might align with
growing season or some other event of interest.

``` r
# Specifying the date range
daterange=c("2019-07-01","2020-05-22")
```

Find available GEDI data:

``` r
# Get path to GEDI data
gLevel1B<-gedifinder(product="GEDI01_B",ul_lat, ul_lon, lr_lat, lr_lon,version="001",daterange=daterange)
gLevel2A<-gedifinder(product="GEDI02_A",ul_lat, ul_lon, lr_lat, lr_lon,version="001",daterange=daterange)
gLevel2B<-gedifinder(product="GEDI02_B",ul_lat, ul_lon, lr_lat, lr_lon,version="001",daterange=daterange)
```

## Downloading GEDI data

Using rGEDI (this is slow!)

``` r
# Set output dir for downloading the files
outdir=getwd()

# Downloading GEDI data
gediDownload(filepath=gLevel1B,outdir=outdir)
gediDownload(filepath=gLevel2A,outdir=outdir)
gediDownload(filepath=gLevel2B,outdir=outdir)
```

The method I have adopted uses Cygwin64 terminal to run a modified GEDI
download bash script:

``` r

gedi_ls<-gLevel2B

#this defines the ftp path, so we can sub in our local directory path
gedi_path<-as.vector(substr(gedi_ls,1,
                            nchar("https://e4ftl01.cr.usgs.gov/GEDI/GEDI02_B.001/2019.04.18/")))

#we modify the list so we know what we have already downloaded
gedi_ls<-as.vector(substr(gedi_ls, 
                          nchar("https://e4ftl01.cr.usgs.gov/GEDI/GEDI02_B.001/2019.04.18/")+1, 
                          200))

#remove GEDI files that are already downloaded to the local drive
gedi_ls_sub<-gedi_ls[!(gedi_ls %in% list.files("D:/GEDI/", pattern="h5"))]

#add the correct file path back to the files for downloading
gedi_ls_sub<-paste0(gedi_path[!(gedi_ls %in% list.files("D:/GEDI/", pattern="h5"))],gedi_ls_sub)

#write this file list
writeLines(gedi_ls_sub, "D:/GEDI/GEDI_ls_sub.txt")
```

Now (yes this is silly...) I paste the massive file list into the
download bash script. Here's where using Nathan's method would be a
winner if you could run it easily from R. I have tried automating the
creation of this script without success. Maybe we can integrate these
two methods for speed!

# Getting started preprocessing GEDI data

## Setup file

Now we will modify a simple setup file so we can blaze thorugh a bunch
of GEDI data:

``` r
library(h5)
```

```         
## Warning: package 'h5' was built under R version 3.6.1
```

``` r
source("R/GEDI_processing_FUN.R")
beams<-c("0000","0001","0010","0011",
         "0101","0110", "1000", "1011")

GEDI.dir<-"D:/GEDI/"
output.dir<-"D:/GEDI/output/"

file_ls<-subset.files(GEDI.dir, output.dir)
```

Here we want to modify the directory where we find the GEDI data
(`GEDI.dir`) and where the new text files will be written
(`output.dir`). `subset.files()` simply makes sure you are not
re-processing h5 GEDI files that already exist.

## GEDI data structure

First, let's take a look

``` r
print(file_ls[1:4])
```

```         
## [1] "D:/GEDI/GEDI02_B_2019146011955_O02550_T02940_02_001_01.h5"
## [2] "D:/GEDI/GEDI02_B_2019162192037_O02810_T02818_02_001_01.h5"
## [3] "D:/GEDI/GEDI02_B_2019175054509_O03003_T00196_02_001_01.h5"
## [4] "D:/GEDI/GEDI02_B_2019179101505_O03068_T05540_02_001_01.h5"
```

For an example of the GEDI data structure I'm going to choose a single
orbit that overlaps with Gabon and Pongara National Park:

``` r
file_ls<-"E:/Gabon/GEDI/current_release_correct/GEDI02_B_2019111040155_O02008_T04616_02_001_01.h5"
```

Normally, this portion of the script is included in a loop or similar,
but this is the basic method of reading h5 GEDI data in:

``` r
i=1

file <- h5::h5file(name = file_ls[i],
                 mode = "a")
```

Now that we have a GEDI file in R we need to look at the stucture. It
has a very complex hierarchical structure:

``` r
  data.sets<-h5::list.datasets(file)

data.sets[1:180]
```

```         
##   [1] "/BEAM0000/algorithmrun_flag"                       
##   [2] "/BEAM0000/ancillary/dz"                            
##   [3] "/BEAM0000/ancillary/l2a_alg_count"                 
##   [4] "/BEAM0000/ancillary/maxheight_cuttoff"             
##   [5] "/BEAM0000/ancillary/rg_eg_constraint_center_buffer"
##   [6] "/BEAM0000/ancillary/rg_eg_mpfit_max_func_evals"    
##   [7] "/BEAM0000/ancillary/rg_eg_mpfit_maxiters"          
##   [8] "/BEAM0000/ancillary/rg_eg_mpfit_tolerance"         
##   [9] "/BEAM0000/ancillary/signal_search_buff"            
##  [10] "/BEAM0000/ancillary/tx_noise_stddev_multiplier"    
##  [11] "/BEAM0000/beam"                                    
##  [12] "/BEAM0000/channel"                                 
##  [13] "/BEAM0000/cover"                                   
##  [14] "/BEAM0000/cover_z"                                 
##  [15] "/BEAM0000/delta_time"                              
##  [16] "/BEAM0000/fhd_normal"                              
##  [17] "/BEAM0000/geolocation/degrade_flag"                
##  [18] "/BEAM0000/geolocation/delta_time"                  
##  [19] "/BEAM0000/geolocation/digital_elevation_model"     
##  [20] "/BEAM0000/geolocation/elev_highestreturn"          
##  [21] "/BEAM0000/geolocation/elev_lowestmode"             
##  [22] "/BEAM0000/geolocation/elevation_bin0"              
##  [23] "/BEAM0000/geolocation/elevation_bin0_error"        
##  [24] "/BEAM0000/geolocation/elevation_lastbin"           
##  [25] "/BEAM0000/geolocation/elevation_lastbin_error"     
##  [26] "/BEAM0000/geolocation/height_bin0"                 
##  [27] "/BEAM0000/geolocation/height_lastbin"              
##  [28] "/BEAM0000/geolocation/lat_highestreturn"           
##  [29] "/BEAM0000/geolocation/lat_lowestmode"              
##  [30] "/BEAM0000/geolocation/latitude_bin0"               
##  [31] "/BEAM0000/geolocation/latitude_bin0_error"         
##  [32] "/BEAM0000/geolocation/latitude_lastbin"            
##  [33] "/BEAM0000/geolocation/latitude_lastbin_error"      
##  [34] "/BEAM0000/geolocation/local_beam_azimuth"          
##  [35] "/BEAM0000/geolocation/local_beam_elevation"        
##  [36] "/BEAM0000/geolocation/lon_highestreturn"           
##  [37] "/BEAM0000/geolocation/lon_lowestmode"              
##  [38] "/BEAM0000/geolocation/longitude_bin0"              
##  [39] "/BEAM0000/geolocation/longitude_bin0_error"        
##  [40] "/BEAM0000/geolocation/longitude_lastbin"           
##  [41] "/BEAM0000/geolocation/longitude_lastbin_error"     
##  [42] "/BEAM0000/geolocation/shot_number"                 
##  [43] "/BEAM0000/geolocation/solar_azimuth"               
##  [44] "/BEAM0000/geolocation/solar_elevation"             
##  [45] "/BEAM0000/l2a_quality_flag"                        
##  [46] "/BEAM0000/l2b_quality_flag"                        
##  [47] "/BEAM0000/land_cover_data/landsat_treecover"       
##  [48] "/BEAM0000/land_cover_data/modis_nonvegetated"      
##  [49] "/BEAM0000/land_cover_data/modis_nonvegetated_sd"   
##  [50] "/BEAM0000/land_cover_data/modis_treecover"         
##  [51] "/BEAM0000/land_cover_data/modis_treecover_sd"      
##  [52] "/BEAM0000/master_frac"                             
##  [53] "/BEAM0000/master_int"                              
##  [54] "/BEAM0000/num_detectedmodes"                       
##  [55] "/BEAM0000/omega"                                   
##  [56] "/BEAM0000/pai"                                     
##  [57] "/BEAM0000/pai_z"                                   
##  [58] "/BEAM0000/pavd_z"                                  
##  [59] "/BEAM0000/pgap_theta"                              
##  [60] "/BEAM0000/pgap_theta_error"                        
##  [61] "/BEAM0000/pgap_theta_z"                            
##  [62] "/BEAM0000/rg"                                      
##  [63] "/BEAM0000/rh100"                                   
##  [64] "/BEAM0000/rhog"                                    
##  [65] "/BEAM0000/rhog_error"                              
##  [66] "/BEAM0000/rhov"                                    
##  [67] "/BEAM0000/rhov_error"                              
##  [68] "/BEAM0000/rossg"                                   
##  [69] "/BEAM0000/rv"                                      
##  [70] "/BEAM0000/rx_processing/algorithmrun_flag_a1"      
##  [71] "/BEAM0000/rx_processing/algorithmrun_flag_a2"      
##  [72] "/BEAM0000/rx_processing/algorithmrun_flag_a3"      
##  [73] "/BEAM0000/rx_processing/algorithmrun_flag_a4"      
##  [74] "/BEAM0000/rx_processing/algorithmrun_flag_a5"      
##  [75] "/BEAM0000/rx_processing/algorithmrun_flag_a6"      
##  [76] "/BEAM0000/rx_processing/pgap_theta_a1"             
##  [77] "/BEAM0000/rx_processing/pgap_theta_a2"             
##  [78] "/BEAM0000/rx_processing/pgap_theta_a3"             
##  [79] "/BEAM0000/rx_processing/pgap_theta_a4"             
##  [80] "/BEAM0000/rx_processing/pgap_theta_a5"             
##  [81] "/BEAM0000/rx_processing/pgap_theta_a6"             
##  [82] "/BEAM0000/rx_processing/pgap_theta_error_a1"       
##  [83] "/BEAM0000/rx_processing/pgap_theta_error_a2"       
##  [84] "/BEAM0000/rx_processing/pgap_theta_error_a3"       
##  [85] "/BEAM0000/rx_processing/pgap_theta_error_a4"       
##  [86] "/BEAM0000/rx_processing/pgap_theta_error_a5"       
##  [87] "/BEAM0000/rx_processing/pgap_theta_error_a6"       
##  [88] "/BEAM0000/rx_processing/rg_a1"                     
##  [89] "/BEAM0000/rx_processing/rg_a2"                     
##  [90] "/BEAM0000/rx_processing/rg_a3"                     
##  [91] "/BEAM0000/rx_processing/rg_a4"                     
##  [92] "/BEAM0000/rx_processing/rg_a5"                     
##  [93] "/BEAM0000/rx_processing/rg_a6"                     
##  [94] "/BEAM0000/rx_processing/rg_eg_amplitude_a1"        
##  [95] "/BEAM0000/rx_processing/rg_eg_amplitude_a2"        
##  [96] "/BEAM0000/rx_processing/rg_eg_amplitude_a3"        
##  [97] "/BEAM0000/rx_processing/rg_eg_amplitude_a4"        
##  [98] "/BEAM0000/rx_processing/rg_eg_amplitude_a5"        
##  [99] "/BEAM0000/rx_processing/rg_eg_amplitude_a6"        
## [100] "/BEAM0000/rx_processing/rg_eg_amplitude_error_a1"  
## [101] "/BEAM0000/rx_processing/rg_eg_amplitude_error_a2"  
## [102] "/BEAM0000/rx_processing/rg_eg_amplitude_error_a3"  
## [103] "/BEAM0000/rx_processing/rg_eg_amplitude_error_a4"  
## [104] "/BEAM0000/rx_processing/rg_eg_amplitude_error_a5"  
## [105] "/BEAM0000/rx_processing/rg_eg_amplitude_error_a6"  
## [106] "/BEAM0000/rx_processing/rg_eg_center_a1"           
## [107] "/BEAM0000/rx_processing/rg_eg_center_a2"           
## [108] "/BEAM0000/rx_processing/rg_eg_center_a3"           
## [109] "/BEAM0000/rx_processing/rg_eg_center_a4"           
## [110] "/BEAM0000/rx_processing/rg_eg_center_a5"           
## [111] "/BEAM0000/rx_processing/rg_eg_center_a6"           
## [112] "/BEAM0000/rx_processing/rg_eg_center_error_a1"     
## [113] "/BEAM0000/rx_processing/rg_eg_center_error_a2"     
## [114] "/BEAM0000/rx_processing/rg_eg_center_error_a3"     
## [115] "/BEAM0000/rx_processing/rg_eg_center_error_a4"     
## [116] "/BEAM0000/rx_processing/rg_eg_center_error_a5"     
## [117] "/BEAM0000/rx_processing/rg_eg_center_error_a6"     
## [118] "/BEAM0000/rx_processing/rg_eg_chisq_a1"            
## [119] "/BEAM0000/rx_processing/rg_eg_chisq_a2"            
## [120] "/BEAM0000/rx_processing/rg_eg_chisq_a3"            
## [121] "/BEAM0000/rx_processing/rg_eg_chisq_a4"            
## [122] "/BEAM0000/rx_processing/rg_eg_chisq_a5"            
## [123] "/BEAM0000/rx_processing/rg_eg_chisq_a6"            
## [124] "/BEAM0000/rx_processing/rg_eg_flag_a1"             
## [125] "/BEAM0000/rx_processing/rg_eg_flag_a2"             
## [126] "/BEAM0000/rx_processing/rg_eg_flag_a3"             
## [127] "/BEAM0000/rx_processing/rg_eg_flag_a4"             
## [128] "/BEAM0000/rx_processing/rg_eg_flag_a5"             
## [129] "/BEAM0000/rx_processing/rg_eg_flag_a6"             
## [130] "/BEAM0000/rx_processing/rg_eg_gamma_a1"            
## [131] "/BEAM0000/rx_processing/rg_eg_gamma_a2"            
## [132] "/BEAM0000/rx_processing/rg_eg_gamma_a3"            
## [133] "/BEAM0000/rx_processing/rg_eg_gamma_a4"            
## [134] "/BEAM0000/rx_processing/rg_eg_gamma_a5"            
## [135] "/BEAM0000/rx_processing/rg_eg_gamma_a6"            
## [136] "/BEAM0000/rx_processing/rg_eg_gamma_error_a1"      
## [137] "/BEAM0000/rx_processing/rg_eg_gamma_error_a2"      
## [138] "/BEAM0000/rx_processing/rg_eg_gamma_error_a3"      
## [139] "/BEAM0000/rx_processing/rg_eg_gamma_error_a4"      
## [140] "/BEAM0000/rx_processing/rg_eg_gamma_error_a5"      
## [141] "/BEAM0000/rx_processing/rg_eg_gamma_error_a6"      
## [142] "/BEAM0000/rx_processing/rg_eg_niter_a1"            
## [143] "/BEAM0000/rx_processing/rg_eg_niter_a2"            
## [144] "/BEAM0000/rx_processing/rg_eg_niter_a3"            
## [145] "/BEAM0000/rx_processing/rg_eg_niter_a4"            
## [146] "/BEAM0000/rx_processing/rg_eg_niter_a5"            
## [147] "/BEAM0000/rx_processing/rg_eg_niter_a6"            
## [148] "/BEAM0000/rx_processing/rg_eg_sigma_a1"            
## [149] "/BEAM0000/rx_processing/rg_eg_sigma_a2"            
## [150] "/BEAM0000/rx_processing/rg_eg_sigma_a3"            
## [151] "/BEAM0000/rx_processing/rg_eg_sigma_a4"            
## [152] "/BEAM0000/rx_processing/rg_eg_sigma_a5"            
## [153] "/BEAM0000/rx_processing/rg_eg_sigma_a6"            
## [154] "/BEAM0000/rx_processing/rg_eg_sigma_error_a1"      
## [155] "/BEAM0000/rx_processing/rg_eg_sigma_error_a2"      
## [156] "/BEAM0000/rx_processing/rg_eg_sigma_error_a3"      
## [157] "/BEAM0000/rx_processing/rg_eg_sigma_error_a4"      
## [158] "/BEAM0000/rx_processing/rg_eg_sigma_error_a5"      
## [159] "/BEAM0000/rx_processing/rg_eg_sigma_error_a6"      
## [160] "/BEAM0000/rx_processing/rg_error_a1"               
## [161] "/BEAM0000/rx_processing/rg_error_a2"               
## [162] "/BEAM0000/rx_processing/rg_error_a3"               
## [163] "/BEAM0000/rx_processing/rg_error_a4"               
## [164] "/BEAM0000/rx_processing/rg_error_a5"               
## [165] "/BEAM0000/rx_processing/rg_error_a6"               
## [166] "/BEAM0000/rx_processing/rv_a1"                     
## [167] "/BEAM0000/rx_processing/rv_a2"                     
## [168] "/BEAM0000/rx_processing/rv_a3"                     
## [169] "/BEAM0000/rx_processing/rv_a4"                     
## [170] "/BEAM0000/rx_processing/rv_a5"                     
## [171] "/BEAM0000/rx_processing/rv_a6"                     
## [172] "/BEAM0000/rx_processing/rx_energy_a1"              
## [173] "/BEAM0000/rx_processing/rx_energy_a2"              
## [174] "/BEAM0000/rx_processing/rx_energy_a3"              
## [175] "/BEAM0000/rx_processing/rx_energy_a4"              
## [176] "/BEAM0000/rx_processing/rx_energy_a5"              
## [177] "/BEAM0000/rx_processing/rx_energy_a6"              
## [178] "/BEAM0000/rx_processing/shot_number"               
## [179] "/BEAM0000/rx_range_highestreturn"                  
## [180] "/BEAM0000/rx_sample_count"
```

**WOW!** That is a lot of info to grab (180 categories X 8 beams!). So
you need to choose what is important to you for processing.

Now, within each beam we can subset by the specific data we need.

# GEDI Preprocessing

## Preprocessing Function

Ideally we want to do all of this automatically for lots of GEDI data,
so we can use this simple function:

``` r
preprocess.GEDI(input.file=x, 
                      filename=filename)
```

`preprocess.GEDI` takes h5 files, selects useful metrics, and outputs a
simple txt file for further processing. This greatly simplifies the
pipeline. We can modify this function to include specific information
from the h5 file, but for now it includes:

"lon" "lat" "elev_lowestmode" "elev_TDX" "rh100" "cover" "pai" "fhd"
"num_detectedmodes" "modis_treecover" "q1" "surface_flag" "sensitivity"
"shot_number" "delta_time" "beam"

I find these variables give me everything I need from the L2B data
product.

## Process in parallel

Now we can use the `parallel` package to process everything with
lightning speeds!

``` r
library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)

parLapply(cl, file_ls, function(x){
  
  source("R/GEDI_setup.R")
  
  filename<-gsub(".h5",".txt", 
                 gsub(GEDI.dir,output.dir, x) )
  
  try(preprocess.GEDI(input.file=x, 
                      filename=filename), silent=TRUE)
})

stopCluster(cl)
```

Check the output directory to see the GEDI files being written to your
drive.

## Read a processed GEDI .txt file

Now we can easily use the preprocessed GEDI data to do some science.

``` r
#Get our processed GEDI files
GEDI.processed.files<-list.files(output.dir, pattern="txt", full.names = TRUE)

library(data.table)
#read a single orbit file in
GEDI<-data.table::fread(GEDI.processed.files[1])

#subset to "quality 1" shots (the good ones!)
GEDI<-GEDI[GEDI$q1==1,]

#load some some packages for viz
library(ggplot2)
```

```         
## Registered S3 methods overwritten by 'ggplot2':
##   method         from 
##   [.quosures     rlang
##   c.quosures     rlang
##   print.quosures rlang
```

``` r
library(reshape)
```

```         
## 
## Attaching package: 'reshape'
```

```         
## The following object is masked from 'package:data.table':
## 
##     melt
```

``` r
GEDI.long<-melt(GEDI, id.vars="shot_number", measure.vars = colnames(GEDI)[-c(11,12,14, 16)])

ggplot(GEDI.long,
       aes(x=value))+
  geom_density(fill="red", color="black")+
  facet_wrap(~variable, scales="free")+
  theme_bw()
```

![](GEDI_intro_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

With this we can see the distribution of values for each variable. Now
we can visualize some GEDI transects:

Let's take a look at the profiles with elevation and canopy height:

``` r
ggplot(GEDI.sub[order(GEDI.sub$delta_time),],
       aes(x=delta_time,y=elev_lowestmode))+
  
  geom_ribbon(aes(ymin=elev_lowestmode,ymax=rh100+elev_lowestmode),
              fill="forestgreen")+
  geom_path(aes(y=rh100+elev_lowestmode), color="forestgreen")+
  geom_path(color="black", size=0.1)+
  theme_bw()+theme(panel.grid = element_blank())
```

![](GEDI_intro_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

``` r
ggplot(GEDI.sub[order(GEDI.sub$delta_time),],
       aes(x=delta_time,y=rh100))+
  geom_ribbon(aes(ymin=0,ymax=rh100),
              fill="forestgreen", alpha=0.5)+
  geom_path(color="forestgreen")+
  theme_bw()+theme(panel.grid = element_blank())
```

![](GEDI_intro_files/figure-html/unnamed-chunk-8-2.png)<!-- -->

Where are we??? Let's look in a map:

```{r echo=TRUE}
library(leaflet)

pal <- colorNumeric(
  palette = "Blues",
  domain = GEDI.sub$rh100)

m1<-leaflet() %>%
  addCircleMarkers(GEDI.sub$lon,
                   GEDI.sub$lat,
                   radius = 1,
                   opacity = 1,
                   color = pal(GEDI.sub$rh100))  %>%
  addScaleBar(options = list(imperial = FALSE)) %>%
  addProviderTiles(providers$Esri.WorldImagery)  %>%
  addLegend(colors = "blue", labels= c("GEDI Shots"),
            title ="GEDI Level2B") 
m1
```

# Gridding GEDI data

Now that we have a bunch of GEDI data processed we can get to work on
gridding some GEDI variables. Let's focus on elevation and canopy height
(rh100) for this tutorial.

First we decide where we want to grid the GEDI data...

In honor of the acceped paper, let's check out Bangladesh!

```         
## Loading required package: sp
```

```         
## rgdal: version: 1.4-4, (SVN revision 833)
##  Geospatial Data Abstraction Library extensions to R successfully loaded
##  Loaded GDAL runtime: GDAL 2.2.3, released 2017/11/20
##  Path to GDAL shared files: C:/Users/aestoval/Documents/R/R-3.6.0/library/rgdal/gdal
##  GDAL binary built with GEOS: TRUE 
##  Loaded PROJ.4 runtime: Rel. 4.9.3, 15 August 2016, [PJ_VERSION: 493]
##  Path to PROJ.4 shared files: C:/Users/aestoval/Documents/R/R-3.6.0/library/rgdal/proj
##  Linking to sp version: 1.3-1
```

```         
## 
## Attaching package: 'raster'
```

```         
## The following object is masked from 'package:data.table':
## 
##     shift
```

![](GEDI_intro_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

Here we just got a list of GEDI files available in Bangladesh and subset
based on the files we have available locally.

Now, want to put all of the GEDI shots within Bangladesh in one file
before we grid. First, we decide on variables and output file name:

``` r
#What are the variables we can choose from?
vars<-colnames(data.table::fread(GEDI_files[1], nrow=1))[-(1:2)]
print(vars)
```

```         
##  [1] "elev_lowestmode"   "elev_TDX"          "rh100"            
##  [4] "cover"             "pai"               "fhd"              
##  [7] "num_detectedmodes" "modis_treecover"   "q1"               
## [10] "surface_flag"      "sensitivity"       "shot_number"      
## [13] "delta_time"        "beam"
```

``` r
#choose the variables you want
output_variables<-vars[c(1,3)]

#Name the output file based on the project and variables included
Q1_output_file<-paste("D:/GEDI/output/rasters/Bangladesh_",
                      paste(output_variables, collapse = "_"),".txt", sep="")
```

Now, we iteratively open GEDI txt files, clip the orbits, and save as a
single new file:

``` r


library(sp)

for (i in 1:length(gedi_ls_sub)) {
  GEDI_df<-na.omit(data.table::fread(gedi_ls_sub[i], 
                             sep=" ", 
                             select = c('lon','lat',output_variables,'q1'), 
                             showProgress = FALSE))
  GEDI_sp<-SpatialPointsDataFrame(GEDI_df[,1:2],data=GEDI_df)
  GEDI_sp<-crop(GEDI_sp, country)
  
  if(is.null(GEDI_sp)) next else {
    GEDI_df<-GEDI_sp@data
  
  data.table::fwrite(GEDI_df[GEDI_df$q1==1,1:(ncol(GEDI_df)-1)], 
                     file = Q1_output_file,
                     sep = " ",
                     na = "NA",
                     append=TRUE)
  }
  
  print(paste0((i/length(gedi_ls_sub))*100,' %'))
} 

```

Now we can read our newly created file and take a look:

``` r
GEDI<-fread(Q1_output_file)

hist(GEDI$elev_lowestmode, breaks=100)
```

![](GEDI_intro_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

``` r
hist(GEDI$rh100, breaks=100)
```

![](GEDI_intro_files/figure-html/unnamed-chunk-12-2.png)<!-- -->

Finally, we can easily convert the GEDI shots into spatial points and
rasterize them to a grid resolution of our choice:

``` r
GEDI_sp<-SpatialPointsDataFrame(coords=cbind(GEDI[,1],GEDI[,2]),
                                data=GEDI)


library(raster)
r<-raster(ext=extent(country), resolution=0.01)
GEDI_r<-rasterize(GEDI_sp, r, field="elev_lowestmode",fun=function(x, na.rm=TRUE) mean(x))
GEDI_rh100<-rasterize(GEDI_sp, r, field="rh100",fun=function(x, na.rm=TRUE) mean(x))


plot(GEDI_r, col=viridis::inferno(250))
plot(country, add=TRUE)
```

![](GEDI_intro_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

``` r
GEDI_rh100[GEDI_rh100>20]<-20
plot(GEDI_rh100,col=viridis::viridis(250))
plot(country, add=TRUE)
```

![](GEDI_intro_files/figure-html/unnamed-chunk-13-2.png)<!-- -->
