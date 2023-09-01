getLevel2A<-function(L2A.h5){
  level2a <- L2A.h5
  groups_id <- grep("BEAM\\d{4}$", gsub("/", "", hdf5r::list.groups(level2a, 
                                                                    recursive = F)), value = T)
  rh.dt <- data.table::data.table()
  # pb <- utils::txtProgressBar(min = 0, max = length(groups_id), 
  #                             style = 3)
  i.s = 0
  L2A.all<-do.call(rbind,lapply(groups_id, function(x) {{
    level2a_i <- level2a[[x]]
    if (any(hdf5r::list.datasets(level2a_i) == "shot_number")) {
      if (length(level2a_i[["rh"]]$dims) == 2) {
        rh = t(level2a_i[["rh"]][, ])
      }
      else {
        rh = t(level2a_i[["rh"]][])
      }
      rhs <- data.table::data.table(beam <- rep(x, length(level2a_i[["shot_number"]][])), 
                                    shot_number = level2a_i[["shot_number"]][], degrade_flag = level2a_i[["degrade_flag"]][], 
                                    quality_flag = level2a_i[["quality_flag"]][], 
                                    delta_time = level2a_i[["delta_time"]][], sensitivity = level2a_i[["sensitivity"]][], 
                                    solar_elevation = level2a_i[["solar_elevation"]][], 
                                    lat_lowestmode = level2a_i[["lat_lowestmode"]][], 
                                    lon_lowestmode = level2a_i[["lon_lowestmode"]][], 
                                    elev_highestreturn = level2a_i[["elev_highestreturn"]][], 
                                    elev_lowestmode = level2a_i[["elev_lowestmode"]][], 
                                    selected_algorithm = level2a_i[["selected_algorithm"]][],
                                    rh)
      rh.dt <- rbind(rh.dt, rhs)
    }
  }
    colnames(rh.dt) <- c("beam", "shot_number", "degrade_flag", 
                         "quality_flag", "delta_time", "sensitivity", "solar_elevation", 
                         "lat_lowestmode", "lon_lowestmode", "elev_highestreturn", 
                         "elev_lowestmode", "selected_algorithm", paste0("rh", seq(0, 100)))
    # close(pb)
    return(rh.dt)
  }))
}


readLevel2A <-function(level2Apath) {
  level2a_h5 <- hdf5r::H5File$new(level2Apath, mode = 'r+')
  return(level2a_h5)
}

readLevel2B <-function(level2Bpath) {
  level2b_h5 <- hdf5r::H5File$new(level2Bpath, mode = 'r+')
  return(level2b_h5)
}

readLevel4A <-function(level4Apath) {
  level4a_h5 <- hdf5r::H5File$new(level4Apath, mode = 'r+')
  return(level4a_h5)
}

source("R/L2BmetricsFUN.R")
source("R/L4AmetricsFUN.R")

getLevel2BPAVDProfile<-function(level2b){
  level2b<-level2b
  groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                     hdf5r::list.groups(level2b, recursive = F)), value = T)
  m.dt<-data.table::data.table()
  i.s=0
  for ( i in groups_id){
    i.s<-i.s+1
    level2b_i<-level2b[[i]]
    m<-data.table::data.table(
      beam<-rep(i,length(level2b_i[["shot_number"]][])),
      shot_number=level2b_i[["shot_number"]][],
      algorithmrun_flag=level2b_i[["algorithmrun_flag"]][],
      l2b_quality_flag=level2b_i[["l2b_quality_flag"]][],
      delta_time=level2b_i[["geolocation/delta_time"]][],
      lat_lowestmode=level2b_i[["geolocation/lat_lowestmode"]][],
      lon_lowestmode=level2b_i[["geolocation/lon_lowestmode"]][],
      elev_highestreturn=level2b_i[["geolocation/elev_highestreturn"]][],
      elev_lowestmode=level2b_i[["geolocation/elev_lowestmode"]][],
      height_lastbin=level2b_i[["geolocation/height_lastbin"]][],
      height_bin0=level2b_i[["geolocation/height_bin0"]][],
      pavd_z=t(level2b_i[["pavd_z"]][,1:level2b_i[["pavd_z"]]$dims[2]]))
    m.dt<-rbind(m.dt,m)
  }
  colnames(m.dt)<-c("beam","shot_number","algorithmrun_flag",
                    "l2b_quality_flag","delta_time","lat_lowestmode",
                    "lon_lowestmode","elev_highestreturn",
                    "elev_lowestmode","height_lastbin",
                    "height_bin0",paste0("pavd_z",seq(0,30*5,5)[-31],"_",seq(5,30*5,5),"m"))
  return(m.dt)
}

getLevel2BPAVDlong<-function(L2B.pavd, L2B.m){
  L2B.pavd.t<-as.data.table(t(L2B.pavd[,c(paste0("pavd_z",paste0(seq(0,30,by=5),"_",seq(5,35, by=5),"m") ))]))
colnames(L2B.pavd.t)<-c(as.character(L2B.pavd$shot_number))
library(reshape2)
L2B.pavd.t$Z<-c(seq(5,35, by=5))
# L2B.pavd.t$shot_number<-NA
L2B.pavd.long<-data.table::melt(L2B.pavd.t, id.vars=c("Z"))
colnames(L2B.pavd.long)[2:3]<-c("shot_number","pavd_z")
L2B.pavd.long$shot_number<-as.character(factor(L2B.pavd.long$shot_number))
L2B.m$shot_number<-as.character(L2B.m$shot_number)
L2B.pavd.long.m<-merge(L2B.pavd.long,L2B.m, by="shot_number")
return(L2B.pavd.long.m)
}


