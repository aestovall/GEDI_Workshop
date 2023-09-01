# var.map.agb
# saveRDS(L4A.vars, file="R/l4a_vars.R")
var.map.agb<-readRDS("R/l4a_vars.R")

getLevel4A<-function(level4A, cols=c(
  "beam",
  "shot_number",
  "lon_lowestmode",
  "lat_lowestmode",
  "delta_time",
  "sensitivity",
  "solar_elevation",
  "selected_algorithm",
  "elev_lowestmode",
  "agbd_pi_lower",
  "agbd_pi_upper",
  "agbd",
  "agbd_se",
  "agbd_t",
  "agbd_t_se",
  "l2_quality_flag" ,
  "degrade_flag" ,
  "l4_quality_flag"
)){
  level4a<-level4A
  groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                     hdf5r::list.groups(level4a, recursive = F)), value = T)
  m.dt<-data.table::data.table()
  i.s=0
  
  #i = groups_id[1]
  L2B.all<-do.call(rbind,lapply(groups_id, function(x) {{
    i.s<-i.s+1
    level4a_i<-level4a[[x]]
    m<-data.table::data.table()
    #col = cols[1]
    for (col in cols) {
      h5.address = var.map.agb[[col]]
      if (is.null(h5.address)) {
        if (col %in% names(level4a_i)) {
          h5.address = col
        } else {
          if (i.s == 1) warning(sprintf("The column '%s' is not available in the GEDI2B product for BEAM %s!", col, i))
          next
        }
      }
      m[,eval(col):=level4a_i[[h5.address]][]]
    }
    m.dt<-rbind(m.dt,m)
  }
    return(m.dt)}))
}
