# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Scripts performs area weighted regridding of gridded data.
# All configuration is done through ./skript_area-weighted-regridding_CONF.R, 
# which must be located in the same directory as this script.
#
# User can specify: 
#   1) several files on the input
#   2) variables to be processed for each file (not possible for geotiff) 
#      and how they should be named in output file
#   3) minimum coverage of output gridcell by input data
#   4) what columns from the output grid file should be transfered to output
#   5) input geotiff should be written as polygon shapefile
#
# Input files are recognised by file extension (case insensitive):
#   *.shp for shapefiles
#   *.tiff, *.tif for geotiffs
#   *.nc ... following assumptions must be met: 
#            1) NetCDF contains dimensions named lat/latitude lon/longitude 
#               dimensions, where geo. coordinates of regular grid are stored.
#            2) ...
#
# Input shapes can be:
#   POLYGON shapefile (not necessarily regular) or 
#   POINT shapefile - it is assumed that points resresent centers 
#     of regular grid and regular polygon grid is created based on 
#     input file extent and estimated resolution
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rm(list = ls())
library(raster)
library(sf)
library(tidyverse)
library(ncdf4)
library(rstudioapi)

scrDir <- dirname(rstudioapi::getSourceEditorContext()$path)

# default settings of non-compulsory variables
out_file                <- NULL
out_dir                 <- NULL
grd_out_cols_to_keep    <- NULL
min_coverage            <- 100

write_inp_netcdf_as_shp <- F
write_inp_raster_as_shp <- F
common_crs              <- NULL
# end of default settings

inp_nc_coord_dec <- 3

source(paste0(scrDir,"/","regrid_aw_CONF.R"))

setwd(wd)



# ======= START PROCESSING =======

tbeg <- Sys.time()

# name of output and log files
# ----------------------------
if (is.null(out_file)) {
  ifelse (length(inp_files)==1, 
          out_file <- gsub(".[a-zA-Z]*$","",basename(inp_files[[1]])), 
          out_file <- "lot_of_stuff")
  out_file <- paste0(out_file,"_ON_",basename(grd_out_file))
}
if (!is.null(out_dir)) {
  if (out_dir!="") {
    out_file <- paste0(out_dir,"/",out_file)  
  }
}

out_log <- gsub(".[a-zA-Z]*$",".log",out_file)

outd <- dirname(out_file)
if (!dir.exists(outd) & !outd%in%c("","./")) {
  dir.create(outd)
} else {
  if (file.exists(out_log)) file.remove(out_log) # Delete log file, if it exists  
}


# read output grid and its projection if common_crs not defined
# -------------------------------------------------------------
grd_out <- st_read(dsn = grd_out_file) 
out_crs <- st_crs(grd_out)
if (is.null(common_crs)) common_crs <- out_crs

grd_out <- st_transform(x = grd_out, crs = common_crs) %>%
  mutate(areaOrg = st_area(.) %>% as.numeric())
grd_out <- grd_out[c(grd_out_id,"areaOrg",grd_out_cols_to_keep)]

i <- 1
for (i in 1:length(inp_files)) {
  tbegf <- Sys.time()
  
  grd_inp_file   <- inp_files[[i]]
  cols_to_regrid <- colnames_inp[[i]]
  cols_out_name  <- colnames_out[[i]]
  
  message <- paste0("\nProcessing ",grd_inp_file)
  cat(paste0(message,"\n"))
  write(message, file=out_log, append=TRUE)
  
  # read input grided data. Filetype determined by extension
  ext <- unlist(strsplit(basename(grd_inp_file), split="\\."))
  tmp <- length(ext)
  ext <- ext[tmp]
  
  if (tolower(ext)=="shp") {
    grd_inp <- st_read(dsn = grd_inp_file)
    
    geom    <- unique(st_geometry_type(grd_inp))
    if (length(geom)==1) {
      if (geom=="POINT") {
        cat("\n   ... Converting input point shp to regular grid\n")
        t1 <- Sys.time()
        # step and number of cells of regular grid that would cover input data
        crd_lim <- st_bbox(grd_inp)
        crd     <- st_coordinates(grd_inp)
        step  <- c(NA,NA)
        ncell <- c(NA,NA)
        for (i in 1:2) {
          sorted <- sort(crd[,i])
          diffs  <- sorted[2:length(sorted)]-sorted[1:length(sorted)-1] 
          step[i]<- min(diffs[diffs > 0])  
          ncell[i] <- (crd_lim[i+2]-crd_lim[i])/step[i]+1
        }
        
        grd_poly <- st_make_grid(
          grd_inp,
          cellsize = step,
          offset   = crd_lim[c("xmin", "ymin")]-step/2,
          n = ncell,
          crs = st_crs(grd_inp),
          what = "polygons",
          square = TRUE,
          flat_topped = FALSE
        ) %>% st_sf()
        
        grd_inp <- st_join(x = grd_poly, y = grd_inp, left = T, join = st_contains)
        rm(grd_poly)
        
        message <- paste0("\n   ... transforming points to polygons took ", round(difftime(time1 = Sys.time(), time2 = t1, units = "min"),1)," min")
        write(message, file=out_log, append=TRUE)
      }
    } else {
      stop("Multiple geometries in input shapefile")
    }  
    
  } else if (tolower(ext) %in% c("tif","tiff")) {
    cat("\n   ... processing raster input\n")
    t1 <- Sys.time()
    raster_inp <- raster(grd_inp_file)
    grd_inp    <- rasterToPolygons(raster_inp) 
    rm(raster_inp)
    cols_to_regrid <- names(grd_inp)
    if (write_inp_raster_as_shp) {
      raster::shapefile(grd_inp,
                        paste0(gsub(".[a-zA-Z]*$","",grd_inp_file),".shp"),
                        overwrite=T)
    }
    grd_inp <- st_as_sf(grd_inp)
    message <- paste0("\n   ... preprocessing raster took ", round(difftime(time1 = Sys.time(), time2 = t1, units = "min"),1)," min")
    write(message, file=out_log, append=TRUE)
    
  } else if (tolower(ext) %in% c("nc")) {
    ncdata <- nc_open(filename = grd_inp_file, write=F)
    latdim <- names(ncdata[["dim"]])[tolower(names(ncdata[["dim"]]))%in%c("lat","latitude" )]
    londim <- names(ncdata[["dim"]])[tolower(names(ncdata[["dim"]]))%in%c("lon","longitude")]
    if (length(latdim)==0 | length(londim)==0) {
      stop("V netCDF nenalezeda dimenze lon nebo lat")
    }
    
    lats <- ncvar_get(ncdata,latdim) %>% round(.,inp_nc_coord_dec)
    sorted   <- sort(lats)
    nSN      <- length(sorted)
    diffs    <- sorted[2:length(sorted)]-sorted[1:length(sorted)-1] 
    lat_step <- min(diffs[diffs > 0])
    lons <- ncvar_get(ncdata,londim) %>% round(.,inp_nc_coord_dec)
    lons[lons>180] <- lons[lons>180]-360 # ensure lons are within -180, 180
    sorted   <- sort(lons)
    nWE      <- length(sorted)
    diffs    <- sorted[2:length(sorted)]-sorted[1:length(sorted)-1] 
    lon_step <- min(diffs[diffs > 0])
    
    vars <- as.data.frame(matrix(nrow = nSN*nWE, 
                                 ncol = 2+length(cols_to_regrid), 
                                 data = NA)) 
    colnames(vars) <- c("lon","lat",cols_to_regrid)
    
    v <- cols_to_regrid[1]
    for (v in cols_to_regrid) {
      nc      <- as.data.frame(ncvar_get(ncdata,v))
      vars[v] <- reshape(data = nc, direction = "long", 
                         varying = list(colnames(nc)), 
                         v.names = v)[v]
    }
    nc_close(ncdata)
    
    for (j in c(1:nSN)){
      vars$lon[(1+(j-1)*nWE) : (j*nWE)] <- lons
      vars$lat[(1+(j-1)*nWE) : (j*nWE)] <- lats[j]
    }
    vars_sf <- vars %>% st_as_sf(.,coords=c("lon","lat"), crs=4326)
    
    e <- as(raster::extent(min(lons)-lon_step/2, max(lons)+lon_step/2,
                           min(lats)-lat_step/2, max(lats)+lat_step/2),
            "SpatialPolygons") %>% st_as_sf()
    grd_poly <- st_make_grid(
      e,
      cellsize = c(lon_step, lat_step),
      offset   = c(min(lons)-lon_step/2, min(lats)-lat_step/2),
      n = c(nWE, nSN),
      crs = 4326,
      what = "polygons",
      square = TRUE,
      flat_topped = FALSE
    ) %>% st_sf()
    
    grd_inp <- st_join(x = grd_poly, y = vars_sf, left = T, join = st_contains)
    rm(grd_poly)
    
    if (write_inp_netcdf_as_shp) {
      st_write(obj = grd_inp, dsn = paste0(gsub(".[a-zA-Z]*$","",grd_inp_file),".shp", delete_dsn = T, delete_layer = T) )
    }
    
  } else {
    stop("unknown file extension")
  }
  
  # transform it into common projection
  grd_inp <- st_transform(x = grd_inp, crs = common_crs)
  
  if (is.null(cols_out_name)) cols_out_name <- cols_to_regrid
  # write to log
  write(paste0("\n   ... variables to regrid: ",paste0(cols_to_regrid,collapse = ", ")), file=out_log, append=TRUE)
  write(paste0(  "   ... their output names : ",paste0(cols_out_name ,collapse = ", ")), file=out_log, append=TRUE)
  
  
  # calculate intersect and its area
  t1 <- Sys.time()
  int <- st_intersection(grd_inp[cols_to_regrid], grd_out[c(grd_out_id)]) %>% 
    mutate(areaInt = st_area(.) %>% as.numeric())
  
  message <- paste0("\n   ... intersecting took ", round(difftime(time1 = Sys.time(), time2 = t1, units = "min"),1)," min")
  write(message, file=out_log, append=TRUE)
  
  
  # calculate area weighted average and merge with output grid
  t1 <- Sys.time()
  tmp     <- as.data.frame(int)
  rm(int)
  areaAgg <- aggregate(formula(paste0(". ~ ",grd_out_id)), 
                       data=tmp[c(grd_out_id,"areaInt")], FUN="sum")
  colnames(areaAgg)[which(colnames(areaAgg)=="areaInt")] <- "areaAgg"
  tmp <- merge(tmp, areaAgg, by=grd_out_id)
  tmp[cols_to_regrid] <- tmp[cols_to_regrid]*tmp$areaInt/tmp$areaAgg
  tmpAgg <- aggregate(formula(paste0(". ~ ",grd_out_id)), 
                      data=tmp[c(grd_out_id,"areaInt",cols_to_regrid)], FUN="sum")
  rm(tmp)
  
  # rename output columns
  colnames(tmpAgg)[which(colnames(tmpAgg)%in%cols_to_regrid)] <- cols_out_name
  
  grd_out <- merge(grd_out, tmpAgg[c(grd_out_id,cols_out_name,"areaInt")], by=grd_out_id, all=T)
  
  # delete values in cell which do not have minimum coverage by the input grid
  grd_out[which(is.na(grd_out$areaInt)),"areaInt"] <- 0
  mask <- (round(grd_out$areaInt/grd_out$areaOrg*100,0) < min_coverage)
  grd_out[which(mask),cols_out_name] <- NA
  grd_out$areaInt <- NULL
  
  
  
  message <- paste0("\n   ... area averaging took ", round(difftime(time1 = Sys.time(), time2 = t1, units = "min"),1)," min")
  write(message, file=out_log, append=TRUE)
  
  message <- paste0("\n   ... TOTAL processing of ", grd_inp_file," took ",round(difftime(time1 = Sys.time(), time2 = tbegf, units = "min"),1)," min")
  write(message, file=out_log, append=TRUE)
  
}

# write output file
if (common_crs!=out_crs) grd_out <- st_transform(x = grd_out, crs = out_crs)

# list.files(gsub(".[a-zA-Z]*$",".*",out_file))
# list.files("ETC_stats_EMEP01_rv4_35_2019met_2018emis_O3_PMx_ON_EEA_10kmgrid_2c.*")
# # file.remove(out_file)
st_write(obj = grd_out, dsn = out_file, delete_dsn = T, delete_layer = T)

tend <- Sys.time()

message <- paste0("\nProcessing took ", round(difftime(time1 = tend, time2 = tbeg, units = "mins"),1)," min")

cat(message)  
write(message, file=out_log, append=TRUE)

