
wd <- "d:/work/Documents/CHMU/projekty/ETC-ATNI/2019"

# grd_out_file <- "../gridy/EEA_10kmgrid_2c.shp" # "gridy/EEA_1kmgrid.shp" #"gridy/EEA_10kmgrid_2c.shp"
# grd_out_id   <- "CellCode"
# grd_out_cols_to_keep <- c("EofOrigin","NofOrigin","POINT_X","POINT_Y","LAT","LON")

grd_out_file <- "../gridy/EEA_2kmgrid.shp"
grd_out_id   <- "CELL_CODE"
grd_out_cols_to_keep <- NULL #c("POINT_X","POINT_Y","LAT","LON")

# grd_out_file <- "../gridy/EEA_1kmgrid.shp"
# grd_out_id   <- "CellCode"
# grd_out_cols_to_keep <- NULL #c("EofOrigin","NofOrigin","POINT_X","POINT_Y","LAT","LON","CODE","COUNTRY","POPULATION","weight_urb","weight_tr")

min_coverage <- 100
out_dir   <- "pregrid/"
out_file  <- NULL # if null name will be derived from output file name

# input file(s)
inp_nc_coord_dec <- 3
write_inp_raster_as_shp <- T
write_inp_netcdf_as_shp <- F

inp_files       <- list("modely/EMEP01_rv4_35_m19_e18_O3_PMx_stats.nc")
colnames_inp    <- list(c("PM10_yAvg","PM10_36Dmax","PM25_yAvg","O3_8HDmx26","AOT40_c","AOT40_f","SOMO35","SOMO10")) # ignored for geotif
colnames_out    <- list(c("PM10_Y"   ,"PM10_D36"   ,"PM25_Y"   ,"O3_8HDmx26","AOT40_c","AOT40_f","SOMO35","SOMO10"))

# inp_files[2]       <- list("modely/MACC-RAQ_2020_ENS_FORECAST_PMx_NO2_O3.nc")
# colnames_inp[2]    <- list(c("PM10_yAvg","PM10_36Dmax","PM25_yAvg","NO2_yAvg","O3_8HDmx26","AOT40_c","AOT40_f","SOMO35","SOMO10")) # ignored for geotif
# colnames_out[2]    <- list(c("PM10_Y"   ,"PM10_D36"   ,"PM25_Y"   ,"NO2_yAvg","O3_8HDmx26","AOT40_c","AOT40_f","SOMO35","SOMO10"))


# common crs for intersects
common_crs <- NULL # if NULL, projection of output grid will be used
