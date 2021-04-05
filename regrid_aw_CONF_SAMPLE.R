# ------------------------
# === OBECNE NASTAVENI ===
# ------------------------

# pracovni adresar
wd <- "d:/work/Documents/CHMU/projekty/ETC-ATNI/2019"

# vystupni adresar, pokud ma byt jiny, nez pracovni
out_dir        <- "pregrid" 

# nazev vystupniho souboru
#    * Pokud je NULL vytvori se automaticky nasledovne:
#        - vstupni soubor je jeden: vystup  se bude jmenovat <nazev vstupniho souboru>_ON_<nazev vystupniho gridu>
#        - vice vstupnich souboru : vystup  se bude jmenovat lot_of_stuff_ON_<nazev vystupniho gridu>
#    * Pokud je nejaky retezec    : vystup  se bude jmenovat <out_file_base>_ON_<nazev vystupniho gridu>
out_file_base  <- "OMI"  #"EMEP01_rv4_35_m19_e18"  #"MACC-RAQ_2020_ENS_FC"

# minimalni pokryti jedne bunky vystupniho gridu vstupnimi daty, aby se spocital prumer
min_coverage   <- 100        

# pocet desetinnych mist, na ktere jsou zaokrouhlovany souradnice stredu gridu, aby se predeslo numerickym chybam
coord_digits   <- 3 

# common crs for intersects - nutne deginovat pouze, pokud je vystupni grid v geograficke projekci
common_crs <- NULL # if NULL, projection of output grid will be used

# ---------------------
# === VYSTUPNI GRID ===
# ---------------------
# 
grd_out_file <- "../gridy/EEA_10kmgrid_2c.shp" # "gridy/EEA_1kmgrid.shp" #"gridy/EEA_10kmgrid_2c.shp"
grd_out_id   <- "CellCode"
# sloupce z puvodniho shapefile, ktere se preulozi do vystupu:
grd_out_cols_to_keep <- NULL #c("EofOrigin","NofOrigin","POINT_X","POINT_Y","LAT","LON")

# grd_out_file <- "../gridy/EEA_2kmgrid.shp"
# grd_out_id   <- "CELL_CODE"
# # sloupce z puvodniho shapefile, ktere se preulozi do vystupu:
# grd_out_cols_to_keep <- NULL #c("POINT_X","POINT_Y","LAT","LON")

# grd_out_file <- "../gridy/EEA_1kmgrid.shp"
# grd_out_id   <- "CellCode"
# # sloupce z puvodniho shapefile, ktere se preulozi do vystupu:
# grd_out_cols_to_keep <- NULL #c("EofOrigin","NofOrigin","POINT_X","POINT_Y","LAT","LON","CODE","COUNTRY","POPULATION","weight_urb","weight_tr")

# --------------
# === VSTUPY ===
# --------------
write_inp_points_as_shp <- F
write_inp_raster_as_shp <- F
write_inp_netcdf_as_shp <- F

# # *** EMEP ***
# inp_files       <- list("modely/EMEP01_rv4_35_m19_e18_O3_PMx_stats.nc")
# colnames_inp    <- list(c("PM10_yAvg","PM10_36Dmax","PM25_yAvg","O3_8HDmx26","AOT40_c","AOT40_f","SOMO35","SOMO10")) # ignored for geotif
# colnames_out    <- list(c("PM10_Y"   ,"PM10_D36"   ,"PM25_Y"   ,"O3_8HDmx26","AOT40_c","AOT40_f","SOMO35","SOMO10"))
# 
# inp_files[2]    <- list("modely/EMEP01_rv4_35_m19_e18_year_pridan_NOx.nc")
# colnames_inp[2] <- list(c("SURF_ug_NOX","SURF_ug_NO2")) # ignored for geotif
# colnames_out[2] <- list(c("NOx_Y"      ,"NO2_Y"   ))
# # 
# inp_files    <- list("modely/EMEP01_rv4_35_m19_e18_year_pridan_NOx.nc")
# colnames_inp <- list(c("SURF_ug_NOX","SURF_ug_NO2")) # ignored for geotif
# colnames_out <- list(c("NOx_Y"      ,"NO2_Y"   ))

# *** MACC-RAQ ***
# inp_files       <- list("modely/MACC-RAQ_2020_ENS_FORECAST_PMx_NO2_O3_stats.nc")
# colnames_inp    <- list(c("PM10_yAvg","PM10_36Dmax","PM25_yAvg","NO2_yAvg","O3_8HDmx26","AOT40_c","AOT40_f","SOMO35","SOMO10")) # ignored for geotif
# colnames_out    <- list(c("PM10_Y"   ,"PM10_D36"   ,"PM25_Y"   ,"NO2_Y"   ,"O3_8HDmx26","AOT40_c","AOT40_f","SOMO35","SOMO10"))

# *** OMI ***
inp_files       <- list("modely/OMI_NO2_2019.tif")
colnames_inp    <- list(c("")) # ignored for geotif
colnames_out    <- list(c("NO2_Y_19"))
# #
inp_files[2]       <- list("modely/OMI_NO2_2020.tif")
colnames_inp[2]    <- list(c("")) # ignored for geotif
colnames_out[2]    <- list(c("NO2_Y_20"))

# *** SENTINEL ***
# inp_files       <- list("modely/s5p_no2_aa_2019_europe.tif")
# colnames_inp    <- list(c("")) # ignored for geotif
# colnames_out    <- list(c("NO2_Y"))
# 
# inp_files[2]       <- list("modely/s5p_no2_aa_2020_europe.tif")
# colnames_inp[2]    <- list(c("")) # ignored for geotif
# colnames_out[2]    <- list(c("NO2_Y"))


