##############################################################################
#########   load packages
##############################################################################

rm(list = ls())
# ENTER COUNTRY OF INTEREST -----------------------------------------------
# Please capitalize the first letter of the country name and replace " " in the country name to "_" if there is.
country <- 'Malawi'
survey_year <- '2015'
survey <- paste0(country,'_',survey_year)

# Load libraries and info ----------------------------------------------------------
library(dplyr)
library(sf)
library(terra)
library(rdhs)


#library(labelled)
#library(haven)
#library(surveyPrev)

################################################################
#########   load survey Info
################################################################
## set directory

# path to home directory
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

home.dir <- paste(code.path.splitted[1: (length(code.path.splitted)-3)], collapse = "/")

data.dir<-paste0(home.dir,'/Data/',survey)

setwd(paste(data.dir))


info.name <- paste0(survey, "_general_info.Rdata")

load(file = paste0(info.name, sep=''))

##############################################################################
#########   load helper functions 
##############################################################################

setwd(paste(code.dir))
source('DataProcessing_helper.R')


##############################################################################
#########   create folders for prepared data
##############################################################################

setwd(paste(data.dir))

# create a folder to store prepared data
#if(!dir.exists(paths = paste0('xxx'))){ 
#  dir.create(path = paste0('xxx'))
#}


##############################################################################
###### RDHS Configration 
##############################################################################


## login
# set API to get DHS data -- you will need to change this to your information!
#rdhs::set_rdhs_config(email = "xxxxx",
#                      project = "xxxxxxxxxxxxxx")


##############################################################################
#########   load survey meta data
##############################################################################

setwd(paste0(home.dir,'/Data'))

if(file.exists('DHS_meta.rda')){
  load('DHS_meta.rda')
  
}else{
  
  ### run for the first time use saved object later
  DHS.country.meta <- rdhs::dhs_countries()
  DHS.survey.meta <- rdhs::dhs_surveys()
  DHS.dataset.meta <- rdhs::dhs_datasets()
  
  save(DHS.country.meta,
       DHS.survey.meta,
       DHS.dataset.meta,
       file='DHS_meta.rda')
  
}


##############################################################################
#########   process polygon files
##############################################################################

setwd(paste0(data.dir,'/shapeFiles_gadm'))

if(file.exists('country_shp_analysis.rds') && file.exists('poly_adm0.rds') &&
   file.exists('poly_adm_strata.rds') && file.exists('strata_adm_info.rds')){
  
  country_shp_analysis <- readRDS('country_shp_analysis.rds')
  poly.adm0 <- readRDS('poly_adm0.rds') ### national boundary 
  poly.adm.strata <- readRDS('poly_adm_strata.rds') ### strata admin level boundary
  strata_admin_info <- readRDS('strata_adm_info.rds') ### strata admin level info
  
  
}else{
  
  ### run for the first time and use stored results later
  
  ### automatically download GADM shapefile
  ### store as a list, with names 'National','Admin-1',etc...
  country_shp_analysis <- get_country_GADM(country=country,resolution=1)
  country_shp_analysis <- lapply(country_shp_analysis, function(x) {
    sf::st_set_crs(x, 4326)
  })
  
  
  ### prepare admin info for the stratification level
  poly.adm.strata <- country_shp_analysis[[paste0('Admin-',strat_level_GADM)]]
  poly.adm.strata$strata.adm.num <- c(1:dim(poly.adm.strata)[1])
  poly.adm.strata$strata.adm.char <- paste0("strata_adm_", 1:dim(poly.adm.strata)[1])
  
  strata_admin_info <- surveyPrev::adminInfo(
    poly.adm =poly.adm.strata,
    admin = 1,by.adm=paste0("NAME_",strat_level_GADM)) 
  
  strata_admin_info$data$strata.adm.char <- paste0("strata_adm_", 1:dim(strata_admin_info$data)[1])
  strata_admin_info$data$strata.adm.num <- c(1:dim(strata_admin_info$data)[1])
  
  poly.adm0 <- country_shp_analysis[['National']]
  
  
  ### save results
  setwd(paste0(data.dir,'/shapeFiles_gadm'))
  
  saveRDS(country_shp_analysis,'country_shp_analysis.rds')
  saveRDS(poly.adm0,'poly_adm0.rds')
  saveRDS(poly.adm.strata,'poly_adm_strata.rds')
  saveRDS(strata_admin_info,'strata_adm_info.rds')
  
}

##############################################################################
#########   process worldpop data for the year of census and calibrated year
##############################################################################


setwd(paste0(data.dir,'/worldpop'))

# census year
### download population
pop.abbrev <- tolower(gadm.abbrev)
pop_file <- paste0(pop.abbrev,'_ppp_',frame_year,'_1km_Aggregated_UNadj.tif')


if(!file.exists(pop_file)){
  
  url <- paste0("https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj/", 
                frame_year, "/", toupper(pop.abbrev),"/",      
                pop.abbrev,'_ppp_',frame_year,'_1km_Aggregated_UNadj.tif')
  
  download.file(url, pop_file, method = "libcurl",mode="wb")
}

### load population
worldpop <- terra::rast(paste0(country.abbrev, '_ppp_', frame_year, '_1km_Aggregated_UNadj.tif'))




# calibrated year
### download population
pop.abbrev <- tolower(gadm.abbrev)
pop_file <- paste0(pop.abbrev,'_ppp_',calibrated_year,'_1km_Aggregated_UNadj.tif')


if(!file.exists(pop_file)){
  
  url <- paste0("https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj/", 
                calibrated_year, "/", toupper(pop.abbrev),"/",      
                pop.abbrev,'_ppp_',calibrated_year,'_1km_Aggregated_UNadj.tif')
  
  download.file(url, pop_file, method = "libcurl",mode="wb")
}

### load population
calibrated_worldpop <- terra::rast(paste0(country.abbrev, '_ppp_', calibrated_year, '_1km_Aggregated_UNadj.tif'))


##############################################################################
#########   process nighttime light data for the year of census
##############################################################################


setwd(paste0(data.dir,'/prepared_dat'))


### load night time light

country_light_file <- paste0(gadm.abbrev,'_night_time_light.tif')

if(file.exists(country_light_file)){
  
  light_ras_buffered <- terra::rast(country_light_file)
  
}else{
  
  setwd(paste0(home.dir,'/Data/Global_rater/night_time_light/'))
  
  
  night_file <- paste0('night_light_',frame_year,'.tif')
  
  light_ras <- terra::rast(night_file)
  
  buffered_adm0 <- sf::st_buffer(poly.adm0, dist = 0.25)
  # Malawi nightime light
  light_ras_buffered <- terra::mask(terra::crop(light_ras, buffered_adm0), buffered_adm0)
  
  setwd(paste0(data.dir,'/prepared_dat'))
  terra::writeRaster(light_ras_buffered,country_light_file)
  
}



#plot(light_ras_buffered)

##############################################################################
#########   set up national grid for prediction
##############################################################################

# national grid without travel time
setwd(paste0(data.dir,'/prepared_dat'))

if(file.exists('natl_grid.rds')){
  
  natl_grid_df <- readRDS('natl_grid.rds')
  
}else{
  ### set up national prediction grid
  # coordinate for pixel pop density, and it's value
  coords <- as.data.frame(terra::xyFromCell(worldpop, 1:ncell(worldpop)))
  pop_density <- terra::values(worldpop)
  pop_density[is.nan(pop_density)] <- NA
  
  # Merge coordinate and population density data
  worldpop_data <- cbind(coords, pop_density)
  colnames(worldpop_data) <- c("LONGNUM", "LATNUM", "Population_Density")
  
  # calibrated pop
  calibrated_pop <- terra::extract(calibrated_worldpop, worldpop_data[c('LONGNUM','LATNUM')])[,2]
  
  # Merge night time light 
  light_values <- terra::extract(light_ras_buffered, worldpop_data[c('LONGNUM','LATNUM')])[,2]
  
  # Match stratification admin
  natl_grid_sf <- st_as_sf(worldpop_data[c('LONGNUM','LATNUM')], coords = c("LONGNUM", "LATNUM"), crs = st_crs(poly.adm.strata))
  
  strata_adm_match <- st_join(natl_grid_sf,poly.adm.strata, join = st_intersects)
  
  # Prepare national grid
  natl_grid_df <- data.frame(
    x = worldpop_data$LONGNUM,
    y = worldpop_data$LATNUM,
    pop_den = worldpop_data$Population_Density,
    calibrated_pop_den = calibrated_pop,
    light = light_values,
    strata.adm.num = strata_adm_match$strata.adm.num,
    strata.adm.char = strata_adm_match$strata.adm.char,
    strata.adm.name = strata_adm_match[[paste0("NAME_",strat_level_GADM)]]
  )
  
  natl_grid_df$pixel_id <- c(1:dim(natl_grid_df)[1])
  
  setwd(paste0(data.dir,'/prepared_dat'))
  saveRDS(natl_grid_df,file='natl_grid.rds')
}


##############################################################################
#########   process DHS surveys (cluster) training data, jittering
##############################################################################
# training data with pop + light, jittering
setwd(paste0(data.dir,'/prepared_dat'))

if(file.exists('train_cluster_dat_cr.rds')){
  
  train_dat <- readRDS('train_cluster_dat_cr.rds')
  
}else{
  
  train_dat <- data.frame()
  
  # iterate through surveys based on the same frame
  for( k in 1:length(svy_training_year)){
    
    # get cluster grid
    tmp.geo <- surveyPrev::getDHSgeo(country = country, year = svy_training_year[k])
    tmp.cluster.info <- surveyPrev::clusterInfo(geo=tmp.geo, 
                                                poly.adm1=poly.adm.strata,
                                                poly.adm2=poly.adm.strata, 
                                                by.adm1 = paste0("NAME_",strat_level_GADM),
                                                by.adm2 = paste0("NAME_",strat_level_GADM))
    
    tmp_cluster_grid <- tmp.cluster.info$data
    
    # # get pop and light
    # tmp_cluster_grid$pop_den <- terra::extract(worldpop, 
    #                                            tmp_cluster_grid[c('LONGNUM','LATNUM')])[,2]
    # tmp_cluster_grid$light <- terra::extract(light_ras_buffered, 
    #                                          tmp_cluster_grid[c('LONGNUM','LATNUM')])[,2]
    
    # get U/R 
    tmp_cluster_grid <- merge(tmp_cluster_grid, as.data.frame(tmp.geo)[,c('DHSCLUST','URBAN_RURA')],
                              by.x='cluster',by.y='DHSCLUST',all.x=T)
    
    tmp_cluster_grid <- tmp_cluster_grid %>%
      mutate(urban = standardize_urban(URBAN_RURA))
    
    ##########
    urban_clus<-tmp_cluster_grid[tmp_cluster_grid$urban=='urban',]
    rural_clus<-tmp_cluster_grid[tmp_cluster_grid$urban=='rural',]
    
    urban_clus$x_adj<-NA
    urban_clus$y_adj<-NA
    rural_clus$x_adj<-rural_clus$LONGNUM
    rural_clus$y_adj<-rural_clus$LATNUM
    
    
    library(sqldf)
    # points have to stay within the same admin2 region 
    for( i in 1:dim(urban_clus)[1]){
      
      print(i)
      temp_frame<-constr_prior(urban_clus[i,],2000,1,poly.adm.strata,worldpop, strat_level_GADM)
      p_mode = sqldf("SELECT * FROM temp_frame GROUP BY x,y ORDER BY SUM(unn_w) DESC LIMIT 1")
      urban_clus[i,]$x_adj<-p_mode$x
      urban_clus[i,]$y_adj<-p_mode$y
      
    }
    
    prep_dat<-rbind(urban_clus,rural_clus)
    #xy <- as.matrix(prep_dat[c('x_adj','y_adj')])
    
    
    crc_dat<-prep_dat
    
    # set up corrected xy for clusters
    xy_crc <- as.matrix(crc_dat[c('x_adj','y_adj')])
    crc_dat$x<-crc_dat$x_adj # x_adj and y_adj: corrected coordinates
    crc_dat$y<-crc_dat$y_adj
    crc_dat$LONGNUM<-crc_dat$x_adj # x_adj and y_adj: corrected coordinates
    crc_dat$LATNUM<-crc_dat$y_adj
    
    
    # extract covariates
    crc_dat$pop_den<-as.vector(terra::extract(worldpop, xy_crc)[,1])
    crc_dat$light <- as.vector(terra::extract(light_ras_buffered, xy_crc)[,1])
    
    # Match stratification admin
    cluster_grid_sf <- st_as_sf(crc_dat[c('LONGNUM','LATNUM')], 
                                coords = c("LONGNUM", "LATNUM"),
                                crs = st_crs(poly.adm.strata))
    
    cluster_strata_adm_match <- st_join(cluster_grid_sf,poly.adm.strata, join = st_intersects)
    
    crc_dat$strata.adm.num = cluster_strata_adm_match$strata.adm.num
    crc_dat$strata.adm.char = cluster_strata_adm_match$strata.adm.char
    crc_dat$strata.adm.name = cluster_strata_adm_match[[paste0("NAME_",strat_level_GADM)]]
    
    tmp_cluster_grid<- crc_dat[,c('x','y','pop_den','light','urban',
                                  'strata.adm.num','strata.adm.char','strata.adm.name')]
    
    tmp_cluster_grid$survey <- svy_training_year[k]
    
    train_dat <- rbind(train_dat,tmp_cluster_grid)
  }
  
  setwd(paste0(data.dir,'/prepared_dat'))
  saveRDS(train_dat, file='train_cluster_dat_cr.rds')
  
}



##############################################################################
#########   next step
##############################################################################

### go to next file
setwd(paste(code.dir))

rstudioapi::navigateToFile("step2_UR_surface.R")
