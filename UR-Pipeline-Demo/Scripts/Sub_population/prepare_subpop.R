##############################################################################
#########   load packages
##############################################################################

rm(list = ls())
# ENTER COUNTRY OF INTEREST -----------------------------------------------
# Please capitalize the first letter of the country name and replace " " in the country name to "_" if there is.
country <- 'Malawi'
gadm.abbrev <- 'mwi'
survey_year <- '2015'   ### if the survey spans two years, use the first year
survey <- paste0(country,'_',survey_year)

# Load libraries and info ----------------------------------------------------------
library(dplyr)
library(sf)
library(terra)
library(rdhs)
library(caret)

################################################################
#########   load survey Info
################################################################
## set directory

# path to home directory
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

home.dir <- paste(code.path.splitted[1: (length(code.path.splitted)-3)], collapse = "/")

data.dir<-paste0(home.dir,'/Data/',survey)


##############################################################################
#########   download sub population at survey year
##############################################################################
setwd(paste0(home.dir,'/Data/Global_rater/subpop'))
pop.abbrev <- tolower(gadm.abbrev)

if (survey_year < 2021){
  ages = c(0, 1, seq(5, 75, by = 5))
  
  for (age in ages) {
    ### download sub pop
    female_pop_file <- paste0(pop.abbrev,'_f_',age,'_',survey_year,'_1km.tif')
    if(!file.exists(female_pop_file)){
      
      if (survey_year < 2000) {
        url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020_1km/unconstrained/", 
                      2000, "/", toupper(pop.abbrev),"/",      
                      pop.abbrev,'_f_',age, '_', 2000,'_1km.tif')
      }else{
        url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020_1km/unconstrained/", 
                      survey_year, "/", toupper(pop.abbrev),"/",      
                      pop.abbrev,'_f_',age, '_', survey_year,'_1km.tif')
      }
      download.file(url, female_pop_file, method = "libcurl",mode="wb")
    }
    
    male_pop_file <- paste0(pop.abbrev,'_m_',age,'_',survey_year,'_1km.tif')
    if(!file.exists(male_pop_file)){
      
      if (survey_year < 2000) {
        url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020_1km/unconstrained/", 
                      2000, "/", toupper(pop.abbrev),"/",      
                      pop.abbrev,'_m_',age, '_', 2000,'_1km.tif')
      }else{
        url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020_1km/unconstrained/", 
                      survey_year, "/", toupper(pop.abbrev),"/",      
                      pop.abbrev,'_m_',age, '_', survey_year,'_1km.tif')
      }
      
      download.file(url, male_pop_file, method = "libcurl",mode="wb")
    }
  }
}else{
  
  # age 0-1 female
  female_pop_file <- paste0(pop.abbrev,'_f_',0,'_',survey_year,'_1km.tif')
  if(!file.exists(female_pop_file)){
    if (survey_year >= 2021 && survey_year <= 2022){
      url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2021_2022_1km_UNadj/unconstrained/", 
                  survey_year, "/single_age/", toupper(pop.abbrev),"/",      
                  pop.abbrev,'_f_',0, '_', survey_year,'_1km_UNadj.tif')
    }else{
      url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2021_2022_1km_UNadj/unconstrained/", 
                    2022, "/single_age/", toupper(pop.abbrev),"/",      
                    pop.abbrev,'_f_',0, '_', survey_year,'_1km_UNadj.tif')
    }
    download.file(url, female_pop_file, method = "libcurl",mode="wb")
  }
  
  # age 0-1 male
  male_pop_file <- paste0(pop.abbrev,'_m_',0,'_',survey_year,'_1km.tif')
  if(!file.exists(male_pop_file)){
    if (survey_year >= 2021 && survey_year <= 2022){
     url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2021_2022_1km_UNadj/unconstrained/", 
                  survey_year, "/single_age/", toupper(pop.abbrev),"/",      
                  pop.abbrev,'_m_',0, '_', survey_year,'_1km_UNadj.tif')
    }else{
      url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2021_2022_1km_UNadj/unconstrained/", 
                    2022, "/single_age/", toupper(pop.abbrev),"/",      
                    pop.abbrev,'_m_',0, '_', survey_year,'_1km_UNadj.tif')
    }
    download.file(url, male_pop_file, method = "libcurl",mode="wb")
    
  }
  
  
  
  
  # other ages
  ages = seq(5, 75, by = 5)
  for (age in ages) {
    ### download sub pop
    female_pop_file <- paste0(pop.abbrev,'_f_',age,'_',survey_year,'_1km.tif')
    if(!file.exists(female_pop_file)){
      
      if (survey_year >= 2021 && survey_year <= 2022){
        url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2021_2022_1km_UNadj/unconstrained/", 
                      survey_year, "/five_year_age_groups/", toupper(pop.abbrev),"/",      
                      pop.abbrev,'_f_',age, '_',(age+4) ,'_',survey_year,'_1km_UNadj.tif')
      }else{
        url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2021_2022_1km_UNadj/unconstrained/", 
                      2022, "/five_year_age_groups/", toupper(pop.abbrev),"/",      
                      pop.abbrev,'_f_',age, '_',(age+4) ,'_',2022,'_1km_UNadj.tif')
      }
      download.file(url, female_pop_file, method = "libcurl",mode="wb")
    }
    
    male_pop_file <- paste0(pop.abbrev,'_m_',age,'_',survey_year,'_1km.tif')
    if(!file.exists(male_pop_file)){
      
      if (survey_year >= 2021 && survey_year <= 2022){
        url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2021_2022_1km_UNadj/unconstrained/", 
                      survey_year, "/five_year_age_groups/", toupper(pop.abbrev),"/",      
                      pop.abbrev,'_m_',age, '_',(age+4) ,'_',survey_year,'_1km_UNadj.tif')
      }else{
        url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2021_2022_1km_UNadj/unconstrained/", 
                      2022, "/five_year_age_groups/", toupper(pop.abbrev),"/",      
                      pop.abbrev,'_m_',age, '_',(age+4) ,'_',2022,'_1km_UNadj.tif')
      }
      
      download.file(url, male_pop_file, method = "libcurl",mode="wb")
    }
  }
}

##############################################################################
#########   aggregated child under 1 at survey year
##############################################################################
setwd(paste0(data.dir,'/worldpop'))
kid_0_1_filename <- paste0(pop.abbrev, "_k0_1_", survey_year, "_1km.tif")

if(!file.exists(kid_0_1_filename)) {
  setwd(paste0(home.dir,'/Data/Global_rater/subpop'))
  age = 0
  female_pop_file <- paste0(pop.abbrev,'_f_',age,'_',survey_year,'_1km.tif')
  male_pop_file <- paste0(pop.abbrev,'_m_',age,'_',survey_year,'_1km.tif')
  
  # Read the rasters
  female_raster <- terra::rast(female_pop_file)
  male_raster <- terra::rast(male_pop_file)
  
  kid_0_1_raster <- female_raster + male_raster
  
  setwd(paste0(data.dir,'/worldpop'))
  terra::writeRaster(kid_0_1_raster, kid_0_1_filename)
}


##############################################################################
#########   aggregated child under 5 at survey year
##############################################################################
setwd(paste0(data.dir,'/worldpop'))
kid_0_5_filename <- paste0(pop.abbrev, "_k0_5_", survey_year, "_1km.tif")

if (survey_year < 2021) {
  if(!file.exists(kid_0_5_filename)) {
    setwd(paste0(home.dir,'/Data/Global_rater/subpop'))
    ages = c(0,1)
    kid_0_5_raster <- NULL
    for (age in ages){
      female_pop_file <- paste0(pop.abbrev,'_f_',age,'_',survey_year,'_1km.tif')
      male_pop_file <- paste0(pop.abbrev,'_m_',age,'_',survey_year,'_1km.tif')
      # Read the rasters
      female_raster <- terra::rast(female_pop_file)
      male_raster <- terra::rast(male_pop_file)
      
      kid_raster <- female_raster + male_raster
      
      if (is.null(kid_0_5_raster)) {
        kid_0_5_raster <- kid_raster
      } else {
        kid_0_5_raster <- kid_0_5_raster + kid_raster
      }
    }
    setwd(paste0(data.dir,'/worldpop'))
    terra::writeRaster(kid_0_5_raster, kid_0_5_filename)
  }
} else{
  
  if(!file.exists(kid_0_5_filename)) {
    setwd(paste0(home.dir,'/Data/Global_rater/subpop'))
    female_pop_file <- paste0(pop.abbrev,'_f_',0,'_',4,'_',survey_year,'_1km.tif')
    if(!file.exists(female_pop_file)){
      if (survey_year >= 2021 && survey_year <= 2022){
        url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2021_2022_1km_UNadj/unconstrained/", 
                      survey_year, "/five_year_age_groups/", toupper(pop.abbrev),"/",      
                      pop.abbrev,'_f_',0, '_',4 ,'_',survey_year,'_1km_UNadj.tif')
      }else{
        url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2021_2022_1km_UNadj/unconstrained/", 
                      2022, "/five_year_age_groups/", toupper(pop.abbrev),"/",      
                      pop.abbrev,'_f_',0, '_',4 ,'_',2022,'_1km_UNadj.tif')
      }
      download.file(url, female_pop_file, method = "libcurl",mode="wb")
    }
    male_pop_file <- paste0(pop.abbrev,'_m_',0,'_',4,'_',survey_year,'_1km.tif')
    if(!file.exists(male_pop_file)){
      if (survey_year >= 2021 && survey_year <= 2022){
        url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2021_2022_1km_UNadj/unconstrained/", 
                      survey_year, "/five_year_age_groups/", toupper(pop.abbrev),"/",      
                      pop.abbrev,'_m_',0, '_',4 ,'_',survey_year,'_1km_UNadj.tif')
      }else{
        url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2021_2022_1km_UNadj/unconstrained/", 
                      2022, "/five_year_age_groups/", toupper(pop.abbrev),"/",      
                      pop.abbrev,'_m_',0, '_',4 ,'_',2022,'_1km_UNadj.tif')
      }
      
      download.file(url, male_pop_file, method = "libcurl",mode="wb")
    }
    
    
    # Read the rasters
    female_raster <- terra::rast(female_pop_file)
    male_raster <- terra::rast(male_pop_file)
    
    kid_0_5_raster <- female_raster + male_raster
    
    setwd(paste0(data.dir,'/worldpop'))
    terra::writeRaster(kid_0_5_raster, kid_0_5_filename)
  }
}

##############################################################################
#########   aggregated female 15_49 at survey year
##############################################################################
setwd(paste0(data.dir,'/worldpop'))
female_15_49_filename <- paste0(pop.abbrev, "_f15_49_", survey_year, "_1km.tif")

if(!file.exists(female_15_49_filename)) {
  setwd(paste0(home.dir,'/Data/Global_rater/subpop'))
  ages = seq(15, 45, by = 5)
  female_15_49_raster <- NULL
  for (age in ages){
    female_pop_file <- paste0(pop.abbrev,'_f_',age,'_',survey_year,'_1km.tif')
    female_raster <- terra::rast(female_pop_file)
    
    if (is.null(female_15_49_raster)) {
      female_15_49_raster <- female_raster
    } else {
      female_15_49_raster <- female_15_49_raster + female_raster
    }
  }
  setwd(paste0(data.dir,'/worldpop'))
  terra::writeRaster(female_15_49_raster, female_15_49_filename)
}


##############################################################################
#########   aggregated male 15_49 at survey year
##############################################################################
setwd(paste0(data.dir,'/worldpop'))
male_15_49_filename <- paste0(pop.abbrev, "_m15_49_", survey_year, "_1km.tif")

if(!file.exists(male_15_49_filename)) {
  setwd(paste0(home.dir,'/Data/Global_rater/subpop'))
  ages = seq(15, 45, by = 5)
  male_15_49_raster <- NULL
  for (age in ages){
    male_pop_file <- paste0(pop.abbrev,'_m_',age,'_',survey_year,'_1km.tif')
    male_raster <- terra::rast(male_pop_file)
    
    if (is.null(male_15_49_raster)) {
      male_15_49_raster <- male_raster
    } else {
      male_15_49_raster <- male_15_49_raster + male_raster
    }
  }
  setwd(paste0(data.dir,'/worldpop'))
  terra::writeRaster(male_15_49_raster, male_15_49_filename)
}

##############################################################################
#########   prepare total population at survey year
##############################################################################
setwd(paste0(data.dir,'/worldpop'))
pop.abbrev <- tolower(gadm.abbrev)

### download total pop
total_pop_file <- paste0(pop.abbrev,'_ppp_',survey_year,'_1km_Aggregated_UNadj.tif')
if(!file.exists(total_pop_file)){
  
  if (survey_year< 2000){
    url <- paste0("https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj/", 
                  2000, "/", toupper(pop.abbrev),"/",      
                  pop.abbrev,'_ppp_',2000,'_1km_Aggregated_UNadj.tif')
  }else if (survey_year >= 2000 && survey_year <= 2020){
    url <- paste0("https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj/", 
                  survey_year, "/", toupper(pop.abbrev),"/",      
                  pop.abbrev,'_ppp_',survey_year,'_1km_Aggregated_UNadj.tif')
  } else if (survey_year >= 2021 && survey_year <= 2022) {
    url <- paste0("https://data.worldpop.org/GIS/Population/Global_2021_2022_1km_UNadj/unconstrained/", 
                  survey_year, "/", toupper(pop.abbrev),"/",      
                  pop.abbrev,'_ppp_',survey_year,'_1km_UNadj.tif')
  }else{
    url <- paste0("https://data.worldpop.org/GIS/Population/Global_2021_2022_1km_UNadj/unconstrained/", 
                  2022, "/", toupper(pop.abbrev),"/",      
                  pop.abbrev,'_ppp_',2022,'_1km_UNadj.tif')
  }
  
  download.file(url, total_pop_file, method = "libcurl",mode="wb")
}

