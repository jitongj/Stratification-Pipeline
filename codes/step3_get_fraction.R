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

setwd(paste(data.dir))


info.name <- paste0(survey, "_general_info.Rdata")

load(file = paste0(info.name, sep=''))

##############################################################################
#########   create folders to store results
##############################################################################

setwd(paste(res.dir))

if(!dir.exists(paths = paste0(res.dir,'/UR_Fraction'))){
  dir.create(path = paste0(res.dir,'/UR_Fraction'))
}

##############################################################################
#########   load helper functions 
##############################################################################

setwd(paste(code.dir))
source('DataProcessing_helper.R')


##############################################################################
#########   load prepared data 
##############################################################################

### load admin info

setwd(paste0(data.dir,'/shapeFiles_gadm'))

country_shp_analysis <- readRDS('country_shp_analysis.rds')
poly.adm0 <- readRDS('poly_adm0.rds') ### national boundary 
#poly.adm.strata <- readRDS('poly_adm_strata.rds') ### strata admin level boundary
#strata_admin_info <- readRDS('strata_adm_info.rds') ### strata admin level info


### load worldpop
setwd(paste0(data.dir,'/worldpop'))
pop.abbrev <- tolower(gadm.abbrev)
pop_file <- paste0(pop.abbrev,'_ppp_',frame_year,'_1km_Aggregated_UNadj.tif')
worldpop <- terra::rast(paste0(country.abbrev, '_ppp_', frame_year, '_1km_Aggregated_UNadj.tif'))


### load indicator surface 
setwd(paste0(res.dir,'/UR_Classification'))
urb_ind_ras <-terra::rast(paste0(country.abbrev,'_','ind_prob.tif'))
urb_prob_ras = urb_ind_ras


# load known fractions for checking
setwd(paste0(data.dir))
ref.tab <- readRDS(paste0(country.abbrev,'_ref_tab.rds'))

##############################################################################
#########   sanity check with sampling frame
##############################################################################

### visualize predicted probability surface
plot(urb_prob_ras)

### produce fraction for total population at strata level
pixel_grid <- as.data.frame(terra::xyFromCell(urb_prob_ras, 1:ncell(urb_prob_ras)))


admin.level.tmp = strat_level_GADM

pixel_adm_grid_tmp <- get_pixel_adm_grid(pixel_grid = pixel_grid,
                                         admin.level = admin.level.tmp,
                                         admin.poly = country_shp_analysis[[admin.level.tmp+1]])


adm_urb_frac_tmp <- get_adm_urb_frac(pixel_adm_grid=pixel_adm_grid_tmp,
                                     pred_surf=urb_prob_ras,
                                     pop_surf=worldpop)

### compare with known truth

matched.order <- match(split_underscore(adm_urb_frac_tmp$adm.name.full),
                       ref.tab$strata.adm.name)

par(pty="s")
plot(ref.tab$urb_frac[matched.order],adm_urb_frac_tmp$urb_frac, col = "blue",
     xlab='sampling_frame_frac',ylab='estimated_frac')
abline(a = 0, b = 1)




##############################################################################
#########   prepare total population and sub population at survey year
##############################################################################
setwd(paste0(data.dir,'/worldpop'))
pop.abbrev <- tolower(gadm.abbrev)

### load total population
total_pop <- terra::rast(paste0(country.abbrev, '_ppp_', survey_year, '_1km_Aggregated_UNadj.tif'))

### load female 15-49 population
f15_49 <- terra::rast(paste0(country.abbrev, '_f15_49_', survey_year, '_1km.tif'))

### load child 0-1 population
k0_1 <- terra::rast(paste0(country.abbrev, '_k0_1_', survey_year, '_1km.tif'))

### load child 0-5 population
k0_5 <- terra::rast(paste0(country.abbrev, '_k0_5_', survey_year, '_1km.tif'))

### load male 15-49 population
m15_49 <- terra::rast(paste0(country.abbrev, '_m15_49_', survey_year, '_1km.tif'))


##############################################################################
#########   create folders for fraction of subpopulation 
##############################################################################
setwd(paste0(res.dir, '/UR_fraction'))

admin.num = length(country_shp_analysis)-1 # exclude the national level

for (admin.level.tmp in 0:admin.num) {
  dir_name <- paste0("Admin_", admin.level.tmp)
  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
  }
}


##############################################################################
#########   fraction for all subpopulation 
##############################################################################

for (admin.level.tmp in 0:admin.num) {
  pixel_adm_grid_tmp <- get_pixel_adm_grid(pixel_grid = pixel_grid,
                                           admin.level = admin.level.tmp,
                                           admin.poly = country_shp_analysis[[admin.level.tmp+1]])
  
  adm_urb_frac_total <- get_adm_urb_frac(pixel_adm_grid=pixel_adm_grid_tmp,
                                         pred_surf=urb_prob_ras,
                                         pop_surf=total_pop)
  
  adm_urb_frac_k0_1 <- get_adm_urb_frac(pixel_adm_grid=pixel_adm_grid_tmp,
                                        pred_surf=urb_prob_ras,
                                        pop_surf=k0_1)
  
  adm_urb_frac_k0_5 <- get_adm_urb_frac(pixel_adm_grid=pixel_adm_grid_tmp,
                                        pred_surf=urb_prob_ras,
                                        pop_surf=k0_5)
  
  adm_urb_frac_f15_49 <- get_adm_urb_frac(pixel_adm_grid=pixel_adm_grid_tmp,
                                          pred_surf=urb_prob_ras,
                                          pop_surf=f15_49)
  
  adm_urb_frac_m15_49 <- get_adm_urb_frac(pixel_adm_grid=pixel_adm_grid_tmp,
                                          pred_surf=urb_prob_ras,
                                          pop_surf=m15_49)
  
  ### save fractions under for this admin
  setwd(paste0(res.dir,'/UR_fraction/Admin_',admin.level.tmp))
  
  paste0(country.abbrev,'_',survey_year,'_admin',admin.level.tmp,'_','total_pop_urb_frac.rds')
  
  saveRDS(adm_urb_frac_total,paste0(country.abbrev,'_',survey_year,'_admin',admin.level.tmp,'_','total_pop_urb_frac.rds'))
  saveRDS(adm_urb_frac_k0_1,paste0(country.abbrev,'_',survey_year,'_admin',admin.level.tmp,'_','k0_1_urb_frac.rds'))
  saveRDS(adm_urb_frac_k0_5,paste0(country.abbrev,'_',survey_year,'_admin',admin.level.tmp,'_','k0_5_urb_frac.rds'))
  saveRDS(adm_urb_frac_f15_49,paste0(country.abbrev,'_',survey_year,'_admin',admin.level.tmp,'_','f15_49_urb_frac.rds'))
  saveRDS(adm_urb_frac_m15_49,paste0(country.abbrev,'_',survey_year,'_admin',admin.level.tmp,'_','m15_49_urb_frac.rds'))
}



##############################################################################
#########   additional checking scripts example, not necessarily useful
##############################################################################
### checking certain areas, no population within 

if(FALSE){
  all_adm_full_name_tmp <- paste0( country_shp_analysis[[admin.level.tmp+1]][[paste0('NAME_',admin.level.tmp-1)]],
                                   '_',
                                   country_shp_analysis[[admin.level.tmp+1]][[paste0('NAME_',admin.level.tmp)]])
  
  all_adm_full_name_tmp[which(!all_adm_full_name_tmp %in%adm_urb_frac_tmp$adm.name.full)]
  
  
  tmp_adm2 <- country_shp_analysis[[3]]
  tmp_adm2_98 <- tmp_adm2[98,]
  
  pol_vect <- vect(tmp_adm2_98)
  
  # Crop and mask
  r_crop <- crop(worldpop, pol_vect)
  r_mask <- mask(r_crop, pol_vect)
  
  # Plotting
  plot(r_mask)
  plot(st_geometry(tmp_adm2))
  
  lines(pol_vect, col = "red")  
  plot(st_geometry(tmp_adm2_98))
}


##############################################################################
#########   fraction for all subpopulation 
##############################################################################

# if(admin.level.tmp==1){
#   
#   adm_info_tmp <- surveyPrev::adminInfo(
#     poly.adm =country_shp_analysis[['Admin-1']],
#     admin = 1, 
#     by.adm=paste0("NAME_",1)) 
#   
# }else{
#   
#   adm_info_tmp <- surveyPrev::adminInfo(
#     poly.adm = country_shp_analysis[[admin.level.tmp+1]],
#     admin =2,
#     by.adm=paste0("NAME_",admin.level.tmp),
#     by.adm.upper = paste0("NAME_",admin.level.tmp-1)) 
#   
# }
