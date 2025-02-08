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
library(purrr)

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

# create a folder to store prepared data
if(!dir.exists(paths = paste0(res.dir,'/GBT_model'))){
  dir.create(path = paste0(res.dir,'/GBT_model'))
}

if(!dir.exists(paths = paste0(res.dir,'/UR_classification'))){
  dir.create(path = paste0(res.dir,'/UR_classification'))
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
poly.adm.strata <- readRDS('poly_adm_strata.rds') ### strata admin level boundary
strata_admin_info <- readRDS('strata_adm_info.rds') ### strata admin level info


### load worldpop
setwd(paste0(data.dir,'/worldpop'))
pop.abbrev <- tolower(gadm.abbrev)
pop_file <- paste0(pop.abbrev,'_ppp_',frame_year,'_1km_Aggregated_UNadj.tif')

worldpop <- terra::rast(paste0(country.abbrev, '_ppp_', frame_year, '_1km_Aggregated_UNadj.tif'))

### load country prediction grid
setwd(paste0(data.dir,'/prepared_dat'))
natl_grid_df <- readRDS('natl_grid.rds')


### load jittering training data
setwd(paste0(data.dir,'/prepared_dat'))
train_dat <- readRDS('train_cluster_dat_cr.rds')



##############################################################################
#########   train GBT
##############################################################################

#train GBT with no travel
setwd(paste0(res.dir,'/GBT_model'))

if(file.exists('fitted_GBT_natl.rds')){
  fit_natl <- readRDS('fitted_GBT_natl.rds')
}else{
  
  ### set up mean structure
  cov_list<-c('pop_den','light','strata.adm.char')
  mean_formula<-as.formula(paste('urban~pop_den+light+as.factor(strata.adm.char)'))
  
  
  ### do not train with missing
  train_dat_complete <- train_dat[complete.cases(train_dat),]
  
  ### set up grid 
  gbtGrid <- expand.grid(interaction.depth=c(3,5,7,9), 
                         n.trees = (20:40)*50,
                         shrinkage=c(0.01,0.005,0.001),
                         n.minobsinnode=10)
  
  
  set.seed(2024)
  
  ### fit model
  fit_natl<-caret::train(mean_formula,
                         data = train_dat_complete,
                         method = "gbm",
                         na.action = na.pass,
                         trControl = caret::trainControl(method="cv", 
                                                         number=10),
                         verbose = 1,
                         tuneGrid=gbtGrid)
  
  #summary(fit_natl)
  #fit_natl$bestTune
  
  ### confusion matrix for training set
  train_pred <- predict(fit_natl, newdata = train_dat_complete)
  train_confusion <- confusionMatrix(train_pred, as.factor(train_dat_complete$urban))
  
  ### save fitted model and results 
  setwd(paste0(res.dir,'/GBT_model'))
  saveRDS(fit_natl,'fitted_GBT_natl.rds')
  saveRDS(train_confusion,'confusion_train.rds')
}


##############################################################################
#########  prediction on national grid
##############################################################################

complete_id <- complete.cases(natl_grid_df)

complete_natl_grid <- natl_grid_df[complete_id,]

pred_prob <- predict(fit_natl, complete_natl_grid,type = "prob",na.action = na.pass)

complete_natl_grid$pred_p <- pred_prob$urban

natl_pred_vec <- complete_natl_grid$pred_p

##############################################################################
#########  probability surface calibration
##############################################################################

# load strata fraction for benchmark
setwd(paste0(data.dir))
ref.tab <- readRDS(paste0(country.abbrev,'_ref_tab.rds'))

# need to standarize the following chunk for all surveys 
### colnames <-(strata.adm.num, strata.adm.char, strata.adm.name, urb_frac)
### values look like (1, strata_adm_1, Balaka, 0.071)

if(FALSE){
  load('strata_adm_urb_ref_tab.rda')
  
  ref.tab$strata.adm.num = c(1:dim(ref.tab)[1])
  ref.tab$strata.adm.char = paste0('strata_adm_',ref.tab$strata.adm.num)
  ref.tab$strata.adm.name = ref.tab$matched_name
  
  
  ref.tab <- ref.tab[,c('strata.adm.num', 'strata.adm.char', 'strata.adm.name', 'urb_frac')]
  
  setwd(paste0(data.dir))
  saveRDS(ref.tab,paste0(country.abbrev,'_ref_tab.rds'))
}


### add known admin fractions to national grid
grid_with_frac = merge(complete_natl_grid, 
                       ref.tab[c('strata.adm.char','urb_frac')],
                       by='strata.adm.char',
                       all.x=T)

if(sum(is.na(grid_with_frac$urb_frac))>0){
  message('Something is wrong, at least one admin region does not have a known fraction.')
}

### find admin specific tau
tau_results <- calculate_optimal_tau(grid_with_frac)

setwd(paste0(res.dir,'/UR_classification'))
saveRDS(tau_results,'tau_results.rds')

# calculate calibrated prob for each pixel
grid_with_frac_calibrated = merge(grid_with_frac, tau_results[c("strata.adm.char", "tau")], by = "strata.adm.char")
grid_with_frac_calibrated$cali_prob = expit(logit(grid_with_frac_calibrated$pred_p)+grid_with_frac_calibrated$tau)

# produce classification surface
natl_grid_prob <- merge(natl_grid_df,grid_with_frac_calibrated[,c('pixel_id','cali_prob')],
                        all.x=T)
natl_grid_prob <- natl_grid_prob[order(natl_grid_prob$pixel_id),]

setwd(paste0(res.dir,'/UR_classification'))
saveRDS(natl_grid_prob,'grid_with_calibrated_prob.rds')

# prepare surface

urb_prob_ras <- worldpop
values(urb_prob_ras) <- natl_grid_prob$cali_prob
plot(urb_prob_ras)

setwd(paste0(res.dir,'/UR_classification'))
terra::writeRaster(urb_prob_ras,paste0(country.abbrev,'_','urb_prob.tif'), overwrite = TRUE)

##############################################################################
#########  probability surface into binary surface (calibrated)
##############################################################################

setwd(paste0(res.dir,'/UR_classification'))
natl_grid_prob <- readRDS(file='grid_with_calibrated_prob.rds')

# load strata fraction for benchmark
setwd(paste0(data.dir))
ref.tab <- readRDS(paste0(country.abbrev,'_ref_tab.rds'))

### add known admin fractions to national grid
grid_prob_frac = merge(natl_grid_prob, 
                       ref.tab[c('strata.adm.char','urb_frac')],
                       by='strata.adm.char',
                       all.x=T)

adm_dat <- split(grid_prob_frac , f = grid_prob_frac$strata.adm.char)

thresh_urb<-function(adm_grid){
  
  # sort grid population
  vals <- adm_grid$cali_prob
  vals[is.na(vals)] <- 0
  
  pop_den <- adm_grid$pop_den
  pop_den[is.na(pop_den)] <- 0
  
  sort.obj <- sort.int(vals, decreasing = TRUE, index.return = TRUE, method = 'shell')
  svals <-  pop_den[sort.obj$ix]
  svals.int <- sort.obj$ix
  
  cutoff <- adm_grid$urb_frac[1]
  
  # determine population threshold and urban rural
  csvals <- cumsum(svals)/sum(svals)
  is.urb <- csvals < cutoff
  org.isurb <- is.urb[Matrix::invPerm(svals.int)]
  threshold <- min(vals[org.isurb == 1]) #cutoff
  
  # prepare return object (grid with urban/rural)
  adm_grid$threshold <- threshold
  adm_grid$urban <- as.numeric(org.isurb)
  #adm_grid[is.na(adm_grid$pop_den),]$urban<-NA
  
  return(adm_grid)
  
}

urb_list<-lapply(adm_dat, FUN=thresh_urb)

urb_class <- do.call("rbind", urb_list)

natl_ind_grid <- merge(natl_grid_prob,urb_class[,c('pixel_id','urban')],
                       by='pixel_id',all.x=T)
natl_ind_grid =natl_ind_grid %>% mutate(urban = case_when (is.na(pop_den) ~ NA,TRUE~urban))
urb_ind_ras <- worldpop

values(urb_ind_ras) <- natl_ind_grid$urban
plot(urb_ind_ras)

setwd(paste0(res.dir,'/UR_classification'))
terra::writeRaster(urb_ind_ras, paste0(country.abbrev, '_', 'ind_prob.tif'), overwrite = TRUE)


##############################################################################
#########  probability surface with pop density threshold
##############################################################################

if (FALSE){
  setwd(paste0(res.dir,'/UR_classification'))
  natl_grid_prob <- readRDS(file='grid_with_calibrated_prob.rds')
  
  # load strata fraction for benchmark
  setwd(paste0(data.dir))
  ref.tab <- readRDS(paste0(country.abbrev,'_ref_tab.rds'))
  
  ### add known admin fractions to national grid
  grid_prob_frac = merge(natl_grid_prob, 
                         ref.tab[c('strata.adm.char','urb_frac')],
                         by='strata.adm.char',
                         all.x=T)
  
  adm_dat <- split(grid_prob_frac , f = grid_prob_frac$strata.adm.char)
  
  pop_thresh_urb<-function(adm_grid){
    
    # sort grid population
    vals <- adm_grid$pop_den
    vals[is.na(vals)] <- 0
    
    pop_den <- adm_grid$pop_den
    pop_den[is.na(pop_den)] <- 0
    
    sort.obj <- sort.int(vals, decreasing = TRUE, index.return = TRUE, method = 'shell')
    svals <-  pop_den[sort.obj$ix]
    svals.int <- sort.obj$ix
    
    cutoff <- adm_grid$urb_frac[1]
    
    # determine population threshold and urban rural
    csvals <- cumsum(svals)/sum(svals)
    is.urb <- csvals < cutoff
    org.isurb <- is.urb[Matrix::invPerm(svals.int)]
    threshold <- min(vals[org.isurb == 1]) #cutoff
    
    # prepare return object (grid with urban/rural)
    adm_grid$threshold <- threshold
    adm_grid$urban <- as.numeric(org.isurb)
    #adm_grid[is.na(adm_grid$pop_den),]$urban<-NA
    
    return(adm_grid)
    
  }
  
  urb_list<-lapply(adm_dat, FUN=pop_thresh_urb)
  
  urb_class <- do.call("rbind", urb_list)
  
  natl_pop_grid <- merge(natl_grid_prob,urb_class[,c('pixel_id','urban')],
                         by='pixel_id',all.x=T)
  natl_pop_grid =natl_pop_grid %>% mutate(urban = case_when (is.na(pop_den) ~ NA,TRUE~urban))
  urb_pop_ras <- worldpop
  
  values(urb_pop_ras) <- natl_pop_grid$urban
  plot(urb_pop_ras)
  
  setwd(paste0(res.dir,'/UR_classification'))
  terra::writeRaster(urb_pop_ras,paste0(country.abbrev,'_','pop_prob.tif'))
  
}

##############################################################################
#########   next step
##############################################################################

### go to next file
setwd(paste(code.dir))

rstudioapi::navigateToFile("step3_get_fraction.R")

