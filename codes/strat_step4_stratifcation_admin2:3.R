##############################################################################
#########   load packages
##############################################################################
rm(list = ls())
# ENTER COUNTRY OF INTEREST -----------------------------------------------
# Please capitalize the first letter of the country name and replace " " in the country name to "_" if there is.
country <- 'Malawi'
country.abbrev <- 'mwi'
gadm.abbrev <-'mwi'
survey_year <- '2015'
year <- as.numeric(survey_year)
survey <- paste0(country,'_',survey_year)
admin.level.tmp <- 2
sub_pop <- "f15_49"

# indicator
indicator <- "HA_HIVP_B_HIV"

# Load libraries and info ----------------------------------------------------------
library(sf)
library(terra)
library(rdhs)
library(surveyPrev)
library(survey)
library(INLA)
library(ggplot2)
library(patchwork)
library(dplyr)

################################################################
#########   load survey Info
################################################################
## set directory

# path to home directory
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

home.dir <- paste(code.path.splitted[1: (length(code.path.splitted)-3)], collapse = "/")

data.dir<-paste0(home.dir,'/Data/',survey)
res.dir<-paste0(home.dir,'/Results/',survey)
code.dir <- paste0(home.dir,'/Scripts/',survey)


##############################################################################
#########   load helper functions 
##############################################################################

setwd(paste(code.dir))
source('DataProcessing_helper.R')
##############################################################################
###### RDHS Configration 
##############################################################################

## login
# set API to get DHS data -- you will need to change this to your information!
#rdhs::set_rdhs_config(email = "xxxxx",
#                      project = "xxxxxxxxxxxxxx")


################################################################
#########   load DHS indicator
################################################################

dhsData <- getDHSdata(country = country, indicator = indicator, year = year)
data <- getDHSindicator(dhsData, indicator = indicator)
head(data)


##############################################################################
#########   Prepare  fraction
##############################################################################
setwd(paste0(data.dir,'/shapeFiles_gadm'))

country_shp_analysis <- readRDS('country_shp_analysis.rds')


if (TRUE) {
  admin.level.tmp = admin.level.tmp
  #sub_pop <- find_population_by_indicator(indicator, indicator_subpop_table)
  # load urban fraction for the subpop at survey year
  setwd(paste0(res.dir,'/UR_fraction/Admin_',admin.level.tmp))
  adm_urb_frac_name <- paste0(country.abbrev,'_',survey_year,'_admin',admin.level.tmp,'_',sub_pop,'_urb_frac.rds')
  adm_urb_frac <- readRDS(adm_urb_frac_name)
}

##############################################################################
#########   load prepared cluster data and admin information
##############################################################################

tmp.geo <- surveyPrev::getDHSgeo(country = country, year = year)

if(admin.level.tmp==1) {
  adm_info_tmp <- surveyPrev::adminInfo(
    poly.adm = country_shp_analysis[[admin.level.tmp+1]],
    admin =1,
    by.adm=paste0("NAME_",admin.level.tmp)) 
  
  tmp.cluster.info <- surveyPrev::clusterInfo(geo=tmp.geo, 
                                              poly.adm1=country_shp_analysis[[paste0('Admin-',admin.level.tmp)]],
                                              poly.adm2=country_shp_analysis[[paste0('Admin-',admin.level.tmp)]], 
                                              by.adm1 = paste0("NAME_",admin.level.tmp),
                                              by.adm2 = paste0("NAME_",admin.level.tmp))
  
  matched.order <- match(adm_info_tmp[["data"]]$admin1.name,
                         adm_urb_frac$adm.name.full)
  
}else{
  adm_info_tmp <- surveyPrev::adminInfo(
    poly.adm = country_shp_analysis[[admin.level.tmp+1]],
    admin =2,
    by.adm=paste0("NAME_",admin.level.tmp),
    by.adm.upper = paste0("NAME_",1)) # or admin.level.tmp-1
  
  ##### CHECK HERE!
  tmp.cluster.info <- surveyPrev::clusterInfo(geo=tmp.geo, 
                                              poly.adm1=country_shp_analysis[[paste0('Admin-',admin.level.tmp-1)]],
                                              poly.adm2=country_shp_analysis[[paste0('Admin-',admin.level.tmp)]], 
                                              by.adm1 = paste0("NAME_1"), ##### CHECK HERE!!!!!
                                              by.adm2 = paste0("NAME_2")) ##### CHECK HERE!!!!!
  matched.order <- match(adm_info_tmp[["data"]]$admin2.name.full,
                         adm_urb_frac$adm.name.full)
  
}



# asign urban fracton to the admin
adm_info_tmp[["data"]]$urban <- adm_urb_frac$urb_frac[matched.order]
# Find the elements of NA in the urban column and replace them with 0
adm_info_tmp[["data"]]$urban[is.na(adm_info_tmp[["data"]]$urban)] <- 0


##############################################################################
#########   stratified model
##############################################################################

if(admin.level.tmp==1) {
  cl_res_ad1_ur <- clusterModel(data=data,
                               cluster.info=tmp.cluster.info,
                               admin.info = adm_info_tmp,
                               model = "bym2",
                               stratification = TRUE,
                               admin = 1, 
                               CI = 0.95)
}else{
  cl_res_ad2_ur <- clusterModel(data=data,
                               cluster.info=tmp.cluster.info,
                               admin.info = adm_info_tmp,
                               model = "bym2",
                               stratification = TRUE,
                               admin = 2, 
                               CI = 0.95)
}

##############################################################################
#########   plot 
##############################################################################
setwd(paste0(res.dir, "/UR_Stratification/Indicators/", indicator))
if(admin.level.tmp==1) {
  head(cl_res_ad1_ur$res.admin1)
  
  cl_res_ad1 <- clusterModel(data=data,
                               cluster.info=tmp.cluster.info,
                               admin.info = adm_info_tmp,
                               model = "bym2",
                               stratification = FALSE,
                               admin = 1, 
                               CI = 0.95)
  
  s1 <- scatterPlot(
    res1=cl_res_ad1$res.admin1,
    res2=cl_res_ad1_ur$res.admin1[cl_res_ad1_ur$res.admin1$type=="full",],
    value1="mean",
    value2="mean",
    by.res1="admin1.name",
    by.res2="admin1.name",
    title="Stratified vs Unstratified model estimate",
    label1="Unstratified model estimates",
    label2="Stratified model estimates")
  s2 <- scatterPlot(
    res1=cl_res_ad1$res.admin1,
    res2=cl_res_ad1_ur$res.admin1[cl_res_ad1_ur$res.admin1$type=="full",],
    value1="sd",
    value2="sd",
    by.res1="admin1.name",
    by.res2="admin1.name",
    title="Stratified vs Unstratified model SD",
    label1="Unstratified model estimates",
    label2="Stratified model estimates")
  
  combined_plot <- s1 + s2
  
    # Display the combined plot
    print(combined_plot)
    ggsave(paste0(country.abbrev,"_",survey_year,"_admin",admin.level.tmp,"_visulization_",indicator,".png"), combined_plot, width = 10, height = 5, units = "in")
    
}else{
  head(cl_res_ad2_ur$res.admin2)
  cl_res_ad2 <- clusterModel(data=data,
                               cluster.info=tmp.cluster.info,
                               admin.info = adm_info_tmp,
                               model = "bym2",
                               stratification = FALSE,
                               admin = 2, 
                               CI = 0.95)
  s1 <- scatterPlot(
    res1=cl_res_ad2$res.admin2,
    res2=cl_res_ad2_ur$res.admin2[cl_res_ad2_ur$res.admin2$type=="full",],
    value1="mean",
    value2="mean",
    by.res1="admin2.name.full",
    by.res2="admin2.name.full",
    title="Stratified vs Unstratified model estimate",
    label1="Unstratified model estimates",
    label2="Stratified model estimates")
  s2 <- scatterPlot(
    res1=cl_res_ad2$res.admin2,
    res2=cl_res_ad2_ur$res.admin2[cl_res_ad2_ur$res.admin2$type=="full",],
    value1="sd",
    value2="sd",
    by.res1="admin2.name.full",
    by.res2="admin2.name.full",
    title="Stratified vs Unstratified model SD",
    label1="Unstratified model estimates",
    label2="Stratified model estimates")
  
  combined_plot <- s1 + s2
  # Display the combined plot
  print(combined_plot)
  ggsave(paste0(country.abbrev,"_",survey_year,"_admin",admin.level.tmp,"_visulization_",indicator,".png"), combined_plot, width = 10, height = 5, units = "in")
  
}

