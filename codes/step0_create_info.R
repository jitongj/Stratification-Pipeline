##################################################################
# This script is used to generate the info file for a given survey
##################################################################
rm(list = ls())

################################################################
#########   set parameters
################################################################

# Files info (For those lines with ### xxx ### above, please fill in as commented)
country <- 'Malawi'

### please fill in the country abbreviation in all upper case of gadm files ### (e.g. fill in SEN for gadm36_SEN_3.shp)
gadm.abbrev <- "MWI"
DHS.abbrev <- 'MWI'

### please fill in the following information ####

survey_year <- 2015 ### if spans two years, use the first one
frame_year <- 2008 ### year of the sampling frame creation
calibrated_year <- 2008

### Which GADM level corresponds to the stratification level used in this survey? 
### i.e. Which GADM level corresponds to DHS defined Admin-1?
strat_level_GADM <- 1

### whether multiple surveys share the same sampling frame
### in that case, use multiple surveys as training data
svy_training_year <- c(2010,2015)




survey <- paste0(country,'_',survey_year)
#survey_year_span <- 0 ### for single year survey, this should be 0, for surveys like Rwanda 2019-2020, this is 1


## setting rest of parameters using info from above
country.abbrev <- tolower(gadm.abbrev)           # lower the country gadm abbreviation 



################################################################
#########   set parameters
################################################################

## set directory

# path to home directory
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

home.dir <- paste(code.path.splitted[1: (length(code.path.splitted)-3)], collapse = "/")

data.dir<-paste0(home.dir,'/Data/',survey)
res.dir<-paste0(home.dir,'/Results/',survey)
code.dir <- paste0(home.dir,'/Scripts/',survey)




setwd(paste(data.dir))

info.name <- paste0(survey, "_general_info.Rdata")

save.image(file = paste0(info.name, sep=''))

##############################################################################
#########   next step
##############################################################################

### go to next file
setwd(paste(code.dir))

rstudioapi::navigateToFile("step1_prepare_dat.R")
