################################################################
#########   Specify survey 
################################################################
rm(list = ls()) # clear the R environment and prepare for the pipeline

### Please type the name of the survey you would like to analyze ### 
# Please capitalize the first letter of the survey name and replace " " in the survey name to "_" if there is.

country <- 'Malawi'
survey_year <- '2015'   ### if the survey spans two years, use the first year


survey <- paste0(country,'_',survey_year)
################################################################
#########   Main folder path creation 
################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set the directory, which is the user-specified folder


# create a folder called Results to store all the results for each survey, including the fitted R models, figures and tables in .csv etc. No actions if the folder is already there.
if(!dir.exists(paths = paste0('Results'))){ 
  dir.create(path = paste0('Results'))
}

# create a folder called Data to store all required data for each survey including DHS survey data, survey shape files etc. No actions if the folder is already there.
if(!dir.exists(paths = paste0('Data'))){ 
  dir.create(path = paste0('Data'))
}


# create a folder called Data to store all scripts
if(!dir.exists(paths = paste0('Scripts'))){ 
  dir.create(path = paste0('Scripts'))
}

# create a text file to store all the countries analyzed for the first run. If created, this file is loaded and the name of the countries analyzed are retrieved.
if(!file.exists("surveys_implemented.txt")){ 
  surveys <- c()
}else{
  surveys <- scan("surveys_implemented.txt", character(), quote = "")
}

if(length(surveys)==0 ||!survey %in% surveys){
  surveys <- c(surveys, survey)
  write(surveys, file = "surveys_implemented.txt") # save the newly added survey into the text file.
}


################################################################
#########   survey specific folder (Data)
################################################################

# create a folder for the user-specified survey in the fold of Data, which contains the following sub-folders:
# DHS_Data: This folder contains all the GPS locations where the DHS survey is carried out.
# shapeFiles_gadm: This folder contains all the shape files of the survey
# worldpop: This folder contains all the population data
# prepared_dat : This folder contains all the prepared data, including survey specific covariate surface, national prediction grid and training data
if(!dir.exists(paths = paste0('Data/',survey))){
  
  
  ### create folder structure within Data subfolder
  dir.create(path = paste0('Data/',survey))
  
  
  if(!dir.exists(paths = paste0('Data/',survey, '/DHS_Data'))){
    dir.create(path = paste0('Data/',survey, '/DHS_Data'))
  }
  if(!dir.exists(paths = paste0('Data/',survey, '/shapeFiles_gadm'))){
    dir.create(path = paste0('Data/',survey, '/shapeFiles_gadm'))
  }
  if(!dir.exists(paths = paste0('Data/',survey, '/worldpop'))){
    dir.create(path = paste0('Data/',survey, '/worldpop'))
  }
  if(!dir.exists(paths = paste0('Data/',survey, '/prepared_dat'))){
    dir.create(path = paste0('Data/',survey, '/prepared_dat'))
  }
  
}



################################################################
#########   survey specific folder (Scripts)
################################################################

# create a sub-folder for the survey in the folder of Results to store the scripts
if(!dir.exists(paths = paste0('Scripts/',survey))){
  
  dir.create(path = paste0('Scripts/',survey))
  
  ### when the pipeline is complete, copy scripts from a folder called 'template'
  
}



################################################################
#########   survey specific folder (Results)
################################################################

# create a sub-folder for the survey in the folder of Results to store the results.
if(!dir.exists(paths = paste0('Results/',survey))){
  
  dir.create(path = paste0('Results/',survey))
  
  ### as constructing the pipeline, add specific result folders
  
  
}


