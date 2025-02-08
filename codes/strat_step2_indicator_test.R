##############################################################################
#########   load packagexs
##############################################################################
rm(list = ls())
# ENTER COUNTRY OF INTEREST -----------------------------------------------
# Please capitalize the first letter of the country name and replace " " in the country name to "_" if there is.
country <- 'Malawi'
country.abbrev <- 'mwi'
survey_year <- '2015'
year <- as.numeric(survey_year)
survey <- paste0(country,'_',survey_year)
alpha <- 0.05
# indicator
indicator <- "CH_DIAT_C_ORT"
indicator <- "HA_HIVP_B_HIV"
# Load libraries and info ----------------------------------------------------------
library(dplyr)
library(sf)
library(terra)
library(rdhs)
library(surveyPrev)
library(survey)


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


################################################################
#########   create folder for this indicator
################################################################
setwd(paste0(res.dir, '/UR_stratification/Indicators'))

if (!dir.exists(indicator)) {
  dir.create(indicator)
}

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

################################################################
#########   likelihood ratio test
################################################################
setwd(paste(code.dir))
setwd(paste(code.dir))
source('DataProcessing_helper.R')

indictor_test <- lr_test(data)


# if the p-value in lr_results is both significant (i.e., < 0.05), then it means we need to do the stratification.
# Thus, go to the step3.

##############################################################################
#########  store results
##############################################################################
setwd(paste0(res.dir, '/UR_stratification/Indicators/', indicator))

# Check if the file exists
sink(paste0(country.abbrev,"_",survey_year,"_","likelihood_ratio_test_",indicator,".txt"))

# Output the formulas of the models
cat("Formula for model1:",indictor_test$formula_model1, "\n")

cat("Formula for model2:",indictor_test$formula_model2, "\n")


cat("Formula for model3:",indictor_test$formula_model3, "\n")

# Print OR and 95% CI for UR
cat("UR's OR in model2：", round(indictor_test$or, 3),"\n")
cat("95% CI：(", round(indictor_test$or_ci_lower, 3), ",", round(indictor_test$or_ci_upper, 3), ")\n")


# Output p-values from ANOVA results
cat("P-value for ANOVA model1 vs model2:", indictor_test$p_value_1_vs_2,"\n")
cat("P-value for ANOVA model1 vs model3:", indictor_test$p_value_1_vs_3,"\n")

sink()
##############################################################################
#########   next step
##############################################################################

if(indictor_test$p_value_1_vs_2 < alpha) {
  cat(indicator,"has urban/rural association. \n")
}else{
  cat(indicator,"does not have urban/rural association. \n")
}

if(TRUE) {
  setwd(paste(code.dir))
  rstudioapi::navigateToFile("strat_step3_check_procedure.R")
}

