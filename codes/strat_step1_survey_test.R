##############################################################################
#########   load packages and set parameters
##############################################################################
rm(list = ls())
# ENTER COUNTRY OF INTEREST -----------------------------------------------
# Please capitalize the first letter of the country name and replace " " in the country name to "_" if there is.
country <- 'Malawi'
country.abbrev <- 'mwi'
survey_year <- '2015'
survey <- paste0(country,'_',survey_year)
alpha <- 0.05

# Load libraries and info ----------------------------------------------------------
library(dplyr)

# path to home directory
code.path <- rstudioapi::getActiveDocumentContext()$path
code.path.splitted <- strsplit(code.path, "/")[[1]]

home.dir <- paste(code.path.splitted[1: (length(code.path.splitted)-3)], collapse = "/")


data.dir<-paste0(home.dir,'/Data/',survey)
res.dir<-paste0(home.dir,'/Results/',survey)
code.dir <- paste0(home.dir,'/Scripts/',survey)


##############################################################################
#########   create folders to store results
##############################################################################

setwd(paste(res.dir))

if(!dir.exists(paths = paste0(res.dir,'/UR_stratification'))){
  dir.create(path = paste0(res.dir,'/UR_stratification'))
}
##############################################################################
#########   load helper functions 
##############################################################################
setwd(paste(code.dir))
source('DataProcessing_helper.R')

##############################################################################
#########   Cochran-Mantel-Haensze Test for survey
##############################################################################
# load strata fraction for sampling frame
setwd(paste0(data.dir))
library(purrr)
# load strata fraction for sampling frame
frame_ea <- readRDS(paste0(country.abbrev,'_frame_ea.rds'))
sample_ea <- readRDS(paste0(country.abbrev,'_sample_ea.rds'))

### make sure frame_ea and sample_ea are matched!!!
order_index <- match(frame_ea$strata, sample_ea$strata)
sample_ea <- sample_ea[order_index, ]

frame_ea$ur_frac <- frame_ea$urban/frame_ea$total
sample_ea$ur_frac <- sample_ea$urban/sample_ea$total

# Identify admin region where frame_ea$ur_frac is 0 or 1
region_to_remove <- which(frame_ea$ur_frac == 0 | frame_ea$ur_frac == 1)
if (length(region_to_remove) > 0) {
  cat("Exclude admin region", frame_ea$strata[region_to_remove], "from chi-squared test.\n")
  
  # Remove the identified regions from both data frames
  frame_ea <- frame_ea[-region_to_remove, ]
  sample_ea <- sample_ea[-region_to_remove, ]
}

n_regions <- dim(frame_ea)[1]
sample_sizes <-sample_ea$total
sample_urban_counts <- round(sample_ea$ur_frac * sample_sizes) # urban count in sample
sample_rural_counts <- sample_sizes - sample_urban_counts # rural count in sample

calculate_chi_square <- function(urban_counts, rural_counts, expected_urban_fraction, sizes) {
  expected_urban <- sizes * expected_urban_fraction
  expected_rural <- sizes - expected_urban
  urban_chi <- (urban_counts - expected_urban)^2 / expected_urban
  rural_chi <- (rural_counts - expected_rural)^2 / expected_rural
  return(sum(urban_chi + rural_chi))
}

# Observed Chi-square
observed_chi_square <- calculate_chi_square(sample_urban_counts, sample_rural_counts, frame_ea$ur_frac, sample_sizes)

set.seed(2024)
# MC
n_simulations <- 10000
simulated_chi_square <- numeric(n_simulations)
for (i in 1:n_simulations) {
  simulated_urban_counts <- rbinom(n_regions, sample_sizes, frame_ea$ur_frac)
  simulated_rural_counts <- sample_sizes - simulated_urban_counts
  simulated_chi_square[i] <- calculate_chi_square(simulated_urban_counts, simulated_rural_counts, frame_ea$ur_frac, sample_sizes)
}
p_value <- sum(simulated_chi_square >= observed_chi_square) / n_simulations


##############################################################################
#########  store results
##############################################################################
setwd(paste0(res.dir, '/UR_stratification'))

# save p-value
sink(paste0(country.abbrev,"_",survey_year,"_","survey_testing.txt"))
cat("Observed Chi-square:", observed_chi_square, "\n")
cat("p-value:", p_value, "\n")
if (p_value < alpha) {
  cat("The test statistic is significant. The country",country,"has sampling issue. \n")
}else{
  cat("The test statistic is not significant. The country",country,"does not have sampling issue. \n")
}

sink()

library(ggplot2)
library(gridExtra)
plot <- create_chi_square_plot(
  simulated_chi_square = simulated_chi_square,
  observed_chi_square = observed_chi_square,
  n_regions = n_regions,
  p_value = p_value,
  country.abbrev = country.abbrev,
  survey_year = survey_year
)



##############################################################################
#########   next step
##############################################################################
# if (overall_p_value < alpha) {
#   setwd(paste0(res.dir, '/UR_stratification'))
#   
#   # create folder for Indicators if needed
#   if (!dir.exists('Indicators')) {
#     dir.create('Indicators')
#   }
#   
#   ### go to next file
#   setwd(paste(code.dir))
#   
#   rstudioapi::navigateToFile("strat_step2_indicator_test.R")
# }


if (TRUE) {
  setwd(paste0(res.dir, '/UR_stratification'))
  
  # create folder for Indicators if needed
  if (!dir.exists('Indicators')) {
    dir.create('Indicators')
  }
  
  ### go to next file
  setwd(paste(code.dir))
  
  rstudioapi::navigateToFile("strat_step2_indicator_test.R")
}


