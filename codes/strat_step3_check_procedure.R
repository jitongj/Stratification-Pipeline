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
strat_level_GADM <- 1

# indicator
indicator <- "CH_DIAT_C_ORT"
indicator <-"HA_HIVP_B_HIV"
sub_pop <- "k0_5" # the corresponding population group to the indicator
sub_pop <- "f15_49"
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

setwd(paste(data.dir))

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
# delete the row here valu=NA
data <- data %>% filter(!is.na(value))
##############################################################################
#########   load prepared data 
##############################################################################
### load admin info
setwd(paste0(data.dir,'/shapeFiles_gadm'))
country_shp_analysis <- readRDS('country_shp_analysis.rds')
poly.adm.strata <- readRDS('poly_adm_strata.rds')


tmp.geo <- surveyPrev::getDHSgeo(country = country, year = year)
tmp.cluster.info <- surveyPrev::clusterInfo(geo=tmp.geo, 
                                            poly.adm1=country_shp_analysis[[paste0('Admin-',1)]],
                                            poly.adm2=country_shp_analysis[[paste0('Admin-',2)]], 
                                            by.adm1 = paste0("NAME_1"),
                                            by.adm2 = paste0("NAME_2"))

cluster_data_tmp <- data %>%
  group_by(cluster) %>%
  dplyr::summarise(
    strata = unique(strata),
    Ntrials = n(),
    value = sum(value)
  )

cluster_data <-merge(cluster_data_tmp, tmp.cluster.info[["data"]], by="cluster")
cluster_data <- as.data.frame(cluster_data)
names(cluster_data)[names(cluster_data) == "strata"] <- "UR"
cluster_data <- cluster_data %>%
  mutate(UR = ifelse(UR == "urban", 1, 0))


##############################################################################
#########   Prepare  worldpop of survey year
##############################################################################

setwd(paste0(data.dir,'/worldpop'))
### download population
pop.abbrev <- tolower(gadm.abbrev)
pop_file <- paste0(pop.abbrev,'_ppp_',year,'_1km_Aggregated_UNadj.tif')


if(!file.exists(pop_file)){
  
  url <- paste0("https://data.worldpop.org/GIS/Population/Global_2000_2020_1km_UNadj/", 
                year, "/", toupper(pop.abbrev),"/",      
                pop.abbrev,'_ppp_',year,'_1km_Aggregated_UNadj.tif')
  
  download.file(url, pop_file, method = "libcurl",mode="wb")
}

### load population
worldpop <- terra::rast(paste0(country.abbrev, '_ppp_', year, '_1km_Aggregated_UNadj.tif'))

##############################################################################
#########   Unstratified prevalence model: only fixed admin
##############################################################################
logit <- function(p) log(p / (1 - p))
expit <- function(x) 1 / (1 + exp(-x))

unstrat_fixedAdmin = function(myData){
  new_data <- myData[, c("value", "Ntrials", "admin1.name", "UR")]
  
  # Create the interaction term in the formula
  formula = value ~ -1 + admin1.name
  
  # Fit
  res.inla = inla(formula = formula,
                  data = new_data,
                  family = "betabinomial",
                  Ntrials = new_data$Ntrials,
                  control.compute = list(config = TRUE, dic = TRUE, cpo = TRUE, waic = TRUE),
                  control.predictor = list(compute = TRUE, link = 1))
  
  return(res.inla)
}

unstrat_fixedAdmin_results = unstrat_fixedAdmin(cluster_data)

unstrat_pi <- data.frame(
  admin = sub("admin1.name", "", rownames(unstrat_fixedAdmin_results[["summary.fixed"]])),
  unstrat_pi = expit(unstrat_fixedAdmin_results$summary.fixed$mean)
)

##############################################################################
#########   Stratified prevalence model: fixed admin + UR
##############################################################################
admin.num = dim(country_shp_analysis[["Admin-1"]])[1]
  
strat_fixedAdmin_fixedUR = function(myData){
  new_data <- myData[, c("value", "Ntrials", "admin1.name", "UR")]
  
  # Create the interaction term in the formula
  formula = value ~ -1 + admin1.name + UR
  
  # Fit
  res.inla = inla(formula = formula,
                  data = new_data,
                  family = "betabinomial",
                  Ntrials = new_data$Ntrials,
                  control.compute = list(config = TRUE, dic = TRUE, cpo = TRUE, waic = TRUE),
                  control.predictor = list(compute = TRUE, link = 1))
  
  return(res.inla)
}

strat_fixedAdmin_fixedUR_results = strat_fixedAdmin_fixedUR(cluster_data)

strat_pi_fixedAdmin_fixedUR <- data.frame(
  admin = sub("admin1.name", "", rownames(strat_fixedAdmin_fixedUR_results[["summary.fixed"]])[1:admin.num]),
  pi_u = expit(strat_fixedAdmin_fixedUR_results$summary.fixed$mean[1:admin.num] # admin effect
               + strat_fixedAdmin_fixedUR_results$summary.fixed$mean[admin.num+1]),  # UR effect
  pi_r = expit(strat_fixedAdmin_fixedUR_results$summary.fixed$mean[1:admin.num]) # admin effect
)

##############################################################################
#########   Prepare urban fraction
##############################################################################

if (TRUE) {
  admin.level.tmp = strat_level_GADM
  # load urban fraction for the subpop at survey year
  setwd(paste0(res.dir,'/UR_fraction/Admin_',admin.level.tmp))
  adm_urb_frac_name <- paste0(country.abbrev,'_',survey_year,'_admin',admin.level.tmp,'_',sub_pop,'_urb_frac.rds')
  adm_urb_frac <- readRDS(adm_urb_frac_name)
}



##############################################################################
#########   Aggregate stratified prevalence model:
##############################################################################
# load fractions: Replace later for our calculated fraction
setwd(paste0(data.dir))



matched.order <- match(strat_pi_fixedAdmin_fixedUR$admin,
                       adm_urb_frac$adm.name.full)

strat_pi_fixedAdmin_fixedUR$urb_frac <- adm_urb_frac$urb_frac[matched.order]

strat_pi_fixedAdmin_fixedUR <- strat_pi_fixedAdmin_fixedUR %>%
  mutate(
    strat_pi = pi_u * urb_frac + pi_r * (1 - urb_frac)
  )


##############################################################################
#########   Direct estimate model
##############################################################################
options(survey.lonely.psu = "adjust")
res_ad <- directEST(data = data,
                     cluster.info = tmp.cluster.info,
                     admin = 1)
res_ad$res.admin1 <- res_ad$res.admin1 %>%
  dplyr::rename(admin= admin1.name)


##############################################################################
#########   Prepare admin population as weight
##############################################################################

## NOTE: should later put this part in the fraction section
setwd(paste0(res.dir,'/UR_fraction'))

if(file.exists('admin_pop_weight.rds')){
  admin_pop_df <- readRDS('admin_pop_weight.rds')
}else{
  pixel_grid_pop <- as.data.frame(terra::xyFromCell(worldpop, 1:ncell(worldpop)))
  
  admin.level.tmp = strat_level_GADM
  
  pixel_adm_grid_pop <- get_pixel_adm_grid(pixel_grid = pixel_grid_pop,
                                           admin.level = admin.level.tmp,
                                           admin.poly = country_shp_analysis[[admin.level.tmp+1]])
  
  pixel_adm_grid_pop$pop_den <- terra::extract(worldpop, pixel_adm_grid_pop[,c('x','y')])[,2]
  
  pixel_adm_grid_pop <- pixel_adm_grid_pop[complete.cases(pixel_adm_grid_pop),]
  
  admin_pop_df <- pixel_adm_grid_pop %>%
    dplyr::group_by(adm.name.full) %>%
    dplyr::summarise(total_pop = sum(pop_den, na.rm = TRUE))
  
  
  admin_pop_df <- admin_pop_df %>%
    dplyr::rename(admin = adm.name.full)
  
  saveRDS(admin_pop_df, "admin_pop_weight.rds")
}
##############################################################################
#########   Numerical Comparison: Unstratified vs Direct, Stratified vs Direct
##############################################################################
setwd(paste0(res.dir, "/UR_Stratification/Indicators/", indicator))


unstrat_vs_direct <- merge(unstrat_pi,res_ad$res.admin1[, c("admin", "direct.est")], by="admin")
unstrat_vs_direct <- merge(unstrat_vs_direct,admin_pop_df, by="admin")

strat_vs_direct <- merge(strat_pi_fixedAdmin_fixedUR,res_ad$res.admin1[, c("admin", "direct.est")], by="admin")
strat_vs_direct <- merge(strat_vs_direct,admin_pop_df, by="admin")


unstrat_weighted_bias <- unstrat_vs_direct %>%
  summarise(
    weighted_bias = sum(abs(unstrat_pi - direct.est) * total_pop) / sum(total_pop)
  ) %>%
  pull(weighted_bias)

strat_weighted_bias <- strat_vs_direct %>%
  summarise(
    weighted_bias = sum(abs(strat_pi - direct.est) * total_pop) / sum(total_pop)
  ) %>%
  pull(weighted_bias)



sink(paste0(country.abbrev,"_",survey_year,"_","check_UR_procedure_",indicator,".txt"))

# Your if-else statement
if (unstrat_weighted_bias > strat_weighted_bias) {
  cat("UR stratification procedure works. \n")
  cat("Weighted bias of unstratified model is: ", unstrat_weighted_bias, "\n")
  cat("Weighted bias of stratified model is: ", strat_weighted_bias, "\n")
} else {
  cat("Warning: something wrong with the stratification procedure. \n")
  cat("Weighted bias of unstratified model is: ", unstrat_weighted_bias, "\n")
  cat("Weighted bias of stratified model is: ", strat_weighted_bias, "\n")
}

sink()

##############################################################################
#########   Visualized comparison: Unstratified vs Direct, Stratified vs Direct
##############################################################################
setwd(paste0(res.dir, "/UR_Stratification/Indicators/", indicator))
# First plot: res_ad1$res.admin1$direct.est vs unstrat_pi$pi
plot1 <- ggplot() +
  geom_point(aes(x = res_ad$res.admin1$direct.est, y = unstrat_pi$unstrat_pi), color = "blue") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(x = "Direct Estimate", y = "Unstratified Estimate", title = "Unstratified Estimate vs Direct Estimate") +
  theme_minimal()

# Second plot: res_ad1$res.admin1$direct.est vs strat_pi_fixedAdmin_fixedUR$pi
plot2 <- ggplot() +
  geom_point(aes(x = res_ad$res.admin1$direct.est, y = strat_pi_fixedAdmin_fixedUR$strat_pi), color = "blue") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(x = "Direct Estimate", y = "Stratified Estimate", title = "Stratified Estimate vs Direct Estimate") +
  theme_minimal()

# Combine the two plots into a 1x2 grid
combined_plot <- plot1 + plot2

# Display the combined plot
print(combined_plot)
ggsave(paste0(country.abbrev,"_",survey_year,"_","visualized_comparison_",indicator,".png"), combined_plot, width = 10, height = 5, units = "in")


##############################################################################
#########   next step
##############################################################################

# if(unstrat_weighted_bias > strat_weighted_bias) {
#   setwd(paste(code.dir))
#   rstudioapi::navigateToFile("strat_step4_stratifcation_admin2:3.R")
# }

setwd(paste(code.dir))
rstudioapi::navigateToFile("strat_step4_stratifcation_admin2_3.R")
