###############################################################
###  load GADM files
###############################################################

get_country_GADM <- function(country,resolution=1) {
  
  country_iso3 <- DHS.country.meta[DHS.country.meta$CountryName==country,'ISO3_CountryCode']
  
  gadm_list <- list()
  levels <- 0
  repeat {
    tmp.gadm <- geodata::gadm(country = country_iso3, resolution=resolution,
                              level = levels,
                              path = tempdir())
    if (is.null(tmp.gadm)) {
      break
    } else {
      tmp.gadm <- sf::st_as_sf(tmp.gadm)
      tmp.gadm <- sf::st_set_crs(tmp.gadm, 4326)
      
      n_region <- dim(tmp.gadm)[1]
      #message(paste0('n region: ',n_region))
      if(n_region >1000){break}
      
      
      if(levels==0){      gadm_list[['National']]  <- tmp.gadm
      }else{
        gadm_list[[paste0('Admin-',levels)]]  <- tmp.gadm}
      levels <- levels + 1
    }
  }
  
  
  return(gadm_list)
}


###############################################################
###  standarize urban/rural naming
###############################################################

standardize_urban <- function(x) {
  x <- tolower(x)  # Convert to lower case to handle case sensitivity
  x <- ifelse(x %in% c("u", "urban"), "urban", x)
  x <- ifelse(x %in% c("r", "rural"), "rural", x)
  return(x)
}

###############################################################
###  helper to split underscore
###############################################################

split_underscore <- function(v) {
  # Check if the first element contains an underscore
  if (grepl("_", v[1])) {
    # Split each element at the underscore and keep the part after the underscore
    split_parts <- sapply(v, function(x) {
      parts <- strsplit(x, "_", fixed = TRUE)[[1]]
      if (length(parts) > 1) parts[2] else NA  # Assume there's something after the underscore
    })
    return(split_parts)
  } else {
    return(v) 
  }
}




###############################################################
###  helper function for calibration 
###############################################################

expit <- function(x) 1 / (1 + exp(-x))
logit <- function(x) log(x / (1 - x))


v_g_star <- function(v_g, tau) {
  expit(logit(v_g) + tau)
}


r_i_star <- function(v_g, N_g, tau) {
  sum(v_g_star(v_g, tau) * N_g) / sum(N_g)
}


L_MSE <- function(tau, v_g, N_g, r_i) {
  (r_i_star(v_g, N_g, tau) - r_i)^2
}


find_optimal_tau <- function(v_g, N_g, r_i) {
  optim(par = 0, fn = L_MSE, v_g = v_g, N_g = N_g, r_i = r_i, method = "BFGS")$par
}


calculate_optimal_tau <- function(data) {
  data %>%
    dplyr::group_by(strata.adm.char) %>%
    dplyr::summarise(
      tau = find_optimal_tau(v_g = pred_p, N_g = calibrated_pop_den, r_i = unique(urb_frac)),
      r_i = unique(urb_frac)
    ) %>%
    dplyr::ungroup()
}


###############################################################
###  get subnational urban fractions for subpopulation 
###############################################################

### Function Assign Admin to national grid

get_pixel_adm_grid <- function(pixel_grid,
                               admin.level,
                               admin.poly){
  
  # Match stratification admin
  pixel_grid_sf <- st_as_sf(pixel_grid, coords = c("x", "y"),
                            crs = st_crs(admin.poly))
  
  adm_match <- st_join(pixel_grid_sf,admin.poly, join = st_intersects)
  
  pixel_grid_df <- data.frame(
    x = pixel_grid$x,
    y = pixel_grid$y
  )
  
  if(admin.level==1){
    pixel_grid_df$adm.name =adm_match[[paste0('NAME_1')]]
    pixel_grid_df$adm.name.full = pixel_grid_df$adm.name
  }else{
    
    pixel_grid_df$adm.name = adm_match[[paste0('NAME_',admin.level)]]
    pixel_grid_df$upper.adm.name = adm_match[[paste0('NAME_',admin.level-1)]]
    pixel_grid_df$adm.name.full = paste0(adm_match[[paste0('NAME_',admin.level-1)]],
                                         '_',
                                         adm_match[[paste0('NAME_',admin.level)]])
  }
  
  
  
  return(pixel_grid_df[complete.cases(pixel_grid_df),])
}





get_adm_urb_frac <- function(pixel_adm_grid,
                             pred_surf,
                             pop_surf){
  
  pixel_adm_grid$pred_prob <- terra::extract(pred_surf, pixel_adm_grid[,c('x','y')])[,2]
  pixel_adm_grid$pop_den <- terra::extract(pop_surf, pixel_adm_grid[,c('x','y')])[,2]
  
  pixel_adm_grid <- pixel_adm_grid[complete.cases(pixel_adm_grid),]
  
  adm_urb_frac <- pixel_adm_grid %>%
    dplyr::group_by(adm.name.full) %>%
    dplyr::summarise(urb_frac = sum(pop_den * pred_prob, na.rm = TRUE) / 
                sum(pop_den, na.rm = TRUE))
  
  return(adm_urb_frac)
}




if(FALSE){
  pred_surf <- urb_prob_ras
  pop_surf <- worldpop
  
  pixel_adm_grid$pred_prob <- terra::extract(pred_surf, pixel_adm_grid[,c('x','y')])[,2]
  pixel_adm_grid$pop_den <- terra::extract(pop_surf, pixel_adm_grid[,c('x','y')])[,2]
  
  pixel_adm_grid <- pixel_adm_grid[complete.cases(pixel_adm_grid),]
  
  adm_urb_frac <- pixel_adm_grid %>%
    group_by(adm.name.full) %>%
    summarise(urb_frac = sum(pop_den * pred_prob, na.rm = TRUE) / 
                sum(pop_den, na.rm = TRUE))
  
  
  
#pixel_grid <- as.data.frame(terra::xyFromCell(worldpop, 1:ncell(worldpop)))
pixel_adm_grid

admin.level <- 1
admin.poly <- country_shp_analysis[[2]]



if(is.null(pixel_grid)){
  
  pixel_grid <- as.data.frame(terra::xyFromCell(pred_surf, 1:ncell(pred_surf)))
}

# Match stratification admin
pixel_grid_sf <- st_as_sf(pixel_grid[c('x','y')], coords = c("x", "y"),
                         crs = st_crs(admin.poly))

adm_match <- st_join(pixel_grid_sf,admin.poly, join = st_intersects)

pixel_grid_df <- data.frame(
  x = pixel_grid$x,
  y = pixel_grid$y,
  pred_prob =  terra::extract(pred_surf, pixel_grid)[,2],
  pop_den = terra::extract(pop_surf, pixel_grid)[,2]
)

if(admin.level==1){
  pixel_grid_df$adm.name =adm_match[[paste0('NAME_1')]]
}else{
  
  pixel_grid_df$adm.name = adm_match[[paste0('NAME_',admin.level)]]
  pixel_grid_df$upper.adm.name = adm_match[[paste0('NAME_',admin.level-1)]]
  pixel_grid_df$adm.name.full = paste0(adm_match[[paste0('NAME_',admin.level-1)]],
                                       '_',
                         adm_match[[paste0('NAME_',admin.level)]])
}


pixel_grid_df <- pixel_grid_df[complete.cases(pixel_grid_df),]


urb_frac <- merge_data %>%
  group_by(admin1.num) %>%
  summarise(frac = sum(f_15_49_pop * cali.prob, na.rm = TRUE) / 
              sum(f_15_49_pop, na.rm = TRUE))


st33rata_admin_info <- surveyPrev::adminInfo(
  poly.adm =country_shp_analysis[[3]],
  admin =2,by.adm=paste0("NAME_",2),by.adm.upper = paste0("NAME_",1)) 





}


##############################################################################
#########   Jittering
##############################################################################
constr_prior <- function(obs,jitter_r,prob_r,poly_admin,pop_ras, strat_level_GADM){
  library(sp)
  
  # for the cluster, find its coordinates and admin1 area
  admin1_index<-obs$admin1.name# actual admin1, eg malawi: 28
  # sp_xy<-SpatialPoints(as.data.frame(obs)[,c('LONGNUM','LATNUM')],
  #                      proj4string = CRS(proj4string(poly_admin)))
  
  # change from sf to sp
  poly_admin_sp <- as(poly_admin, "Spatial")
  
  sp_xy <- SpatialPoints(as.data.frame(obs)[, c('LONGNUM', 'LATNUM')],
                         proj4string = CRS(proj4string(poly_admin_sp)))
  
  
  #pt<-as.data.frame(obs)[,c('LONGNUM','LATNUM')]
  #colnames(pt)<-c('x','y')
  
  # generate gitter 
  #jitter_r<-2000
  cluster_buffer<-buffer(sp_xy, width=jitter_r)
  
  # # extract pixels within the buffer
  # temp_pop_cir<-mask(crop(pop_ras,cluster_buffer),
  #                    cluster_buffer)
  
  # change SpatialPolygons to SpatVector
  cluster_buffer_vect <- vect(cluster_buffer)
  temp_pop_cir <- mask(crop(pop_ras, cluster_buffer_vect), cluster_buffer_vect)
  
  
  # put admin area restriction
  admin_poly<-poly_admin[poly_admin[[paste0("NAME_",strat_level_GADM)]]==admin1_index,]
  temp_pop_admin<-mask(crop(temp_pop_cir,admin_poly),
                       admin_poly)
  
  
  # check whether need to adjust for constraint 
  cir_val<-values(temp_pop_cir)
  admin_val<-values(temp_pop_admin)
  
  admin_adj<-length(which(!is.na(cir_val)))!=length(which(!is.na(admin_val)))
  
  if(admin_adj){
    #normc<-admin1_normc(pt,jitter_r,admin_poly,ntrial=1000)
    normc<-1
  }else{  normc<-1}
  
  
  
  ## prepare sample frame
  
  temp_val<-values(temp_pop_admin)
  pop_index<-which(!is.na(temp_val))
  
  
  #temp_frame<-as.data.frame(coordinates(temp_pop_admin))
  temp_frame <- as.data.frame(terra::xyFromCell(temp_pop_admin, 1:ncell(temp_pop_admin)))
  
  pixel_candidate<-temp_frame[pop_index,]
  pixel_candidate$pop_den<-temp_val[pop_index]
  
  pixel_candidate$center_x<-obs$LONGNUM
  pixel_candidate$center_y<-obs$LATNUM
  library(geosphere)
  pixel_candidate$dist<-diag(distm(pixel_candidate[,c('x','y')], 
                                   pixel_candidate[,c('center_x','center_y')]))
  
  pixel_candidate$unn_w<-pixel_candidate$pop_den*
    1/(2*pi * 2 * pixel_candidate$dist)*normc*prob_r
  
  pixel_candidate$normc<-normc
  return(pixel_candidate[,c("x","y",'normc','unn_w')])
  #return(pixel_candidate)
  
}

##################################################################
####### Cochran-Mantel-Haensze Test
##################################################################

perform_cmh_test <- function(frame_ea, sample_ea, significance_level = 0.05) {
  
  # Create contingency tables for each strata
  contingency_tables <- lapply(1:nrow(frame_ea), function(i) {
    matrix(c(frame_ea$urban[i], frame_ea$rural[i],
             sample_ea$urban[i], sample_ea$rural[i]),
           nrow = 2, byrow = TRUE,
           dimnames = list(Group = c("Frame", "Sample"),
                           Area = c("Urban", "Rural")))
  })
  
  # Perform Cochran-Mantel-Haenszel test
  cmh_result <- mantelhaen.test(array(unlist(contingency_tables), 
                                      dim = c(2, 2, nrow(frame_ea))))
  
  # Extract p-value
  p_value <- cmh_result$p.value
  
  # Determine significance
  is_significant <- ifelse(p_value < significance_level, 
                           "The result is statistically significant.", 
                           "The result is not statistically significant.")
  
  # Print p-value and significance
  cat("P-value:", p_value, "\n")
  cat(is_significant, "\n")
  
  # Return p-value and significance message as a list
  return(list(p_value = p_value, significance = is_significant))
}

# Create the plot
create_chi_square_plot <- function(simulated_chi_square, observed_chi_square, n_regions, p_value, country.abbrev, survey_year) {
  # Create a data frame for the histogram
  df_hist <- data.frame(chi_square = simulated_chi_square)
  
  # Create sequence for the theoretical density curve
  x_seq <- seq(0, max(simulated_chi_square), length.out = 200)
  df_density <- data.frame(
    x = x_seq,
    y = dchisq(x_seq, df = n_regions)
  )
  
  # Create the plot
  p <- ggplot(df_hist, aes(x = chi_square)) +
    # Add histogram
    geom_histogram(aes(y = ..density..), bins = 30, 
                   fill = "lightgray", color = "white") +
    # Add theoretical density curve
    geom_line(data = df_density, aes(x = x, y = y, color = "Chi-Square Density"),
              linewidth = 1) +
    # Add vertical line for observed value
    geom_vline(aes(xintercept = observed_chi_square, 
                   color = "Observed Chi-Square"),
               linewidth = 1, linetype = "dashed") +
    # Customize colors
    scale_color_manual(name = "",
                       values = c("Chi-Square Density" = "blue",
                                  "Observed Chi-Square" = "red"),
                       guide = guide_legend(direction = "horizontal",
                                            override.aes = list(
                                              linetype = c("solid", "dashed")))) +
    # Add labels and title
    labs(x = "Chi-Square Value",
         y = "Density",
         title = "Histogram of Simulated Chi-Square Values with Density Curve",
         subtitle = paste("p-value =", p_value)) +  # Remove rounding
    # Add theme elements
    theme_light() +
    theme(
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      legend.position = "top",
      legend.box = "horizontal",
      legend.spacing.x = unit(0.5, 'cm'),  # Add some spacing between legend items
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  
  # Save the plot
  ggsave(
    filename = paste0(country.abbrev, "_", survey_year, "_chi_square_histogram.png"),
    plot = p,
    width = 10,
    height = 7.5,
    dpi = 150
  )
  
  return(p)
}

# Usage example:
# plot <- create_chi_square_plot(
#   simulated_chi_square = simulated_chi_square,
#   observed_chi_square = observed_chi_square,
#   n_regions = n_regions,
#   p_value = p_value,
#   country.abbrev = country.abbrev,
#   survey_year = survey_year
# )
##################################################################
####### Find indicator's corresponding population
##################################################################

# This is the table of indicators with their corresponding targeted population
# Use this for choosing the urban fraction
# has not finished!!!!
indicator_subpop_table <- data.frame(
  indicator = c("WS_SRCE_P_BAS", "HA_HIVP_B_HIV", "CN_NUTS_C_WH2"),
  population = c("total_pop", "f15_49", "k0_5"),
  stringsAsFactors = FALSE
)


# Define the function to find the population based on the input indicator
find_population_by_indicator <- function(indicator, data) {
  # Search for the matching indicator in the data frame
  result <- data[data$indicator == indicator, "population"]
  
  # Check if a matching result is found
  if (length(result) == 0) {
    cat("No matching population found for the given indicator.\n")
    return(NULL)
  } else {
    cat("Population corresponding to the indicator", indicator, "is:", result, "\n")
    return(result)
  }
}

#population <- find_population_by_indicator("HA_HIVP_B_HIV", indicator_pop_table)


################################################################
#########   likelihood ratio test funrction for indicator
################################################################

lr_test <- function(tmp.analysis.dat) {
  
  # Remove NAs
  if (sum(is.na(tmp.analysis.dat$value)) > 0) {
    tmp.analysis.dat <- tmp.analysis.dat[rowSums(is.na(tmp.analysis.dat)) == 0, ]
  }
  
  tmp.analysis.dat <- dplyr::mutate(tmp.analysis.dat, strata = dplyr::case_when(
    strata == "urban" ~ 1,
    strata == "rural" ~ 0,
    TRUE ~ NA_real_  # Sets to NA for any other case
  ))
  
  # Rename strata to UR and process admin column
  names(tmp.analysis.dat)[names(tmp.analysis.dat) == "strata"] <- "UR"
  tmp.analysis.dat$value <- as.integer(tmp.analysis.dat$value)
  tmp.analysis.dat$admin <- gsub(" - urban| - rural", "", tmp.analysis.dat$v022)
  tmp.analysis.dat$admin <- gsub(" urban| rural", "", tmp.analysis.dat$admin)
  
  tmp.analysis.dat$admin_UR <- paste(tmp.analysis.dat$admin, tmp.analysis.dat$UR, sep = "_")
  admin_UR_counts <- table(tmp.analysis.dat$admin_UR)
  valid_admin_UR <- names(admin_UR_counts[admin_UR_counts > 1])
  tmp.analysis.dat <- tmp.analysis.dat[tmp.analysis.dat$admin_UR %in% valid_admin_UR, ]
  
  
  tmp.dhs.design <- survey::svydesign(id = ~cluster, weights = ~weight, strata = ~v022,
                                      nest = TRUE, survey.lonely.psu = "average", data = tmp.analysis.dat)
  
  # Fit models
  model1 <- survey::svyglm(value ~ admin, design = tmp.dhs.design, family = quasibinomial)
  model2 <- survey::svyglm(value ~ admin + UR, design = tmp.dhs.design, family = quasibinomial)
  model3 <- survey::svyglm(value ~ admin + admin:UR, design = tmp.dhs.design, family = quasibinomial)
  
  #formulas of the models
  formula_model1 <- paste(deparse(formula(model1)), collapse = " ")
  formula_model2 <- paste(deparse(formula(model2)), collapse = " ")
  formula_model3 <- paste(deparse(formula(model3)), collapse = " ")
  
  # Compute OR for UR in model2
  ur_coef <- coef(model2)["UR"]
  ur_se <- summary(model2)$coefficients["UR", "Std. Error"]
  
  ci_lower <- ur_coef - 1.96 * ur_se
  ci_upper <- ur_coef + 1.96 * ur_se
  
  or <- exp(ur_coef)
  or_ci_lower <- exp(ci_lower)
  or_ci_upper <- exp(ci_upper)
  
  # Perform ANOVA between models and extract p-values
  anova_result_1_vs_2 <- anova(model1, model2)
  anova_result_1_vs_3 <- anova(model1, model3)
  
  p_value_1_vs_2 <- anova_result_1_vs_2$p
  p_value_1_vs_3 <- anova_result_1_vs_3$p
  
  
  return(list(
    formula_model1 = formula_model1,
    formula_model2 = formula_model2,
    formula_model3 = formula_model3,
    or = or,
    or_ci_lower = or_ci_lower,
    or_ci_upper = or_ci_upper,
    p_value_1_vs_2 = p_value_1_vs_2,
    p_value_1_vs_3 = p_value_1_vs_3
  ))
}

