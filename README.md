# Urban-Rural Stratification Pipeline

This repository contains scripts for conducting urban-rural stratification using DHS survey data. The pipeline is designed to process subpopulation data, generate urban-rural fractions, and evaluate stratification needs.
Please download the the pipline folder from this link and the code folder in this repository.

## Repository Structure

```
Urban-Rural-Pipeline/
├── create_folder.R
├── Data/
│   ├── DHS_meta.rda
│   ├── Global_rater/
│   │   ├── night_time_light/ (TIF files)
│   │   ├── subpop/
├── Results/
├── Scripts/
│   ├── Sub_population/
│   │   ├── prepare_subpop.R
```

## How to Run the Pipeline

### **1. Setup the Folder Structure**

This script aim to generate the relevant folder structure.

Modify the following variables in `create_folder.R` according to your chosen survey and then run the script:
```r
country <- 'Malawi'
survey_year <- '2015'   # Use the first year if the survey spans two years
```

Then, inside the `surveys_implemented.txt`, you will see the survey entry be added, for example:
```
Malawi_2015
```

---

### **2. Prepare Subpopulation Data**
Modify the following variables in `prepare_subpop.R` according to your chosen survey under pathway `Scripts/Sub_population/` and then run the script

```r
country <- 'Malawi'
gadm.abbrev <- 'mwi'
survey_year <- '2015'
```
Check the `Data/Global_rater/subpop/` folder for the generated population TIF files.

---

### **3. Move Scripts to the Country-Specific Folder**
Copy all scripts from `code` folder into the respective country-survey folder under `Scripts`, e.g.`Scripts/Malawi_2015/`

---

### We have prepared all structure so from the following we are going to start the steps for stratification
---

### **4. Run Initial Data Preparation**
Under `Scripts/Malawi_2015/`,modify the following variables in `step0_create_info.R` according to your chosen survey and then run the script:

```r
country <- 'Malawi'
gadm.abbrev <- 'MWI'
DHS.abbrev <- 'MWI'
survey_year <- 2015
frame_year <- 2008
calibrated_year <- 2008
strat_level_GADM <- 1
svy_training_year <- c(2010,2015)
```
Here `frame_year` is the census year which the DHS survey this based on,  `calibrated_year` is the. `svy_training_year` is referring to the DHS surveys who were conducted based on the same census survey. 
For example, Malawi 2015 survey and 2010 survey are all based on Malawi 2008 census survey. We combine them together to increase the training sample size for classification modeling.

---

### **5. Generate Spatial Data**
Run:
```r
source("Scripts/step1_prepare_dat.R")
```
Modify:
```r
country <- 'Malawi'
survey_year <- '2015'
```
This step generates required datasets and polygon information.

---

### **6. Create Reference DataFrames**
Generate the following `.rds` files:
- `mwi_ref_tab.rds` (contains urban fractions)
- `mwi_frame_ea.rds` (EA clusters from census)
- `mwi_sample_ea.rds` (EA clusters from DHS survey)

Ensure they are stored under:
```
Data/Malawi_2015/
```

---

### **7. Generate Urban-Rural Surface**
Run:
```r
source("Scripts/step2_UR_surface.R")
```
Modify:
```r
country <- 'Malawi'
survey_year <- '2015'
```

---

### **8. Compute Urban-Rural Fractions**
Run:
```r
source("Scripts/step3_get_fraction.R")
```
Modify:
```r
country <- 'Malawi'
survey_year <- '2015'
```
The results will be saved under:
```
Results/Malawi_2015/
```

---

### **9. Test for Oversampling Bias**
Run:
```r
source("Scripts/strat_step1_survey_test.R")
```
Modify:
```r
country <- 'Malawi'
country.abbrev <- 'mwi'
survey_year <- '2015'
```
Check the output in:
```
Results/Malawi_2015/UR_stratification/
```
This will produce a histogram and statistical tests for oversampling bias.

---

### **10. Indicator Association Testing**
Run:
```r
source("Scripts/strat_step2_indicator_test.R")
```
Modify:
```r
country <- 'Malawi'
country.abbrev <- 'mwi'
survey_year <- '2015'
indicator <- "CH_DIAT_C_ORT"
```
The output will indicate whether stratification is recommended.

---

### **11. Validate Stratification Procedure**
Run:
```r
source("Scripts/strat_step3_check_procedure.R")
```
Modify:
```r
country <- 'Malawi'
country.abbrev <- 'mwi'
gadm.abbrev <- 'mwi'
survey_year <- '2015'
indicator <- "CH_DIAT_C_ORT"
sub_pop <- "k0_5"
```
Valid `sub_pop` values:
- "k0_5"
- "total_pop"
- "k0_1"
- "m15_49"
- "f15_19"

---

### **12. Stratification at Admin Levels**
Run:
```r
source("Scripts/strat_step4_stratifcation_admin2:3.R")
```
Modify:
```r
country <- 'Malawi'
country.abbrev <- 'mwi'
gadm.abbrev <- 'mwi'
survey_year <- '2015'
admin.level.tmp <- 2
sub_pop <- "f15_49"
indicator <- "HA_HIVP_B_HIV"
```
Adjust `admin.level.tmp` for different administrative levels.

---

## Results & Interpretation
- The pipeline produces urban-rural fraction estimates for subpopulations.
- Statistical tests check for oversampling bias.
- Stratification recommendations depend on bias detection and indicator association.
- Final stratification procedures refine estimates at various admin levels.

For further details, check the results stored under:
```
Results/Malawi_2015/
```

---

## Contributors
- **Project Lead:** [Your Name]
- **Institution:** [Your Organization]
- **Contact:** [Your Email]

---

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments
- DHS Program for survey data.
- Contributors to the R packages used in this pipeline.
