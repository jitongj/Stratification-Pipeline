# Urban-Rural Stratification Pipeline

In this document, we give an overview of the urban-rural stratification pipeline and suggest the order at which scripts are run. Note that our pipeline is not aimed to replicate the estimates presented in DHS final report. Due to distinctions in methods, our estimates will not be exactly the same as DHS estimates, but they will be consistent.

This repository contains scripts for conducting urban-rural stratification using DHS survey data. The pipeline is designed to process subpopulation data, generate urban-rural fractions, and evaluate stratification needs. Please download the the pipline folder from [UR-Pipeline-Demo](https://uwnetid-my.sharepoint.com/:f:/g/personal/jitongj_uw_edu/EplTkYWv0eFCjjwmH7PW-qwBiWzINyB8VIr6VYBMcHPdQQ?e=7wLIAi) and the **codes** folder in this repository.

The detailed methodology of the stratification modelcan be found in the referenced document here: [Yunhan Wu and Jon Wakefield](https://academic.oup.com/jrsssa/article/187/3/811/7578683), and the **Stratification_Pipeline.pdf** in this repository.

We also recommend you to try on our Shiny App: [SAE4Health](https://sae4health.stat.uw.edu)


## UR-Pipeline-Demo Structure

```
UR-Pipeline-Demo/
├── create_folder.R
├── Data/
│   ├── DHS_meta.rda
│   ├── Global_rater/
│   │   ├── night_time_light/ (TIF files)
├── Scripts/
│   ├── Sub_population/
│   │   ├── prepare_subpop.R
```

## Packages
We use R to conduct all the statistical analysis. Packages and dependencies could be installed as following through the pipeline. Most packages are found under CRAN. For **surveyPrev**, we recommend you download from [surveyPrev GitHub repository](https://github.com/richardli/surveyPrev).


## How to Run the Pipeline

### **1. Setup the Folder Structure**

This script aims to generate the relevant folder structure.

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
### With the structure in place, we will now proceed with the stratification steps.
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
Here `frame_year` is the census year which the DHS survey this based on,  `calibrated_year` is the year which DHS survey's u/r population fraction is based on (See section 6 for more details). `svy_training_year` is referring to the DHS surveys who were conducted based on the same census survey. 
For example, Malawi 2015 survey and 2010 survey are all based on Malawi 2008 census survey, so We combine them together to increase the training sample size for classification modeling.

---

### **5. Generate Spatial Data**
Modify the following variables in `step1_prepare_dat.R` according to your chosen survey and then run the script:

```r
country <- 'Malawi'
survey_year <- '2015'
```
This step generates required datasets and polygon information.

---

### **6. Create Reference DataFrames**
Generate the following `.rds` files, naming them with the `country.abbrev` prefix. For example:

- `mwi_ref_tab.rds` (contains urban fractions)
- `mwi_frame_ea.rds` (EA clusters in census)
- `mwi_sample_ea.rds` (EA clusters in DHS survey)

(change `mwi` for your country's abbrev)

Each of these files should be structured as follows:

For `mwi_ref_tab.rds`, create a data frame in the following format:

![image](https://github.com/user-attachments/assets/dfb0b11f-5ad7-4a4e-b516-983ba3de2ba3)


- `urb_frac` should be derived from the urban-rural population table in the DHS survey, typically located in Appendix A of the survey reports. For greater accuracy, use the actual population numbers instead of the urban percentage provided in the table. For example.
  ![image](https://github.com/user-attachments/assets/154dc67d-3a95-49b7-adb0-5a2e47ea65c1)

- If the DHS survey does not provide this information, refer to the corresponding census survey or estimate it based on urban and rural household numbers (you can calculate the u/r population based on the number of average u/r population per household in the census survey)

- Ensure region names's spelling and order match exactly with those in `country_shp_analysis.rds`.
- The year on which this fraction is based will serve as our calibration year, as it will be used to adjust the fractions calculated in subsequent steps.

For `mwi_frame_ea.rds` and `mwi_sample_ea.rds`, create the data frame in the following format:

![image](https://github.com/user-attachments/assets/9e934b66-5ad6-454b-804b-6a636f7e4423)


- Ensure region names's spelling and order match exactly with those in `country_shp_analysis.rds`. ((eg. Nigeria 2018 survey report Table A.2 and Table A.3 in Appendix A)
  

Ensure they are stored under the country_survey folder under `Data`, for example:
```
Data/Malawi_2015/
```

---

### **7. Generate Urban-Rural Surface**
Modify the following variables in `step2_UR_surface.R` according to your chosen survey and then run the script:

```r
country <- 'Malawi'
survey_year <- '2015'
```
This will generate the country's urban/rural indicator tif map file and urban/rural probability tif map file.

---

### **8. Compute Urban-Rural Fractions for different admin level and subpopluation**
Modify the following variables in `step3_get_fraction.R` according to your chosen survey and then run the script:

```r
country <- 'Malawi'
survey_year <- '2015'
```
The results will be saved under `Results/Malawi_2015/UR_fraction` with separate folders for different admin level:
<img width="324" alt="image" src="https://github.com/user-attachments/assets/9eacabb1-63a6-4fd4-973e-a43936a1e5aa" />


Each folder will contain 5 kinds of population:

- k0_1: childern aged 0-1
- k0_5: childern aged 0-5
- f15_49: women aged 15-49 
- m15_49: men aged 15-49
- total_pop: total population 



---

### Now we start to determine if we need to do the stratification and how to do it

The basic idea is that if the survey has oversampling issue and the indicator has urban-rural association, then it will introduce bias when we apply the cluster-level model.
Therefore we need the stratification model to address the bias.

---
### **9. Test for Oversampling Bias**

Firt we use chi-square test with Monte Carlo method to test if there is oversampling issue in the survey.

Modify the following variables in `strat_step1_survey_test.R` according to your chosen survey and then run the script:
```r
country <- 'Malawi'
country.abbrev <- 'mwi'
survey_year <- '2015'
```
Check the output in `Results/Malawi_2015/UR_stratification/`

This will produce a histogram and statistical tests for oversampling bias.

---

### **10. Indicator's Urban-Rural Association Testing**

Next, we conduct the Likelihood Ratio test to see if the indicator we are interested has u/r association.

Modify the following variables in `strat_step2_indicator_test.R` according to your chosen survey and then run the script:

```r
country <- 'Malawi'
country.abbrev <- 'mwi'
survey_year <- '2015'
indicator <- "HA_HIVP_B_HIV"
```
Check the output in `Results/Malawi_2015/UR_stratification/Indicators/HA_HIVP_B_HIV/`

Meanwhile, the console output will show you whether the indicator has urban/rural association.

---

### **11. Validate Stratification Procedure**

Modify the following variables in `strat_step3_check_procedure.R` according to your chosen survey and then run the script:

```r
country <- 'Malawi'
country.abbrev <- 'mwi'
gadm.abbrev <- 'mwi'
survey_year <- '2015'
indicator <- "HA_HIVP_B_HIV"
sub_pop <- "f15_49"
```

Make sure `sub_pop` is consistent with the definition of your indicator. For example, the `HA_HIVP_B_HIV` is referring to the female aged 15-49 according to the DHS report.

Valid `sub_pop` values:
- "k0_5"
- "total_pop"
- "k0_1"
- "m15_49"
- "f15_19"


This will generate a plot comparing the stratified and unstratified models against direct estimates at the admin1 level. Additionally, it will calculate the weighted absolute bias to quantify the differences between the (un)stratified models and the direct estimate.

---

### **12. Stratification at Admin Levels**
Modify the following variables in `strat_step4_stratifcation_admin2:3.R` according to your chosen survey and then run the script:

```r
country <- 'Malawi'
country.abbrev <- 'mwi'
gadm.abbrev <- 'mwi'
survey_year <- '2015'
admin.level.tmp <- 2
sub_pop <- "f15_49"
indicator <- "HA_HIVP_B_HIV"
```
Adjust `admin.level.tmp` for different administrative levels which you want to model at, such as, `1`. But please ensure it does not exceed the country's maximum admin level.

This will generate a plot comparing the stratified against unstratified models with mean and standard deviation.

