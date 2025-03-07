# Urban-Rural Stratification Pipeline

In this document, we give an overview of the urban-rural stratification pipeline and suggest the order at which scripts are run. Note that our pipeline is not aimed to replicate the estimates presented in DHS final report. Due to distinctions in methods, our estimates will not be exactly the same as DHS estimates, but they will be consistent.

This repository contains scripts for conducting urban-rural stratification using DHS survey data. The pipeline is designed to process subpopulation data, generate urban-rural fractions, and evaluate stratification needs. Please download the the pipline folder from [UR-Pipeline-Demo](https://uwnetid-my.sharepoint.com/:f:/g/personal/jitongj_uw_edu/EplTkYWv0eFCjjwmH7PW-qwBiWzINyB8VIr6VYBMcHPdQQ?e=7wLIAi) and the **codes** folder in this repository.

**Note for Mac user**: After the **UR-Pipeline-Demo** is fully downloaded, you may encounter an error saying: `Unable to expand "UR-Pipeline-Demo.zip". (Error 79 - Inappropriate file type or format.)` This is due to the unzip limit in Mac system for lage file. To solve this there are two ways:

  1. First solution:
  Instead of downloading the entire folder at once,  download only the .tif files you need for your specific survey years (the ones you set up        later in the steps section) and then download the remaining folders separately. Then reconstruct the folder structure.

  2. Second solution:
     - Go to the terminal, then cd to the zip file location (To be safe, do not put it in the download path).
     - Enter the code: `zip -FF "UR-Pipeline-Demo.zip" --out fixed.zip`
     - Enter the code: `unzip fixed.zip -d UR-Pipeline-Demo-Extracted`
     - Rename the `UR-Pipeline-Demo-Extracted` to `UR-Pipeline-Demo`

**Note for Windows user**: you can also try the **first solution** if you want to save time from downloading UR-Pipeline-Demo


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
strat_level_GADM <- 1 # NO NEED TO CHANGE
svy_training_year <- c(2010,2015)
```
Here `frame_year` is the census year which the DHS survey this based on,  `calibrated_year` is the year which DHS survey's u/r population fraction is based on (See section 6 for more details). `strat_level_GADM` is the admin level which the sampling fram is based on, typical equal to `1`. `svy_training_year` is referring to the DHS surveys who were conducted based on the same census survey. 
For example, Malawi 2015 survey and 2010 survey are all based on Malawi 2008 census survey, so we combine them together to increase the training sample size for classification modeling.

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

#### **6.1 Create New Script**
Now in this step you need to **open a new empty script with cleaned environment** and load the data `country_shp_analysis.rds` under `Data/Malawi_2015/shapeFiles_gadm/`.

The `country_shp_analysis.rds` file contains a list of administrative levels, each of which includes detailed information on the areas within that level. For Malawi, since it has 2 admin levels beyond National level, you will see the following:


<img width="654" alt="Screenshot 2025-02-12 at 10 44 11" src="https://github.com/user-attachments/assets/d63533ae-4cce-4195-8afc-9f63a719326a" />

Generally we look at the **Admin-1**, so by `View(country_shp_analysis[["Admin-1"]])`, you will see:

<img width="905" alt="Screenshot 2025-02-12 at 10 49 20" src="https://github.com/user-attachments/assets/ae2a8fae-96cd-452f-a61d-577d9ea31af4" />

In `country_shp_analysis[["Admin-1"]]`, focus on the `NAME_1` column, which lists the names of all areas in Admin-1. Pay close attention to the spelling and order of these names, as you will later use them as a baseline to verify your own created data.

Next, generate the following three data frames and save them as `.rds` files, using `country.abbrev` as the prefix. For example:

- `mwi_ref_tab` -> `mwi_ref_tab.rds` (contains urban fractions)
- `mwi_frame_ea` -> `mwi_frame_ea.rds` (EA clusters in census)
- `mwi_sample_ea` -> `mwi_sample_ea.rds` (EA clusters in DHS sampling survey)

(change `mwi` for your country's abbrev)


#### **6.2 Create `mwi_ref_tab`**
`mwi_ref_tab` has four columns:

![image](https://github.com/user-attachments/assets/dfb0b11f-5ad7-4a4e-b516-983ba3de2ba3)


- `strata.adm.name`: Obtained from `country_shp_analysis[["Admin-1"]]$NAME_1`. In other words, the areas's name and area's order in `mwi_ref_tab$strata.adm.name` should **exactly match** with `country_shp_analysis[["Admin-1"]]$NAME_1`.
- `strata.adm.num`: the corresponding order after you obtain `mwi_ref_tab$strata.adm.name`.
- `strata.adm.char`: it is `paste0("strata_adm_", mwi_ref_tab$strata.adm.num)`.
- `urb_frac`: it is derived from the urban-rural population table in the DHS survey, typically located in Appendix A of the survey reports. To be noticed, some country may not directly have this table in the DHS survey. If it happens, please refer to the corresponding census survey or estimate it based on urban and rural household numbers (you can calculate the u/r population based on the number of average u/r population per household in the census survey).  For greater accuracy, use the actual population numbers to calculate the fraction instead of the urban percentage provided in the table. **We highly recommend using your favorite LLM model to work this out.** For example, take a screenshot of the following figure and use the prompt 'Based on this table, can you extract all the regions inside, with their corresponding population of urban,rural, total? Give me the r code for creating the data frame using data.frame function.'
  ![image](https://github.com/user-attachments/assets/154dc67d-3a95-49b7-adb0-5a2e47ea65c1)
    - The year on which this fraction is based will serve as our calibration year, as it will be used to adjust the fractions calculated in subsequent steps.


#### **6.3 Create `mwi_frame_ea` and `mwi_sample_ea`**
For `mwi_frame_ea` and `mwi_sample_ea`, write codes to create the data frame in the following format:

![image](https://github.com/user-attachments/assets/9e934b66-5ad6-454b-804b-6a636f7e4423)

`mwi_frame_ea` and `mwi_sample_ea` both have 3 columns:

- `strata`: Obtained from `country_shp_analysis[["Admin-1"]]$NAME_1`. In other words, the areas's name and area's order in `strata` should **exactly match** with `country_shp_analysis[["Admin-1"]]$NAME_1`.
- `urban`, `rural`, `total`: You can find the number of EA cluster for census or DHS sampling survey in the Appendix A (eg. Nigeria 2018 survey report Table A.2 and Table A.3 in Appendix A)


Make sure that `mwi_ref_tab.rds`, `mwi_frame_ea.rds`, and `mwi_sample_ea.rds` are saved in the country_survey folder within `Data`. For example:

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
Modify the following variables in `strat_step4_stratifcation_admin2_3.R` according to your chosen survey and then run the script:

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

