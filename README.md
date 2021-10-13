# HAV-outbreak-Louisville

The R code in this repository reproduces the results in the paper:

Emmanuelle A. Dankwa, Christl A. Donnelly, Andrew F. Brouwer, Rui Zhao, Martha Montgomery, Mark K. Weng, Natasha K. Martin. 
Estimating vaccination threshold and impact in the 2017-2019 hepatitis A virus outbreak among persons experiencing homelessness or who use drugs in Louisville, Kentucky, United States

This repository has three folders for the three main analyses sections: 

* Displaying surveillance data: `surveillance` folder
* Comparison of relative risks: `risk_analysis` folder
* Hepatitis A virus (HAV) transmission model: `model` folder. 

These have been structured such that each analysis can be run independently of the others. The contents of each folder and guidance on reproducing the results are outlined below. 


# Displaying surveillance data

The `surveillance` folder contains: 

* `surveillance_main.Rmd`: This is the main script for this section. It reproduces the results in Figure 2 and Table 2.

* `surveillance_data`: folder containing the following data files;

    + `dat2.RDS`: A line list of all cases with the following nine variables - 
         * `MMWRYear`: year of detection
         * `MMWRWeek`: week of detection
         * `IVDU`: intravenous drug use status (levels: "Yes", "No", "Unknown")
         * `Non.IV.DU`: non-intravenous drug use status  (levels: "Yes", "No", "Unknown")
         * `Homelessness`: housing status (levels: "Not Homeless", "Unstable Housing", "Homeless - Shelter/Streets", "Contact with Homeless-Social".)
         * `Hospitalization`: hospitalization status (levels: "Yes", "No", "Unknown")
         * `Mortality`: mortality (levels: "Yes", "No")
         * `Sex`: sex (levels: "Male", "Female")
         * `Age_group`: age group (levels: "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70+")
         
    + `dat2_yrwk.RDS`: case counts by year, week and risk group 
    + `target.RDS`: a line list of cases among persons experiencing homelessness or who use drugs (PEH/PWUD) with the same variables as in `dat2.RDS`
    + `target_yrwk.RDS`: case counts by year and week among PEH/PWUD
    + `vacc_counts_target.RDS`: weekly vaccinations administered to PEH/PWUD
    
* `surveillance_functions.R`: functions for producing the results in this section
    


# Comparison of relative risks 

The `risk_analysis` folder contains: 

* `risk_analysis_main.Rmd`: This is the main script for this section. It reproduces the results in Figure 3 (main text), Supplementary tables S1-S3 and Supplementary Figure S1.

* `risk_analysis_data`: This folder contains the data files `dat2.RDS`, `dat2_yrwk.RDS`, `target.RDS`, and `target_yrwk.RDS`. File descriptions are as in the previous section.  

* `risk_analysis_functions.R`: functions for producing the results in this section
 



# HAV transmission model

The `model` folder contains:


 * `model_main.Rmd`: This is the main script for this section, with the code for 
 
    + parameter estimation, 
    + estimating the basic reproduction number, critical vaccination coverage level and the herd immunity threshold, and 
    + sensitivity analysis. 

    
*  `model_data`: folder containing the data required to run the model:

    + `casecountperwk`: observed weekly case count among persons who experience homelessness or who use drugs (PEH/PWUD).
    
    + `vaccdata`: weekly hepatitis A vaccination counts among PEH/PWUD.
    
* `model_functions.R`: functions for producing the results in this section 
        
 
* `results`:  All relevant outputs of `model_main.Rmd` will be stored in this folder. 
 
 


# Running the scripts

To run the scripts, 

1) Download this repository to a desired location on your computer (by clicking the green "Code" button above and then "Download Zip" or other preferred option). 

2) Open the project in the downloaded folder (by double-clicking the file "HAV-outbreak-Louisville.Rproj") 

3) As each section is structured as a stand-alone, you would need to first set the directory to the location of the folder for the section you desire to replicate. Example, for the transmission model, run `setwd(~\HAV-outbreak-model\model)`.

4) Run the main Rmd file in the folder (by clicking the "Knit" button). NOTE: As some computations may require long runtimes, the user is advised to select relevant code chunks and run these separately rather than compiling the file in one go.



# Queries

Queries and suggestions are always welcome. Please email Emmanuelle at dankwa@stats.ox.ac.uk with such. 


# License 

MIT License
