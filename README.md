# Absolute Risk Modeling and Conditional Average Derivative Estimator (CADE)
This repository contains the R code for calculating absolute risk of disease incidence and mortality and conditional average derivative estimator in the All of Us cohort. All data were taken from the All of Us research cohort version 7. <a href="https://www.medrxiv.org/content/10.64898/2026.06.03.26354842v2">Read the paper on medRxiv.</a>

<b>The file structure is as follows:</b>
1. Data Set-up contains the code for setting up the datasets for fitting the absolute risk model.
2. Model Training and Calibration in All of Us uses the training datasets to calculate the baseline risk of disease and set up model calibration calculations.
3. Projecting Absolute Risk of Cancer Incidence calculates the absolute risk of cancer incidence for the validation dataset.
4. Projecting Absolute Risk of Cancer Mortality calculates the 10-year disease mortality risk.
5. Code for Tables Figures and CADE contains all the code for generating the tables and figures, and an example on calculating the CADE for disease incidence risk.
