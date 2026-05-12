library(tidyverse)
library(data.table)
library(bigrquery)

# Polygenic risk score PGS000004 for breast cancer calculated using PLINK from the ACAF files.
prs = fread("/home/jupyter/plink_prs/brca313_harmonise_prs.sscore")

brca_female_ancestry = readRDS("/home/jupyter/outcome_datasets/brca_female_ancestry.rds")

prs = prs %>% select(IID, SCORE1_SUM)
colnames(prs) = c("IID", "PRS")

brca_female_outcome = merge(brca_female_ancestry, prs, by.x = "person_id", by.y = "IID")

# PCs computed from the HapMap3 SNPs shared between 1000G, UKBB, and AoU Array Files
pcs <- fread("/home/jupyter/hapmap3_pca/pca_aou_1000G_weight.sscore")

pc_mat = pcs[,c(-1:-4)]

tot_alleles = as.vector(pcs[,4])

pc_mat_average = pc_mat / tot_alleles

pca_df = as.data.frame(cbind(pcs[,'IID'], pc_mat_average))

colnames(pca_df) = c("IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

brca_female_outcome = merge(brca_female_outcome, pca_df, by.x = "person_id", by.y = "IID")

# Reading in the person dataset

library(bigrquery)

# Read the data directly from Cloud Storage into memory.
# NOTE: Alternatively you can `gsutil -m cp {person_33309842_path}` to copy these files
#       to the Jupyter disk.
read_bq_export_from_workspace_bucket <- function(export_path) {
  col_types <- cols(gender = col_character(), race = col_character(), ethnicity = col_character(), sex_at_birth = col_character())
  bind_rows(
    map(system2('gsutil', args = c('ls', export_path), stdout = TRUE, stderr = TRUE),
        function(csv) {
          message(str_glue('Loading {csv}.'))
          chunk <- read_csv(pipe(str_glue('gsutil cat {csv}')), col_types = col_types, show_col_types = FALSE)
          if (is.null(col_types)) {
            col_types <- spec(chunk)
          }
          chunk
        }))
}
dataset_33309842_person_df <- read_bq_export_from_workspace_bucket("gs://fc-secure-9fd19e9b-dbea-4d98-8704-72be2dff6c18/bq_exports/fumm95@researchallofus.org/20240404/person_33309842/person_33309842_*.csv")

x# Reading in the survey dataset

# Read the data directly from Cloud Storage into memory.
# NOTE: Alternatively you can `gsutil -m cp {survey_33309842_path}` to copy these files
#       to the Jupyter disk.
read_bq_export_from_workspace_bucket <- function(export_path) {
  col_types <- cols(survey = col_character(), question = col_character(), answer = col_character(), survey_version_name = col_character())
  bind_rows(
    map(system2('gsutil', args = c('ls', export_path), stdout = TRUE, stderr = TRUE),
        function(csv) {
          message(str_glue('Loading {csv}.'))
          chunk <- read_csv(pipe(str_glue('gsutil cat {csv}')), col_types = col_types, show_col_types = FALSE)
          if (is.null(col_types)) {
            col_types <- spec(chunk)
          }
          chunk
        }))
}
dataset_33309842_survey_df <- read_bq_export_from_workspace_bucket("gs://fc-secure-9fd19e9b-dbea-4d98-8704-72be2dff6c18/bq_exports/fumm95@researchallofus.org/20240404/survey_33309842/survey_33309842_*.csv")

# Calculating the age of all participants

date_of_birth_df = dataset_33309842_person_df %>% select(person_id, date_of_birth)

survey_unique = dataset_33309842_survey_df %>%
  arrange(person_id) %>%
  filter(survey == "The Basics") %>% 
  distinct(person_id, .keep_all=T)

basic_survey_date = survey_unique %>% select(person_id, survey_datetime)

age_df = merge(date_of_birth_df, basic_survey_date, by = "person_id")

age_df = age_df %>% mutate(date_of_birth = as.Date(date_of_birth)) %>% 
  mutate(survey_datetime = as.Date(survey_datetime)) %>%
  mutate(age = as.numeric((survey_datetime - date_of_birth)/365))

# Reading in Potential Confounding Information
education_df = dataset_33309842_survey_df %>% filter(question_concept_id == 1585940)

education_df_small = education_df %>% select(person_id, answer)

colnames(education_df_small) = c("person_id", "education_answer")

income_df = dataset_33309842_survey_df %>% filter(question_concept_id == 1585375)

income_df_small = income_df %>% select(person_id, answer)

colnames(income_df_small) = c("person_id", "income_answer")

confounders_df = merge(education_df_small, income_df_small, by = "person_id")

# Adding confounders to the outcome dataframe

brca_female_confounders = merge(brca_female_outcome, confounders_df, by = "person_id")

# Add age to the confounders dataframe

keep_races = c("Asian", "Black or African American", "Hispanic or Latino", "Non Hispanic White")

brca_female_final = merge(brca_female_confounders, age_df, by = "person_id")

brca_female_inc = brca_female_final %>% filter(hisp_race %in% keep_races)

# Disease Incidence and Mortality Rates

bc_incidence = read.csv("/home/jupyter/seerstat_data/breast_cancer_incidence_race_age.csv")
colnames(bc_incidence) = c("Age", "Race", "Rate", "Count", "Population")

factor_race = as.factor(bc_incidence$Race)

factor_race = factor_race %>% recode_factor("0" = "White",
                                            "1" = "Black",
                                            "2" = "Native",
                                            "3" = "Asian",
                                            "4" = "Hispanic",
                                            "5" = "Unknown")

bc_incidence = bc_incidence %>% 
  mutate(Race = factor_race) %>%
  mutate(Rate = as.numeric(Rate)) %>%
  mutate(Rate = Rate / 100000)


bc_mortality = read.csv("/home/jupyter/seerstat_data/breast_cancer_mortality_race_age.csv")
colnames(bc_mortality) = c("Age", "Race", "Rate", "Count", "Population")

factor_race = as.factor(bc_mortality$Race)

factor_race = factor_race %>% recode_factor("0" = "White",
                                            "1" = "Black",
                                            "2" = "Native",
                                            "3" = "Asian",
                                            "4" = "Hispanic",
                                            "5" = "Unknown")

bc_mortality = bc_mortality %>% 
  mutate(Race = factor_race) %>%
  mutate(Rate = as.numeric(Rate)) %>%
  mutate(Rate = Rate / 100000)





