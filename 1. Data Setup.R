################################
####### 1. Data Set Up #########
################################


# Note that we use breast cancer data obtained from All of Us to illustrate the code.

library(tidyverse)
library(data.table)
library(bigrquery)

# Create main dataset with cases from observation and condition

observation_ids = unique(dataset_33309842_observation_df$person_id)
condition_ids = unique(dataset_33309842_condition_df$person_id)
brca_outcome <- dataset_33309842_person_df
brca_outcome$observation_case = as.numeric(brca_outcome$person_id %in% observation_ids)
brca_outcome$condition_case = as.numeric(brca_outcome$person_id %in% condition_ids)
brca_female = brca_outcome %>% filter(sex_at_birth == "Female")
brca_female$ehr_case = as.numeric((brca_female$condition_case + brca_female$observation_case) > 0)

# Approximate date of entry using date of completing The Basics survey

survey_unique = dataset_33309842_survey_df %>%
  arrange(person_id) %>%
  filter(survey == "The Basics") %>% 
  distinct(person_id, .keep_all=T)

basic_survey_date = survey_unique %>% select(person_id, survey_datetime)



# Identify baseline vs incident cases using survey data for observation

observation_ordered <- dataset_33309842_observation_df %>% arrange(person_id)
first_obs <- observation_ordered %>% distinct(person_id, .keep_all = T)
first_obs <- first_obs %>% mutate(observation_datetime = as.Date(observation_datetime))
obs_outcome = first_obs %>% select(person_id, observation_datetime)
observation_survey = merge(basic_survey_date, obs_outcome, by = "person_id") %>% 
  mutate(survey_datetime = as.Date(survey_datetime))
observation_survey$baseline_observation = as.numeric(observation_survey$observation_datetime < observation_survey$survey_datetime)
observation_ordered$observation_datetime = as.Date(observation_ordered$observation_datetime)
observation_female = observation_ordered %>% filter(person_id %in% brca_female$person_id)

age_df = brca_female %>% select(person_id, date_of_birth) %>% mutate(date_of_birth = as.Date(date_of_birth))
observation_female_age = merge(observation_female, age_df)
observation_female_age = observation_female_age %>% mutate(age_at_observation = as.numeric((observation_datetime - date_of_birth)/365))


# Identify baseline vs incident cases using survey data for condition

condition_ordered <- dataset_33309842_condition_df %>% arrange(person_id)
first_cond <- condition_ordered %>% distinct(person_id, .keep_all = T)
first_cond <- first_cond %>% 
  mutate(condition_start_datetime = as.Date(condition_start_datetime)) %>% 
  mutate(condition_end_datetime = as.Date(condition_end_datetime))
cond_outcome = first_cond %>% select(person_id, condition_start_datetime)
condition_survey = merge(basic_survey_date, cond_outcome, by = "person_id") %>% 
  mutate(survey_datetime = as.Date(survey_datetime))
condition_survey$baseline_condition = as.numeric(condition_survey$condition_start_datetime < condition_survey$survey_datetime)

# Adding observation and condition dates to main dataset 

obs_survey_base = observation_survey %>% select(person_id, baseline_observation)
brca_female_obs_base = left_join(brca_female, obs_survey_base, by = "person_id")
brca_female_obs_base = brca_female_obs_base %>% mutate(baseline_observation = replace_na(baseline_observation, 0))
cond_survey_base = condition_survey %>% select(person_id, baseline_condition)
brca_female_base = left_join(brca_female_obs_base, cond_survey_base, by = "person_id")
brca_female_base = brca_female_base %>% mutate(baseline_condition = replace_na(baseline_condition, 0))
brca_female_base$ehr_baseline = as.numeric((brca_female_base$baseline_observation + brca_female_base$baseline_condition) > 0)

# Dividing out hispanic ancestry to create new ancestry groups

brca_female_ancestry = brca_female_base %>% mutate(hisp_race = ifelse(ethnicity == "Hispanic or Latino", "Hispanic or Latino", race)) %>%
  mutate(hisp_race = recode(hisp_race, "White" = "Non Hispanic White"))

# Collecting in medical history survey data

brca_self = dataset_33309842_survey_df %>% filter(answer_concept_id == 43528500)
history_brca = dataset_33309842_survey_df %>% filter(question_concept_id == 836772)
self_id = unique(brca_self$person_id)
all_history_brca = unique(history_brca$person_id)

# Adding medical history to main database

brca_female_ancestry = brca_female_ancestry %>% mutate(personal_history = as.numeric(person_id %in% self_id)) %>%
  mutate(total_case = as.numeric(personal_history + ehr_case > 0)) %>%
  mutate(total_baseline = as.numeric(personal_history + ehr_baseline > 0))

saveRDS(brca_female_ancestry, "/home/jupyter/outcome_datasets/brca_female_ancestry.rds")




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

# Reading in the person dataset for All of Us

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

x# Reading in the survey dataset for All of Us

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

# Defining Potential Confounding Information
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

# Reading in Disease Incidence and Mortality Rates

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