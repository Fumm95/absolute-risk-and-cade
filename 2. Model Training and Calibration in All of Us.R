################################################################
####### 2. Model Training and Calibration in All of Us #########
################################################################

# Setting up the training datasets

## Filter to male only

ehr_male <- dataset_33309842_person_df %>% filter(sex_at_birth == "Male")

## Dividing out hispanic ancestry to create new ancestry groups

ehr_male_ancestry = ehr_male %>% mutate(hisp_race = ifelse(ethnicity == "Hispanic or Latino", "Hispanic or Latino", race)) %>%
  mutate(hisp_race = recode(hisp_race, "White" = "Non Hispanic White"))

brca_male_prs <- merge(ehr_male_ancestry, prs, by.x = "person_id", by.y = "IID")

brca_male_total = merge(brca_male_prs, pca_df, by.x = "person_id", by.y = "IID")

## Creating Training Datasets

### EUR Female

brca_white_inc = brca_female_inc %>% filter(hisp_race == "Non Hispanic White")

brca_train_ids <- sample(c(1:nrow(brca_white_inc)), size = 20000)

brca_white_train <- brca_white_inc[brca_train_ids,]
brca_white_test <- brca_white_inc[-brca_train_ids,]

### AFR Female

brca_black_inc = brca_female_inc %>% filter(hisp_race == "Black or African American")

black_train_ids <- sample(c(1:nrow(brca_black_inc)), size = 8000)

brca_black_train <- brca_black_inc[black_train_ids,]
brca_black_test <- brca_black_inc[-black_train_ids,]

### Hisp Female

brca_hisp_inc = brca_female_inc %>% filter(hisp_race == "Hispanic or Latino")

hisp_train_ids <- sample(c(1:nrow(brca_hisp_inc)), size = 8000)

brca_hisp_train <- brca_hisp_inc[hisp_train_ids,]
brca_hisp_test <- brca_hisp_inc[-hisp_train_ids,]

### ASN Female

brca_asn_inc = brca_female_inc %>% filter(hisp_race == "Asian")

asn_train_ids <- sample(c(1:nrow(brca_asn_inc)), size = 1000)

brca_asn_train <- brca_asn_inc[asn_train_ids,]
brca_asn_test <- brca_asn_inc[-asn_train_ids,]

### EUR Male

brca_male_eur = brca_male_total %>% filter(hisp_race == "Non Hispanic White")

brca_male_eur$stdprs <- brca_male_eur$PRS / sd(brca_male_eur$PRS)

male_train_ids <- sample(c(1:nrow(brca_male_eur)), size = 20000)

male_white_train <- brca_male_eur[male_train_ids,]
male_white_test <- brca_male_eur[-male_train_ids,]

### AFR Male

brca_male_afr = brca_male_total %>% filter(hisp_race == "Black or African American")

brca_male_afr$stdprs <- brca_male_afr$PRS / sd(brca_male_afr$PRS)

male_train_ids <- sample(c(1:nrow(brca_male_afr)), size = 8000)

male_black_train <- brca_male_afr[male_train_ids,]
male_black_test <- brca_male_afr[-male_train_ids,]

### Hisp Male

brca_male_hisp = brca_male_total %>% filter(hisp_race == "Hispanic or Latino")

brca_male_hisp$stdprs <- brca_male_hisp$PRS / sd(brca_male_hisp$PRS)

male_train_ids <- sample(c(1:nrow(brca_male_hisp)), size = 6000)

male_hisp_train <- brca_male_hisp[male_train_ids,]
male_hisp_test <- brca_male_hisp[-male_train_ids,]

### ASN Male

brca_male_asn = brca_male_total %>% filter(hisp_race == "Asian")

brca_male_asn$stdprs <- brca_male_asn$PRS / sd(brca_male_asn$PRS)

male_train_ids <- sample(c(1:nrow(brca_male_asn)), size = 1000)

male_asian_train <- brca_male_asn[male_train_ids,]
male_asian_test <- brca_male_asn[-male_train_ids,]

## Setting up female training data

brca_female_train = rbind(brca_white_train, brca_black_train, brca_hisp_train, brca_asn_train)

brca_female_train = brca_female_train %>% mutate(hisp_race = factor(hisp_race, levels = c("Non Hispanic White", "Asian",
                                                                                          "Black or African American", "Hispanic or Latino")))

brca_female_train$stdprs = brca_female_train$PRS / sd(brca_female_train$PRS)

standardize_df = brca_female_train %>% select(stdprs, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)

standardized_df = scale(standardize_df)


## Modeling for PC effect sizes in training dataset

fit.alpha = glm(total_case ~ (hisp_race) * stdprs + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +  education_answer + income_answer,
                data = brca_female_train, family = binomial(link = "logit"))


# Modeling PRS Distribution in Male Training Data

## European

###  Use Male Train for finding relationship between PRS and PCs

prsfiteur <- lm(stdprs ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = male_white_train)

residvareur <- as.vector(prsfiteur$residuals^2)

varfiteur <- glm(residvareur ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = male_white_train, family = Gamma(link = "log"))

### Use Male Test to calculate alphas for calibration

male_white_test$meanscale <- cbind(1, as.matrix(male_white_test[, c("PC1", "PC2", "PC3","PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(prsfiteur)[c(1:11)]

male_white_test$varscale <- exp(cbind(1, as.matrix(male_white_test[, c("PC1", "PC2", "PC3","PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(varfiteur)[c(1:11)])

male_white_test$alphapc <- cbind(as.matrix(male_white_test[, c("PC1", "PC4", "PC6")])) %*% coefficients(fit.alpha)[c(6, 9, 11)]

eurmean <- mean(male_white_test$meanscale)
eurvar <- mean(male_white_test$varscale)
eurpcs <- mean(male_white_test$alphapc)

## African

### Use Male Train for finding relationship between PRS and PCs

prsfitafr <- lm(stdprs ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = male_black_train)

residvarafr <- as.vector(prsfitafr$residuals^2)

varfitafr <- glm(residvarafr ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = male_black_train, family = Gamma(link = "log"))

### Use Male Test to calculate alphas for calibration

male_black_test$meanscale <- cbind(1, as.matrix(male_black_test[, c("PC1", "PC2", "PC3","PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(prsfitafr)[c(1:11)]

male_black_test$varscale <- exp(cbind(1, as.matrix(male_black_test[, c("PC1", "PC2", "PC3","PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(varfitafr)[c(1:11)])

male_black_test$alphapc <- cbind(as.matrix(male_black_test[, c("PC1", "PC4", "PC6")])) %*% coefficients(fit.alpha)[c(6, 9, 11)]

afrmean <- mean(male_black_test$meanscale)
afrvar <- mean(male_black_test$varscale)
afrpcs <- mean(male_black_test$alphapc)


## Hispanic

### Use Male Train for finding relationship between PRS and PCs

prsfithisp <- lm(stdprs ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = male_hisp_train)

residvarhisp <- as.vector(prsfithisp$residuals^2)

varfithisp <- glm(residvarhisp ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = male_hisp_train, family = Gamma(link = "log"))

### Use Male Test to calculate alphas for calibration

male_hisp_test$meanscale <- cbind(1, as.matrix(male_hisp_test[, c("PC1", "PC2", "PC3","PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(prsfithisp)[c(1:11)]

male_hisp_test$varscale <- exp(cbind(1, as.matrix(male_hisp_test[, c("PC1", "PC2", "PC3","PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(varfithisp)[c(1:11)])

male_hisp_test$alphapc <- cbind(as.matrix(male_hisp_test[, c("PC1", "PC4", "PC6")])) %*% coefficients(fit.alpha)[c(6, 9, 11)]

hispmean <- mean(male_hisp_test$meanscale)
hispvar <- mean(male_hisp_test$varscale)
hisppcs <- mean(male_hisp_test$alphapc)


## Asian

### Use Male Train for finding relationship between PRS and PCs

prsfitasn <- lm(stdprs ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = male_asian_train)

residvarasn <- as.vector(prsfitasn$residuals^2)

varfitasn <- glm(residvarasn ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = male_asian_train, family = Gamma(link = "log"))

### Use Male Test to calculate alphas for calibration

male_asian_test$meanscale <- cbind(1, as.matrix(male_asian_test[, c("PC1", "PC2", "PC3","PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(prsfitasn)[c(1:11)]

male_asian_test$varscale <- exp(cbind(1, as.matrix(male_asian_test[, c("PC1", "PC2", "PC3","PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(varfitasn)[c(1:11)])

male_asian_test$alphapc <- cbind(as.matrix(male_asian_test[, c("PC1", "PC4", "PC6")])) %*% coefficients(fit.alpha)[c(6, 9, 11)]

asnmean <- mean(male_asian_test$meanscale)
asnvar <- mean(male_asian_test$varscale)
asnpcs <- mean(male_asian_test$alphapc)
