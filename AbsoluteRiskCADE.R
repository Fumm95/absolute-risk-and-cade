

# Filter to male only

ehr_male <- dataset_33309842_person_df %>% filter(sex_at_birth == "Male")

# Dividing out hispanic ancestry to create new ancestry groups

ehr_male_ancestry = ehr_male %>% mutate(hisp_race = ifelse(ethnicity == "Hispanic or Latino", "Hispanic or Latino", race)) %>%
  mutate(hisp_race = recode(hisp_race, "White" = "Non Hispanic White"))

brca_male_prs <- merge(ehr_male_ancestry, prs, by.x = "person_id", by.y = "IID")

brca_male_total = merge(brca_male_prs, pca_df, by.x = "person_id", by.y = "IID")

# EUR Female

brca_white_inc = brca_female_inc %>% filter(hisp_race == "Non Hispanic White")

brca_train_ids <- sample(c(1:nrow(brca_white_inc)), size = 20000)

brca_white_train <- brca_white_inc[brca_train_ids,]
brca_white_test <- brca_white_inc[-brca_train_ids,]

# AFR Female

brca_black_inc = brca_female_inc %>% filter(hisp_race == "Black or African American")

black_train_ids <- sample(c(1:nrow(brca_black_inc)), size = 8000)

brca_black_train <- brca_black_inc[black_train_ids,]
brca_black_test <- brca_black_inc[-black_train_ids,]

# Hisp Female

brca_hisp_inc = brca_female_inc %>% filter(hisp_race == "Hispanic or Latino")

hisp_train_ids <- sample(c(1:nrow(brca_hisp_inc)), size = 8000)

brca_hisp_train <- brca_hisp_inc[hisp_train_ids,]
brca_hisp_test <- brca_hisp_inc[-hisp_train_ids,]

# ASN Female

brca_asn_inc = brca_female_inc %>% filter(hisp_race == "Asian")

asn_train_ids <- sample(c(1:nrow(brca_asn_inc)), size = 1000)

brca_asn_train <- brca_asn_inc[asn_train_ids,]
brca_asn_test <- brca_asn_inc[-asn_train_ids,]

# EUR Male

brca_male_eur = brca_male_total %>% filter(hisp_race == "Non Hispanic White")

brca_male_eur$stdprs <- brca_male_eur$PRS / sd(brca_male_eur$PRS)

male_train_ids <- sample(c(1:nrow(brca_male_eur)), size = 20000)

male_white_train <- brca_male_eur[male_train_ids,]
male_white_test <- brca_male_eur[-male_train_ids,]

# AFR Male

brca_male_afr = brca_male_total %>% filter(hisp_race == "Black or African American")

brca_male_afr$stdprs <- brca_male_afr$PRS / sd(brca_male_afr$PRS)

male_train_ids <- sample(c(1:nrow(brca_male_afr)), size = 8000)

male_black_train <- brca_male_afr[male_train_ids,]
male_black_test <- brca_male_afr[-male_train_ids,]

# Hisp Male

brca_male_hisp = brca_male_total %>% filter(hisp_race == "Hispanic or Latino")

brca_male_hisp$stdprs <- brca_male_hisp$PRS / sd(brca_male_hisp$PRS)

male_train_ids <- sample(c(1:nrow(brca_male_hisp)), size = 6000)

male_hisp_train <- brca_male_hisp[male_train_ids,]
male_hisp_test <- brca_male_hisp[-male_train_ids,]

# ASN Male

brca_male_asn = brca_male_total %>% filter(hisp_race == "Asian")

brca_male_asn$stdprs <- brca_male_asn$PRS / sd(brca_male_asn$PRS)

male_train_ids <- sample(c(1:nrow(brca_male_asn)), size = 1000)

male_asian_train <- brca_male_asn[male_train_ids,]
male_asian_test <- brca_male_asn[-male_train_ids,]

brca_female_train = rbind(brca_white_train, brca_black_train, brca_hisp_train, brca_asn_train)

brca_female_train = brca_female_train %>% mutate(hisp_race = factor(hisp_race, levels = c("Non Hispanic White", "Asian",
                                                                                          "Black or African American", "Hispanic or Latino")))

brca_female_train$stdprs = brca_female_train$PRS / sd(brca_female_train$PRS)

standardize_df = brca_female_train %>% select(stdprs, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)

standardized_df = scale(standardize_df)

fit.alpha = glm(total_case ~ (hisp_race) * stdprs + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +  education_answer + income_answer,
                data = brca_female_train, family = binomial(link = "logit"))


# Subgroup PRS Mean and Variance

## Use Male Train for finding relationship between PRS and PCs

prsfiteur <- lm(stdprs ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = male_white_train)

residvareur <- as.vector(prsfiteur$residuals^2)

varfiteur <- glm(residvareur ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = male_white_train, family = Gamma(link = "log"))

## Use Male Test for calculating alphas

male_white_test$meanscale <- cbind(1, as.matrix(male_white_test[, c("PC1", "PC2", "PC3","PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(prsfiteur)[c(1:11)]

male_white_test$varscale <- exp(cbind(1, as.matrix(male_white_test[, c("PC1", "PC2", "PC3","PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(varfiteur)[c(1:11)])

male_white_test$alphapc <- cbind(as.matrix(male_white_test[, c("PC1", "PC4", "PC6")])) %*% coefficients(fit.alpha)[c(6, 9, 11)]

eurmean <- mean(male_white_test$meanscale)
eurvar <- mean(male_white_test$varscale)
eurpcs <- mean(male_white_test$alphapc)


## Use Male Train for finding relationship between PRS and PCs

prsfitafr <- lm(stdprs ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = male_black_train)

residvarafr <- as.vector(prsfitafr$residuals^2)

varfitafr <- glm(residvarafr ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = male_black_train, family = Gamma(link = "log"))

## Use Male Test for calculating alphas

male_black_test$meanscale <- cbind(1, as.matrix(male_black_test[, c("PC1", "PC2", "PC3","PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(prsfitafr)[c(1:11)]

male_black_test$varscale <- exp(cbind(1, as.matrix(male_black_test[, c("PC1", "PC2", "PC3","PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(varfitafr)[c(1:11)])

male_black_test$alphapc <- cbind(as.matrix(male_black_test[, c("PC1", "PC4", "PC6")])) %*% coefficients(fit.alpha)[c(6, 9, 11)]

afrmean <- mean(male_black_test$meanscale)
afrvar <- mean(male_black_test$varscale)
afrpcs <- mean(male_black_test$alphapc)



## Use Male Train for finding relationship between PRS and PCs

prsfithisp <- lm(stdprs ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = male_hisp_train)

residvarhisp <- as.vector(prsfithisp$residuals^2)

varfithisp <- glm(residvarhisp ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = male_hisp_train, family = Gamma(link = "log"))

## Use Male Test for calculating alphas

male_hisp_test$meanscale <- cbind(1, as.matrix(male_hisp_test[, c("PC1", "PC2", "PC3","PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(prsfithisp)[c(1:11)]

male_hisp_test$varscale <- exp(cbind(1, as.matrix(male_hisp_test[, c("PC1", "PC2", "PC3","PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(varfithisp)[c(1:11)])

male_hisp_test$alphapc <- cbind(as.matrix(male_hisp_test[, c("PC1", "PC4", "PC6")])) %*% coefficients(fit.alpha)[c(6, 9, 11)]

hispmean <- mean(male_hisp_test$meanscale)
hispvar <- mean(male_hisp_test$varscale)
hisppcs <- mean(male_hisp_test$alphapc)


## Use Male Train for finding relationship between PRS and PCs

prsfitasn <- lm(stdprs ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = male_asian_train)

residvarasn <- as.vector(prsfitasn$residuals^2)

varfitasn <- glm(residvarasn ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = male_asian_train, family = Gamma(link = "log"))

## Use Male Test for calculating alphas

male_asian_test$meanscale <- cbind(1, as.matrix(male_asian_test[, c("PC1", "PC2", "PC3","PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(prsfitasn)[c(1:11)]

male_asian_test$varscale <- exp(cbind(1, as.matrix(male_asian_test[, c("PC1", "PC2", "PC3","PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(varfitasn)[c(1:11)])

male_asian_test$alphapc <- cbind(as.matrix(male_asian_test[, c("PC1", "PC4", "PC6")])) %*% coefficients(fit.alpha)[c(6, 9, 11)]

asnmean <- mean(male_asian_test$meanscale)
asnvar <- mean(male_asian_test$varscale)
asnpcs <- mean(male_asian_test$alphapc)

eurbeta <- log(1.61)

hispbeta <- log(1.3) #Liu et al.

afrbeta <- log(1.27) #Du et al.

asnbeta <- log(1.52) #Ho et al.

afrbeta <- afrbeta * sd(brca_male_afr$stdprs) / sd(brca_male_eur$stdprs)
hispbeta <- hispbeta * sd(brca_male_hisp$stdprs) / sd(brca_male_eur$stdprs)
asnbeta <- asnbeta * sd(brca_male_asn$stdprs) / sd(brca_male_eur$stdprs)

#EUR

bc_white_inc = bc_incidence %>% filter(Race == "White")

brca_white_test$bc_inc <- bc_white_inc$Rate[floor(brca_white_test$age/5) + 1]
brca_white_test$alpha <- log(brca_white_test$bc_inc) - eurbeta * eurmean - 0.5 * eurbeta^2 * eurvar - eurpcs
brca_white_test$alphapc <- cbind(as.matrix(brca_white_test[, c("PC1", "PC4", "PC6")])) %*% coefficients(fit.alpha)[c(6, 9, 11)]

brca_white_test$meanscale <- cbind(1, as.matrix(brca_white_test[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(prsfiteur)[c(1:11)]
brca_white_test$varscale <- exp(cbind(1, as.matrix(brca_white_test[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(varfiteur)[c(1:11)])

brca_white_test$stdprs <- brca_white_test$PRS / sd(brca_white_test$PRS)

brca_white_test$abs_prob <- exp(brca_white_test$alpha + brca_white_test$alphapc + eurbeta * brca_white_test$stdprs)
brca_white_test$exp_prob <- exp(brca_white_test$alpha + brca_white_test$alphapc + eurbeta * brca_white_test$stdprs + 0.5 * eurbeta^2 * brca_white_test$varscale)

brca_white_test$scaleprs <- (brca_white_test$stdprs - brca_white_test$meanscale) / sqrt(brca_white_test$varscale)
brca_white_test$arprs <- brca_white_test$stdprs - brca_white_test$meanscale - 0.5 * brca_white_test$varscale

#AFR

bc_black_inc = bc_incidence %>% filter(Race == "Black")

brca_black_test$bc_inc <- bc_black_inc$Rate[floor(brca_black_test$age/5) + 1]
brca_black_test$alpha <- log(brca_black_test$bc_inc) - afrbeta * afrmean - 0.5 * afrbeta^2 * afrvar - afrpcs
brca_black_test$alphapc <- cbind(as.matrix(brca_black_test[, c("PC1", "PC4", "PC6")])) %*% coefficients(fit.alpha)[c(6, 9, 11)]

brca_black_test$meanscale <- cbind(1, as.matrix(brca_black_test[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(prsfitafr)[c(1:11)]
brca_black_test$varscale <- exp(cbind(1, as.matrix(brca_black_test[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(varfitafr)[c(1:11)])

brca_black_test$stdprs <- brca_black_test$PRS / sd(brca_black_test$PRS)

brca_black_test$abs_prob <- exp(brca_black_test$alpha + brca_black_test$alphapc + afrbeta * brca_black_test$stdprs)
brca_black_test$exp_prob <- exp(brca_black_test$alpha + brca_black_test$alphapc + afrbeta * brca_black_test$stdprs + 0.5 * afrbeta^2 * brca_black_test$varscale)

brca_black_test$scaleprs <- (brca_black_test$stdprs - brca_black_test$meanscale) / sqrt(brca_black_test$varscale)
brca_black_test$arprs <- brca_black_test$stdprs - brca_black_test$meanscale - 0.5 * brca_black_test$varscale

#AMR

bc_hisp_inc = bc_incidence %>% filter(Race == "Hispanic")

brca_hisp_test$bc_inc <- bc_hisp_inc$Rate[floor(brca_hisp_test$age/5) + 1]
brca_hisp_test$alpha <- log(brca_hisp_test$bc_inc) - hispbeta * hispmean - 0.5 * hispbeta^2 * hispvar - hisppcs
brca_hisp_test$alphapc <- cbind(as.matrix(brca_hisp_test[, c("PC1", "PC4", "PC6")])) %*% coefficients(fit.alpha)[c(6, 9, 11)]

brca_hisp_test$meanscale <- cbind(1, as.matrix(brca_hisp_test[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(prsfithisp)[c(1:11)]
brca_hisp_test$varscale <- exp(cbind(1, as.matrix(brca_hisp_test[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(varfithisp)[c(1:11)])

brca_hisp_test$stdprs <- brca_hisp_test$PRS / sd(brca_hisp_test$PRS)

brca_hisp_test$abs_prob <- exp(brca_hisp_test$alpha + brca_hisp_test$alphapc + hispbeta * brca_hisp_test$stdprs)
brca_hisp_test$exp_prob <- exp(brca_hisp_test$alpha + brca_hisp_test$alphapc + hispbeta * brca_hisp_test$stdprs + 0.5 * hispbeta^2 * brca_hisp_test$varscale)

brca_hisp_test$scaleprs <- (brca_hisp_test$stdprs - brca_hisp_test$meanscale) / sqrt(brca_hisp_test$varscale)
brca_hisp_test$arprs <- brca_hisp_test$stdprs - brca_hisp_test$meanscale - 0.5 * brca_hisp_test$varscale

#ASN

bc_asn_inc = bc_incidence %>% filter(Race == "Asian")

brca_asn_test$bc_inc <- bc_asn_inc$Rate[floor(brca_asn_test$age/5) + 1]
brca_asn_test$alpha <- log(brca_asn_test$bc_inc) - asnbeta * asnmean - 0.5 * asnbeta^2 * asnvar - asnpcs
brca_asn_test$alphapc <- cbind(as.matrix(brca_asn_test[, c("PC1", "PC4", "PC6")])) %*% coefficients(fit.alpha)[c(6, 9, 11)]

brca_asn_test$meanscale <- cbind(1, as.matrix(brca_asn_test[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(prsfitasn)[c(1:11)]
brca_asn_test$varscale <- exp(cbind(1, as.matrix(brca_asn_test[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(varfitasn)[c(1:11)])

brca_asn_test$stdprs <- brca_asn_test$PRS / sd(brca_asn_test$PRS)

brca_asn_test$abs_prob <- exp(brca_asn_test$alpha + brca_asn_test$alphapc + asnbeta * brca_asn_test$stdprs)
brca_asn_test$exp_prob <- exp(brca_asn_test$alpha + brca_asn_test$alphapc + asnbeta * brca_asn_test$stdprs + 0.5 * asnbeta^2 * brca_asn_test$varscale)

brca_asn_test$scaleprs <- (brca_asn_test$stdprs - brca_asn_test$meanscale) / sqrt(brca_asn_test$varscale)
brca_asn_test$arprs <- brca_asn_test$stdprs - brca_asn_test$meanscale - 0.5 * brca_asn_test$varscale

aoutest <- rbind(brca_white_test, brca_black_test, brca_hisp_test, brca_asn_test)


# Calculating CADE Example

aoutest_white = aoutest %>% filter(hisp_race == "Non Hispanic White" & age >= 55 & age < 60)
cade_eur = eurbeta * mean(aoutest_white$abs_prob, na.rm = T)

aoutest_black = aoutest %>% filter(hisp_race == "Black or African American" & age >= 55 & age < 60)
cade_afr = afrbeta * mean(aoutest_black$abs_prob, na.rm = T)

aoutest_hisp = aoutest %>% filter(hisp_race == "Hispanic or Latino" & age >= 55 & age < 60)
cade_hisp = hispbeta * mean(aoutest_hisp$abs_prob, na.rm = T)

aoutest_asian = aoutest %>% filter(hisp_race == "Asian" & age >= 55 & age < 60)
cade_asn = asnbeta * mean(aoutest_asian$abs_prob, na.rm = T)

