###################################################################
####### 3. Calculating Absolute Risk of Disease Incidence #########
###################################################################

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
brca_white_test$scaleprs <- (brca_white_test$stdprs - brca_white_test$meanscale) / sqrt(brca_white_test$varscale)


#AFR

bc_black_inc = bc_incidence %>% filter(Race == "Black")

brca_black_test$bc_inc <- bc_black_inc$Rate[floor(brca_black_test$age/5) + 1]
brca_black_test$alpha <- log(brca_black_test$bc_inc) - afrbeta * afrmean - 0.5 * afrbeta^2 * afrvar - afrpcs
brca_black_test$alphapc <- cbind(as.matrix(brca_black_test[, c("PC1", "PC4", "PC6")])) %*% coefficients(fit.alpha)[c(6, 9, 11)]

brca_black_test$meanscale <- cbind(1, as.matrix(brca_black_test[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(prsfitafr)[c(1:11)]
brca_black_test$varscale <- exp(cbind(1, as.matrix(brca_black_test[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(varfitafr)[c(1:11)])

brca_black_test$stdprs <- brca_black_test$PRS / sd(brca_black_test$PRS)
brca_black_test$abs_prob <- exp(brca_black_test$alpha + brca_black_test$alphapc + afrbeta * brca_black_test$stdprs)
brca_black_test$scaleprs <- (brca_black_test$stdprs - brca_black_test$meanscale) / sqrt(brca_black_test$varscale)


#AMR

bc_hisp_inc = bc_incidence %>% filter(Race == "Hispanic")

brca_hisp_test$bc_inc <- bc_hisp_inc$Rate[floor(brca_hisp_test$age/5) + 1]
brca_hisp_test$alpha <- log(brca_hisp_test$bc_inc) - hispbeta * hispmean - 0.5 * hispbeta^2 * hispvar - hisppcs
brca_hisp_test$alphapc <- cbind(as.matrix(brca_hisp_test[, c("PC1", "PC4", "PC6")])) %*% coefficients(fit.alpha)[c(6, 9, 11)]

brca_hisp_test$meanscale <- cbind(1, as.matrix(brca_hisp_test[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(prsfithisp)[c(1:11)]
brca_hisp_test$varscale <- exp(cbind(1, as.matrix(brca_hisp_test[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(varfithisp)[c(1:11)])

brca_hisp_test$stdprs <- brca_hisp_test$PRS / sd(brca_hisp_test$PRS)
brca_hisp_test$abs_prob <- exp(brca_hisp_test$alpha + brca_hisp_test$alphapc + hispbeta * brca_hisp_test$stdprs)
brca_hisp_test$scaleprs <- (brca_hisp_test$stdprs - brca_hisp_test$meanscale) / sqrt(brca_hisp_test$varscale)


#ASN

bc_asn_inc = bc_incidence %>% filter(Race == "Asian")

brca_asn_test$bc_inc <- bc_asn_inc$Rate[floor(brca_asn_test$age/5) + 1]
brca_asn_test$alpha <- log(brca_asn_test$bc_inc) - asnbeta * asnmean - 0.5 * asnbeta^2 * asnvar - asnpcs
brca_asn_test$alphapc <- cbind(as.matrix(brca_asn_test[, c("PC1", "PC4", "PC6")])) %*% coefficients(fit.alpha)[c(6, 9, 11)]

brca_asn_test$meanscale <- cbind(1, as.matrix(brca_asn_test[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(prsfitasn)[c(1:11)]
brca_asn_test$varscale <- exp(cbind(1, as.matrix(brca_asn_test[, c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])) %*% coefficients(varfitasn)[c(1:11)])

brca_asn_test$stdprs <- brca_asn_test$PRS / sd(brca_asn_test$PRS)
brca_asn_test$abs_prob <- exp(brca_asn_test$alpha + brca_asn_test$alphapc + asnbeta * brca_asn_test$stdprs)
brca_asn_test$scaleprs <- (brca_asn_test$stdprs - brca_asn_test$meanscale) / sqrt(brca_asn_test$varscale)


aoutest <- rbind(brca_white_test, brca_black_test, brca_hisp_test, brca_asn_test)