###################################################################
####### 4. Calculating Absolute Risk of Disease Mortality #########
###################################################################

# Calculating Mortality Age Groups and Yearly Mortality Risk

## EUR
bc_white_mort = bc_mortality %>% filter(Race == "White")

brca_white_test$bc_inc5 <- bc_white_inc$Rate[floor(brca_white_test$age/5) + 2]
brca_white_test$alpha5 <- log(brca_white_test$bc_inc5) - eurbeta * eurmean - 0.5 * eurbeta^2 * eurvar - eurpcs
brca_white_test$abs_prob5 <- exp(brca_white_test$alpha5 + brca_white_test$alphapc + eurbeta * brca_white_test$stdprs)

brca_white_test$bc_mortnow <- bc_white_mort$Rate[floor(brca_white_test$age/5) + 1]
brca_white_test$bc_mort5 <- bc_white_mort$Rate[floor(brca_white_test$age/5) + 2]

## AFR
bc_black_mort = bc_mortality %>% filter(Race == "Black")

brca_black_test$bc_inc5 <- bc_black_inc$Rate[floor(brca_black_test$age/5) + 2]
brca_black_test$alpha5 <- log(brca_black_test$bc_inc5) - afrbeta * afrmean - 0.5 * afrbeta^2 * afrvar - afrpcs
brca_black_test$abs_prob5 <- exp(brca_black_test$alpha5 + brca_black_test$alphapc + afrbeta * brca_black_test$stdprs)

brca_black_test$bc_mortnow <- bc_black_mort$Rate[floor(brca_black_test$age/5) + 1]
brca_black_test$bc_mort5 <- bc_black_mort$Rate[floor(brca_black_test$age/5) + 2]


## AMR
bc_hisp_mort = bc_mortality %>% filter(Race == "Hispanic")

brca_hisp_test$bc_inc5 <- bc_hisp_inc$Rate[floor(brca_hisp_test$age/5) + 2]
brca_hisp_test$alpha5 <- log(brca_hisp_test$bc_inc5) - hispbeta * hispmean - 0.5 * hispbeta^2 * hispvar - hisppcs
brca_hisp_test$abs_prob5 <- exp(brca_hisp_test$alpha5 + brca_hisp_test$alphapc + hispbeta * brca_hisp_test$stdprs)

brca_hisp_test$bc_mortnow <- bc_hisp_mort$Rate[floor(brca_hisp_test$age/5) + 1]
brca_hisp_test$bc_mort5 <- bc_hisp_mort$Rate[floor(brca_hisp_test$age/5) + 2]


## ASN
bc_asn_mort = bc_mortality %>% filter(Race == "Asian")

brca_asn_test$bc_inc5 <- bc_asn_inc$Rate[floor(brca_asn_test$age/5) + 2]
brca_asn_test$alpha5 <- log(brca_asn_test$bc_inc5) - asnbeta * asnmean - 0.5 * asnbeta^2 * asnvar - asnpcs
brca_asn_test$abs_prob5 <- exp(brca_asn_test$alpha5 + brca_asn_test$alphapc + asnbeta * brca_asn_test$stdprs)

brca_asn_test$bc_mortnow <- bc_asn_mort$Rate[floor(brca_asn_test$age/5) + 1]
brca_asn_test$bc_mort5 <- bc_asn_mort$Rate[floor(brca_asn_test$age/5) + 2]



# Calculating 10 Year Risk of Disease Mortality

aoutest <- rbind(brca_white_test, brca_black_test, brca_hisp_test, brca_asn_test)

calcprob = function(df){
  final_prob = rep(0, nrow(df))
  for(i in c(1:9)){
    tot_inc = (df$abs_prob * min(i, 5) + df$abs_prob5 * (max(i, 5) - 5))
    tot_mort = (df$bc_mortnow * (5 - min(i, 5)) + df$bc_mort5 * (10 - max(i, 5)))
    final_prob = final_prob + tot_inc * tot_mort
  }
  
  return(final_prob)
}

aoutest$mort_prob = calcprob(aoutest)