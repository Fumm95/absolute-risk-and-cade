###########################################################################
####### 5. Absolute Risk Tables and Figures, and Calculating CADE #########
###########################################################################

# Calculating Performance of Ranking Methods

## Top 5%, 10%, and 20% of population

aoutesttop5 = floor(nrow(aoutest) * 0.05)
aoutesttop10 = floor(nrow(aoutest) * 0.1)
aoutesttop20 = floor(nrow(aoutest) * 0.2)

rank_prs <- aoutest[order(-aoutest$stdprs),]
rank_scprs <- aoutest[order(-aoutest$scaleprs),]
rank_absprob <- aoutest[order(-aoutest$abs_prob),]

### Cases Identified
top5 <- c(sum(rank_prs[c(1:aoutesttop5),]$total_case), sum(rank_scprs[c(1:aoutesttop5),]$total_case), sum(rank_absprob[c(1:aoutesttop5),]$total_case))
top10 <- c(sum(rank_prs[c(1:aoutesttop10),]$total_case), sum(rank_scprs[c(1:aoutesttop10),]$total_case), sum(rank_absprob[c(1:aoutesttop10),]$total_case))
top20 <- c(sum(rank_prs[c(1:aoutesttop20),]$total_case), sum(rank_scprs[c(1:aoutesttop20),]$total_case), sum(rank_absprob[c(1:aoutesttop20),]$total_case))

aoutesttop5
aoutesttop10
aoutesttop20

top_mat <- cbind(top5, top10, top20)
top_mat


## Bottom 5%, 10%, and 20% of population

aoutestbottom5 = floor(nrow(aoutest) * 0.05)
aoutestbottom10 = floor(nrow(aoutest) * 0.1)
aoutestbottom20 = floor(nrow(aoutest) * 0.2)

rank_prs <- aoutest[order(aoutest$stdprs),]
rank_scprs <- aoutest[order(aoutest$scaleprs),]
rank_absprob <- aoutest[order(aoutest$abs_prob),]

### Cases Identified
bottom5 <- c(sum(rank_prs[c(1:aoutestbottom5),]$total_case), sum(rank_scprs[c(1:aoutestbottom5),]$total_case), sum(rank_absprob[c(1:aoutestbottom5),]$total_case))
bottom10 <- c(sum(rank_prs[c(1:aoutestbottom10),]$total_case), sum(rank_scprs[c(1:aoutestbottom10),]$total_case), sum(rank_absprob[c(1:aoutestbottom10),]$total_case))
bottom20 <- c(sum(rank_prs[c(1:aoutestbottom20),]$total_case), sum(rank_scprs[c(1:aoutestbottom20),]$total_case), sum(rank_absprob[c(1:aoutestbottom20),]$total_case))

aoutestbottom5
aoutestbottom10
aoutestbottom20

bottom_mat <- cbind(bottom5, bottom10, bottom20)
bottom_mat



# Figures for Absolute Risk by Age Group (Ex. Ages 45-50)

aoutest_age = aoutest %>% filter(age >= 45 & age < 50) %>%
  filter(total_case == 0)

aoutest_ageeur = aoutest_age %>% filter(hisp_race == "Non Hispanic White")

aoutest_ageafr = aoutest_age %>% filter(hisp_race == "Black or African American")

aoutest_agehisp = aoutest_age %>% filter(hisp_race == "Hispanic or Latino")

aoutest_ageasn = aoutest_age %>% filter(hisp_race == "Asian")

nrow(aoutest_ageeur)
nrow(aoutest_ageafr)
nrow(aoutest_agehisp)
nrow(aoutest_ageasn)

# Age 45-50

options(repr.plot.width=10, repr.plot.height=5)

aoutest_ageeur_small = aoutest_ageeur[sample(c(1:nrow(aoutest_ageeur)), size = 280),]
aoutest_ageafr_small = aoutest_ageafr[sample(c(1:nrow(aoutest_ageafr)), size = 140),]
aoutest_agehisp_small = aoutest_agehisp[sample(c(1:nrow(aoutest_agehisp)), size = 160),]
aoutest_ageasn_small = aoutest_ageasn[sample(c(1:nrow(aoutest_ageasn)), size = 20),]

aoutest_agesmall = rbind(aoutest_ageeur_small, aoutest_ageafr_small, aoutest_agehisp_small, aoutest_ageasn_small)

aoutest_agesmall = aoutest_agesmall[sample(nrow(aoutest_agesmall)),]


myplot <- aoutest_agesmall %>% 
  mutate(hisp_race = factor(hisp_race, levels = c("Black or African American", "Asian",
                                                  "Non Hispanic White", "Hispanic or Latino"))) %>% 
  ggplot(aes(x = scaleprs, y = 5 * abs_prob)) +
  geom_point(aes(fill = hisp_race), shape = 21, size = 4, stroke = 0.5) +
  scale_x_continuous(limits = c(-3.75, 3.75)) +
  scale_y_continuous(n.breaks = 5) +
  #ggtitle("5 Year Disease Incidence for Age 45-50") +
  #labs(fill = "Self-reported Race") +
  xlab("Standardized PRS") +
  ylab("5 Year Risk of Disease Incidence") +
  theme_bw(base_size = 20) +
  theme(legend.position = "none")




# Comparing Distribution of Top 10% (Ex. Age 45-50)

brca_full_ordered = aoutest %>% filter(age >= 45 & age < 50) %>% filter(total_case == 0)

brca_full_ordered$scprs_order <- rank(-brca_full_ordered$scaleprs)
brca_full_ordered$mortprob_order <- rank(-brca_full_ordered$mort_prob)
brca_full_ordered$absprob_order <- rank(-brca_full_ordered$abs_prob)

brca_full_ordered$scprs_top1000 <- as.numeric(brca_full_ordered$scprs_order <= floor(nrow(brca_full_ordered) * 0.1))
brca_full_ordered$mortprob_top1000 <- as.numeric(brca_full_ordered$mortprob_order <= floor(nrow(brca_full_ordered) * 0.1))
brca_full_ordered$absprob_top1000 <- as.numeric(brca_full_ordered$absprob_order <= floor(nrow(brca_full_ordered) * 0.1))

# Standardized PRS

distr_plot_final10 = cbind(prop.table(table(brca_full_ordered$hisp_race, brca_full_ordered$scprs_top1000), margin = 2)[,2], rep("Standardized PRS", 4), rep("Age 45-50", 4))
colnames(distr_plot_final10) = c("Proportion", "Ranking Method", "Age Group")
distr_plot_final10 = as.data.frame(distr_plot_final10)
distr_plot_final10$Race = unlist(rownames(distr_plot_final10))

# Absolute Risk of Incidence

distr_plot_df= cbind(prop.table(table(brca_full_ordered$hisp_race, brca_full_ordered$absprob_top1000), margin = 2)[,2], rep("Disease Incidence", 4), rep("Age 45-50", 4))
colnames(distr_plot_df) = c("Proportion", "Ranking Method", "Age Group")
distr_plot_df = as.data.frame(distr_plot_df)
distr_plot_df$Race = unlist(rownames(distr_plot_df))

distr_plot_final10 = rbind(distr_plot_final10, distr_plot_df)


# Absolute Risk of Incidence & Mortality

distr_plot_df= cbind(prop.table(table(brca_full_ordered$hisp_race, brca_full_ordered$mortprob_top1000), margin = 2)[,2], rep("Disease Mortality", 4), rep("Age 45-50", 4))
colnames(distr_plot_df) = c("Proportion", "Ranking Method", "Age Group")
distr_plot_df = as.data.frame(distr_plot_df)
distr_plot_df$Race = unlist(rownames(distr_plot_df))

distr_plot_final10 = rbind(distr_plot_final10, distr_plot_df)

distr_plot_final10 %>% 
  mutate(Proportion = as.numeric(Proportion)) %>%
  mutate(`Age Group` = factor(`Age Group`, levels = c("Age 45-50", "Age 55-60", "Age 65-70"))) %>%
  mutate(`Ranking Method` = factor(`Ranking Method`, levels = c("Standardized PRS", "Disease Incidence", "Disease Mortality"))) %>%
  mutate(hisp_race = factor(Race, levels = c("Black or African American", "Asian",
                                             "Non Hispanic White", "Hispanic or Latino"))) %>%
  filter(`Age Group` == "Age 45-50") %>%
  ggplot(aes(fill=hisp_race, y=Proportion, x = `Ranking Method`)) +
  geom_bar(position="fill", stat="identity", color = "black") +
  scale_x_discrete(labels = c("Standardized\nPRS", "Disease\nIncidence", "Disease\nMortality")) +
  # coord_flip() +
  labs(fill = "Self-reported Race and Ethnicity") +
  theme_bw(base_size = 35) +
  theme(legend.position = "none", 
        axis.text.x = element_text(face = "bold", color = "black"),
        axis.text.y = element_text(face = "bold", color = "black"),
        axis.title.y = element_text(face = "bold"), 
        axis.title.x = element_blank())


# Calculating CADE (Ex. Age 45-50)

aoutest_white = aoutest %>% filter(hisp_race == "Non Hispanic White" & age >= 45 & age < 50)
cade_eur = eurbeta * mean(aoutest_white$abs_prob, na.rm = T)

aoutest_black = aoutest %>% filter(hisp_race == "Black or African American" & age >= 45 & age < 50)
cade_afr = afrbeta * mean(aoutest_black$abs_prob, na.rm = T)

aoutest_hisp = aoutest %>% filter(hisp_race == "Hispanic or Latino" & age >= 45 & age < 50)
cade_hisp = hispbeta * mean(aoutest_hisp$abs_prob, na.rm = T)

aoutest_asian = aoutest %>% filter(hisp_race == "Asian" & age >= 45 & age < 50)
cade_asn = asnbeta * mean(aoutest_asian$abs_prob, na.rm = T)



