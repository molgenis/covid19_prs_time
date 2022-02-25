#!/usr/bin/env Rscript

## Load libraries
library(foreign)
library(MASS)
library(readr)

## Main
workdir <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/"
phenoPath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/combined_questionnaires_v22_14-05-2021_genome_fitered_participants_filtered_correct_recoded/questionnaire_subset_participants_filtered_recoded_answers_14-05-2021"
prsLabelFile <- "prsLables.txt"
prsGsaFile <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/input_PGS_data_ugli_v5/PGS_combined_ugli_14-05-2021.txt"

pheno <- read_delim(paste0(phenoPath, ".txt"), delim = "\t", quote = "", guess_max = 100000)
dim(pheno)
pheno2 <- as.data.frame(pheno)
row.names(pheno2) <- pheno2[,1]

colnames(pheno2)[1] <- "PROJECT_PSEUDO_ID"

prsGsa <- read.delim(prsGsaFile, stringsAsFactors = F)

merged <- merge(pheno2, prsGsa, by = "PROJECT_PSEUDO_ID")

#age_recent + age2_recent + chronic_recent + household_recent + have_childs_at_home_recent + gender_recent
covs <-c("age_recent", "age2_recent", "chronic_recent", "household_recent", "have_childs_at_home_recent", "gender_recent")


## test 1 ##
# PRS: BMI
# vraag: covt06_society_adu_q_1_c
# Hercodering: geen

testDf1 <- merged[c("covt06_society_adu_q_1_c", "EduYears", covs)]
testDf1$covt06_society_adu_q_1_c_ordered <- ordered(testDf1$covt06_society_adu_q_1_c, levels = 1:5)
testDf1 <- testDf1[complete.cases(testDf1),]
dim(testDf1)
olmod1 <- polr(covt06_society_adu_q_1_c_ordered ~ EduYears + age_recent + age2_recent + chronic_recent + household_recent + have_childs_at_home_recent + gender_recent, data = testDf1, method = "logistic", Hess = TRUE)
solmod1 <- summary(olmod1)
solmod1
print(2*pt(solmod1$coefficients[1,3], solmod1$df.residual, lower=FALSE))


## test 2 ##
# PRS: Schizophrenia
# vraag: covt17_covidvaccine_adu_q_2_e1 (ik vind het corona vaccin eng)
# Hercodering: {“5.0”: 1.0, “4.0": 2.0, “3.0”: 3.0, “2.0": 4.0, “1.0”: 5.0, “6.0": NaN}

recoding2 <- c(`5.0` = 1.0, `4.0` = 2.0, `3.0` = 3.0, `2.0` = 4.0, `1.0` = 5.0, `6.0` = NA_real_)

testDf2 <- merged[c("covt17_covidvaccine_adu_q_2_e1", "Schizophrenia", covs)]
testDf2$covt17_covidvaccine_adu_q_2_e1_ordered <- ordered(as.numeric(names(recoding2[testDf2$covt17_covidvaccine_adu_q_2_e1])), levels = 1:5)
testDf2 <- testDf2[complete.cases(testDf2),]
dim(testDf2)
olmod2 <- polr(covt17_covidvaccine_adu_q_2_e1_ordered ~ Schizophrenia + age_recent + age2_recent + chronic_recent + household_recent + have_childs_at_home_recent + gender_recent, data = testDf2, method = "logistic", Hess = TRUE)
solmod2 <- summary(olmod2)
solmod2
print(2*pt(solmod2$coefficients[1,3], solmod2$df.residual, lower=FALSE))

## test 3 ##
# PRS: OCD
# vraag: covt01_scl90som11_adu_q_1 (je lichamelijk ergens slap voelen / in welke mate had u de afgelopen 7 dagen last van:)
# Hercodering: geen

testDf3 <- merged[c("covt01_scl90som11_adu_q_1", "OCD", covs)]
testDf3$covt01_scl90som11_adu_q_1_ordered <- ordered(testDf3$covt01_scl90som11_adu_q_1, levels = 1:5)
testDf3 <- testDf3[complete.cases(testDf3),]
dim(testDf3)
olmod3 <- polr(covt01_scl90som11_adu_q_1_ordered ~ OCD + age_recent + age2_recent + chronic_recent + household_recent + have_childs_at_home_recent + gender_recent, data = testDf3, method = "logistic", Hess = TRUE)
solmod3 <- summary(olmod3)
solmod3
print(2*pt(solmod3$coefficients[1,3], solmod3$df.residual, lower=FALSE))

## test 4 ##

testDf4 <- merged[c("covt10_washhands_adu_q_1_a", "Anxiety.tension")]
testDf4$covt10_washhands_adu_q_1_a_ordered <- ordered(testDf4$covt10_washhands_adu_q_1_a, levels = 1:4)
testDf4 <- testDf4[complete.cases(testDf4),]
dim(testDf4)
olmod4 <- polr(covt10_washhands_adu_q_1_a_ordered ~ Anxiety.tension + age_recent + age2_recent + chronic_recent + household_recent + have_childs_at_home_recent + gender_recent, data = testDf4, method = "logistic", Hess = TRUE)
olmod4 <- polr(covt10_washhands_adu_q_1_a_ordered ~ Anxiety.tension, data = testDf4, method = "logistic", Hess = TRUE)

solmod4 <- summary(olmod4)
solmod4
print(2*pt(solmod4$coefficients[1,3], solmod4$df.residual, lower=FALSE))

testDf5 <- merged[c("covt15_controlsituation_adu_q_1", "General.risky.behavior", covs)]
testDf5$covt15_controlsituation_adu_q_1_ordered <- ordered(testDf5$covt15_controlsituation_adu_q_1, levels = 1:5)
testDf5 <- testDf5[complete.cases(testDf5),]

olmod4 <- polr(covt15_controlsituation_adu_q_1_ordered ~ General.risky.behavior + age_recent + age2_recent + chronic_recent + household_recent + have_childs_at_home_recent + gender_recent, data = testDf5, method = "logistic", Hess = TRUE)
olmod4 <- polr(covt15_controlsituation_adu_q_1_ordered ~ General.risky.behavior, data = testDf5, method = "logistic", Hess = TRUE)

testDf6 <- merged[,c("covt07_fatigue_adu_q_2_a", "Life.satisfaction")]
testDf6$covt07_fatigue_adu_q_2_a_ordered <- ordered(testDf6$covt07_fatigue_adu_q_2_a, levels = 1:7)
testDf6 <- testDf6[complete.cases(testDf6),]

olmod6 <- polr(covt07_fatigue_adu_q_2_a_ordered ~ Life.satisfaction, data = testDf6, method = "logistic", Hess = TRUE)


testDf6 <- merged[,c("covt01_fatigue_adu_q_1_a", "Life.satisfaction")]
table(testDf6$covt01_fatigue_adu_q_1_a)
testDf6$covt01_fatigue_adu_q_1_a_ordered <- ordered(testDf6$covt01_fatigue_adu_q_1_a, levels = 1:7)
table(testDf6$covt01_fatigue_adu_q_1_a_ordered)
testDf6 <- testDf6[complete.cases(testDf6),]

olmod6 <- polr(covt01_fatigue_adu_q_1_a_ordered ~ Life.satisfaction, data = testDf6, method = "logistic", Hess = TRUE)




pdf("controlSituation.pdf")
boxplot(General.risky.behavior ~ covt15_controlsituation_adu_q_1, data = testDf5[testDf5$gender_recent == 1,])
dev.off()
  
# Test
dat <- read.dta("https://stats.idre.ucla.edu/stat/data/ologit.dta")
solmod <- 
