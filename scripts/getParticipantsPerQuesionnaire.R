#!/usr/bin/env Rscript

# Script that processes pheno files to extract how many participants are used in both baseline analyses and longitudinal analyses

## Load libraries

library(tidyverse)

## Main

phenoBaselinePath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/combined_questionnaires_v22_14-05-2021_genome_fitered_participants_filtered_correct_recoded/questionnaire_subset_participants_filtered_recoded_answers_14-05-2021"
phenoLongitudinalPath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/combined_questionnaires_v22_14-05-2021_genome_fitered_participants_filtered_correct_recoded/questionnaire_subset_participants_filtered_recoded_answers_longitudinal_filtered_15-05-2021"
qOverviewFile <- "quest_overview_nl_new_quest17_codes_updated_14-days_v2.txt"

## Baseline samples
# Read pheno samples
pheno <- read_delim(paste0(phenoBaselinePath, ".txt"), delim = "\t", quote = "", guess_max = 100000)
pheno3 <- as.data.frame(pheno)
row.names(pheno3) <- pheno3[,1]

colnames(pheno3)[1] <- "PROJECT_PSEUDO_ID"

#na col to use in reshapre for missing questions 
pheno3$naCOl <- NA

##Load and format questions meta data
qOverview <- as.matrix(read.delim(qOverviewFile, stringsAsFactors = F, row.names = 1))
vls <- colnames(qOverview)[-c(20,21)]

qOverview[2,]

qOverview[,-c(20,21)] <- apply(qOverview[,-c(20,21)], 1:2, function(x){
  if(x==""){
    return ("")
  } else if(!x %in% colnames(pheno3)) {
    return("")
  } else if(grepl("responsedate_adu", x)){
    #always include data question
    return(x)
  } else {
    return(x)#no masking
    # mask questions with >75% missing
    #if(sum(is.na(pheno3[,x])) >= (totalPart/4)){
    #  return (x)
    #} else {
    #return ("")
    #}
  }
})


qPassQc <- sapply(rownames(qOverview), function(q){sum(qOverview[q,1:19]!="")>=2})
table(qPassQc)

#qOverview2 <- qOverview[qPassQc,-c(20,21)]
qOverview2 <- qOverview[,-c(20,21)]
#head(qOverview2[,1])
dim(qOverview2)

qNameMap <- data.frame(orginal = row.names(qOverview2), new = make.names(row.names(qOverview2), unique = T), stringsAsFactors = F)
row.names(qNameMap) <- row.names(qOverview2)

qList <- lapply(qNameMap[,1], function(q){
  qs <- qOverview2[q,]
  qs[qs==""]="naCOl"
  qs[!qs %in% colnames(pheno3)] <- "naCOl"
  return(qs)
})
names(qList) <- qNameMap[,2]
#qList[[qNameMap["responsdatum covid-vragenlijst",2]]]
#qList[[qNameMap["als u moet kiezen, denkt u zelf dat u een coronavirus/covid-19 infectie hebt (gehad)?",2]]]

## remove quesed body weight
qList[[qNameMap["op welk moment hebt u dit lichaamsgewicht gemeten?",2]]]
qList[[qNameMap["wat is uw lichaamsgewicht (in kg)?",2]]]
qList[[qNameMap["BMI",2]]]

for(vl in colnames(qOverview)[-c(20,21)]){
  #1 = morning; 2 = afternoon/evening; 3 = guess
  weightMeasureType <- pheno[,qList[[qNameMap["op welk moment hebt u dit lichaamsgewicht gemeten?",2]]][vl]]
  weightMask <- is.na(weightMeasureType) | weightMeasureType == 3
  
  pheno[weightMask,qList[[qNameMap["wat is uw lichaamsgewicht (in kg)?",2]]][vl]] <- NA
  pheno[weightMask,qList[[qNameMap["BMI",2]]][vl]] <- NA
}

## Reshape to long format and clean some variables
if(any(names(qList) %in% colnames(pheno3))){
  stop("Column name clash after reshape")
}

vragenLong <- reshape(pheno3, direction = "long", idvar = "PROJECT_PSEUDO_ID", varying = qList, v.names = names(qList), times = vls, timevar = "vl")
vragenLong$vl2 <- as.numeric(factor(vragenLong$vl, levels = vls, ordered = T))
vragenLong$vl3 <- factor(vragenLong$vl, levels = vls, ordered = F)

summarized <- vragenLong %>% select(PROJECT_PSEUDO_ID, vl2, responsdatum.covid.vragenlijst) %>% filter(!is.na(responsdatum.covid.vragenlijst)) %>% group_by(vl2) %>%
  summarise(nSamples = length(unique(PROJECT_PSEUDO_ID)))

# Table 1, g1

summaryPerSample <- vragenLong %>% group_by(PROJECT_PSEUDO_ID) %>% summarise(
  nQuestionnaires = n_distinct(responsdatum.covid.vragenlijst),
  BMI = ifelse(length(na.omit(BMI)) <= 0, NA_real_, head(na.omit(BMI),1)),
  gender = head(na.omit(gender_recent),1),
  age = head(na.omit(age_recent),1),
  chronic = head(na.omit(chronic_recent),1),
  household = head(na.omit(household_recent),1),
  children = head(na.omit(have_childs_at_home_recent),1)
)

print(summaryPerSample %>%
        ungroup() %>%
        summarise(
          nSamples = n(),
          averageQuestionnaireCount = mean(nQuestionnaires),
          BMI_mean = mean(BMI, na.rm=T),
          BMI_sd = sd(BMI, na.rm=T),
          age_min = min(age),
          age_mean = mean(age),
          age_sd = sd(age),
          age_max = max(age),
          male_N =sum(gender == 1),
          male_percentage = sum(gender == 1) / n() * 100,
          female_N = sum(gender == 0),
          female_percentage = sum(gender == 0) / n() * 100,
          chronic_N = sum(chronic == 1),
          chronic_percentage = sum(chronic == 1) / n() * 100,
          household_percentage = sum(household == 1) / n() * 100,
          household_N = sum(household == 1),
          children_percentage = sum(children == 1) / n() * 100,
          children_N = sum(children == 1),
        ), width = Inf)
## Longitudinal samples

# read pheno
pheno <- read_delim(paste0(phenoLongitudinalPath, ".txt"), delim = "\t", quote = "", guess_max = 100000)
pheno3 <- as.data.frame(pheno)
row.names(pheno3) <- pheno3[,1]

colnames(pheno3)[1] <- "PROJECT_PSEUDO_ID"

#na col to use in reshapre for missing questions 
pheno3$naCOl <- NA

##Load and format questions meta data
qOverview <- as.matrix(read.delim(qOverviewFile, stringsAsFactors = F, row.names = 1))
vls <- colnames(qOverview)[-c(20,21)]

qOverview[2,]

qOverview[,-c(20,21)] <- apply(qOverview[,-c(20,21)], 1:2, function(x){
  if(x==""){
    return ("")
  } else if(!x %in% colnames(pheno3)) {
    return("")
  } else if(grepl("responsedate_adu", x)){
    #always include data question
    return(x)
  } else {
    return(x)#no masking
    # mask questions with >75% missing
    #if(sum(is.na(pheno3[,x])) >= (totalPart/4)){
    #  return (x)
    #} else {
    #return ("")
    #}
  }
})


qPassQc <- sapply(rownames(qOverview), function(q){sum(qOverview[q,1:19]!="")>=2})
table(qPassQc)

#qOverview2 <- qOverview[qPassQc,-c(20,21)]
qOverview2 <- qOverview[,-c(20,21)]
#head(qOverview2[,1])
dim(qOverview2)

qNameMap <- data.frame(orginal = row.names(qOverview2), new = make.names(row.names(qOverview2), unique = T), stringsAsFactors = F)
row.names(qNameMap) <- row.names(qOverview2)

qList <- lapply(qNameMap[,1], function(q){
  qs <- qOverview2[q,]
  qs[qs==""]="naCOl"
  qs[!qs %in% colnames(pheno3)] <- "naCOl"
  return(qs)
})
names(qList) <- qNameMap[,2]
#qList[[qNameMap["responsdatum covid-vragenlijst",2]]]
#qList[[qNameMap["als u moet kiezen, denkt u zelf dat u een coronavirus/covid-19 infectie hebt (gehad)?",2]]]

## remove quesed body weight
qList[[qNameMap["op welk moment hebt u dit lichaamsgewicht gemeten?",2]]]
qList[[qNameMap["wat is uw lichaamsgewicht (in kg)?",2]]]
qList[[qNameMap["BMI",2]]]

for(vl in colnames(qOverview)[-c(20,21)]){
  #1 = morning; 2 = afternoon/evening; 3 = guess
  weightMeasureType <- pheno[,qList[[qNameMap["op welk moment hebt u dit lichaamsgewicht gemeten?",2]]][vl]]
  weightMask <- is.na(weightMeasureType) | weightMeasureType == 3
  
  pheno[weightMask,qList[[qNameMap["wat is uw lichaamsgewicht (in kg)?",2]]][vl]] <- NA
  pheno[weightMask,qList[[qNameMap["BMI",2]]][vl]] <- NA
}

## Reshape to long format and clean some variables
if(any(names(qList) %in% colnames(pheno3))){
  stop("Column name clash after reshape")
}


vragenLong <- reshape(pheno3, direction = "long", idvar = "PROJECT_PSEUDO_ID", varying = qList, v.names = names(qList), times = vls, timevar = "vl")
vragenLong$vl2 <- as.numeric(factor(vragenLong$vl, levels = vls, ordered = T))
vragenLong$vl3 <- factor(vragenLong$vl, levels = vls, ordered = F)

summarized <- vragenLong %>% select(PROJECT_PSEUDO_ID, vl2, responsdatum.covid.vragenlijst) %>% filter(!is.na(responsdatum.covid.vragenlijst)) %>% group_by(vl2) %>%
  summarise(nSamples = length(unique(PROJECT_PSEUDO_ID)))

# Table 1, g2

summaryPerSample <- vragenLong %>% group_by(PROJECT_PSEUDO_ID) %>% summarise(
  nQuestionnaires = n_distinct(responsdatum.covid.vragenlijst),
  BMI = ifelse(length(na.omit(BMI)) <= 0, NA_real_, head(na.omit(BMI),1)),
  gender = head(na.omit(gender_recent),1),
  age = head(na.omit(age_recent),1),
  chronic = head(na.omit(chronic_recent),1),
  household = head(na.omit(household_recent),1),
  children = head(na.omit(have_childs_at_home_recent),1)
)

print(summaryPerSample %>%
  ungroup() %>%
  summarise(
    nSamples = n(),
    averageQuestionnaireCount = mean(nQuestionnaires),
    BMI_mean = mean(BMI, na.rm=T),
    BMI_sd = sd(BMI, na.rm=T),
    age_min = min(age),
    age_mean = mean(age),
    age_sd = sd(age),
    age_max = max(age),
    male_N =sum(gender == 1),
    male_percentage = sum(gender == 1) / n() * 100,
    female_N = sum(gender == 0),
    female_percentage = sum(gender == 0) / n() * 100,
    chronic_N = sum(chronic == 1),
    chronic_percentage = sum(chronic == 1) / n() * 100,
    household_percentage = sum(household == 1) / n() * 100,
    household_N = sum(household == 1),
    children_percentage = sum(children == 1) / n() * 100,
    children_N = sum(children == 1),
  ), width = Inf)