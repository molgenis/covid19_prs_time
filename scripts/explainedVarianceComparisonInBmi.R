#!/usr/bin/env Rscript

### Load libraries

library(data.table)
library(rjson)
library(tidyverse)

### Source data

workdir <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/"
preparedDataFile <- "longitudinal.RData"

setwd(workdir)
load(preparedDataFile)

### Constants

## The weeks since the first questionnaire + 1
weeksSinceStart <- c("1" = 1, 
                     "2" = 2, 
                     "3" = 3, 
                     "4" = 4, 
                     "5" = 5, 
                     "6" = 6, 
                     "7" = 8, 
                     "8" = 10, 
                     "9" = 12, 
                     "10" = 15, 
                     "11" = 16, 
                     "12" = 20, 
                     "13" = 24, 
                     "14" = 29, 
                     "15" = 32, 
                     "16" = 34, 
                     "17" = 36, 
                     "18" = 38, 
                     "19" = 42)

chosenWeeks <- c("4.0", "9.0", "14.0", "17.0")
confounders <- c("gender_recent", "age_recent", "age2_recent", "chronic_recent", "household_recent", "have_childs_at_home_recent")
bmiFormula <- as.formula("BMI~gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent+BMI_gwas")
bmiDiffFormula <- as.formula("BmiDiff~gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent+BMI_gwas")
bmiFormulaExt <- as.formula("BMI~gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent+BMI_gwas*validationSet")
bmiGwasFormulaExt <- as.formula("BMI~gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent+BMI_gwas*validationSet")


# Set plotting defaults
old <- theme_set(theme_classic(base_family="Helvetica"))
theme_update(line = element_line(
  colour = "black", size = 1 / (ggplot2::.pt * 72.27/96), 
  linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  axis.line = element_line(
    colour = "#595A5C", size = 1 / (ggplot2::.pt * 72.27/96), 
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  axis.ticks = element_line(
    colour = "#595A5C", size = 1 / (ggplot2::.pt * 72.27/96), 
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  plot.subtitle = element_text(color="gray50"),
  plot.title = element_text(face = "bold"),
)

### Functions
compareCoefficient <- function(b1, b2, se1, se2) {
  zObserved <- (b1 - b2)/sqrt(se1^2 + se2^2)
  return(2*pnorm(-abs(zObserved)))
}

fitBmi <- function(x, coeffLab = "BMI_gwas") {
  validationData <- x[x$validationSet == "validationSamples", c(confounders, "BMI", "BMI_gwas")]
  linearModelValidation = lm(
    bmiFormula,
    data = validationData)

  droppedData <- x[x$validationSet == "incompleteSamples", c(confounders, "BMI", "BMI_gwas")]
  linearModelDropouts = lm(
    bmiFormula,
    data = droppedData)

  valuesLm1 <- summary(linearModelValidation)$coefficients[coeffLab,]
  valuesLm2 <- summary(linearModelDropouts)$coefficients[coeffLab,]

  b1 <- unname(valuesLm1["Estimate"])
  b2 <- unname(valuesLm2["Estimate"])

  se1 <- unname(valuesLm1["Std. Error"])
  se2 <- unname(valuesLm2["Std. Error"])
  
  p1 <- unname(valuesLm1["Pr(>|t|)"])
  p2 <- unname(valuesLm2["Pr(>|t|)"])

  return(c("p.value" = compareCoefficient(b1, b2, se1, se2),
              "nSamples_ValidationSet" = nrow(validationData),
              "nSamples_DropoutSet" = nrow(droppedData),
              "beta_ValidationSet" = b1,
              "beta_DropoutSet" = b2,
           "pVal_ValidationSet" = p1,
           "pVal_DropoutSet" = p2,
              "se_ValidationSet" = se1,
              "se_DropoutSet" = se1))
}

### Main
bmiDataSetLong <- vragenLong %>% 
  select(PROJECT_PSEUDO_ID, BMI_gwas, BMI, vl2, !!!confounders) %>% 
  as_tibble() %>%
  mutate(
    validationSet = factor(case_when(
      PROJECT_PSEUDO_ID %in% validationSamples ~ "validationSamples", 
      TRUE ~ "incompleteSamples")),
    weekSinceStart = weeksSinceStart[vl2]) %>%
  filter(vl2 %in% c("4", "9", "14", "19")) %>%
  filter(!is.na(BMI) & !is.na(BMI_gwas))

bmiExplainedVariance <- bmiDataSetLong %>%
  group_by(weekSinceStart, validationSet) %>%
  summarise(
    nSamples = n(),
    rho = cor(BMI_gwas, BMI, use = "complete.obs"),
  ) %>% mutate(
    fisherZ = (1/2) * log((1+rho) / (1-rho)),
    se = 1/sqrt(nSamples - 3)
  ) %>% 
  pivot_wider(names_from = validationSet,
              values_from  = c(nSamples, rho, fisherZ, se)) %>%
  mutate(
    zObserved = (fisherZ_validationSamples - fisherZ_incompleteSamples) / sqrt((1/(nSamples_validationSamples-3)) + (1/(nSamples_incompleteSamples-3)))
  )

bmiCoefficients <- data.frame(t(sapply(unique(bmiDataSetLong$weekSinceStart), 
       function(week) {
         curData <- bmiDataSetLong[bmiDataSetLong$weekSinceStart == week,]
         return(fitBmi(curData))
         }, simplify = T)))

rownames(bmiCoefficients) <- unique(bmiDataSetLong$weekSinceStart)

bmiData4 <- vragenLong %>% 
  select(PROJECT_PSEUDO_ID, BMI_gwas, BMI, vl2, responsdatum.covid.vragenlijst, !!!confounders) %>% 
  as_tibble() %>%
  group_by(PROJECT_PSEUDO_ID) %>%
  mutate(
    validationSet = factor(case_when(
      any(vl2 == 19 & !is.na(responsdatum.covid.vragenlijst))  ~ "validationSamples", 
      TRUE ~ "incompleteSamples"))) %>%
  ungroup() %>%
  filter(vl2 %in% c("4")) %>%
  filter(!is.na(BMI) & !is.na(BMI_gwas))

fitBmi(bmiData4)

linearModelExt = lm(bmiFormulaExt, data = bmiDataSetLong %>% filter(vl2 == "4"))
t.test(BMI_gwas ~ validationSet, data = bmiDataSetLong %>% filter(vl2 == "4"))

#Verschillen tussen week 1,2,3,4,5

diff(c(1,2,3))

bmiDataSetLong2 <- bmiDataSetLong %>%
  group_by(PROJECT_PSEUDO_ID) %>%
  filter(validationSet == "incompleteSamples") %>%
  mutate(BmiDiff = c(NA_real_, diff(BMI))) %>%
  filter(!is.na(BmiDiff))

lm(bmiDiffFormula, data = bmiDataSetLong2)