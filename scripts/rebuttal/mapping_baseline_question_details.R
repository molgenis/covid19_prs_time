#!/usr/bin/env Rscript

## ----
## Author:  C.A. (Robert) Warmerdam
## Email:   c.a.warmerdam@umcg.nl
##
## Copyright (c) C.A. Warmerdam, 2021
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## A copy of the GNU General Public License can be found in the LICENSE file in the
## root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
## ----

# Load libraries
library(data.table)
library(rjson)
library(tidyverse)

# Declare constants
custom_mapping <- c("Mini combined, depressief" = "Depressive score, calculated with the M.I.N.I. questionnaire items",
  "Positive tested cumsum" = "Ever self reported positive SARS-CoV-2 PCR test",
  "BMI" = "BMI calculated from self reported weight and height",
  "ik voelde me fit / in welke mate had u de afgelopen 7 dagen last van:" = "i felt good / to what extent have you experienced following in the last 7 days:",
  "vermijden van diensten van de contactberoepen (kapper/masseur/fysiotherapeut/rij-instructeurs ect.) / mocht de overheid weer besluiten terug te gaan naar een intelligente lockdown, ben ik bereid me te houden aan de volgende maatregelen:" = "avoid services from contact professionals (hairdresser/masseur/physical therapist/driving instructor etc.) / should the government decide to go back into an intelligent lockdown, i am prepared to comply with the following measures:")

questionnaire_mapping <- setNames(
  object = 1:19, 
  nm = paste0("X", formatC(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,15.5,16,16.5,17), digits = 1, format = "f")))

# Declare function definitions

load_pheno_file <- function(phenoPath) {
  if(!file.exists(paste0(phenoPath, ".rds"))){
    #Run once to conver pheno data to RDS format
    pheno <- read_delim(paste0(phenoPath, ".txt"), delim = "\t", quote = "", guess_max = 100000)
    dim(pheno)
    pheno2 <- as.data.frame(pheno)
    row.names(pheno2) <- pheno2[,1]
    
    colnames(pheno2)[1] <- "PROJECT_PSEUDO_ID"
    
    saveRDS(pheno2, paste0(phenoPath, ".rds"))
  } else {
    pheno2 <- readRDS(paste0(phenoPath, ".rds"))
  }
  
}

get_total_sample_size <- function(my_json_string_fixed) {
  list_sample_size = fromJSON(json_str = my_json_string_fixed)
  return(sum(unlist(list_sample_size)[names(list_sample_size) != "NaN"]))
}

# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv=NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  
  # Process input
  phenoGenotypedPath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/combined_questionnaires_v22_14-05-2021_genome_fitered_participants_filtered_correct_recoded/questionnaire_subset_participants_filtered_recoded_answers_14-05-2021"
  phenoGeno <- load_pheno_file(phenoGenotypedPath)
  input_question_overview <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/combined_questionnaires_v22_14-05-2021_genome_fitered_participants_filtered_correct_recoded/questionnaire_subset_recode_info_overview_14-05-2021.txt" # datafrane with the question ids over time
  workdir <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/"
  qOverviewFile <- "quest_overview_nl_new_quest17_codes_updated_14-days_v2.txt"
  questions_baseline <- fread("baseline_questions_with_label.txt", header = T, sep = "\t",
                              col.names = c("Outcome_label", "DEFINITION_EN"))
  
  qOverview <- as.matrix(read.delim(qOverviewFile, stringsAsFactors = F, row.names = 1))
  vls <- colnames(qOverview)[-c(20,21)]
  
  question_overview <- fread(input_question_overview) %>%
    select(-V1)
  
  question_overview_tibble <- as_tibble(qOverview, rownames = "question_nl") %>% 
    pivot_longer(cols=starts_with("X"), names_to = "vl", values_to = "question_id") %>%
    filter(question_id != "") %>%
    full_join(question_overview, by = c("question_id" = "question_id"))
  
  question_overview_filtered <- question_overview_tibble %>% filter(
    (is_first_time_asked == T & question_nl != "Positive tested cumsum") 
    | question_id == "covt17_positive_tested_cumsum")
  
  configurationFilePath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/barometer_statistics/coronabarometer-statistics/config.v3.json"
  configuration <- fromJSON(file = configurationFilePath)
  covidVariablesPerWeek <- configuration$covidResFilePerWeekVariables
  covidResFileTypeVersionPerWeek <- configuration$fileTypeVersionPerWeek
  
  variablesList <- list()
  
  for (weekLabel in names(covidVariablesPerWeek)) {
    print(paste0("Loading: ", weekLabel, ", ", covidVariablesPerWeek[weekLabel]))
    
    if (covidResFileTypeVersionPerWeek[[weekLabel]] %in% c("v3", "v4")) {
      fileList <- lapply(covidVariablesPerWeek[[weekLabel]], function(covidResFile) {
        fread(as.character(covidResFile))
      })
      
      variablesList[[weekLabel]] <- Reduce(merge,fileList)
      
    } else if (covidResFileTypeVersionPerWeek[[weekLabel]] == "v2") {
      fileList <- lapply(covidVariablesPerWeek[[weekLabel]], function(covidResFile) {
        fread(as.character(covidResFile))
      })
      
      variablesList[[weekLabel]] <- Reduce(merge,fileList)
      
    }
  }
  
  variables <- bind_rows(variablesList) %>%
    filter(VARIABLE_NAME %in% question_overview_filtered$question_id)
  
  question_overview_filtered_with_en <- question_overview_filtered %>% 
    left_join(variables, c("question_id" = "VARIABLE_NAME")) %>%
    mutate(DEFINITION_EN = ifelse(question_nl %in% names(custom_mapping), custom_mapping[question_nl], DEFINITION_EN),
           DEFINITION_EN = sub("…", "...", DEFINITION_EN)) %>%
    inner_join(questions_baseline, by = c("DEFINITION_EN")) %>%
    distinct(DEFINITION_EN, vl, .keep_all = TRUE)
  
  questions_missing <- questions_baseline %>% filter(!(DEFINITION_EN %in% question_overview_filtered_with_en$DEFINITION_EN)) %>% pull(DEFINITION_EN)

  question_overview_filtered_with_en %>% filter(!DEFINITION_EN %in% questions_baseline$DEFINITION_EN) %>% pull(question_nl)
  
  question_overview_filtered_with_en_ordered <- tibble(DEFINITION_EN=unique(questions_baseline$DEFINITION_EN)) %>% 
    right_join(question_overview_filtered_with_en) %>%
    group_by(question_id) %>%
    mutate(sample_size = sum(!is.na(phenoGeno[, question_id]))) %>% ungroup() %>%
    mutate(questionnaire_index = questionnaire_mapping[vl])
  
  write.table(
    question_overview_filtered_with_en_ordered %>% select(DEFINITION_EN, sample_size, questionnaire_index), 
    sep = ";", col.names= F, row.names = F, quote = F)
  # Perform method
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}