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

load_files <- function(base_path, files_path = "files.txt") {
  files <- fread(files_path, header = F, sep = "\t") %>%
    filter(V1 != "") %>%
    pull(V1)
  
  my_tibble <- bind_rows(mapply(function(file_path) {
    pgs <- str_split_fixed(file_path, "/", n = Inf)[2]
    
    dat <- fread(file.path(base_path, file_path), header = T, sep = "\t") %>% 
      pivot_longer(-question, names_to = "vl", values_to = "value") %>%
      mutate(pgs = pgs)
    
    return(dat)
  }, files, USE.NAMES = F, SIMPLIFY = F)) %>%
    mutate(vl = paste0("X", vl))
  
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
  matrixPath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/covid_combined_matrix_z_scores_first_values_25-05-2021.txt"
  mat <- fread(matrixPath) %>%
    mutate(V1 = case_when(V1 == "bezoek kreeg van vrienden en/of familie / hoe goed kon u zich aan de 1,5 meter maatregelen houden toen u" ~ "bezoek kreeg van vrienden en/of familie / hoe goed kon u zich aan de 1,5 meter maatregelen houden toen u…",
                          V1 == "thuis bezoek ontving / hoe goed kon u zich aan de 1,5 meter maatregelen houden toen u" ~ "thuis bezoek ontving / hoe goed kon u zich aan de 1,5 meter maatregelen houden toen u…",
                          TRUE ~ V1))
  
  ugliPathPval <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/pgs_correlation_output_linear_model_all_quest_v30_final_20_all_genome_participants_test_with_cof/all_part_lin_new_ugli/BMI/covid_export_correlations_BMI_p_values_allParti_new_lin_ugli_20-05-2021.txt"
  ugliPathCorr <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/pgs_correlation_output_linear_model_all_quest_v30_final_20_all_genome_participants_test_with_cof/all_part_lin_new_ugli/BMI/covid_export_correlations_BMI_correlations_allParti_new_lin_ugli_20-05-2021.txt"

  ugliMat <- fread(ugliPathPval, header = T) %>%
    pivot_longer(cols=-question, names_to = "vl", values_to = "p_val_ugli") %>%
    mutate(vl = paste0("X", vl))
  
  ugliMat <- fread(ugliPathCorr, header = T) %>%
    pivot_longer(cols=c(-question, -`Question answers`), names_to = "vl", values_to = "corr_ugli") %>%
    mutate(vl = paste0("X", vl)) %>%
    inner_join(ugliMat)
  
  cytoPathPval <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/pgs_correlation_output_linear_model_all_quest_v30_final_20_all_genome_participants_test_with_cof/all_part_lin_new_cyto/BMI/covid_export_correlations_BMI_p_values_allParti_new_lin_cyto_23-05-2021.txt"
  cytoPathCorr <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/pgs_correlation_output_linear_model_all_quest_v30_final_20_all_genome_participants_test_with_cof/all_part_lin_new_cyto/BMI/covid_export_correlations_BMI_correlations_allParti_new_lin_cyto_23-05-2021.txt"
  
  cytoMat <- fread(cytoPathPval, header = T) %>%
    pivot_longer(cols=-question, names_to = "vl", values_to = "p_val_cyto") %>%
    mutate(vl = paste0("X", vl))
  
  cytoMat <- fread(cytoPathCorr, header = T) %>%
    pivot_longer(cols=c(-question, -`Question answers`), names_to = "vl", values_to = "corr_cyto") %>%
    mutate(vl = paste0("X", vl)) %>%
    inner_join(cytoMat)
  
  metaPathPval <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/analyses/20211125/regenerate_baseline_output/output/BMI/covid_meta_analysis_pvalues_BMI_20-01-2022.txt"
  metaPathCorr <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/analyses/20211125/regenerate_baseline_output/output/BMI/covid_meta_analysis_correlations_BMI_20-01-2022.txt"
  
  metaMat <- fread(metaPathPval, header = T) %>%
    pivot_longer(cols=-question, names_to = "vl", values_to = "p_val_meta") %>%
    mutate(vl = paste0("X", vl))
  
  metaMat <- fread(metaPathCorr, header = T) %>%
    pivot_longer(cols=-question, names_to = "vl", values_to = "corr_meta") %>%
    mutate(vl = paste0("X", vl)) %>%
    inner_join(metaMat)
  
  combined_table <- ugliMat %>% full_join(cytoMat) %>% full_join(metaMat) 
  
  table_annotated <- question_overview_filtered_with_en %>% inner_join(combined_table, by = c("question_nl" = "question", "vl" = "vl"))
  
  table_annotated_2 <- tibble(DEFINITION_EN=unique(baseline_pval_labels$V1)) %>% 
    left_join(table_annotated)
    
  
  write.table(
    table_annotated_2,
    file.path(base_path, "covid_combined_BMI_results_cyto_ugli_meta_baseline_18_02_2022.txt"),
    sep = ";", col.names= T, row.names = F, quote = F)
  
  input_question_overview <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/combined_questionnaires_v22_14-05-2021_genome_fitered_participants_filtered_correct_recoded/questionnaire_subset_recode_info_overview_14-05-2021.txt" # datafrane with the question ids over time
  workdir <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/"
  qOverviewFile <- "quest_overview_nl_new_quest17_codes_updated_14-days_v2.txt"
  questions_baseline <- fread("baseline_questions_with_label.txt", header = T, sep = "\t",
                              col.names = c("Outcome_label", "DEFINITION_EN"))
  
  qOverview <- as.matrix(read.delim(qOverviewFile, stringsAsFactors = F, row.names = 1))
  vls <- colnames(qOverview)[-c(20,21)]
  
  baseline_pval_labels <- fread("baseline_pval_labels.txt", sep = "\t", header = F) %>%
    mutate(V1 = sub("…", "...", V1))
  
  base_path <- "analyses/20211125/regenerate_baseline_output/"
  
  
  question_overview <- fread(input_question_overview) %>%
    select(-V1)
  
  question_overview_tibble <- as_tibble(qOverview, rownames = "question_nl") %>% 
    pivot_longer(cols=starts_with("X"), names_to = "vl", values_to = "question_id") %>%
    filter(question_id != "") %>%
    full_join(question_overview, by = c("question_id" = "question_id"))
  
  question_overview_filtered <- question_overview_tibble %>% filter(
    (is_first_time_asked == T & question_nl != "Positive tested cumsum") 
    | question_id == "covt17_positive_tested_cumsum",
    question_nl != "hoe waardeert u uw kwaliteit van leven over de afgelopen 7 dagen?")
  
  configurationFilePath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/barometer_statistics/coronabarometer-statistics/config.v3.json"
  configuration <- fromJSON(file = configurationFilePath)
  covidVariablesPerWeek <- configuration$covidResFilePerWeekVariables
  covidResFileTypeVersionPerWeek <- configuration$fileTypeVersionPerWeek
  
  variablesList <- list()
  
  for (weekLabel in names(covidVariablesPerWeek)[1:19]) {
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
    filter(variable_name %in% question_overview_filtered$question_id)
  
  question_overview_filtered_with_en <- question_overview_filtered %>% 
    left_join(variables, c("question_id" = "variable_name")) %>%
    mutate(DEFINITION_EN = ifelse(question_nl %in% names(custom_mapping), custom_mapping[question_nl], definition_en),
           DEFINITION_EN = sub("…", "...", definition_en)) %>%
    select(DEFINITION_EN, vl, question_id, question_nl) %>% 
    filter(question_nl != "hebt u een (pleeg)kind of (pleeg)kinderen van 18 jaar of jonger?_y",
           question_nl != "hebt u een (pleeg)kind of (pleeg)kinderen van 18 jaar of jonger?") %>%
    distinct(question_nl, vl, .keep_all = TRUE)
  
  table_annotated_2 <- read_delim("~/Documents/projects/lifelines_covid_research/analysis/pgs_correlations/analyses/20211125/compare_bmi_baseline_results/covid_combined_BMI_results_cyto_ugli_meta_baseline_18_02_2022.txt", delim = ";")
  # Process output
  
  ggplot(table_annotated_2 %>% filter(p_val_meta < 1.0004001600640257e-05), aes(x = -log10(p_val_ugli), y = -log10(p_val_cyto))) + geom_point()
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}