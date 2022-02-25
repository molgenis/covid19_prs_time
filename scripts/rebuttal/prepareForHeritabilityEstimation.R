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
library(tidyverse)

# Declare constants

# Declare function definitions

# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv=NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  # Process input
  workdir <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/"
  intermediatesdir <-  "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/longiIntermediates/"
  preparedDataFile <- "longitudinal.RData"
  pseudoIdLinkageFile <- "/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v1/phenotype_linkage_file_project_pseudo_id.txt"
  
  setwd(workdir)
  load(preparedDataFile)
  pseudoIdLinkage <- fread(pseudoIdLinkageFile)
  
  # Get a vector of outcomes. Here we use the validation results table.
  outcomes <- fread("interactionSummaryValidation.txt") %>%
    select(q2) %>%
    distinct() %>%
    pull(q2)
  
  # Perform method
  d <- vragenLong %>% as_tibble() %>%
    inner_join(pseudoIdLinkage, by = "PROJECT_PSEUDO_ID") %>%
    select(all_of(c("PSEUDOIDEXT", outcomes, "gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "vl")))
    
  wide <- d %>% 
    pivot_wider(
      id_cols = all_of(c("PSEUDOIDEXT","gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent")),
      names_from = vl,
      names_sep = "_",
      values_from = all_of(outcomes)
    ) %>%
    select(where(~!all(is.na(.x))))
  
  write.table(
    wide, 
    "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/plosgenetics_heritability_estimation/data/20211125/samples_genotyped_only/combined_flat_20211125.txt",
    sep = "\t",
    col.names = T,
    row.names = F,
    quote = F)
  
  column_description_outcomes <- d %>%
    select(all_of(c(outcomes, "vl"))) %>%
    pivot_longer(cols = all_of(outcomes), names_to = "q", values_to = "values") %>%
    select(-values) %>% distinct() %>%
    mutate(column_name = paste(q, vl, sep="_")) %>%
    filter(column_name %in% colnames(wide)) %>%
    inner_join(selectedQ %>% as_tibble(rownames = "q")) %>%
    select(column_name, c(questionnaire_label = vl, question_nl = Question, question_en = question_en, model_type = Type))
  
  write.table(
    column_description_outcomes, 
    "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/plosgenetics_heritability_estimation/data/20211125/samples_genotyped_only/combined_flat_20211125_column_description_outcomes.txt",
    sep = "\t",
    col.names = T,
    row.names = F,
    quote = F)
  
  column_description_participants_covariates <- tibble(column_name = c("PSEUDOIDEXT","gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent"),
                                                       column_description = c("lifelines participant identifier", 
                                                                              "gender or sex of the participant",
                                                                              "age of the participant",
                                                                              "age of the participant squared",
                                                                              "the sum of the number of household members",
                                                                              "wether or not the participant has children living with them in the same home",
                                                                              "wether or not the pariticpant suffers from a chronic illness"))
  
  write.table(
    column_description_participants_covariates, 
    "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/plosgenetics_heritability_estimation/data/20211125/samples_genotyped_only/combined_flat_20211125_column_description_participant_ids_and_covariates.txt",
    sep = "\t",
    col.names = T,
    row.names = F,
    quote = F)
  # Process output

  
  # Wegschrijven colnames en data.
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}