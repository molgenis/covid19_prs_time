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
library(glue)
library(ggtext)

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

# Declare constants
model1_output_excel_path <- "/Users/cawarmerdam/UMCG/FrankeSwertz - Documenten/Manuscripts/Covid19_Prs_Interactions/plos_genetics/Rebuttal/model1_output.xlsx"
model1_output_excel_path <- "/Users/cawarmerdam/UMCG/FrankeSwertz - Documenten/Manuscripts/Covid19_Prs_Interactions/plos_genetics/Rebuttal/model2_output_V2.xlsx"
dates <- tibble(vl = c("X1.0", "X2.0", "X3.0", "X4.0", "X5.0", "X6.0", "X7.0", "X8.0", "X9.0", "X10.0", 
                       "X11.0", "X12.0", "X13.0", "X14.0", "X15.0", "X15.5", "X16.0", "X16.5", "X17.0"),
                date = c("30/03/20", "6/04/20", "13/04/20", "20/04/20", "27/04/20", "4/05/20", "18/05/20",
                         "1/06/20", "15/06/20", "6/07/20", "13/07/20", "10/08/20", "7/09/20", "12/10/20",
                         "2/11/20", "17/11/20", "30/11/20", "14/12/20", "11/01/21"),
                vl2 = c(1:19))

questions <- tibble(
  question_id = c("q1", "q2", "q3", "q4", "q5", "q6", "q7"),
  question_en = c(
    "how much have you been concerned about the corona crisis in the past 7 days?",
    "how would you rate your quality of life over the last 7 days?",
    "i felt good / to what extent have you experienced following in the last 7 days:",
    "i felt tired / to what extent have you experienced following in the last 7 days:",
    "i was easily tired / to what extent have you experienced following in the last 7 days:",
    "i felt physically exhausted / to what extent have you experienced following in the last 7 days:",
    "Self reported positive SARS-CoV-2 PCR test"),
  question_label = c(
    "Concerned about the COVID-19 pandemic",
    "Quality of life",
    "Felt good",
    "Felt tired",
    "Was easily tired",
    "Felt physically exhausted",
    "Ever positive SARS-CoV-2 PCR test"))

# Declare function definitions
read_results_sheet <- function(question_id, path) {
  model1_output_excel <- read_excel(
    path, 
    sheet = question_id) %>%
  mutate(question_id = question_id)
  
  print(model1_output_excel)
  
  return(model1_output_excel)
}

# Linear res
lm_res <- function(x, y) {
  m <- lm(y ~ x)
  print(summary(m))
  
  p_val <- summary(m)$coefficients["x", "Pr(>|t|)"]
  r_square <- summary(m)$r.squared
  
  print(p_val)
  
  label <- glue("*r*<sup>2</sup> = {round(r_square, 2)}, *p*-value = {formatC(p_val,2)}")
  return(label)
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
  results_long <- bind_rows(
    mapply(read_results_sheet, 
           questions$question_id, 
           MoreArgs = list(path = model1_output_excel_path),
           USE.NAMES = F, SIMPLIFY = F)
  ) %>%
  inner_join(questions, by = c("question_id")) %>%
  inner_join(dates, by = c("Phenotype" = "vl")) %>%
  mutate(date = as.Date(date, "%d/%m/%y"))
  
  results_formatted <- results_long %>% select(c(
    "Question label" = "question_label",
    "Full question" = "question_en",
    "Questionnaire" = "vl2", 
    "Heritability estimate (h2)" = "H2", 
    "Standard error of heritability estimate" = "SE_H2", 
    "P-value of heritability estimate" = "P_H2", 
    "Variance explained by shared environment (c2)" = "C2", 
    "Standard error of variance explained by shared environment"= "SE_C2", 
    "P-value of variance explained by shared environment" = "P_C2"))
  
  write.table(
    results_formatted, 
    "/Users/cawarmerdam/Documents/projects/lifelines_covid_research/analysis/plosgenetics_heritability_estimation/heritability_baseline_results_20220117.txt", 
    col.names = T, row.names = F, quote = F, sep = "\t")
  
  results_long_filtered <- results_long %>% filter(question_en != "Self reported positive SARS-CoV-2 PCR test")
  
  # H2
  
  results_h2_summarised <- results_long_filtered %>% group_by(question_label) %>%
    summarise(label = lm_res(date, H2))
  
  # Perform method
  time_series_facets_h2 <- ggplot(results_long_filtered, aes(x = date, y = H2)) +
    geom_smooth(method = "lm", se = TRUE, fill = "grey80",
                size = 1.5 / (ggplot2::.pt * 72.27/96)) +
    geom_errorbar(aes(ymin = H2 - SE_H2, ymax = H2 + SE_H2), 
                  color = "gray50", 
                  size = 0.75 / (ggplot2::.pt * 72.27/96)) +
    geom_point() +
    scale_x_date(date_labels = "%b %Y") +
    ylab("Heritability estimates and standard error of heritabilities") +
    xlab("Dates") +
    facet_wrap(vars(question_label), scales = "free_y", nrow = 3)

  time_series_facets_h2 + 
    geom_richtext(
      inherit.aes = F,
      data = results_h2_summarised,
      aes(label = label),
      hjust = 1.02, vjust = -0.4,
      # remove label background and outline
      fill = NA, label.color = NA,
      # remove label padding, since we have removed the label outline
      label.padding = grid::unit(rep(0, 4), "pt"),
      size = 7 / (ggplot2::.pt * 72.27/96),
      x = Inf,
      y = -Inf
    )
  
  ggsave("/Users/cawarmerdam/Documents/projects/lifelines_covid_research/analysis/plosgenetics_heritability_estimation/model_2_heritability.pdf",
         width = 6, height = 8, units = "in")
  ggsave("/Users/cawarmerdam/Documents/projects/lifelines_covid_research/analysis/plosgenetics_heritability_estimation/model_2_heritability.png",
         width = 6, height = 8, units = "in")
  
  # C2
  
  results_c2_summarised <- results_long_filtered %>% group_by(question_label) %>%
    summarise(label = lm_res(date, C2))
  
  # Perform method
  time_series_facets_c2 <- ggplot(results_long_filtered, aes(x = date, y = C2)) +
    geom_smooth(method = "lm", se = TRUE, fill = "grey80",
                size = 1.5 / (ggplot2::.pt * 72.27/96)) +
    geom_errorbar(aes(ymin = C2 - SE_C2, ymax = C2 + SE_C2), 
                  color = "gray50", 
                  size = 0.75 / (ggplot2::.pt * 72.27/96)) +
    geom_point() +
    scale_x_date(date_labels = "%b %Y") +
    ylab("Variance explained by shared environment and standard error") +
    xlab("Dates") +
    facet_wrap(vars(question_label), scales = "free_y", nrow = 3)
  
  time_series_facets_c2 + 
    geom_richtext(
      inherit.aes = F,
      data = results_c2_summarised,
      aes(label = label),
      hjust = 1.02, vjust = -0.4,
      # remove label background and outline
      fill = NA, label.color = NA,
      # remove label padding, since we have removed the label outline
      label.padding = grid::unit(rep(0, 4), "pt"),
      size = 7 / (ggplot2::.pt * 72.27/96),
      x = Inf,
      y = -Inf
    )
  
  ggsave("/Users/cawarmerdam/Documents/projects/lifelines_covid_research/analysis/plosgenetics_heritability_estimation/model_2_environment.pdf",
         width = 6, height = 8, units = "in")
  ggsave("/Users/cawarmerdam/Documents/projects/lifelines_covid_research/analysis/plosgenetics_heritability_estimation/model_2_environment.png",
         width = 6, height = 8, units = "in")
  
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}