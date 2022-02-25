#!/usr/bin/env Rscript

## Load libraries

library(data.table)
library(rjson)
library(tidyverse)
library(ggpubr)
library(ggsignif)
library(viridis)

## Constants

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

## Function definitions
to_matrix <- function(tbl, id_cols=NULL, names_from=name, values_from=value, id_cols_sep = "_", ...){

  long <- tbl %>% 
    pivot_wider(id_cols = {{id_cols}}, names_from = {{names_from}}, values_from = {{values_from}}, ...)
  
  if (is.null(id_cols)) {
    # id_cols_char <- colnames(tbl %>% select(-{{names_from}}, -{{values_from}}))
    # colnames_vect <- colnames(long)[(length(id_cols_char)+1):ncol(long)]
    # print(length(colnames_vect))
    # rownames_col <- do.call(paste, c(long %>% select(all_of(id_cols_char)), sep=id_cols_sep))
    # print(length(rownames_col))
    # print(long %>% select(-any_of(id_cols_char)) %>% as.data.frame())
    # 
    # return(matrix(long %>% select(-any_of(id_cols_char)) %>% as.data.frame(), dimnames = list(
    #   make.names(rownames_col),
    #   make.names(colnames_vect)
    # )))
    stop()
  } else {
    mat <- long %>% select(-any_of({{id_cols}})) %>% as.matrix()
    rownames(mat) <- do.call(paste, c(long %>% select(all_of(id_cols)), sep=id_cols_sep))
    colnames(mat) <- colnames(long %>% select(-any_of(id_cols)))
    return(mat)
  }
}

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

to_pretty_string <- function(named_vector) {
  paste0(sapply(names(named_vector), function(this_name) {
    val <- named_vector[this_name]
    val_as_character <- ifelse(val < 10, "<10", as.character(val))
    paste0(this_name, " = ", val_as_character)
    }), collapse = ", ")
}

get_outcome_comparison <- function(model_type, dat, group_col, value_col, question_label, answer_labels, comparison, make_plot = T) {
  res <- list()
  
  res[["method"]] <- NA_character_
  
  res[[paste0("frequencies_", comparison[1])]] <- NA_character_
  res[[paste0("frequencies_", comparison[2])]] <- NA_character_
  res[["fisher_exact_p"]] <- NA_complex_
  res[["fisher_exact_estimate"]] <- NA_complex_
  
  res[[paste0("mean_", comparison[1])]] <- NA_complex_
  res[[paste0("mean_", comparison[2])]] <- NA_complex_
  res[["t_test_p"]] <- NA_complex_
  
  print(model_type)
  # input 'dat' should be 
  if (model_type == "" | model_type == "binomial" | model_type == "ordinal-ordered" | model_type == "ordinal-ordered-turned") {
    # For the binomial and the ordinal data, we can simply make bar charts
    # and do a fisher-exact test to compare both of the groups
    
    dat_summarised <- dat %>% 
      filter(!is.na(get(value_col))) %>%
      group_by(across(all_of(c(group_col, value_col)))) %>% 
      summarize(freq = n()) %>% ungroup() %>%
      group_by(across(all_of(c(group_col)))) %>%
      mutate(percent = freq / sum(freq) * 100) %>%
      ungroup()
    
    matrix <- dat_summarised %>% 
      to_matrix(id_cols = group_col, names_from = value_col, values_from = freq)
    
    #print(matrix)
    fisher.results <- fisher.test(matrix, simulate.p.value=TRUE, B=1e7)
    res[["method"]] <- fisher.results$method
    res[[paste0("frequencies_", comparison[1])]] <- to_pretty_string(matrix[1,])
    res[[paste0("frequencies_", comparison[2])]] <- to_pretty_string(matrix[2,])
    res[["fisher_exact_p"]] <- fisher.results$p.value
    if("estimate" %in% names(fisher.results)) {
      res[["fisher_exact_or"]] <- as.numeric(fisher.results$estimate)
    }
    #res[["p.value"]] <- 0.4
    
    print(res)
    
    #print(dat_summarised)
    #print(list(unname(comparison)))
    
    if (make_plot == T) {
      p <- ggplot(dat_summarised, aes(fill=get(value_col), y=percent, x=get(group_col))) +
        geom_col(size = 1 / (ggplot2::.pt * 72.27/96), colour = "grey20") +
        scale_x_discrete(breaks = unname(comparison), labels = names(comparison), name = NULL) +
        scale_y_continuous(breaks = seq(0, 100, 25)) +
        scale_fill_viridis(discrete = T, labels = function(x) str_wrap(x, width = 12)) +
        labs(title=str_wrap(question_label, width = 32), y="Percentage", x = "Group") + 
        coord_cartesian(ylim = c(1, 120)) +
        guides(fill=guide_legend(title=NULL)) +
        geom_signif(y_position = 105, annotation = formatC(res$fisher_exact_p, digits=1), 
                    comparisons = list(unname(comparison))) +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
              plot.title = element_text(size = 8, face = "plain"))
      res[["plot"]] <- p
    }
    
  } else if (model_type == "gaussian") {
    # For the gaussian data, we can simply do a t-test
    
    print(dat)
    
    print(head(as.numeric(dat %>% filter(!!as.name(group_col) == comparison[1]) %>% pull(!!as.name(value_col)))))
    print(head(as.numeric(dat %>% filter(!!as.name(group_col) == comparison[2]) %>% pull(!!as.name(value_col)))))

    t.results <- t.test(
      x = as.numeric(dat %>% filter(!!as.name(group_col) == comparison[1]) %>% pull(!!as.name(value_col))), 
      y = as.numeric(dat %>% filter(!!as.name(group_col) == comparison[2]) %>% pull(!!as.name(value_col))),
      na.action = na.omit)

    # res[["input"]] <- c(
    #   mean(as.numeric(dat[dat[group_col] == comparison[1], value_col], na.rm = T)),
    #   mean(as.numeric(dat[dat[group_col] == comparison[2], value_col], na.rm = T)))
    # 
    res[["method"]] <- t.results$method
    res[["t_test_p"]] <- t.results$p.value
    res[[paste0("mean_", comparison[1])]] <- as.numeric(t.results$estimate[1])
    res[[paste0("mean_", comparison[2])]] <- as.numeric(t.results$estimate[2])
    
    print(res)
    
    if (make_plot == T) {
      y_scale <- scale_y_continuous(expand = expansion(mult = c(0.02, .1)), name = NULL)
      
      if (!is.null(answer_labels)) {
        y_scale <- scale_y_continuous(
          breaks = 1:length(answer_labels),
          labels = answer_labels,
          expand = expansion(mult = c(0.02, .1)), name = NULL)
      }
      
      p <- ggplot(dat, aes(y=as.numeric(get(value_col)), x=get(group_col))) +
        geom_violin(trim=FALSE, size = 1 / (ggplot2::.pt * 72.27/96), fill = "gray", colour = "grey20", adjust=2) +
        geom_boxplot(width=0.1, outlier.shape = NA, colour = "grey10") +
        scale_x_discrete(breaks = unname(comparison), labels = names(comparison), name = NULL) +
        y_scale +
        labs(title=str_wrap(question_label, width = 32), x = "Group") + 
        geom_signif(annotations = formatC(res$t_test_p, digits=1), 
                    comparisons = list(unname(comparison))) +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
              plot.title = element_text(size = 8, face = "plain"))
      
      res[["plot"]] <- p
    }
  }
  
  return(res)
}

## Main

phenoGenotypedPath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/combined_questionnaires_v22_14-05-2021_genome_fitered_participants_filtered_correct_recoded/questionnaire_subset_participants_filtered_recoded_answers_14-05-2021"
phenoUngenotypedPath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/rebuttal/combined_questionnaires_v10_06-05-2021_nogeno_fitered_participants_filtered_correct_recoded/questionnaire_subset_participants_filtered_recoded_answers_25-10-2021"

invitedSamplesFilePath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/data/raw/participants_program/deelname_covid-19_onderzoeksprogramma.dat"

phenoGeno <- load_pheno_file(phenoGenotypedPath)
phenoNogeno <- load_pheno_file(phenoUngenotypedPath)

#Constants
startdate <- as.Date("30/03/2020","%d/%m/%Y")
confounders <- c("gender_recent", "age_recent", "age2_recent", "chronic_recent", "household_recent", "have_childs_at_home_recent")

# dataframe with the combined questionnaire data
input_question_overview <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/combined_questionnaires_v22_14-05-2021_genome_fitered_participants_filtered_correct_recoded/questionnaire_subset_recode_info_overview_14-05-2021.txt" # datafrane with the question ids over time
model_selection_file <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/Question_selection_PRS_paper-question_selection_20210525.tsv"

workdir <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/"
sample_inclusion_per_questionnaire_genotyped_file <- "inclusionPerVl.txt"
sample_inclusion_per_questionnaire_ungenotyped_file <- "inclusionPerVl_NoGenotypeData.txt"
qOverviewFile <- "quest_overview_nl_new_quest17_codes_updated_14-days_v2.txt"

##Load and format questions meta data
qOverview <- as.matrix(read.delim(qOverviewFile, stringsAsFactors = F, row.names = 1))
vls <- colnames(qOverview)[-c(20,21)]

sample_inclusion_per_questionnaire_genotyped <- fread(sample_inclusion_per_questionnaire_genotyped_file)
sample_inclusion_per_questionnaire_ungenotyped <- fread(sample_inclusion_per_questionnaire_ungenotyped_file)

# Load the questions. We ideally loop through all of these
model_selection <- fread(model_selection_file)
question_overview <- fread(input_question_overview)

question_ids <- apply(model_selection, 1, function(row) {
  question <- as.character(row["Question"])

  if (!question %in% rownames(qOverview)) {
    print(question)
    return(NULL)
  }
  
  ids <- qOverview[question, vls]
  if (question == "Positive tested cumsum") {
    question_row <- question_overview[question_overview$question_id == "covt17_positive_tested_cumsum"]
  } else {
    question_row <- question_overview[question_overview$question_id %in% ids 
                                      & question_overview$is_first_time_asked == TRUE]
  }
  question_row$Question <- question
  return(question_row)
})

question_table <- merge(model_selection, Reduce(rbind, question_ids), by = "Question")

out <- mapply(function(row_index) {
  row <- question_table[row_index,]
  question <- row$Question
  model_type <- row$Type
  question_id <- row$question_id
  question_en <- row$question_en
  answer_options <- row$Answers_options_plots
  
  print(row_index)
  print(question)
  # print(model_type)
  # print(question_id)
  # print(question_en)
  # print(answer_options)
  
  # Grab the data
  dat <- bind_rows(list(
    "genotyped" = phenoGeno %>% select(all_of(c(question_id, "PROJECT_PSEUDO_ID"))),
    "ungenotyped" = phenoNogeno %>% select(all_of(c(question_id, "PROJECT_PSEUDO_ID")))), .id = "group") %>%
    filter(!is.na(get(question_id)))
  
  answer_labels = NULL
   
  if (answer_options != "" & answer_options != "{}") {
    answer_options_list <- unlist(fromJSON(gsub("\"\"", "\"", answer_options)))
    results_as_values <- format(as.double(dat[,question_id]), nsmall = 1, justify = "none", trim = T)
    results_as_labels <- answer_options_list[results_as_values]
    answer_labels = unname(answer_options_list)[order(as.double(names(answer_options_list)))]
    dat[,question_id] <- ordered(results_as_labels, levels = answer_labels)
  }
  
  dat <- dat %>% as_tibble() %>%
    filter(!is.na(get(question_id)))
  
  # print(dat)
  
  res <- get_outcome_comparison(
    model_type = model_type, dat = dat, value_col = question_id, group_col = "group",
    question_label = question_en, answer_labels = answer_labels, 
    comparison = c("Genotyped" = "genotyped", "Not genotyped" = "ungenotyped"))
  
  res[["question_en"]] <- question_en
  #print(res)
  
  return(res)
  
}, seq(nrow(question_table)))

#dev.off()

this_plots <- out["plot",]
multi.page <- ggarrange(plotlist = this_plots,
                        ncol = 2, nrow = 2)



ggexport(multi.page, filename = "genotypedVsUngenotyped_20211117.pdf")

out2 <- out
out2 <- apply(out2, 2, function(i) {
  i[["plot"]] <- NULL
  return(i)})

dd  <-  as.data.frame(t(matrix(unlist(out2), nrow=length(unlist(out2[1])))))
colnames(dd) <- names(out2[[1]])

my_table <- as_tibble(dd) %>% 
  select(question_en, method, frequencies_genotyped, frequencies_ungenotyped, fisher_exact_p, fisher_exact_estimate, mean_genotyped, mean_ungenotyped, t_test_p) %>%
  mutate(method = str_replace_all(method, "[\r\n\t]" , ""))

write.table(my_table, "genotypedVsUngenotyped_20211117.txt", row.names = F, col.names = T, quote = F, sep = "\t")
