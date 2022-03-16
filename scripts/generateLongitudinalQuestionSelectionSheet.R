library(data.table)
library(tidyverse)

setwd("/Users/cawarmerdam/Documents/projects/lifelines_covid_research/analysis/pgs_correlations")
recodingDataFrame <- fread("robert/covid19_prs_time/data/recode_info_overview_14-05-2021.txt", header = T, quote = "")
gwasSelectionDataFrameOld <- fread("robert/covid19_prs_time/data/gwasses_to_perform_filtered.txt", header = T)
gwasSelectionDataFrame <- fread("robert/covid19_prs_time/data/gwasses_to_perform_filtered_new.txt", header = T)
questionOverviewCodesLarge <- fread("robert/covid19_prs_time/data/question_overview_codes_include_19.txt", header = T, quote = "")
recodingDataFrame2 <- fread("robert/covid19_prs_time/data/recode.txt", header = T, quote = "")

  
myFun <- function(x) {
  sorted <- "{}"
  if (!is.na(x)) {
    aAsList <- fromJSON(x)
    names(aAsList) <- as.character(as.integer(names(aAsList)))
    sorted <- toJSON(aAsList[order(as.integer(names(aAsList)),decreasing=FALSE)])
  }
  return(sorted)
}

filteredQuestionOverview <- questionOverviewCodesLarge %>%
  select(V1, `Question en`, `Question answers en`) %>%
  mutate(`Question answers en` = sapply(`Question answers en`, myFun)) %>%
  distinct()

data2 <- gwasSelectionDataFrame %>%
  left_join(recodingDataFrame, by = c("V1" = "question_label_nl")) %>%
  rename("question_label_nl" = "V1") %>%
  group_by(question_label_nl) %>%
  select(question_label_nl, model_type, answer_options_en_after_recode_json) %>%
  mutate(TypeLongnitudinal = case_when(model_type == "gaussian" ~ "gaussian",
                                       TRUE ~ "binomial"),
         Mixed = case_when(question_label_nl == "Positive tested cumsum" ~ F,
                           TRUE ~ T),
         answer_options_en_after_recode_json = sapply(answer_options_en_after_recode_json, myFun),
         `Number of timepoints` = n()) %>%
  distinct() %>%
  left_join(filteredQuestionOverview,
             by = c("question_label_nl" = "V1")) %>%
  mutate(`Answer options` = case_when(!is.na(`Question answers en`) & answer_options_en_after_recode_json == "{}" ~ `Question answers en`,
                                      !is.na(answer_options_en_after_recode_json) & answer_options_en_after_recode_json != "{}" ~ answer_options_en_after_recode_json,
                                      TRUE ~ "{}")) %>%
  select(-`Question answers en`, -answer_options_en_after_recode_json) %>%
  inner_join(recodingDataFrame2 %>% select(Question, `Recode value labels`), by = c("question_label_nl" = "Question"))

Encoding(data[,"question_label_nl", drop=T]) <- "UTF-16"

write.table(data, "robert/covid19_prs_time/data/longitudinalQuestionSelection_20210525.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
