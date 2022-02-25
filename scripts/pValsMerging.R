library(data.table)
library(tidyverse)

pVals <- fread("../covid_combined_matrix_p_values_first_values_25-05-2021.txt")

processedPVals <- masterTable %>% select("Question en", "X1") %>%
  full_join(pVals, by = c("X1" = "V1"))

write.table(processedPVals, "../covid_combined_matrix_p_values_first_values_25-05-2021_with_question_en.txt", col.names = T, row.names = F, quote = F, sep = "\t")
