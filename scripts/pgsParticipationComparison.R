#!/usr/bin/env Rscript

## Load libraries

library(data.table)
library(rjson)
library(tidyverse)

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

## Main

prsGsa <- read_delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/input_PGS_data_ugli_v4/PGS_combined_ugli_07-04-2021.txt", delim = "\t")
prsCyto <- read_delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/input_PGS_data_cyto_v4_duplicate_filtered/PGS_combined_cyto_duplicate_from_ugli_removed_07-04-2021.txt", delim = "\t")
traitTable <- read_delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/trait-gwas-mapping-pruned-07-04-2021.txt", delim = "\t")
includedSamplesFilePath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/analyses/pgs-validation/includedSamples.txt"
invitedSamplesFilePath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/data/raw/participants_program/deelname_covid-19_onderzoeksprogramma.dat"
selectedTraits <- traitTable %>% pull(trait)

# Configuration file that lists COVID-19 related file paths
configuration <- fromJSON(file = "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/barometer_statistics/coronabarometer-statistics/config.v3.json")
pseudoIdLinkageFile <- configuration$pseudoIdLinkageFile
pseudoIdLink <- fread(pseudoIdLinkageFile)

out <- "../output"

isRespondent <- read_csv(
  includedSamplesFilePath,
  col_names = c("PROJECT_PSEUDO_ID")) %>%
  inner_join(pseudoIdLink) %>% select(PSEUDOIDEXT) %>%
  mutate(isRespondent = 1L)

sampleGroups <- fread(invitedSamplesFilePath) %>%
  rename(isInvited = "COVID19_PROGRAMME") %>%
  left_join(isRespondent, by = c("PSEUDOIDEXT")) %>%
  mutate(PSEUDOIDEXT = as.character(PSEUDOIDEXT))

pgsFull <- bind_rows(list("GSA" = prsGsa, "HumanCytoSNP" = prsCyto), .id = "ARRAY") %>%
  inner_join(pseudoIdLink) %>%
  select(-PROJECT_PSEUDO_ID) %>%
  pivot_longer(c(-PSEUDOIDEXT, -ARRAY), names_to = "trait", values_to = "PGS") %>%
  mutate(PSEUDOIDEXT = as.character(PSEUDOIDEXT)) %>%
  left_join(sampleGroups, by = "PSEUDOIDEXT") %>%
  mutate(isRespondentVsNot = case_when(
    isRespondent == 1 ~ "Respondent", 
    TRUE ~ "Not respondent"),
         isRespondentVsInvited = factor(case_when(
           isInvited == 1 & isRespondent == 1 ~ "Respondent", 
           isInvited == 1 ~ "Invited"))) %>%
  filter(!is.na(isRespondentVsInvited))

tTestsRespondentVsInvitedGsa <- bind_rows(mapply(function(thisTrait, pgsTibble) {
  message(thisTrait)
  thisTraitData = pgsTibble %>% filter(trait == thisTrait)
  tTestRespondentVsInvited <- t.test(thisTraitData$PGS ~ thisTraitData$isRespondentVsInvited)
  print(tTestRespondentVsInvited$method)
  return(data.frame(pValueRespondentVsInvited = tTestRespondentVsInvited[["p.value"]],
                    meanIncluded = tTestRespondentVsInvited[["estimate"]][["mean in group Respondent"]],
                    meanInvitedNotIncluded = tTestRespondentVsInvited[["estimate"]][["mean in group Invited"]]))
}, selectedTraits, list(pgsTibble = pgsFull %>% filter(ARRAY == "GSA")), SIMPLIFY = F, USE.NAMES = T), .id = "trait")

tTestsRespondentVsInvitedCyto <- bind_rows(mapply(function(thisTrait, pgsTibble) {
  message(thisTrait)
  thisTraitData = pgsTibble %>% filter(trait == thisTrait)
  tTestRespondentVsInvited <- t.test(thisTraitData$PGS ~ thisTraitData$isRespondentVsInvited)
  print(tTestRespondentVsInvited$method)
  return(data.frame(pValueRespondentVsInvited = tTestRespondentVsInvited[["p.value"]],
                    meanIncluded = tTestRespondentVsInvited[["estimate"]][["mean in group Respondent"]],
                    meanInvitedNotIncluded = tTestRespondentVsInvited[["estimate"]][["mean in group Invited"]]))
}, 
selectedTraits, list(pgsTibble = pgsFull %>% filter(ARRAY == "HumanCytoSNP")), 
SIMPLIFY = F, USE.NAMES = T), .id = "trait")

tTestResults <- bind_rows(list("GSA" = tTestsRespondentVsInvitedGsa, 
                               "HumanCytoSNP" = tTestsRespondentVsInvitedCyto),
                          .id = "array")

# Write summaries for PGSs
write.table(tTestResults, 
            file.path(out, paste0("pgsCoronaSamplesTtestResults_", format(Sys.Date(), "%Y%m%d"), ".txt")),
            row.names = F, col.names = T, quote = F, sep = "\t")

pgsFull <- pgsFull %>% inner_join(traitTable, by = "trait")

# Plot bloxplots for all traits
pdf(file.path(out, paste0("pgsCoronaSamplesTtestResults_", format(Sys.Date(), "%Y%m%d"), ".pdf")), 
    useDingbats = FALSE, width = 8.5, height = 11)

par(xpd = NA)

global_labeller <- labeller(
  traitLabel = label_wrap_gen(24),
  .default = label_both
)

ggplot(pgsFull %>% filter(!is.na(isRespondentVsInvited) & ARRAY == "GSA" & trait %in% selectedTraits), 
       aes(x = isRespondentVsInvited, y = PGS)) +
  geom_violin(trim=T, size = 1 / (ggplot2::.pt * 72.27/96), fill = "gray", colour = "grey20") +
  geom_boxplot(width=0.1, size = 1 / (ggplot2::.pt * 72.27/96), outlier.shape = NA, colour = "grey10") +
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("Respondent", "Invited")),
                     bracket.nudge.y = 1, vjust = -0.5) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_x_discrete(breaks = c("Respondent", "Invited"),
                   labels = c("Included samples", "Invited, not \nincluded samples")) +
  xlab("Included samples compared to samples that \nhave been invited but were not included") +
  ylab("Polygenic scores") +
  labs(title = "A: GSA") +
  facet_wrap(~traitLabel, ncol = 4, scales = "free_y",
             labeller = global_labeller) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggplot(pgsFull %>% filter(!is.na(isRespondentVsInvited) & ARRAY == "HumanCytoSNP" & trait %in% selectedTraits), 
       aes(x = isRespondentVsInvited, y = PGS)) +
  geom_violin(trim=T, size = 1 / (ggplot2::.pt * 72.27/96), fill = "gray", colour = "grey20") +
  geom_boxplot(width=0.1, size = 1 / (ggplot2::.pt * 72.27/96), outlier.shape = NA, colour = "grey10") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("Respondent", "Invited")),
                     bracket.nudge.y = 1, vjust = -0.5) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_x_discrete(breaks = c("Respondent", "Invited"),
                   labels = c("Included samples", "Invited, not \nincluded samples")) +
  xlab("Included samples compared to samples that \nhave been invited but were not included") +
  ylab("Polygenic scores") +
  labs(title = "B: HumanCytoSNP") +
  facet_wrap(~traitLabel, ncol = 4, scales = "free_y",
             labeller = global_labeller) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

dev.off()
