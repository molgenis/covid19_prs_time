#!/usr/bin/env Rscript

## Load libraries

library(data.table)
library(rjson)
library(tidyverse)
library(ggpubr)

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

phenoBaselinePath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/combined_questionnaires_v22_14-05-2021_genome_fitered_participants_filtered_correct_recoded/questionnaire_subset_participants_filtered_recoded_answers_14-05-2021"
phenoLongitudinalPath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/combined_questionnaires_v22_14-05-2021_genome_fitered_participants_filtered_correct_recoded/questionnaire_subset_participants_filtered_recoded_answers_longitudinal_filtered_15-05-2021"

prsGsaFile <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/input_PGS_data_ugli_v4/PGS_combined_ugli_07-04-2021.txt"
prsCytoFile <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/input_PGS_data_cyto_v4_duplicate_filtered/PGS_combined_cyto_duplicate_from_ugli_removed_07-04-2021.txt"
invitedSamplesFilePath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/data/raw/participants_program/deelname_covid-19_onderzoeksprogramma.dat"

# Get traits
traitTable <- read_delim("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/trait-gwas-mapping-pruned-07-04-2021.txt", delim = "\t")
selectedTraits <- traitTable %>% pull(trait)

# Read pheno samples
pheno <- read_delim(paste0(phenoBaselinePath, ".txt"), delim = "\t", quote = "", guess_max = 100000)
baselineSamples <- unique(pheno$X1)
rm(pheno)
gc()

pheno <- read_delim(paste0(phenoLongitudinalPath, ".txt"), delim = "\t", quote = "", guess_max = 100000)
longitudinalSamples <- unique(pheno$X1)
rm(pheno)
gc()

# Read PGS data
prsGsa <- read_delim(prsGsaFile, delim = "\t")
prsCyto <- read_delim(prsCytoFile, delim = "\t")

# Configuration file that lists COVID-19 related file paths
configuration <- fromJSON(file = "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/barometer_statistics/coronabarometer-statistics/config.v3.json")
pseudoIdLinkageFile <- configuration$pseudoIdLinkageFile
pseudoIdLink <- fread(pseudoIdLinkageFile)

out <- "robert"

isBaselineRespondent <- tibble(PROJECT_PSEUDO_ID = baselineSamples) %>%
  inner_join(pseudoIdLink) %>% select(PSEUDOIDEXT) %>%
  mutate(isBaselineRespondent = 1L)

isLongitudinalRespondent <- tibble(PROJECT_PSEUDO_ID = longitudinalSamples) %>%
  inner_join(pseudoIdLink) %>% select(PSEUDOIDEXT) %>%
  mutate(isLongitudinalRespondent = 1L)

sampleGroups <- fread(invitedSamplesFilePath) %>%
  rename(isInvited = "COVID19_PROGRAMME") %>%
  left_join(isLongitudinalRespondent, by = c("PSEUDOIDEXT")) %>%
  left_join(isBaselineRespondent, by = c("PSEUDOIDEXT")) %>%
  mutate(PSEUDOIDEXT = as.character(PSEUDOIDEXT))

pgsFull <- bind_rows(list("GSA" = prsGsa, "HumanCytoSNP" = prsCyto), .id = "ARRAY") %>%
  inner_join(pseudoIdLink) %>%
  select(-PROJECT_PSEUDO_ID) %>%
  pivot_longer(c(-PSEUDOIDEXT, -ARRAY), names_to = "trait", values_to = "PGS") %>%
  mutate(PSEUDOIDEXT = as.character(PSEUDOIDEXT)) %>%
  left_join(sampleGroups, by = "PSEUDOIDEXT") %>%
  mutate(
    isBaselineRespondentVsInvited = factor(case_when(
      isInvited == 1 & isBaselineRespondent == 1 ~ "Baseline", 
      isInvited == 1 ~ "Invited")),
    isLongitudinalRespondentVsInvited = factor(case_when(
      isInvited == 1 & isLongitudinalRespondent == 1 ~ "Longitudinal", 
      isInvited == 1 ~ "Invited")),
    isLongitudinalRespondentVsOnlyBaseline = factor(case_when(
      isInvited == 1 & isLongitudinalRespondent == 1 ~ "Longitudinal",
      isInvited == 1 & isBaselineRespondent == 1 ~ "Baseline"
    ))) %>%
    filter(!is.na(isBaselineRespondentVsInvited))

tTestsBaselineVsInvitedGsa <- bind_rows(mapply(function(thisTrait, pgsTibble) {
  message(thisTrait)
  thisTraitData = pgsTibble %>% filter(trait == thisTrait)
  tTestRespondentVsInvited <- t.test(thisTraitData$PGS ~ thisTraitData$isBaselineRespondentVsInvited)
  print(tTestRespondentVsInvited$method)
  return(data.frame(pValueRespondentVsInvited = tTestRespondentVsInvited[["p.value"]],
                    meanIncluded = tTestRespondentVsInvited[["estimate"]][["mean in group Baseline"]],
                    meanInvitedNotIncluded = tTestRespondentVsInvited[["estimate"]][["mean in group Invited"]]))
}, selectedTraits, list(pgsTibble = pgsFull %>% filter(ARRAY == "GSA")), SIMPLIFY = F, USE.NAMES = T), .id = "trait")

tTestsBaselineVsInvitedCyto <- bind_rows(mapply(function(thisTrait, pgsTibble) {
  message(thisTrait)
  thisTraitData = pgsTibble %>% filter(trait == thisTrait)
  print(table(thisTraitData$isBaselineRespondentVsInvited))
  tTestRespondentVsInvited <- t.test(thisTraitData$PGS ~ thisTraitData$isBaselineRespondentVsInvited)
  print(tTestRespondentVsInvited$method)
  return(data.frame(pValueRespondentVsInvited = tTestRespondentVsInvited[["p.value"]],
                    meanIncluded = tTestRespondentVsInvited[["estimate"]][["mean in group Baseline"]],
                    meanInvitedNotIncluded = tTestRespondentVsInvited[["estimate"]][["mean in group Invited"]]))
}, selectedTraits, list(pgsTibble = pgsFull %>% filter(ARRAY == "HumanCytoSNP")), SIMPLIFY = F, USE.NAMES = T), .id = "trait")

tTestsLongitudinalVsInvitedGsa <- bind_rows(mapply(function(thisTrait, pgsTibble) {
  message(thisTrait)
  thisTraitData = pgsTibble %>% filter(trait == thisTrait) %>% filter(!is.na(isLongitudinalRespondentVsInvited))
  tTestRespondentVsInvited <- t.test(thisTraitData$PGS ~ thisTraitData$isLongitudinalRespondentVsInvited)
  print(tTestRespondentVsInvited$method)
  return(data.frame(pValueRespondentVsInvited = tTestRespondentVsInvited[["p.value"]],
                    meanIncluded = tTestRespondentVsInvited[["estimate"]][["mean in group Longitudinal"]],
                    meanInvitedNotIncluded = tTestRespondentVsInvited[["estimate"]][["mean in group Invited"]]))
}, selectedTraits, list(pgsTibble = pgsFull %>% filter(ARRAY == "GSA")), SIMPLIFY = F, USE.NAMES = T), .id = "trait")

tTestsLongitudinalVsInvitedCyto <- bind_rows(mapply(function(thisTrait, pgsTibble) {
  message(thisTrait)
  thisTraitData = pgsTibble %>% filter(trait == thisTrait) %>% filter(!is.na(isLongitudinalRespondentVsInvited))
  tTestRespondentVsInvited <- t.test(thisTraitData$PGS ~ thisTraitData$isLongitudinalRespondentVsInvited)
  print(tTestRespondentVsInvited$method)
  return(data.frame(pValueRespondentVsInvited = tTestRespondentVsInvited[["p.value"]],
                    meanIncluded = tTestRespondentVsInvited[["estimate"]][["mean in group Longitudinal"]],
                    meanInvitedNotIncluded = tTestRespondentVsInvited[["estimate"]][["mean in group Invited"]]))
}, selectedTraits, list(pgsTibble = pgsFull %>% filter(ARRAY == "HumanCytoSNP")), SIMPLIFY = F, USE.NAMES = T), .id = "trait")

tTestsLongitudinalVsBaselineGsa <- bind_rows(mapply(function(thisTrait, pgsTibble) {
  message(thisTrait)
  thisTraitData = pgsTibble %>% filter(trait == thisTrait) %>% filter(!is.na(isLongitudinalRespondentVsOnlyBaseline))
  tTestRespondentVsBaseline <- t.test(thisTraitData$PGS ~ thisTraitData$isLongitudinalRespondentVsOnlyBaseline)
  print(tTestRespondentVsBaseline$method)
  return(data.frame(pValueRespondentVsInvited = tTestRespondentVsBaseline[["p.value"]],
                    meanIncluded = tTestRespondentVsBaseline[["estimate"]][["mean in group Longitudinal"]],
                    meanInvitedNotIncluded = tTestRespondentVsBaseline[["estimate"]][["mean in group Baseline"]]))
}, selectedTraits, list(pgsTibble = pgsFull %>% filter(ARRAY == "GSA")), SIMPLIFY = F, USE.NAMES = T), .id = "trait")

tTestsLongitudinalVsBaselineCyto <- bind_rows(mapply(function(thisTrait, pgsTibble) {
  message(thisTrait)
  thisTraitData = pgsTibble %>% filter(trait == thisTrait) %>% filter(!is.na(isLongitudinalRespondentVsOnlyBaseline))
  tTestRespondentVsBaseline <- t.test(thisTraitData$PGS ~ thisTraitData$isLongitudinalRespondentVsOnlyBaseline)
  print(tTestRespondentVsBaseline$method)
  return(data.frame(pValueRespondentVsInvited = tTestRespondentVsBaseline[["p.value"]],
                    meanIncluded = tTestRespondentVsBaseline[["estimate"]][["mean in group Longitudinal"]],
                    meanInvitedNotIncluded = tTestRespondentVsBaseline[["estimate"]][["mean in group Baseline"]]))
}, selectedTraits, list(pgsTibble = pgsFull %>% filter(ARRAY == "HumanCytoSNP")), SIMPLIFY = F, USE.NAMES = T), .id = "trait")


tTestResultsBaseline <- bind_rows(list("GSA" = tTestsBaselineVsInvitedGsa, 
                               "HumanCytoSNP" = tTestsBaselineVsInvitedCyto),
                          .id = "array")

# Write summaries for PGSs
write.table(tTestResultsBaseline, 
            file.path(out, paste0("pgsCoronaSamplesTtestResultsBaseline_", format(Sys.Date(), "%Y%m%d"), ".txt")),
            row.names = F, col.names = T, quote = F, sep = "\t")

tTestResultsLongitudinal <- bind_rows(list("GSA" = tTestsLongitudinalVsInvitedGsa, 
                                       "HumanCytoSNP" = tTestsLongitudinalVsInvitedCyto),
                                  .id = "array")

# Write summaries for PGSs
write.table(tTestResultsLongitudinal, 
            file.path(out, paste0("pgsCoronaSamplesTtestResultsLongitudinal_", format(Sys.Date(), "%Y%m%d"), ".txt")),
            row.names = F, col.names = T, quote = F, sep = "\t")

tTestResultsLongitudinalVsBaseline <- bind_rows(list("GSA" = tTestsLongitudinalVsBaselineGsa, 
                                           "HumanCytoSNP" = tTestsLongitudinalVsBaselineCyto),
                                      .id = "array")

# Write summaries for PGSs
write.table(tTestResultsLongitudinalVsBaseline, 
            file.path(out, paste0("pgsCoronaSamplesTtestResultsLongitudinalVsBaseline_", format(Sys.Date(), "%Y%m%d"), ".txt")),
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

ggplot(pgsFull %>% filter(!is.na(isBaselineRespondentVsInvited) & ARRAY == "GSA" & trait %in% selectedTraits), 
       aes(x = isBaselineRespondentVsInvited, y = PGS)) +
  geom_violin(trim=T, size = 1 / (ggplot2::.pt * 72.27/96), fill = "gray", colour = "grey20") +
  geom_boxplot(width=0.1, size = 1 / (ggplot2::.pt * 72.27/96), outlier.shape = NA, colour = "grey10") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("Baseline", "Invited")),
                     bracket.nudge.y = 1, vjust = -0.5) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_x_discrete(breaks = c("Baseline", "Invited"),
                   labels = c("Baseline samples", "Invited, not \nincluded samples")) +
  xlab("Included baseline samples compared to samples that \nhave been invited but were not included") +
  ylab("Polygenic scores") +
  labs(title = "A: Baseline samples versus not included samples (GSA)") +
  facet_wrap(~traitLabel, ncol = 4, scales = "free_y",
             labeller = global_labeller) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggplot(pgsFull %>% filter(!is.na(isBaselineRespondentVsInvited) & ARRAY == "HumanCytoSNP" & trait %in% selectedTraits), 
       aes(x = isBaselineRespondentVsInvited, y = PGS)) +
  geom_violin(trim=T, size = 1 / (ggplot2::.pt * 72.27/96), fill = "gray", colour = "grey20") +
  geom_boxplot(width=0.1, size = 1 / (ggplot2::.pt * 72.27/96), outlier.shape = NA, colour = "grey10") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("Baseline", "Invited")),
                     bracket.nudge.y = 1, vjust = -0.5) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_x_discrete(breaks = c("Baseline", "Invited"),
                   labels = c("Baseline samples", "Invited, not \nincluded samples")) +
  xlab("Included baseline samples compared to samples that \nhave been invited but were not included") +
  ylab("Polygenic scores") +
  labs(title = "B: Baseline samples vs not included samples (HumanCytoSNP)") +
  facet_wrap(~traitLabel, ncol = 4, scales = "free_y",
             labeller = global_labeller) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggplot(pgsFull %>% filter(!is.na(isLongitudinalRespondentVsInvited) & ARRAY == "GSA" & trait %in% selectedTraits), 
       aes(x = isLongitudinalRespondentVsInvited, y = PGS)) +
  geom_violin(trim=T, size = 1 / (ggplot2::.pt * 72.27/96), fill = "gray", colour = "grey20") +
  geom_boxplot(width=0.1, size = 1 / (ggplot2::.pt * 72.27/96), outlier.shape = NA, colour = "grey10") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("Longitudinal", "Invited")),
                     bracket.nudge.y = 1, vjust = -0.5) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_x_discrete(breaks = c("Longitudinal", "Invited"),
                   labels = c("Longitudinal samples", "Invited, not \nincluded samples")) +
  xlab("Included longitudinal samples compared to samples that \nhave been invited but were not included in the longitudinal analysis") +
  ylab("Polygenic scores") +
  labs(title = "C: Longitudinal samples vs not included samples (GSA)") +
  facet_wrap(~traitLabel, ncol = 4, scales = "free_y",
             labeller = global_labeller) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggplot(pgsFull %>% filter(!is.na(isLongitudinalRespondentVsInvited) & ARRAY == "HumanCytoSNP" & trait %in% selectedTraits), 
       aes(x = isLongitudinalRespondentVsInvited, y = PGS)) +
  geom_violin(trim=T, size = 1 / (ggplot2::.pt * 72.27/96), fill = "gray", colour = "grey20") +
  geom_boxplot(width=0.1, size = 1 / (ggplot2::.pt * 72.27/96), outlier.shape = NA, colour = "grey10") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("Longitudinal", "Invited")),
                     bracket.nudge.y = 1, vjust = -0.5) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_x_discrete(breaks = c("Longitudinal", "Invited"),
                   labels = c("Longitudinal samples", "Invited, not \nincluded samples")) +
  xlab("Included longitudinal samples compared to samples that \nhave been invited but were not included in the longitudinal analysis") +
  ylab("Polygenic scores") +
  labs(title = "D: Longitudinal samples vs not included samples (HumanCytoSNP)") +
  facet_wrap(~traitLabel, ncol = 4, scales = "free_y",
             labeller = global_labeller) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggplot(pgsFull %>% filter(!is.na(isLongitudinalRespondentVsOnlyBaseline) & ARRAY == "GSA" & trait %in% selectedTraits), 
       aes(x = isLongitudinalRespondentVsOnlyBaseline, y = PGS)) +
  geom_violin(trim=T, size = 1 / (ggplot2::.pt * 72.27/96), fill = "gray", colour = "grey20") +
  geom_boxplot(width=0.1, size = 1 / (ggplot2::.pt * 72.27/96), outlier.shape = NA, colour = "grey10") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("Longitudinal", "BaselineOnly")),
                     bracket.nudge.y = 1, vjust = -0.5) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_x_discrete(breaks = c("Longitudinal", "BaselineOnly"),
                   labels = c("Longitudinal samples", "Samples only included \nin baseline analyses")) +
  xlab("Included longitudinal samples compared to samples that \nhave only been inlcuded in baseline analyses") +
  ylab("Polygenic scores") +
  labs(title = "E: Longitudinal samples vs baseline samples (GSA)") +
  facet_wrap(~traitLabel, ncol = 4, scales = "free_y",
             labeller = global_labeller) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggplot(pgsFull %>% filter(!is.na(isLongitudinalRespondentVsOnlyBaseline) & ARRAY == "HumanCytoSNP" & trait %in% selectedTraits), 
       aes(x = isLongitudinalRespondentVsOnlyBaseline, y = PGS)) +
  geom_violin(trim=T, size = 1 / (ggplot2::.pt * 72.27/96), fill = "gray", colour = "grey20") +
  geom_boxplot(width=0.1, size = 1 / (ggplot2::.pt * 72.27/96), outlier.shape = NA, colour = "grey10") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("Longitudinal", "BaselineOnly")),
                     bracket.nudge.y = 1, vjust = -0.5) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_x_discrete(breaks = c("Longitudinal", "BaselineOnly"),
                   labels = c("Longitudinal samples", "Samples only included \nin baseline analyses")) +
  xlab("Included longitudinal samples compared to samples that \nhave only been included in baseline analyses") +
  ylab("Polygenic scores") +
  labs(title = "F: Longitudinal samples vs baseline samples (HumanCytoSNP)") +
  facet_wrap(~traitLabel, ncol = 4, scales = "free_y",
             labeller = global_labeller) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

dev.off()
