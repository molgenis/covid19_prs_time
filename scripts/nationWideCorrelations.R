#!/usr/bin/env Rscript

### Load libraries
library(rjson)
library(tidyverse)
library(grid)
library(gridExtra)
library(cowplot)
library(heatmap3)

### Constants

# Quality of life data
workdir <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/"
preparedDataFile <- "longitudinal.RData"

setwd(workdir)
load(preparedDataFile)

# Plotting defaults
old <- theme_set(theme_classic(base_family="Helvetica"))
theme_update(line = element_line(
  colour = "grey70", size = 2 / (ggplot2::.pt * 72.27/96), 
  linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  axis.line = element_line(
    colour = "#595A5C", size = 2 / (ggplot2::.pt * 72.27/96), 
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  axis.ticks = element_line(
    colour = "#595A5C", size = 2 / (ggplot2::.pt * 72.27/96), 
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  plot.subtitle = element_text(color="gray50"),
  plot.title = element_text(face = "bold"),
  plot.margin = unit(c(0,0.96,0,0.24), "cm")
)

dates <- as.Date(c(
  "2020-03-30",
  "2020-04-06",
  "2020-04-13",
  "2020-04-20",
  "2020-04-27",
  "2020-05-04",
  "2020-05-18",
  "2020-06-01",
  "2020-06-15",
  "2020-07-06",
  "2020-07-13",
  "2020-08-10",
  "2020-09-07",
  "2020-10-12",
  "2020-11-02",
  "2020-11-16",
  "2020-11-30",
  "2020-12-14",
  "2021-01-11"))

testingAvailable <- as.Date("01-06-2020", format = "%d-%m-%Y")
secondWaveA <- as.Date("10-09-2020", format = "%d-%m-%Y")
secondWaveB <- as.Date("01-12-2020", format = "%d-%m-%Y")

figureTitle <- "Development of quality of life during the COVID-19 pandemic in the Netherlands"
annotTable <- tibble(startDate = dates)
dateLimits <- c(as.Date("2020-03-16"), max(annotTable$startDate))
questionnaireCol <- "firebrick2"
movingAvgCol <- "darkturquoise"
annotCol <- "grey80"

# KNMI weather data from Eelde
weatherData <- "nation-wide_correlations/data/etmgeg_280.txt"

# LCPS data on icu occupancy
covidIcuOccupancyData <- "nation-wide_correlations/data/covid-19.csv"

# RIVM data per day and per municipality
covidCasesData <- "nation-wide_correlations/data/COVID-19_aantallen_gemeente_per_dag.csv"

# Google mobility data
mobilityData <- "nation-wide_correlations/data/2020_NL_Region_Mobility_Report.csv"

# Oxford government response data
oxfordGovernmentResponseData <- "../risky_behaviour/jobs/OxCGRT/OxCGRT_Netherlands.txt"


### Functions

## Moving average function (incorporating only prior info)
ma <- function(x, n = 7){stats::filter(x, rep(1 / n, n), sides = 1)}

### Main
## Weather
weatherRaw <- read_delim(weatherData, trim_ws = T, delim = ",", skip = 46,
                         col_types = cols_only(YYYYMMDD = col_character(),
                                               SP = col_double(),
                                               SQ = col_double(),
                                               RH = col_double(),
                                               TG = col_double(),
                                               DR = col_double())) %>%
  mutate(DATE = as.Date(YYYYMMDD, format = "%Y%m%d")) %>%
  filter(DATE > as.Date("2020-01-01")) 

weatherClassification <- weatherRaw %>%
  mutate(SUNNY = SP >= 50,
         RAINY = RH > 0.2) %>%
  pivot_longer(c(SUNNY, RAINY)) %>%
  filter(value == T)

weatherProcessed <- weatherRaw %>%
  mutate(across(c(SP, SQ, RH, TG, DR), ma, .names = "ma_{.col}")) %>%
  pivot_longer(ends_with("SP") | ends_with("SQ") | ends_with("RH") | ends_with("TG") | ends_with("DR"),
               names_to = c("type", ".value"),
               names_pattern = "(ma_|)(.+)") %>%
  mutate(valueType = case_when(type == "ma_" ~ "movingAvg", TRUE ~ "raw"))

## Covid IC occupation
icCovid <- read_delim(covidIcuOccupancyData, delim = ",") %>%
  mutate(Date = as.Date(Datum, format = "%d-%m-%Y"))

## Covid infection data
covid <- read_delim(covidCasesData, delim = ";") %>% 
  group_by(Date_of_publication) %>%
  summarize(Total_reported = sum(Total_reported),
            Hospital_admission = sum(Hospital_admission)) %>%
  mutate(ma_Total_reported = 
           ma(Total_reported))

## Mobility data
mobility <- read_csv(mobilityData) %>%
  filter(place_id == "ChIJu-SH28MJxkcRnwq9_851obM") %>%
  mutate(across(ends_with("_percent_change_from_baseline"), ma, .names = "ma_{.col}")) %>%
  pivot_longer(ends_with("_percent_change_from_baseline"),
               names_to = c("type", ".value"),
               names_pattern = "(ma_|)(.+)") %>% 
  mutate(valueType = case_when(type == "ma_" ~ "movingAvg", TRUE ~ "raw"))

## Government response index
gri <- read_delim(oxfordGovernmentResponseData, delim = "\t")

## QoL
qualityOfLifeQuestion <- "hoe.waardeert.u.uw.kwaliteit.van.leven.over.de.afgelopen.14.dagen...include.7.days."
qualityOfLifeTable <- vragenLong %>% 
  select(vl2, !!!qualityOfLifeQuestion, responsdatum.covid.vragenlijst) %>%
  filter(!is.na(responsdatum.covid.vragenlijst) & !is.na(get(qualityOfLifeQuestion))) %>%
  group_by(vl2) %>%
  summarise(avgResponseDate = as.Date(as.character(mean(responsdatum.covid.vragenlijst))), avgQoL = mean(get(qualityOfLifeQuestion)))

testingAvailableLayer <- annotate("segment", x = as.Date(testingAvailable), xend = as.Date(testingAvailable),
                                  y = -Inf, yend = Inf, colour = annotCol, linetype = "dashed",
                                  size = 0.5 / (ggplot2::.pt * 72.27/96))
secondWaveALayer <- annotate("segment", x = as.Date(secondWaveA), xend = as.Date(secondWaveA),
                             y = -Inf, yend = Inf, colour = annotCol, linetype = "dashed", 
                             size = 0.5 / (ggplot2::.pt * 72.27/96))
secondWaveBLayer <- annotate("segment", x = as.Date(secondWaveB), xend = as.Date(secondWaveB),
                             y = -Inf, yend = Inf, colour = annotCol, linetype = "dashed", 
                             size = 0.5 / (ggplot2::.pt * 72.27/96))
summerLayer <- annotate("rect",
                        xmin = as.Date("2020-06-01"), xmax = as.Date("2020-08-31"), 
                        ymin = -Inf, ymax = Inf, fill = annotCol, alpha=.2)

tempPlot <- ggplot(weatherProcessed %>% filter(valueType == "movingAvg"), 
                   aes(x = DATE, y = TG/10)) +
  summerLayer + testingAvailableLayer + secondWaveALayer + secondWaveBLayer +
  geom_line(colour = movingAvgCol) +
  coord_cartesian(xlim = dateLimits) + scale_x_date(date_labels = "%b %Y") +
  ylab(str_wrap(
    "Average temperature in 24 hours (°C)", width = 21)) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_text(angle = 90, size = 9))

sunshinePlot <- ggplot(weatherProcessed %>% filter(valueType == "movingAvg"), 
                       aes(x = DATE, y = SQ/10)) +
  geom_line(colour = movingAvgCol) +
  summerLayer + testingAvailableLayer + secondWaveALayer + secondWaveBLayer +
  coord_cartesian(xlim = dateLimits) + scale_x_date(date_labels = "%b %Y") +
  ylab(str_wrap("Hours of sunshine", width = 21)) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_text(angle = 90, size = 9))

casesPlot <- ggplot(covid, aes(x = Date_of_publication, y = ma_Total_reported)) +
  geom_line(colour = movingAvgCol) +
  summerLayer + testingAvailableLayer + secondWaveALayer + secondWaveBLayer +
  coord_cartesian(xlim = dateLimits) + scale_x_date(date_labels = "%b %Y") +
  ylab(str_wrap("Nation-wide confirmed cases", width = 18)) +
  theme(legend.position="none", 
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.line.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_text(angle = 90, size = 9))

icPlot <- ggplot(icCovid, aes(x = Date, y = IC_Bedden_COVID)) +
  geom_line(colour = "grey70") +
  summerLayer + testingAvailableLayer + secondWaveALayer + secondWaveBLayer +
  coord_cartesian(xlim = dateLimits) + scale_x_date(date_labels = "%b %Y") +
  ylab(str_wrap("Nation-wide COVID-19 ICU occupancy", width = 18)) +
  theme(legend.position="none",
        axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), axis.text.x = element_blank(), 
        axis.title.y = element_text(angle = 90, size = 9))

griPlot <- ggplot(gri, 
                  aes(x = Date, y = StringencyIndex)) +
  geom_line(colour = "grey70") + 
  summerLayer + testingAvailableLayer + secondWaveALayer + secondWaveBLayer +
  coord_cartesian(xlim = dateLimits) + scale_x_date(date_labels = "%b %Y") +
  ylab(str_wrap("Stringency Index", width = 21)) +
  theme(legend.position="none",
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.line.x = element_blank(), axis.text.x = element_blank(),
        axis.title.y = element_text(angle = 90, size = 9))

mobilityWorkPlot <- ggplot(mobility, 
                           aes(x = date, 
                               y = workplaces_percent_change_from_baseline,
                               colour = valueType,
                               size = valueType)) +
  geom_line() +
  scale_colour_manual(
    breaks = c("raw", "movingAvg"),
    values = c("grey70", movingAvgCol),
    labels = c("Raw", "Moving average over the past 7 days")) +
  scale_size_manual(breaks = c("raw",
                               "movingAvg"),
                    values = c(0.5 / (ggplot2::.pt * 72.27/96),
                               1 / (ggplot2::.pt * 72.27/96))) +
  guides(size = "none") +
  summerLayer + testingAvailableLayer + secondWaveALayer + secondWaveBLayer +
  coord_cartesian(xlim = dateLimits) + scale_x_date(date_labels = "%b %Y") +
  ylab(str_wrap("Work-related mobility relative to baseline (%)", width = 21)) +
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), axis.text.x = element_blank(),
        axis.title.y = element_text(angle = 90, size = 9),
        legend.direction = "horizontal")

mobilityRecrPlot <- ggplot(mobility, 
                           aes(x = date, 
                               y = retail_and_recreation_percent_change_from_baseline,
                               colour = valueType,
                               size = valueType)) +
  geom_line() +
  scale_colour_manual(breaks = c("raw",
                                 "movingAvg"),
                      values = c("grey70",
                                 movingAvgCol),
                      labels = c("-", 
                                 str_wrap("Moving average over the past 7 days", width = 16)),
                      name = "Legend") +
  scale_size_manual(breaks = c("raw",
                               "movingAvg"),
                    values = c(0.5 / (ggplot2::.pt * 72.27/96),
                               1 / (ggplot2::.pt * 72.27/96))) +
  summerLayer + testingAvailableLayer + secondWaveALayer + secondWaveBLayer +
  coord_cartesian(xlim = dateLimits) + scale_x_date(date_labels = "%b %Y") +
  ylab(str_wrap("Retail and recreation-related mobility relative to baseline (%)", width = 21)) +
  theme(legend.position="none",
        axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), axis.text.x = element_blank(),
        axis.title.y = element_text(angle = 90, size = 9))

legend <- get_legend(mobilityWorkPlot)

mobilityWorkPlot <- mobilityWorkPlot + theme(legend.position="none")

qolPlot <- ggplot(qualityOfLifeTable, 
                  aes(x = avgResponseDate, 
                      y = avgQoL)) +
  geom_point(colour = "grey70") +
  geom_line(colour = "grey70") +
  summerLayer + testingAvailableLayer + secondWaveALayer + secondWaveBLayer +
  coord_cartesian(xlim = dateLimits) + 
  scale_x_date(date_labels = "%b %Y") +
  ylab(str_wrap("Mean perceived quality of life", width = 21)) +
  xlab("Date") +
  theme(legend.position="none", 
        axis.title.y = element_text(angle = 90, size = 9))

myTitle <- ggdraw() +
  draw_label(
    str_wrap(figureTitle, width = 48), fontface = "bold"
  )

prow <- plot_grid(myTitle,
                  plot_grid(ggplot() + theme_void(), casesPlot, icPlot, griPlot, 
                            mobilityWorkPlot, mobilityRecrPlot, 
                            tempPlot, sunshinePlot, qolPlot,
                            align = 'v', nrow = 10,
                            rel_heights = c(0.06, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.24, 0.32)),
                  legend, nrow = 3, rel_heights = c(0.1, 2.06, 0.04))

pdf("nation-wide_correlations/qualityOfLifeOverPandemic_20210414.pdf", useDingbats = FALSE,
    width = 5.4, height = 11)
par(xpd = NA)
ggdraw(prow)
dev.off()

## Correlations
pandemicFiguresCombined <- qualityOfLifeTable %>% 
  left_join(covid, by = c("avgResponseDate" = "Date_of_publication")) %>%
  left_join(gri, by = c("avgResponseDate" = "Date")) %>%
  left_join(mobility %>% filter(valueType == "movingAvg"), 
            by = c("avgResponseDate" = "date")) %>%
  left_join(icCovid, by = c("avgResponseDate" = "Date")) %>%
  left_join(weatherProcessed %>% filter(valueType == "movingAvg"), 
            by = c("avgResponseDate" = "DATE")) %>%
  select(avgQoL, ma_Total_reported, StringencyIndex,
         retail_and_recreation_percent_change_from_baseline,
         workplaces_percent_change_from_baseline,
         IC_Bedden_COVID, SQ, TG) %>%
  rename("Mean perceived quality of life" = "avgQoL",
         "Average nation-wide confirmed \ncases over the previous 7 days" = "ma_Total_reported",
         "Stringency Index In the Netherlands" = "StringencyIndex",
         "Retail and recreation-related \nmobility over previous 7 days" = "retail_and_recreation_percent_change_from_baseline",
         "Work-related mobility over \nprevious 7 days" = "workplaces_percent_change_from_baseline",
         "Nation-wide COVID-19 ICU occupancy" = "IC_Bedden_COVID",
         "Sunshine in hours over the \nprevious 7 days" = "SQ",
         "Average temperature over the \nprevious 7 days (°C)" = "TG")

sapply(colnames(pandemicFiguresCombined), function (name) {cor.test(pandemicFiguresCombined[,"Mean perceived quality of life", drop = T], pandemicFiguresCombined[,name, drop = T])})

pandemicCorrelation <- cor(pandemicFiguresCombined)
write.table(pandemicCorrelation, "nation-wide_correlations/correlationMatrix_20210414.txt",
            sep = "\t", quote = F, col.names = T, row.names = T)

pdf("nation-wide_correlations/correlationHeatmap_20210414.pdf", width = 8, height = 8, useDingbats = FALSE)
par(xpd = NA)

heatmap3(pandemicCorrelation, balanceColor = T, margins = c(21,21), scale = "none", 
         cexCol = 0.96, cexRow = 0.96)
dev.off()
