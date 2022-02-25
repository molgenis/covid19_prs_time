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
#!/usr/bin/env Rscript

### Load libraries
library(rjson)
library(tidyverse)
library(grid)
library(gridExtra)
library(cowplot)
library(heatmap3)

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

figureTitle <- "Comparison of the number of cases in "
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

# Municipality data
populationData <- "nation-wide_correlations/data/Population_per_Province.csv"

# Google mobility data
mobilityData <- "nation-wide_correlations/data/2020_NL_Region_Mobility_Report.csv"

# Oxford government response data
oxfordGovernmentResponseData <- "../risky_behaviour/jobs/OxCGRT/OxCGRT_Netherlands.txt"


### Functions

## Moving average function (incorporating only prior info)
ma <- function(x, n = 7){stats::filter(x, rep(1 / n, n), sides = 1)}
ma <- function(x, n = 7){stats::filter(x, rep(1 / n, n), sides = 1)}

### Main
## Covid infection data

population <- read_delim(populationData, delim = "\t") %>% 
  select(c("Population" = "aantal", "Province" = "Regio's")) %>%
  filter(!is.na(Population))

covid <- read_delim(covidCasesData, delim = ";") %>% 
  group_by(Date_of_publication, Province) %>%
  summarize(Total_reported = sum(Total_reported),
            Hospital_admission = sum(Hospital_admission)) %>%
  full_join(population) %>%
  mutate(LifelinesRegion = case_when(Province %in% c("Drenthe", "Groningen", "Friesland") ~ "Yes",
                                     TRUE ~ "No")) %>%
  group_by(LifelinesRegion, Date_of_publication) %>%
  summarize(Total_reported = sum(Total_reported),
            Hospital_admission = sum(Hospital_admission),
            Population = sum(Population)) %>%
  mutate(ma_Adjusted = ma(Total_reported / (Population / 100000)),
         ma_Change = (Total_reported / (Population / 100000)) - ma_Adjusted)

casesPlot <- ggplot(covid, aes(x = Date_of_publication, y = ma_Adjusted)) +
  geom_line(aes(colour = LifelinesRegion)) +
  scale_colour_manual(breaks = c("Yes", "No"), values = c("darkturquoise", "grey70"), name = "Region", labels = c("Lifelines", "Other regions in \nthe Netherlands")) +
  coord_cartesian(xlim = dateLimits) + scale_x_date(date_labels = "%b %Y") +
  ggtitle("Comparison of COVID-19 case statistics") +
  ylab(str_wrap("Average number of confirmed cases per 100 000 inhabitants over the past 7 days", width = 40)) +
  xlab("Date") +
  theme(axis.title.y = element_text(angle = 90))

pdf("nation-wide_correlations/comparisonOfCaseStatistic.pdf", useDingbats = FALSE,
    width = 7, height = 5)
par(xpd = NA)
ggdraw(casesPlot)
dev.off()
