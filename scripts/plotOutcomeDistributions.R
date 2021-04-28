#!/usr/bin/env Rscript

### Source data

workdir <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/"
preparedDataFile <- "longitudinal.RData"

setwd(workdir)
load(preparedDataFile)

### Constants

## The weeks since the first questionnaire + 1
weeksSinceStart <- c("1.0" = 1, 
                     "2.0" = 2, 
                     "3.0" = 3, 
                     "4.0" = 4, 
                     "5.0" = 5, 
                     "6.0" = 6, 
                     "7.0" = 8, 
                     "8.0" = 10, 
                     "9.0" = 12, 
                     "10.0" = 15, 
                     "11.0" = 16, 
                     "12.0" = 20, 
                     "13.0" = 24, 
                     "14.0" = 29, 
                     "15.0" = 32, 
                     "15.5" = 34, 
                     "16.0" = 36, 
                     "16.5" = 38, 
                     "17.0" = 42)

plotSubtitle <- "Frequency distributions across all included questionnaire timepoints."
plotContinuousSubtitle <- "Frequency distributions across all included questionnaire timepoints."

## Plotting defaults
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

# Plot answers per week for every question
pdf("answerDistributions_all_20210413.pdf", useDingbats = FALSE,
    width = 8.5, height = 11
)
par(xpd = NA)

for (qIndex in 1:nrow(selectedQ)) {
  q <- selectedQ[qIndex, "qId"]
  qLabel <- selectedQ[qIndex, "English.label"]
  
  d <- vragenLong[!is.na(vragenLong[,q]), c(q, "vl")]
  
  questionnaireNumber <- gsub("X", "", d$vl)
  questionnaireNumberNew <- sapply(questionnaireNumber, function(qN) which(names(weeksSinceStart) == qN))
  
  d$questionnaireLabel <- factor(
    questionnaireNumberNew, 
    levels = questionnaireNumberNew, 
    labels = paste0("Questionnaire ", questionnaireNumberNew," (week ", weeksSinceStart[questionnaireNumber] - 1, ")"))
  
  valueLabelsAsJson <- selectedQ[qIndex, "English.value.labels"]
  print(valueLabelsAsJson)
  if (!is.na(valueLabelsAsJson) 
      && !is.null(valueLabelsAsJson) 
      && valueLabelsAsJson != "") {
    valueLabels <- fromJSON(valueLabelsAsJson)
    
    d[,q] <- factor(as.numeric(d[,q]), 
                    levels = sort(as.numeric(names(valueLabels))), 
                    labels = valueLabels[sort(as.numeric(names(valueLabels)), index.return = T)$ix])
  } else if (q == "hoe.waardeert.u.uw.kwaliteit.van.leven.over.de.afgelopen.14.dagen." |
             q == "hoe.waardeert.u.uw.kwaliteit.van.leven.over.de.afgelopen.14.dagen...include.7.days.") {
    
    valueLabels <- c("1" = "1 (Terrible)", 
                     "2" = "2", "3" = "3", "4" = "4", "5" = "5", "6" = "6", "7" = "7", "8" = "8", "9" = "9",
                     "10" = "10 (Excellent)")
    
    d[,q] <- factor(as.numeric(d[,q]), 
                    levels = sort(as.numeric(names(valueLabels))), 
                    labels = valueLabels[sort(as.numeric(names(valueLabels)), index.return = T)$ix])
    
  } else if (q == "ik.ben.bang.dat.het.fout.zal.gaan.in.de.samenleving..in.de.afgelopen.14.dagen.") {
    
    valueLabels <- c("1" = "1 (Absolutely not)", 
                     "2" = "2", "3" = "3", "4" = "4", "5" = "5", "6" = "6", "7" = "7 (Very much)")
    
    d[,q] <- factor(as.numeric(d[,q]), 
                    levels = sort(as.numeric(names(valueLabels))), 
                    labels = valueLabels[sort(as.numeric(names(valueLabels)), index.return = T)$ix])
  }
  
  p <- NULL
  
  if (length(unique(d[,q])) <= 2) {
    
    p <- ggplot(d, aes(x=get(q), fill=get(q))) +
      stat_count(position = position_dodge2(width = 0.9, preserve = "single"), size = 1 / (ggplot2::.pt * 72.27/96), colour = "grey20") +
      scale_fill_manual(values = c("gray", "gray50")) +
      labs(title=qLabel, subtitle = plotSubtitle, y="Frequency", x = "Values") + 
      guides(fill=guide_legend(title=NULL)) +
      facet_wrap(vars(questionnaireLabel), nrow = 5, ncol = 4, as.table = T)
    
  } else if (length(unique(d[,q])) <= 10) {
    
    p <- ggplot(d, aes(x=get(q))) +
      stat_count(position = position_dodge2(width = 0.9, preserve = "single"), size = 1 / (ggplot2::.pt * 72.27/96), fill = "gray", colour = "grey20") +
      labs(title=qLabel, subtitle = plotContinuousSubtitle, y="Frequency", x = "Values") + 
      guides(fill=guide_legend(title=NULL)) +
      # scale_x_continuous(labels = labs, breaks = brks) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      facet_wrap(vars(questionnaireLabel), nrow = 5, ncol = 4, as.table = T)
    
  } else if (q == "hoeveel.glazen.alcohol.heeft.u.de.afgelopen.14.dagen.gedronken.") {
    
    p <- ggplot(d, aes(y=get(q), x=questionnaireLabel)) +
      geom_violin(trim=TRUE, bw = 3, size = 1 / (ggplot2::.pt * 72.27/96), fill = "gray", colour = "grey20") +
      geom_boxplot(width=0.1, outlier.shape = NA, colour = "grey10") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      labs(title=qLabel, subtitle = "Distributions across all included questionnaire timepoints", 
           y=qLabel, 
           x = "Questionnaire (week since first questionnaire)")
    
  } else {
    
    p <- ggplot(d, aes(y=get(q), x=questionnaireLabel)) +
      geom_violin(trim=FALSE, size = 1 / (ggplot2::.pt * 72.27/96), fill = "gray", colour = "grey20") +
      geom_boxplot(width=0.1, outlier.shape = NA, colour = "grey10") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      labs(title=qLabel, subtitle = "Distributions across all included questionnaire timepoints", 
           y=qLabel, 
           x = "Questionnaire (week since first questionnaire)")
    
  }
  
  print(p)
}

dev.off()
