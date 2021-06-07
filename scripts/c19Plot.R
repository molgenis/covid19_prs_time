workdir <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/"
phenoPath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/combined_questionnaires_v22_14-05-2021_genome_fitered_participants_filtered_correct_recoded/questionnaire_subset_participants_filtered_recoded_answers_longitudinal_filtered_15-05-2021"
prsLabelFile <- "prsLables.txt"
preparedDataFile <- "longitudinal.RData"
intermediatesdir <-  "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/longiIntermediates/"

setwd(workdir)
load(preparedDataFile)

library(survival)
#library("viridis")
#library(scales)

#show_col(viridis_pal(option= "C")(12))
#dev.off()

qName <- "Positive tested cumsum"
q<-qNameMap[qName,2]


vlDates <- read.delim("sendDates.txt", stringsAsFactors = T, row.names = 1)
vlDates[,1] <- as.Date(vlDates[,1], format = "%Y-%m-%d")
vlDates <- vlDates[order(vlDates[,1]),,drop=F]




colHigh = "#6300A7"
colMedium = "#D5546E"
colLow = "#FCD225"
colAxis = "grey70"
colMean = "#808080"


startdate <- as.Date("30/03/2020","%d/%m/%Y")
startday <- 0
endday <- startday + 307
axisAt <- c(startday,startday+100,startday+200, endday)


qInfo <- selectedQ[q,]

qLable <- qInfo[,"label_en"]



fullRes <- readRDS(paste0(intermediatesdir,"/Positive.tested.cumsum_COVID.19.susceptibility.rds"))
array = "Gsa"




fixedModel <- fullRes$fixedModel
randomModel <- fullRes$randomModel
fullModel <- fullRes$fullModel



qPrs <- fullRes$qPrs
usedPrs <- qPrs["prsTrait"]
prsLabel = prsLabels[usedPrs]






covid19Events <- lapply(pheno3$PROJECT_PSEUDO_ID, function(id){
  idVragen <- vragenLong[vragenLong[,"PROJECT_PSEUDO_ID"] == id,c(qNameMap["hebt u een coronavirus/covid-19 infectie (gehad)?",2],"vl", qNameMap["responsdatum covid-vragenlijst",2],usedPrs,"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "array", "days", "days2")]
  
  
  firstVlPos <- which(idVragen[,qNameMap["hebt u een coronavirus/covid-19 infectie (gehad)?",2]] == 1)[1]
  
  
  if(!is.na(firstVlPos)){
    firstPosDate <- idVragen[firstVlPos,qNameMap["responsdatum covid-vragenlijst",2]]
    
    eventInfo <- idVragen[firstVlPos,]
    eventInfo$eventDate <- firstPosDate
    eventInfo$eventDays <- as.numeric(difftime(firstPosDate, startdate ,units="days"))
    
    return(eventInfo)
  } else{
    #now test for censoring 
    
    censorDate <- 0
    eventInfo <- idVragen[max(which(!is.na(idVragen[,"responsdatum.covid.vragenlijst"]))),]
    lastCompletedVl = eventInfo[,"vl"]
    
    if(lastCompletedVl == "X17.0"){
      #completed
      censorDate <- endDate
    } else {
      #censoring
      censorDate <- vlDates[which(rownames(vlDates) == lastCompletedVl)+1,]
    }
    
    eventInfo$eventDate <- censorDate
    eventInfo$eventDays <- as.numeric(difftime(censorDate, startdate ,units="days"))
    eventInfo[,qNameMap["hebt u een coronavirus/covid-19 infectie (gehad)?",2]] <- 0
    
    
    return(eventInfo)
    
    
  }
  
  
  
})


covid19Events2 <- do.call("rbind", covid19Events)

prsRangeKm <- quantile(prs[,usedPrs],probs = c(0,0.1,0.4,0.6,0.9,1))
covid19Events2$prsQuantile <- cut(covid19Events2[,usedPrs],breaks = prsRangeKm, include.lowest = T)
length(table(covid19Events2$prsQuantile))


table(covid19Events2[,"eventDays"], useNA="always")
table(covid19Events2[,q], useNA="always")

survModel <- as.formula(paste(" Surv(eventDays, hebt.u.een.coronavirus.covid.19.infectie..gehad..) ~prsQuantile"))


km_fit <- survfit(survModel, data=covid19Events2)
summary(km_fit, times = c(1,180,max(covid19Events2$eventDays)))          



pdf("c19Figure.pdf", width = 3.5, height = 7, useDingbats = F)
#rpng()
layout(matrix(c(1,2,3,4), nrow=4, byrow = T), heights = c(1,1,1,1))

par(mar = c(3,5,1,1), xpd = NA, las = 1, cex.lab = 0.6, cex.axis = 0.6)


plot.new()
plot.window(xlim =c(0,307), ylim = c(0.8,1))
axis(side = 1, at = axisAt, labels = format.Date( axisAt+startdate, "%d-%b-%Y"), col = colAxis, col.axis = colMean)
axis(side = 2, col = colAxis, col.axis = colMean)
title(ylab = paste0("Mean ",qLable), col.lab = colMean)
lines(km_fit, xlab="Days", main = 'C19 incidence', col = c(colLow, NA,colMedium,NA, colHigh), lwd=2, mark.time=FALSE, ylim = c(0.8, 1)) 
#legend("bottomleft", fill = c(colLow, colMedium, colHigh), legend = paste0(c("Lowest 10% PGS of ", "Median PGS of ", "Highest 10% PGS of "), prsLabel), bty = "n")


##Panel B


interactionTerm <- paste0(usedPrs, ":days")

interactionP <- fullRes$metaRes[interactionTerm,"p"]
interactionZ <- fullRes$metaRes[interactionTerm,"z"]

q <- qPrs["question2"]
qInfo <- selectedQ[q,]

print(q)
print(usedPrs)

#daysSeq <- qInfo[,"firstDay"]:qInfo[,"lastDay"]
daysSeq <- seq(qInfo[,"firstDay"],qInfo[,"lastDay"],10)

qLable <- qInfo[,"label_en"]

prsRange <- quantile(vragenLong[!is.na(vragenLong[,q]),usedPrs],probs = seq(0,1,0.1))
prsRange2 <- prsRange[c(10,6,2)]


res <- fullRes$resPerArray[[array]]

d <<- vragenLong[!is.na(vragenLong[,q]) & vragenLong$array == array,c("PROJECT_PSEUDO_ID", q,usedPrs,"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2", "days3", "vl")]

length(unique(d[,"PROJECT_PSEUDO_ID",]))

predictions <- ggeffect(res, terms = c("days[daysSeq]", paste0(usedPrs,"[prsRange2]")) , type = "fixed" ) 


vlForThisQ <- unique(d$vl)
yRange <- range(predictions$conf.low,predictions$conf.high, na.rm=T)
sapply(vlForThisQ, function(vl){
  dVl <- d[d$vl==vl,]
  meanHigh <- mean(dVl[dVl[,usedPrs] >= prsRange2[1],q])
  meanLow <- mean(dVl[dVl[,usedPrs] <= prsRange2[3],q])
  meanRest <- mean(dVl[dVl[,usedPrs] > prsRange2[3] & dVl[,usedPrs] < prsRange2[1],q])
  yRange <<- range(yRange,meanHigh, meanLow, meanRest)
})




par(mar = c(3,5,1,1), xpd = NA)
plot.new()
plot.window(xlim = c(startday,endday), ylim = yRange)
#plot.window(xlim = c(startday,endday), ylim = c(1,10))
axis(side = 1, at = axisAt, labels = format.Date( axisAt+startdate, "%d-%b-%Y"), col = colAxis, col.axis = colMean)
axis(side = 2, col = colAxis, col.axis = colMean)

arrayLable <- array
if(array == "Gsa"){
  arrayLable <- "Global Screening Array"
} else if(array == "Cyto"){
  arrayLable <- "HumanCytoSNP-12"
}

title(ylab = paste0(qLable), col.lab = colMean, main = arrayLable)

predictionsLow <- predictions[predictions$group == levels(predictions$group)[1], ]
predictionsMedium <- predictions[predictions$group == levels(predictions$group)[2], ]
predictionsHigh <- predictions[predictions$group == levels(predictions$group)[3], ]

#medium first by design
polygon(c(predictionsLow$x, rev(predictionsLow$x)), y = c(predictionsMedium$conf.low, rev(predictionsMedium$conf.high)), col = adjustcolor(colMedium, alpha.f = 0.2), border = NA)
polygon(c(predictionsLow$x, rev(predictionsLow$x)), y = c(predictionsLow$conf.low, rev(predictionsLow$conf.high)), col = adjustcolor(colLow, alpha.f = 0.2), border = NA)
polygon(c(predictionsHigh$x, rev(predictionsHigh$x)), y = c(predictionsHigh$conf.low, rev(predictionsHigh$conf.high)), col = adjustcolor(colHigh, alpha.f = 0.2), border = NA)

vlForThisQ <- unique(d$vl)






for(vl in vlForThisQ){
  
  dVl <- d[d$vl==vl,]
  day <- median(dVl$days)
  meanHigh <- mean(dVl[dVl[,usedPrs] >= prsRange2[1],q])
  meanLow <- mean(dVl[dVl[,usedPrs] <= prsRange2[3],q])
  meanRest <- mean(dVl[dVl[,usedPrs] > prsRange2[3] & dVl[,usedPrs] < prsRange2[1],q])
  
  points(day, meanHigh, pch = 5, col = colHigh)
  points(day, meanRest, pch = 5, col = colMedium)
  points(day, meanLow, pch = 5, col = colLow)
  
}

points(predictionsLow$x, predictionsLow$predicted, col = colLow, type = "l", lwd = 2)
points(predictionsMedium$x, predictionsMedium$predicted, col = colMedium, type = "l", lwd = 2)
points(predictionsHigh$x, predictionsHigh$predicted, col = colHigh, type = "l", lwd = 2)


legend("topleft", lty = c(NA,NA,NA,1,1,1), pch = c(5,5,5, NA,NA,NA), lwd = 2, col = rep(c(colLow, colMedium, colHigh),2), legend = paste0(c("Mean for participants with lowest 10% PGS for ", "Mean for participants with average PGS for ", "Mean for participants with highest 10% PGS for ", "Fit for lowest 10% PGS for ", "Fit for median PGS for ", "Fit for highest 10% PGS for "), prsLabel), bty = "n")




##panel c

fam <- NA

if(qInfo["Type"] == "gaussian"){
  fam <- gaussian()
} else if (qInfo["Type"] == "binomial"){
  fam <- binomial(link = "logit")
} else {
  stop("error")
}

coef=coefficients(summary(res))[,1]

coef[!grepl(usedPrs, names(coef))] <- 0

dummy <- d[daysSeq,]
dummy$days <- daysSeq
dummy$days2 <- dummy$days * dummy$days
for(prsCol in colnames(prs)[-1]){
  dummy[,prsCol] <- mean(prs[,prsCol])
  #dummy[,prsCol] <- 0
}
dummy[,"age_recent"] <- mean(pheno3$age_recent)
#dummy[,"age_recent"] <- 0
dummy[,"age2_recent"] <- dummy[,"age_recent"] * dummy[,"age_recent"]
dummy[,"household_recent"] <- factor(levels(dummy[,"household_recent"])[1], levels = levels(dummy[,"household_recent"]))
dummy[,"have_childs_at_home_recent"] <- factor(levels(dummy[,"have_childs_at_home_recent"])[1], levels = levels(dummy[,"have_childs_at_home_recent"]))
dummy[,"gender_recent"] <- factor(levels(dummy[,"gender_recent"])[1], levels = levels(dummy[,"gender_recent"]))
dummy[,"chronic_recent"] <- factor(levels(dummy[,"chronic_recent"])[1], levels = levels(dummy[,"chronic_recent"]))

str(coef)

dummy[,usedPrs] <- prsRange2[1]
highPrs <- predict_meta(df = dummy, coefficients = coef, family = fam)#, family = binomial(link = "logit")
dummy[,usedPrs] <- prsRange2[2]
medianPrs <- predict_meta(df = dummy, coefficients = coef, family = fam)#, family = binomial(link = "logit")
dummy[,usedPrs] <- prsRange2[3]
lowPrs <- predict_meta(df = dummy, coefficients = coef, family = fam)#, family = binomial(link = "logit")



plot.new()
plot.window(xlim = c(startday,endday), ylim = range(lowPrs, medianPrs, highPrs))
axis(side = 1, at = axisAt, labels = format.Date( axisAt+startdate, "%d-%b-%Y"), col = colAxis, col.axis = colMean)
axis(side = 2, col = colAxis, col.axis = colMean)
title(ylab = paste0("Contribution of ", prsLabel, " PGS\non ", qLable), col.lab = colMean)
points(daysSeq, lowPrs, col = colLow, type = "l", lwd = 2)
points(daysSeq, medianPrs, col = colMedium, type = "l", lwd = 2)
points(daysSeq, highPrs, col = colHigh, type = "l", lwd = 2)


legend("bottomleft", lty = c(1,1,1), lwd = 2, col = rep(c(colLow, colMedium, colHigh),2), legend = paste0(c("Relative contribution of PGS for lowest 10% PGS for ", "Relative contribution of PGS for median PGS for ", "Relative contribution of PGS for highest 10% PGS for "), prsLabel), bty = "n")


#panel D


replication <- read.delim("replication_data_points2.txt",row.names =1 , stringsAsFactors = F, header = F)
replicationSe <- read.delim("replication_data_points2_se.txt",row.names =1 , stringsAsFactors = F, header = F)
str(replication)

betas <- replication["COVID-19 susceptibility / Ever positive SARS-CoV-2 PCR test",!is.na(replication["COVID-19 susceptibility / Ever positive SARS-CoV-2 PCR test",])]
se <- replicationSe["COVID-19 susceptibility / Ever positive SARS-CoV-2 PCR test",!is.na(replicationSe["COVID-19 susceptibility / Ever positive SARS-CoV-2 PCR test",])]
days <- replication[1,!is.na(replication["COVID-19 susceptibility / Ever positive SARS-CoV-2 PCR test",])]

lowCi <- betas - (se * 1.96)
highCi <- betas + (se * 1.96)

  plot.new()
  plot.window(xlim = c(startday,endday), ylim = range(betas, lowCi, highCi))
  axis(side = 1, at = axisAt, labels = format.Date( axisAt+startdate, "%d-%b-%Y"), col = colAxis, col.axis = colMean)
  axis(side = 2, col = colAxis, col.axis = colMean)
  title(ylab = paste0("Regression coeffcients for ", prsLabel, " PGS \n per questionnaire"), col.lab = colMean)
  
  for(a in names(days)){
    lines(c(days[a],days[a]), c(lowCi[a], highCi[a]), col = colAxis)
  }
  
  points(as.numeric(days), as.numeric(betas), col = colMean, type = "p", lwd = 2, pch = 16)
  
  
  
  
  
  days2 <-  as.numeric(days)
  abline(lm(as.numeric(betas) ~ days2), xpd = F, col = colMean)
  
  cor.test(as.numeric(betas), as.numeric(days))
  
  
  dev.off()
  

