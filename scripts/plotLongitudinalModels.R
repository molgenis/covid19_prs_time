fitModel <- function(qPrs, selectedQ, arrayList){
  
  # library(nlme)
  #  library(lme4)
  
  #zScores = tryCatch({
  print(paste(qPrs["question"], qPrs["prsTrait"], sep = " - "))
  
  q <- qPrs["question2"]
  
  intermediateFile <- make.names(paste0(qPrs["question"], "_", qPrs["prsTrait"]))
  intermediatePath <- paste0(intermediatesdir, "/" , intermediateFile, ".rds")
  
  fullRes <- NA
  
  if(file.exists(intermediatePath)){
    #Load exising results
    fullRes <- readRDS(intermediatePath)
  }
  
}
library(ggeffects)

fullRes <- resultList[[4]]

fullRes <- readRDS(paste0(intermediatesdir,"/hoeveel.zorgen.maakte.u.zich.de.afgelopen.7.dagen.over.de.corona.crisis._Schizophrenia.rds"))
fullRes <- readRDS(paste0(intermediatesdir,"/ik.voelde.me.moe...in.welke.mate.had.u.de.afgelopen.7.dagen.last.van._Depression..broad..rds"))
fullRes <- readRDS(paste0(intermediatesdir,"/hoe.waardeert.u.uw.kwaliteit.van.leven.over.de.afgelopen.14.dagen._Life.satisfaction.rds"))
fullRes <- readRDS(paste0(intermediatesdir,"/Positive.tested.cumsum_COVID.19.susceptibility.rds"))

res <- fullRes$resPerArray$Gsa

summary(res)

fixedModel <- fullRes$fixedModel
randomModel <- fullRes$randomModel
fullModel <- fullRes$fullModel


qPrs <- fullRes$qPrs

usedPrs <- qPrs["prsTrait"]

prsRange <- quantile(prs[,usedPrs],probs = seq(0,1,0.1))

array <- "Gsa"
q <- qPrs["question2"]
qInfo <- selectedQ[q,]

#qLable <- qInfo[,"English.label"]

daysSeq <- qInfo[,"firstDay"]:qInfo[,"lastDay"]


 
d <- vragenLong[!is.na(vragenLong[,q]) & vragenLong$array == array,c("PROJECT_PSEUDO_ID", q,usedPrs,"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2", "days3", "vl")]

y <- prsRange[c(10,6,2)]
  
(predictions <- ggeffect(res, terms = c("days[daysSeq]", paste0(usedPrs,"[y]")) , type = "fixed" ) )
rpng(width = 1000, height = 1000)
plot(predictions)
dev.off()


colHigh = "#6300A7"
colMedium = "#D5546E"
colLow = "#FCD225"
colAxis = "grey70"
colMean = "#808080"

startdate <- as.Date("30/03/2020","%d/%m/%Y")

startday <- 0
endday <- startday + 307

axisAt <- c(startday,startday+100,startday+200, endDate)

table(d$vl)

rpng(width = 800, height  = 800)
par(mar = c(3,5,1,0), xpd = NA)
plot.new()
plot.window(xlim = c(startday,endday), ylim = range(predictions$conf.low,predictions$conf.high))
#plot.window(xlim = c(startday,endday), ylim = c(1,10))
axis(side = 1, at = axisAt, labels = format.Date( axisAt+startdate, "%d-%b-%Y"), col = colAxis, col.axis = colMean)
axis(side = 2, col = colAxis, col.axis = colMean)
title(ylab = paste0("Predicted of ", qLable, "\n by full model"), col.lab = colMean)

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
  meanHigh <- mean(dVl[dVl[,usedPrs] >= y[1],q])
  meanLow <- mean(dVl[dVl[,usedPrs] <= y[3],q])
  meanRest <- mean(dVl[dVl[,usedPrs] > y[3] & dVl[,usedPrs] < y[1],q])
  
  sdHigh <- sd(dVl[dVl[,usedPrs] >= y[1],q])
  sdLow <- sd(dVl[dVl[,usedPrs] <= y[3],q])
  sdRest <- sd(dVl[dVl[,usedPrs] > y[3] & dVl[,usedPrs] < y[1],q])
  
  points(day, meanHigh, pch = 5, col = colHigh)
  points(day, meanRest, pch = 5, col = colMedium)
  points(day, meanLow, pch = 5, col = colLow)
  
}

# 
# boxplot(d[d$vl=="X2.0" & d[,usedPrs] <= y[3],q], at = 1, add=TRUE, width = 10, boxwex =5, col = colLow)
# boxplot(d[d$vl=="X2.0" & d[,usedPrs] >= y[1],q], at = 8, add=TRUE, width = 10, boxwex =5, col = colHigh)
# 
# boxplot(d[d$vl=="X7.0" & d[,usedPrs] <= y[3],q], at = 50, add=TRUE, width = 10, boxwex =5, col = colLow)
# boxplot(d[d$vl=="X7.0" & d[,usedPrs] >= y[1],q], at = 58, add=TRUE, width = 10, boxwex =5, col = colHigh)
# 
# boxplot(d[d$vl=="X11.0" & d[,usedPrs] <= y[3],q], at = 110, add=TRUE, width = 10, boxwex =5, col = colLow)
# boxplot(d[d$vl=="X11.0" & d[,usedPrs] >= y[1],q], at = 118, add=TRUE, width = 10, boxwex =5, col = colHigh)
# 
# boxplot(d[d$vl=="X11.0" & d[,usedPrs] <= y[3],q], at = 110, add=TRUE, width = 10, boxwex =5, col = colLow)
# boxplot(d[d$vl=="X11.0" & d[,usedPrs] >= y[1],q], at = 118, add=TRUE, width = 10, boxwex =5, col = colHigh)
# 
# boxplot(d[d$vl=="X11.0" & d[,usedPrs] <= y[3],q], at = 200, add=TRUE, width = 10, boxwex =5, col = colLow)
# boxplot(d[d$vl=="X11.0" & d[,usedPrs] >= y[1],q], at = 208, add=TRUE, width = 10, boxwex =5, col = colHigh)
# 
# 
# boxplot(d[d$vl=="X17.0" & d[,usedPrs] <= y[3],q], at = 299, add=TRUE, width = 10, boxwex =5, col = colLow)
# boxplot(d[d$vl=="X17.0" & d[,usedPrs] >= y[1],q], at = 307, add=TRUE, width = 10, boxwex =5, col = colHigh)

range(d[d$vl=="X7.0","days"])

points(predictionsLow$x, predictionsLow$predicted, col = colLow, type = "l", lwd = 2)
points(predictionsMedium$x, predictionsMedium$predicted, col = colMedium, type = "l", lwd = 2)
points(predictionsHigh$x, predictionsHigh$predicted, col = colHigh, type = "l", lwd = 2)

str(d)



test <- residuals(res)
str(test)

res2 <- res
str(res2, max.level = 1)


dev.off()

table(d$vl)

points(startdate +daysSeq, medianPrs, col = colMedium, type = "l", lwd = 2)
points(startdate +daysSeq, highPrs, col = colHigh, type = "l", lwd = 2)

c(startdate +daysSeq, rev(startdate +daysSeq))

ciExtend <- 0.03
polygon(c(startdate +daysSeq, rev(startdate +daysSeq)), y = c(lowPrs + ciExtend, rev(lowPrs) - ciExtend))
dev.off()


