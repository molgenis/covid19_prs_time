#Plot marginal effects of longitudinalModels

library(ggeffects)

workdir <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/"
phenoPath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/combined_questionnaires_v22_14-05-2021_genome_fitered_participants_filtered_correct_recoded/questionnaire_subset_participants_filtered_recoded_answers_longitudinal_filtered_15-05-2021"
prsLabelFile <- "prsLables.txt"
preparedDataFile <- "longitudinal.RData"
intermediatesdir <-  "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/longiIntermediates/"

setwd(workdir)
load(preparedDataFile)

#Functions:
spread_factor_columns <- function(df) {
  variables <- colnames(df)
  # For all columns, check if it is a factor,
  # if so, spread levels across columns as 0/1.
  # else, copy the column as is.
  out <- do.call("cbind", lapply(variables, function(variable) {
    col <- df[,variable]
    if (is.factor(col)) {
      # For all levels, add a numeric column (0/1) indicating the
      # factor value for the row.
      uniqueLevels <- levels(col)
      colDf <- as.data.frame(lapply(uniqueLevels, function(factorLevel) {
        levelBool <- as.numeric(col == factorLevel)
      }))
      colnames(colDf) <- paste0(variable, uniqueLevels)
      return(colDf)
    } else {
      colDf <- data.frame(col)
      colnames(colDf) <- variable
      return(colDf)
    }
  }))
  return(out)
}
# Function to predict outcome values given a new dataset, 
# named coefficient list (from a meta-analysis),
# the model on which analysis is based
predict_meta <- function(df, coefficients, family) {
  # Spread all factor levels across columns for all columns that are factors.
  dfSpread <- spread_factor_columns(df)
  # Calculate the predicted values per term
  predictionPerTerm <- do.call("cbind", sapply(
    names(coefficients), 
    function(coefficientLabel) {
      
      if (coefficientLabel == "I(days^2)") {
        dfSpread[coefficientLabel] <- dfSpread$days^2
      }
      
      # Coefficient
      coefficient <- coefficients[coefficientLabel]
      # Parse coefficient label,
      # Find out which columns should be multiplied
      variablesToMultiply <- strsplit(coefficientLabel, ":", fixed = T)[[1]]
      if (all(variablesToMultiply == "(Intercept)")) {
        return(data.frame(coefficientLabel = rep(coefficient, nrow(dfSpread))))
      } else if (length(variablesToMultiply) == 2) {
        return(Reduce(`*`, dfSpread[variablesToMultiply]) * coefficient)
      } else if (length(variablesToMultiply == 1)) {
        return(dfSpread[variablesToMultiply] * coefficient)
      } else {
        stop("Something went wrong, or we are going to multiply more than two terms. Check if this works first!")
      }
    }))
  # Calculate the predicted values on a regular linear scale
  predicted <- rowSums(predictionPerTerm)
  # Return values according to the linear inverse of the link function
  return(family$linkinv(predicted))
}

getScientificNUmberExpressionPrefix <- function(prefix, x, suffics, digits = 2){
  if(is.infinite(x)){
    stop("Error inf")
  }
  b <- x
  e <- 0
  if(b >= 10){
    while(b >= 10){
      b <- b /10
      e <- e + 1
    }
  } else if(b < 1 & b != 0) {
    while(b < 1){
      b <- b * 10
      e <- e - 1
    }
  }
  b <- round(b, digits)
  return(bquote(paste(.(prefix), .(b) %*% 10 ^{.(e)},.(suffics))))
}




prsLabels <- as.matrix(read.delim(prsLabelFile, stringsAsFactors = F, row.names = 1))[,1]


resultList <- sapply(list.files(intermediatesdir, pattern = "*.rds", full.names = T), readRDS)



resultsListNames <- sapply(resultList, function(x){
  return(paste0(x$qPrs["question2"],"_", x$qPrs["prsTrait"]))
})
names(resultList) <- resultsListNames



i <- 0 
summary <- lapply(resultList, function(x){
  i <<- i + 1
  print(i)
  
  q <- x[["qPrs"]]["question"]
  q2 <- x[["qPrs"]]["question2"]
  t <- x[["qPrs"]]["prsTrait"]
  if(!is.na(x["metaRes"])){
    r <- nrow(x[["metaRes"]])
    z <- x[["metaRes"]][r,"z"]
    p <- x[["metaRes"]][r,"p"]
  } else{
    z <- NA
    p <- NA
  }
  e <- x[["error"]]
  
  return(data.frame(q,q2,t,p,z,e))
  
  
  
  
})

summary <- do.call(rbind, summary)




#not skip 7 days and not NA result
summary$qEn <- selectedQ$label_en[match(summary$q, selectedQ$Question)]

summary$t2 <- prsLabels[match(summary$t, names(prsLabels))]

summary2 <- summary[!selectedQ$skip_7_days[match(summary$q, selectedQ$Question)] & !is.na(summary[,"z"]), ]

summary2$fdr <- p.adjust(summary2[,"p"], method = "BH")



write.table(summary2, file = "interactionSummary.txt", sep = "\t", quote = F, row.names = F)




sum(summary2$fdr <= 0.05)

summary3 <- summary2[summary2$fdr <= 0.05,]
dim(summary3)

summary3$p

colHigh = "#6300A7"
colMedium = "#D5546E"
colLow = "#FCD225"
colAxis = "grey70"
colMean = "#808080"

startdate <- as.Date("30/03/2020","%d/%m/%Y")
startday <- 0
endday <- startday + 307
axisAt <- c(startday,startday+100,startday+200, endday)

#fullRes <- resultList[[123]]
#
#fullRes <- readRDS(paste0(intermediatesdir,"/hoeveel.zorgen.maakte.u.zich.de.afgelopen.7.dagen.over.de.corona.crisis._Schizophrenia.rds"))
#fullRes <- readRDS(paste0(intermediatesdir,"/ik.voelde.me.moe...in.welke.mate.had.u.de.afgelopen.7.dagen.last.van._Depression..broad..rds"))
#fullRes <- readRDS(paste0(intermediatesdir,"/hoe.waardeert.u.uw.kwaliteit.van.leven.over.de.afgelopen.14.dagen._Life.satisfaction.rds"))
#fullRes <- readRDS(paste0(intermediatesdir,"/Positive.tested.cumsum_COVID.19.susceptibility.rds"))



fullRes <- readRDS(paste0(intermediatesdir,"/hoeveel.zorgen.maakte.u.zich.de.afgelopen.7.dagen.over.de.corona.crisis._Schizophrenia.rds"))
fullRes$metaRes
#lapply(resultList, function(fullRes){



######

i <- 1
i<- 11

for(i in 1:nrow(summary3)){
  
  print(i)
  
  
  
  pdf(paste0(summary3[i, "q2"],"_",summary3[i, "t"],"interactionPlots.pdf"), width = 10, height = 10)
  
  
  fullRes <- resultList[[paste0(summary3[i, "q2"],"_",summary3[i, "t"])]]
  
  
  fixedModel <- fullRes$fixedModel
  randomModel <- fullRes$randomModel
  fullModel <- fullRes$fullModel
  
  
  layout(matrix(c(1,1,3,5,4,6,2,2), ncol = 2, byrow = T), heights = c(0.2,1,1,0.2))
  
  qPrs <- fullRes$qPrs
  usedPrs <- qPrs["prsTrait"]
  prsLabel = prsLabels[usedPrs]
  
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
  
  
  par(mar = c(0,0,0,0), xpd = NA)
  plot.new()
  plot.window(xlim = 0:1, ylim = 0:1)
  text(0.5,0.75,paste0("Model fitted on '", qLable, "' stratified by '", prsLabel, "'"), cex = 1.5 , font = 1.5)
  
  if(interactionP < 0.01){
    text(0.5,0.35, getScientificNUmberExpressionPrefix("Interaction P-value: ", interactionP,paste0(" Z-score: ", format(interactionZ,nsmall = 2, digits = 2))))
  } else {
    text(0.5,0.35,paste0("Interaction P-value: ", format(interactionP,nsmall = 2, digits = 2), " Z-score: ", format(interactionZ,nsmall = 2, digits = 2)), cex = 1 , font = 1)
  }
  
  
  par(mar = c(0,0,0,0), xpd = NA)
  plot.new()
  plot.window(xlim = 0:1, ylim = 0:1)
  legend("center", lty = c(NA,NA,NA,1,1,1), pch = c(5,5,5, NA,NA,NA), lwd = 2, col = rep(c(colLow, colMedium, colHigh),2), legend = paste0(c("Mean for participants with lowest 10% PGS for ", "Mean for participants with average PGS for ", "Mean for participants with highest 10% PGS for ", "Fit for lowest 10% PGS for ", "Fit for median PGS for ", "Fit for highest 10% PGS for "), prsLabel), bty = "n", ncol = 2, text.col = colAxis)
  
  
  
  
  
  yRanges <- sapply(arrayList, function(array){
    
    res <- fullRes$resPerArray[[array]]
    
    d <<- vragenLong[!is.na(vragenLong[,q]) & vragenLong$array == array,c("PROJECT_PSEUDO_ID", q,usedPrs,"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2", "days3", "vl")]
    
    
    
    predictions <- ggeffect(res, terms = c("days[daysSeq]", paste0(usedPrs,"[prsRange2]")) , type = "fixed" ) 
    
    
    vlForThisQ <- unique(d$vl)
    plotYRange <- range(predictions$conf.low,predictions$conf.high, na.rm=T)
    sapply(vlForThisQ, function(vl){
      dVl <- d[d$vl==vl,]
      meanHigh <- mean(dVl[dVl[,usedPrs] >= prsRange2[1],q])
      meanLow <- mean(dVl[dVl[,usedPrs] <= prsRange2[3],q])
      meanRest <- mean(dVl[dVl[,usedPrs] > prsRange2[3] & dVl[,usedPrs] < prsRange2[1],q])
      plotYRange <<- range(plotYRange,meanHigh, meanLow, meanRest)
    })
    return(plotYRange)
    
  })
  
  yRange <- range(yRanges)
  

  
  
  
  for(array in arrayList){
    
    
    
    res <- fullRes$resPerArray[[array]]
    
    d <<- vragenLong[!is.na(vragenLong[,q]) & vragenLong$array == array,c("PROJECT_PSEUDO_ID", q,usedPrs,"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2", "days3", "vl")]
    
    
    
    predictions <- ggeffect(res, terms = c("days[daysSeq]", paste0(usedPrs,"[prsRange2]")) , type = "fixed" ) 
    #  plot(predictions)
    
    
    
  
    ## 
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
    # dev.off()
    
    
    
    
    
    
    
    
    
    
    ##
    
    
    
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
    
    
    par(mar = c(3,5,1,1), xpd = NA)
    plot.new()
    plot.window(xlim = c(startday,endday), ylim = range(lowPrs, medianPrs, highPrs))
    axis(side = 1, at = axisAt, labels = format.Date( axisAt+startdate, "%d-%b-%Y"), col = colAxis, col.axis = colMean)
    axis(side = 2, col = colAxis, col.axis = colMean)
    title(ylab = paste0("Contribution of ", prsLabel, " PGS\non ", qLable), col.lab = colMean)
    points(daysSeq, lowPrs, col = colLow, type = "l", lwd = 2)
    points(daysSeq, medianPrs, col = colMedium, type = "l", lwd = 2)
    points(daysSeq, highPrs, col = colHigh, type = "l", lwd = 2)
    
    
  }
  
  dev.off()
  
  
}

