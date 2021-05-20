# Run longitudinal models

library(nlme)
library(heatmap3)
library(readr)
library(lme4)
#library(future.apply)
library(survival)
library(parallel)

workdir <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/"
preparedDataFile <- "longitudinal.RData"

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


inverseVarianceMeta <- function(resultsPerArray, seCol, valueCol){
  x <- as.data.frame(resultsPerArray[[1]][,FALSE])
  x$sumYDivSe2 <- 0
  x$sum1DivSe2 <- 0
  
  for(array in names(resultsPerArray)){
    se2 <- resultsPerArray[[array]][,seCol] * resultsPerArray[[array]][,seCol]
    x$sumYDivSe2 <- x$sumYDivSe2 + (resultsPerArray[[array]][,valueCol]/ se2)
    x$sum1DivSe2 <- x$sum1DivSe2 + (1/se2)
  }
  
  metaRes <- as.data.frame(resultsPerArray[[1]][,FALSE])
  metaRes$y <- x$sumYDivSe2/x$sum1DivSe2
  metaRes$se <- sqrt(1/x$sum1DivSe2)
  metaRes$z <- metaRes$y/metaRes$se 
  metaRes$p <- 2*pnorm(-abs(metaRes$z))
  return(metaRes)
}

#qPrs <- qVsPrsRecode2[qVsPrsRecode2[,"prsTrait"] == "Life.satisfaction",][1,]
fitModel <- function(qPrs, selectedQ, vragenLong, arrayList){
  
  library(nlme)
  library(lme4)
  
  #zScores = tryCatch({
  print(paste(qPrs["question"], qPrs["prsTrait"], sep = " - "))
  
  q <- qPrs["question2"]
  
  if(!q %in% rownames(selectedQ)){
    print("skip")
    return(list("resPerArray" = NA, "metaRes" = NA, "qPrs" = qPrs, "error" = "Not in selected Q"))
  }
  
  qInfo <- selectedQ[q,]
  usedPrs <- qPrs["prsTrait"]
  # usedPrs <- usedPrs[!usedPrs %in% "Cigarettes.per.day"]
  #usedPrs <- c("BMI_gwas", "Life.satisfaction", "Neuroticism")
  #usedPrs <- "Life.satisfaction"
  #usedPrs <- "Neuroticism"
  #usedPrs <- "BMI_gwas"
  #usedPrs <- "Cigarettes.per.day"
  #usedPrs <- "COVID.19.susceptibility"
  #usedPrs <- "Anxiety.tension"
  #usedPrs <- "COVID.19.severity"
  #usedPrs <- "Worry.vulnerability"
  
  fixedString <- paste(q, "~((gender_recent+age_recent+age2_recent+household_recent+have_childs_at_home_recent+chronic_recent +", paste0(usedPrs, collapse = " + ") ,")*days + I(days^2) ) ")
  randomString <- "1+days+days2|PROJECT_PSEUDO_ID"#
  fixedModel <- as.formula(fixedString)
  randomModel <- as.formula(paste0("~",randomString))
  fullModel <- as.formula(paste0(fixedString, "+ (", randomString, ")"))
  
  resultsPerArray <- tryCatch(
    {
      lapply(arrayList, function(array){
        
        d <- vragenLong[!is.na(vragenLong[,q]) & vragenLong$array == array,c("PROJECT_PSEUDO_ID", q,usedPrs,"gender_recent","age_recent","age2_recent","household_recent","have_childs_at_home_recent","chronic_recent", "days", "days2", "days3", "vl")]
        table(d[,q])
        coef <- 0
        
        if(qInfo["Type"] == "gaussian" & qInfo["Mixed"]){
          res <-  lme(fixed = fixedModel, random=randomModel, data=d,na.action=na.omit, control = lmeControl(opt = "optim"))#
          return(res)
        } else if (qInfo["Type"] == "gaussian" & !qInfo["Mixed"]) {
          stop("Not implement")
        } else if (qInfo["Type"] == "binomial" & qInfo["Mixed"]) {
          if(max(d[,q])==2){
            d[,q] <- d[,q] -1
          }
          if(sum(range(d[,q])==0:1)!=2){
            stop("not binomal")
          }
          glmMerFit <- glmer(fullModel, data = d, family = binomial, nAGQ=0 )
          return(glmMerFit)
        } else if (qInfo["Type"] == "binomial" & !qInfo["Mixed"]) {
          d[,q] <- as.factor(d[,q])
          glmBinomFit <- glm(fixedModel ,family=binomial(link='logit'),data=d)
          return(glmBinomFit)
        }
        
        return(coef)
      })
    },
    error=function(cond){
      message(paste("ERROR:", qPrs["question"], qPrs["prsTrait"],cond))
      return(list("resPerArray" = NA, "metaRes" = NA, "qPrs" = qPrs, "error" = cond$message))
    }
  )#end try catch
  
  
  if(is.na(resultsPerArray[[1]])){
    #contains list with info
    return(resultsPerArray)
  }
  
  
  coefPerArray <- lapply(resultsPerArray, function(res){
    coef <- 0
    if(qInfo["Type"] == "gaussian" & qInfo["Mixed"]){
      print("test1")
      coef <- summary(res)$tTable
    } else if (qInfo["Type"] == "gaussian" & !qInfo["Mixed"]) {
      print("test2")
      stop("Not implement")
    } else if (qInfo["Type"] == "binomial" & qInfo["Mixed"]) {
      print("test3")
      coef <- summary(res)$coefficients
      colnames(coef)[1:2]<-c("Value", "Std.Error")
    } else if (qInfo["Type"] == "binomial" & !qInfo["Mixed"]) {
      print("test4b")
      coef <- summary(res)$coefficients
      colnames(coef)[1:2]<-c("Value", "Std.Error")
    }
    return(coef)
  })
  
  metaRes <- as.matrix(inverseVarianceMeta(coefPerArray, "Std.Error", "Value"))
  
  coefPerArray[[1]]
  coefPerArray[[2]]
  
  fullRes <- list("resPerArray" = resultsPerArray, "metaRes" = metaRes, "qPrs" = qPrs, "error" = NA)
  
  return(fullRes)
  
}



qVsPrs <- read.delim("gwasses_to_perform_filtered.txt", stringsAsFactors = F)
if(!all(qVsPrs$X %in% rownames(qNameMap))){
  stop("Not all Q found")
}
qVsPrs$saveName <- qNameMap[qVsPrs$X,2]
rownames(qVsPrs) <- qVsPrs$saveName

qVsPrsRecode <- apply(qVsPrs, 1,function(x){
  prsTraits <- strsplit(x[2], ";")[[1]]
  combi <- cbind(x[1],x[3],prsTraits)
  return(combi)
})

qVsPrsRecode[[2]][1,3]

qVsPrsRecode2 <- do.call(rbind, qVsPrsRecode)
colnames(qVsPrsRecode2) <- c("question", "question2", "prsTrait")

head(qVsPrsRecode2)

qLoop2 <- qLoop[names(qLoop) %in% qVsPrs$X]

q<-qNameMap["BMI",2]
q<-qNameMap["BMI",2]
q<-qNameMap["hoeveel zorgen maakte u zich de afgelopen 14 dagen over de corona-crisis?",2]
q<-qNameMap["hoe waardeert u uw kwaliteit van leven over de afgelopen 14 dagen? (include 7 days)",2]



clust <- makeCluster(8)
clusterExport(clust, "inverseVarianceMeta")
resultList <- parApply(clust, qVsPrsRecode2[qVsPrsRecode2[,"prsTrait"] == "Life.satisfaction",], 1, fitModel, selectedQ = selectedQ, vragenLong = vragenLong, arrayList = arrayList)
stopCluster(clust)



a <- lapply(resultList, function(x){
  if(!is.na(x[["metaRes"]]))
  {
    print(paste(x[["qPrs"]]["question"], x[["qPrs"]]["prsTrait"], sep = " - "))
    r <- nrow(x[["metaRes"]])
    print(x[["metaRes"]][r,"z"])
  }
})


str(resultList)

qVsPrsRecode2[c(12,13,16),"question2"] %in% rownames(selectedQ)


names(resultList) <- paste(qVsPrsRecode2[,1], qVsPrsRecode2[,3], sep = "-")

resultList[[3]]
str(resultList, max.level = 1)

dim(vragenLongTest)
dim(vragenLongValidation)

plot(ranef(res))
dev.off()

plot(res)
dev.off()

names(resultList)

resultList

zScoreList[[1]]

zScoreList <- lapply(resultList, function(x){return(tail(x[["metaRes"]][,"z", drop = F],n=1))})
pvalueList <- lapply(resultList, function(x){return(x[["metaRes"]][,"p", drop = F])})

# combine into z-score matrix excluding intercept
zscores <- do.call("cbind", zScoreList)
pvalues <- do.call("cbind", pvalueList)
#zscores <- zscores[,colnames(zscores) != "(Intercept)"]

str(zscores)

write.table(zscores, file = "zscoreMatrix.txt", sep = "\t", quote = F, col.names = NA)
write.table(pvalues, file = "pvaluesMatrix.txt", sep = "\t", quote = F, col.names = NA)


