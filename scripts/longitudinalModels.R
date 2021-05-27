# Run longitudinal models

library(nlme)
library(heatmap3)
library(readr)
library(lme4)
library(future.apply)
library(survival)
library(parallel)

workdir <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/"
intermediatesdir <-  "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/longiIntermediates/"
preparedDataFile <- "longitudinal.RData"

setwd(workdir)
load(preparedDataFile)


#Functions:

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


fitModel <- function(qPrs, selectedQ, arrayList){
  
 # library(nlme)
#  library(lme4)
  
  #zScores = tryCatch({
  print(paste(qPrs["question"], qPrs["prsTrait"], sep = " - "))
  
  q <- qPrs["question2"]
  
  intermediateFile <- make.names(paste0(qPrs["question"], "_", qPrs["prsTrait"]))
  intermediatePath <- paste0(intermediatesdir, "/" , intermediateFile, ".rds")
  
  fullRes <- NA
  
  if(FALSE ){
    #file.exists(intermediatePath)
    #Load exising results
    fullRes <- readRDS(intermediatePath)
  } else {
    #calculate new results
  
  if(!q %in% rownames(selectedQ)){
    print("skip")
    return(list("resPerArray" = NA, "metaRes" = NA, "qPrs" = qPrs, "error" = "Not in selected Q"))
  }
  
  qInfo <- selectedQ[q,]
  
  if(!qPrs["Used_in_longitudinal_analysis"]){
    print("skip")
    return(list("resPerArray" = NA, "metaRes" = NA, "qPrs" = qPrs, "error" = "Not longitudinal"))
  }
    
  
  
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
        coef <- NA
        
        if(qInfo["TypeLongitudinal"] == "gaussian" & qInfo["Mixed"]){
          print("test1")
          #res <-  lme(fixed = fixedModel, random=randomModel, data=d,na.action=na.omit, control = lmeControl(opt = "optim"))#
          res <- 1
          return(res)
        } else if (qInfo["TypeLongitudinal"] == "gaussian" & !qInfo["Mixed"]) {
          print("test2")
          stop("Not implement")
        } else if (qInfo["TypeLongitudinal"] == "binomial" & qInfo["Mixed"]) {
          print("test3")
          if(max(d[,q])==2){
            d[,q] <- d[,q] -1
          }
          if(sum(range(d[,q])==0:1)!=2){
            stop("not binomal: " )
          }
          glmMerFit <- glmer(fullModel, data = d, family = binomial, nAGQ=0 )
          return(glmMerFit)
        } else if (qInfo["Type"] == "binomial" & !qInfo["Mixed"]) {
          print("test4")
          d[,q] <- as.factor(d[,q])
          glmBinomFit <- glm(fixedModel ,family=binomial(link='logit'),data=d)
          return(glmBinomFit)
        } else {
          error("no model")
        }
        
        return(coef)
      })
    },
    error=function(cond){
      print("ERROR")
      message(paste("ERROR:", qPrs["question"], qPrs["prsTrait"],cond))
      fullRes <- list("resPerArray" = NA, "metaRes" = NA, "qPrs" = qPrs, "error" = cond$message)
      return(fullRes)
    }
  )#end try catch
  
  
  if(is.na(resultsPerArray[1])){
    #contains list with info
    saveRDS(resultsPerArray, intermediatePath)
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
  
  fullRes <- list("resPerArray" = resultsPerArray, "metaRes" = metaRes, "qPrs" = qPrs, "error" = NA, "fixedModel" = fixedModel, "randomModel" = randomModel, "fullModel" = fullModel)
  saveRDS(fullRes, intermediatePath)
  }
  
  return(fullRes)
  
}

qPrs <- qVsPrsRecode2[10,]

test <- apply(qVsPrsRecode2[10,,drop =F], 1, fitModel, selectedQ = selectedQ, arrayList = arrayList)
test


qVsPrs <- read.delim("gwasses_to_perform_filtered_include_14_days.txt", stringsAsFactors = F)
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


clust <- makeForkCluster(nnodes = 8)
clusterExport("vragenLong")
clusterExport("intermediatesdir")
clusterEvalQ(clust, {
  library(lme4)
  library(nlme)
})
resultList <- parApply(clust, qVsPrsRecode2, 1, fitModel, selectedQ = selectedQ, arrayList = arrayList)
stopCluster(clust)

str(resultList[[1]], max.level = 1)

resultList <- apply(qVsPrsRecode2, 1, fitModel, selectedQ = selectedQ, arrayList = arrayList)


resultList <- parApply(clust, qVsPrsRecode2[qVsPrsRecode2[,"prsTrait"] == "Schizophrenia",], 1, fitModel, selectedQ = selectedQ, arrayList = arrayList)

i <- 0 
a <- lapply(resultList, function(x){
  
  i <<- i + 1
  
  if(!is.na(x["metaRes"]))
  {
    r <- nrow(x[["metaRes"]])
    z <- x[["metaRes"]][r,"z"]
    
    if(!is.null(z)){
      if(abs(z)>= 2.5){
      print(paste(i, x[["qPrs"]]["question"], x[["qPrs"]]["prsTrait"], sep = " - "))
        print(paste("Zscore", z))
      }
      
      }
    }
})



str(summary)

summary(resultList[[9]])

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


