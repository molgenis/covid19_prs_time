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