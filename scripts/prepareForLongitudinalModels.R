# Prepare data for longitudinal models



library(heatmap3)
library(readr)
library(rjson)



## Paths
workdir <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/"
phenoPath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/combined_questionnaires_v22_14-05-2021_genome_fitered_participants_filtered_correct_recoded/questionnaire_subset_participants_filtered_recoded_answers_longitudinal_filtered_15-05-2021"
prsLabelFile <- "prsLables.txt"
inclusionPrVlFile <- "inclusionPerVl.txt"
prsGsaFile <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/input_PGS_data_ugli_v4/PGS_combined_ugli_07-04-2021.txt"
prsCytoFile <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/input_PGS_data_cyto_v4_duplicate_filtered/PGS_combined_cyto_duplicate_from_ugli_removed_07-04-2021.txt"
qOverviewFile <- "quest_overview_nl_new_quest17_codes_updated_14-days_v2.txt"
selectedPrsFile <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/selectedTraits.txt"
validationSamplesFile <- "validationSamples.txt"
longitudinalSelectionRecodingFile <- "question_selection_20210525.tsv"
preparedDataFile <- "longitudinal.RData"

setwd(workdir)

if(FALSE){
  #Run once to conver pheno data to RDS format
  pheno <- read_delim(paste0(phenoPath, ".txt"), delim = "\t", quote = "", guess_max = 100000)
  dim(pheno)
  pheno2 <- as.data.frame(pheno)
  row.names(pheno2) <- pheno2[,1]
  
  colnames(pheno2)[1] <- "PROJECT_PSEUDO_ID"
  
  saveRDS(pheno2, paste0(phenoPath, ".rds"))
  
}

prsLabels <- as.matrix(read.delim(prsLabelFile, stringsAsFactors = F, row.names = 1))[,1]

#Constants
startdate <- as.Date("30/03/2020","%d/%m/%Y")
confounders <- c("gender_recent", "age_recent", "age2_recent", "chronic_recent", "household_recent", "have_childs_at_home_recent")


##load pheno and prs
pheno2 <- readRDS(paste0(phenoPath, ".rds"))

sampleQc <- read.delim(inclusionPrVlFile, stringsAsFactors = F, row.names = 1)

pheno2 <- pheno2[pheno2$PROJECT_PSEUDO_ID %in% row.names(sampleQc),]
dim(pheno2)

if(!all(confounders %in% colnames(pheno2))){
  stop("Not all confounders found")
}

prsGsa <- read.delim(prsGsaFile, stringsAsFactors = F)
prsCyto <- read.delim(prsCytoFile, stringsAsFactors = F)

if(!all(colnames(prsGsa) == colnames(prsCyto))){
  stop("Colnames must be equal")
}
if(!all(!row.names(prsGsa$PROJECT_PSEUDO_ID) %in% row.names(prsCyto$PROJECT_PSEUDO_ID))){
  stop("Overlapping samples")
}



#Scale to zero per array
prsGsa[,-1] <- scale(prsGsa[,-1])
prsCyto[,-1] <- scale(prsCyto[,-1])


prs <- rbind(prsGsa, prsCyto)


selectedTraits <- read.delim(selectedPrsFile, header = F, stringsAsFactors = F, comment.char = "#")[,1]
if(!all(selectedTraits %in% colnames(prs))){
  stop("Not all traits found")
}
prs <- prs[,c("PROJECT_PSEUDO_ID",selectedTraits)]
colnames(prs)[colnames(prs) == "BMI"] <- "BMI_gwas"

if(!all(pheno2$PROJECT_PSEUDO_ID %in% prs$PROJECT_PSEUDO_ID)){
  stop("Not all pheno have genetics")#easly solved but code makes this assumtion
}



pheno2$array <- factor(as.numeric(pheno2$PROJECT_PSEUDO_ID %in% prsGsa$PROJECT_PSEUDO_ID), levels = 0:1, labels = c("Cyto", "Gsa"))

arrayList <- as.list(levels(pheno2$array))
names(arrayList) <- levels(pheno2$array)



pheno3 <- merge(pheno2, prs, by = "PROJECT_PSEUDO_ID")


#center age
meanAge <- mean(pheno3[,"age_recent"])

pheno3[,"age_recent"] <- pheno3[,"age_recent"] - meanAge
pheno3[,"age2_recent"] <- pheno3[,"age_recent"] * pheno3[,"age_recent"]


totalPart <- nrow(pheno3)



#na col to use in reshapre for missing questions 
pheno3$naCOl <- NA


##Load and format questions meta data
qOverview <- as.matrix(read.delim(qOverviewFile, stringsAsFactors = F, row.names = 1))
vls <- colnames(qOverview)[-c(20,21)]

qOverview[2,]

qOverview[,-c(20,21)] <- apply(qOverview[,-c(20,21)], 1:2, function(x){
  if(x==""){
    return ("")
  } else if(!x %in% colnames(pheno3)) {
    return("")
  } else if(grepl("responsedate_adu", x)){
    #always include data question
    return(x)
  } else {
    return(x)#no masking
    # mask questions with >75% missing
    #if(sum(is.na(pheno3[,x])) >= (totalPart/4)){
    #  return (x)
    #} else {
    #return ("")
    #}
  }
})


qPassQc <- sapply(rownames(qOverview), function(q){sum(qOverview[q,1:19]!="")>=2})
table(qPassQc)

#qOverview2 <- qOverview[qPassQc,-c(20,21)]
qOverview2 <- qOverview[,-c(20,21)]
#head(qOverview2[,1])
dim(qOverview2)


qNameMap <- data.frame(orginal = row.names(qOverview2), new = make.names(row.names(qOverview2), unique = T), stringsAsFactors = F)
row.names(qNameMap) <- row.names(qOverview2)

qList <- lapply(qNameMap[,1], function(q){
  qs <- qOverview2[q,]
  qs[qs==""]="naCOl"
  qs[!qs %in% colnames(pheno3)] <- "naCOl"
  return(qs)
})
names(qList) <- qNameMap[,2]
#qList[[qNameMap["responsdatum covid-vragenlijst",2]]]
#qList[[qNameMap["als u moet kiezen, denkt u zelf dat u een coronavirus/covid-19 infectie hebt (gehad)?",2]]]



## remove quesed body weight
qList[[qNameMap["op welk moment hebt u dit lichaamsgewicht gemeten?",2]]]
qList[[qNameMap["wat is uw lichaamsgewicht (in kg)?",2]]]
qList[[qNameMap["BMI",2]]]

for(vl in colnames(qOverview)[-c(20,21)]){
  #1 = morning; 2 = afternoon/evening; 3 = guess
  weightMeasureType <- pheno3[,qList[[qNameMap["op welk moment hebt u dit lichaamsgewicht gemeten?",2]]][vl]]
  weightMask <- is.na(weightMeasureType) | weightMeasureType == 3
  
  
  pheno3[weightMask,qList[[qNameMap["wat is uw lichaamsgewicht (in kg)?",2]]][vl]] <- NA
  pheno3[weightMask,qList[[qNameMap["BMI",2]]][vl]] <- NA
  
  
}

## Reshape to long format and clean some variables

if(any(names(qList) %in% colnames(pheno3))){
  stop("Column name clash after reshape")
}




vragenLong <- reshape(pheno3, direction = "long", idvar = "PROJECT_PSEUDO_ID", varying = qList, v.names = names(qList), times = vls, timevar = "vl")
vragenLong$vl2 <- as.numeric(factor(vragenLong$vl, levels = vls, ordered = T))
vragenLong$vl3 <- factor(vragenLong$vl, levels = vls, ordered = F)
vragenLong$days <- as.numeric(difftime(vragenLong[,qNameMap["responsdatum covid-vragenlijst",2]], startdate ,units="days"))
vragenLong$days2 <- vragenLong$days*vragenLong$days
vragenLong$days3 <- vragenLong$days*vragenLong$days*vragenLong$days
vragenLong$days4 <- vragenLong$days*vragenLong$days*vragenLong$days*vragenLong$days
#vragenLong$days5 <- vragenLong$days*vragenLong$days*vragenLong$days*vragenLong$days*vragenLong$days


vragenLong$gender_recent <- factor(vragenLong$gender_recent, levels = 0:1, labels = c("female","male"))
vragenLong$household_recent <- factor(vragenLong$household_recent, levels = 0:1, labels = c("single-person household","multi-person household"))
vragenLong$have_childs_at_home_recent <- factor(vragenLong$have_childs_at_home_recent, levels = 0:1, labels = c("No childeren at home","Childeren at home"))
vragenLong$chronic_recent <- factor(vragenLong$chronic_recent, levels = 0:1, labels = c("Healthy","Chronic disease"))


hist(vragenLong$days, breaks = 330)
dev.off()

## Recode smoking
table(vragenLong[,"hebt.u.de.afgelopen.14.dagen.gerookt."], useNA = "always")
vragenLong[!is.na(vragenLong[,"hebt.u.de.afgelopen.14.dagen.gerookt."]) & vragenLong[,"hebt.u.de.afgelopen.14.dagen.gerookt."] == 2,"hebt.u.de.afgelopen.14.dagen.gerookt."] <- 0
table(vragenLong[,"hebt.u.de.afgelopen.14.dagen.gerookt."], useNA = "always")

## Read selected questions


selectedQ <- read.table(
  longitudinalSelectionRecodingFile, sep = "\t",
  stringsAsFactors = F, header = T, quote = "")


selectedQ <- selectedQ[selectedQ[,"Question"] %in% qNameMap[,1],]
selectedQ$qId <- qNameMap[selectedQ[,"Question"],2]
rownames(selectedQ) <- selectedQ[,"qId"]


## Add first / last day for selectedQ

str(!is.na(vragenLong[,selectedQ[,"qId"][1]]))
qRange <- t(sapply(selectedQ[,"qId"], function(x){
  range(vragenLong[!is.na(vragenLong[,x]),"days"],na.rm = T)
}))
str(qRange)
colnames(qRange) <- c("firstDay", "lastDay")
selectedQ <- cbind(selectedQ, qRange)

## Correlate PRS
library(heatmap3)
prsCor <- cor(prs[prs[,1] %in% pheno3[,"PROJECT_PSEUDO_ID"],-1])#pheno3$array=="Gsa"
diag(prsCor) <- 0
rownames(prsCor) <- prsLabels[rownames(prsCor)]
colnames(prsCor) <- prsLabels[colnames(prsCor)]
png("pgsCorrelation.png", width = 1000, height = 1000)
heatmap3(prsCor, balanceColor = T, margins = c(21,21), scale = "none")
dev.off()
colnames(prs)
str(cor.test(prs[prs[,1] %in% pheno3[,"PROJECT_PSEUDO_ID"],"Life.satisfaction"], prs[prs[,1] %in% pheno3[,"PROJECT_PSEUDO_ID"],"Neuroticism"]))
str(cor.test(prs[prs[,1] %in% pheno3[,"PROJECT_PSEUDO_ID"],"Life.satisfaction"], prs[prs[,1] %in% pheno3[,"PROJECT_PSEUDO_ID"],"Depression..broad."]))
# 
# longitudinalSelectionRecoding <- read.table(
#   longitudinalSelectionRecodingFile, sep = "\t",
#   stringsAsFactors = F, row.names = 1, header = T, quote = "")


## Convert ordinal to binary
for (qIndex in (1:nrow(selectedQ))) {
  q <- rownames(selectedQ)[qIndex]
  qName <- selectedQ[qIndex, "Question"]
  qInfo <- selectedQ[q,]
  if (!is.na(qInfo["Type"]) && qInfo["Type"] == "ordinal") {
    print(q)
    recodedQId <- paste0(q, "_binary")
    ordinalAnswers <- vragenLong[,q]
    recoded <- rep(NA_integer_, length(ordinalAnswers))
    
    valueLabelsAsJson <- selectedQ[qIndex, "Recode.value.labels"]
    
    print(valueLabelsAsJson)
    
    if (!is.na(valueLabelsAsJson) 
        && !is.null(valueLabelsAsJson) 
        && valueLabelsAsJson != "") {
      
      valueLabels <- unlist(fromJSON(valueLabelsAsJson))
      
      recoded <- valueLabels[ordinalAnswers]
    } else {
      recoded <- ordinalAnswers
    }
    
    print(table(recoded))
    print(table(ordinalAnswers))
    if (sum(table(recoded)) != sum(table(ordinalAnswers))) {
      stop("Sum of answer frequencies not equal")
    }
    vragenLong[,recodedQId] <- recoded
    selectedQ[q, "Type"] <- "binomial"
    selectedQ[q, "qId"] <- recodedQId
    qNameMap[selectedQ[qIndex,"Question"],2] <- recodedQId
    rownames(selectedQ)[qIndex] <- recodedQId
  }
}

qLoop <- as.list(selectedQ[,"qId"])
names(qLoop) <- selectedQ[,"Question"]

qLoop <- qLoop[!names(qLoop)=="hoeveel verschillende mensen, ouder dan 12 jaar, buiten uw eigen huishouden, hebt u in totaal tijdens de kerstvakantie bezocht en/of als bezoek ontvangen?"]
qLoop <- qLoop[!names(qLoop)=="ik ben bereid de coronaregels te overtreden om kerst en/of oud en nieuw te kunnen vieren zoals ik gewend ben"]
qLoop <- qLoop[!names(qLoop)=="ik vind het ongeacht de corona crisis fijn dat mensen meer onderlinge afstand houden."]

#table(vragenLong[,"everC19Pos"], useNA = "always")
table(vragenLong[,"Positive.tested.cumsum" ], useNA = "always")


validationSamples <- read.delim(validationSamplesFile, header = F)[,1]
validationRounds <- c("X4.0","X9.0", "X14.0", "X17.0")

vragenLongValidation <- vragenLong[(vragenLong$vl %in% validationRounds) & (vragenLong$PROJECT_PSEUDO_ID %in% validationSamples),]

vragenLongOther <- vragenLong[(vragenLong$vl %in% validationRounds) & !(vragenLong$PROJECT_PSEUDO_ID %in% validationSamples),]
vragenLongOther2 <- vragenLong[ !(vragenLong$PROJECT_PSEUDO_ID %in% validationSamples),]

vragenLongOther3 <- vragenLong[ (vragenLong$vl %in% validationRounds),]


table(vragenLongValidation$vl)

save(pheno3, vragenLong, qLoop, selectedQ, prs, qNameMap, validationSamples, validationRounds, vragenLongValidation, arrayList, file = preparedDataFile)

