# Calculating PGS correlations and correlations between set of questions on the baseline sample set.

# Load libraries
library(heatmap3)
library(readr)
library(rjson)

# Paths
workdir <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/"
phenoPath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/combined_questionnaires_v22_14-05-2021_genome_fitered_participants_filtered_correct_recoded/questionnaire_subset_participants_filtered_recoded_answers_14-05-2021"
prsLabelFile <- "prsLables.txt"
inclusionPrVlFile <- "inclusionPerVl.txt"
prsGsaFile <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/input_PGS_data_ugli_v4/PGS_combined_ugli_07-04-2021.txt"
prsCytoFile <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/input_PGS_data_cyto_v4_duplicate_filtered/PGS_combined_cyto_duplicate_from_ugli_removed_07-04-2021.txt"
qOverviewFile <- "quest_overview_nl_new_quest17_codes_updated_14-days_v2.txt"
selectedPrsFile <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/selectedTraits.txt"


# Read pheno
pheno <- read_delim(paste0(phenoPath, ".txt"), delim = "\t", quote = "", guess_max = 100000)
pheno2 <- as.data.frame(pheno)

row.names(pheno2) <- pheno2[,1]
colnames(pheno2)[1] <- "PROJECT_PSEUDO_ID"

startdate <- as.Date("30/03/2020","%d/%m/%Y")

# Read the prs labels
prsLabels <- as.matrix(read.delim(prsLabelFile, stringsAsFactors = F, row.names = 1))[,1]
dim(pheno2)

# List the questions for which we want to calculate correlations
qSelection <- c(
  "Quality of life" = "covt02_qualityoflife_adu_q_1",
  "Felt fine" = "covt01_fatigue_adu_q_1_c", 
  "Felt tired" = "covt01_fatigue_adu_q_1_a", 
  "Was easily tired" = "covt01_fatigue_adu_q_1_b", 
  "Felt physically exhausted" = "covt01_fatigue_adu_q_1_d")
qSelectionRev <- names(qSelection)
names(qSelectionRev) <- unname(qSelection)

correlationDf <- pheno2[, qSelection]

## Correlate PRS
library(heatmap3)
tiredCor <- cor(correlationDf, use = "pairwise.complete.obs", method = "spearman")#pheno3$array=="Gsa"
diag(tiredCor) <- 1
rownames(tiredCor) <- qSelectionRev[rownames(tiredCor)]
colnames(tiredCor) <- qSelectionRev[colnames(tiredCor)]
rownames(tiredCor)[rownames(tiredCor) == "Felt fine"] <- "Felt good"
colnames(tiredCor)[colnames(tiredCor) == "Felt fine"] <- "Felt good"
png("tiredCorrelation_20210607.png", width = 800, height = 800)
heatmap3(tiredCor, balanceColor = T, margins = c(21,21), scale = "none")
dev.off()

# Move on to the PRS correlations
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

# Only use the selected traits
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

# Merge the phenotype table and the PGS table
pheno3 <- merge(pheno2, prs, by = "PROJECT_PSEUDO_ID")

dim(pheno3)

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


qOverview2 <- qOverview[,-c(20,21)]
#head(qOverview2[,1])
dim(qOverview2)

# Generate a map of question names
qNameMap <- data.frame(orginal = row.names(qOverview2), new = make.names(row.names(qOverview2), unique = T), stringsAsFactors = F)
row.names(qNameMap) <- row.names(qOverview2)

# Do some processing of the question list to be able to reshape the dataframe
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

## Correlate PRS
library(heatmap3)
prsCor <- cor(prs[prs[,1] %in% pheno3[,"PROJECT_PSEUDO_ID"],-1])#pheno3$array=="Gsa"
diag(prsCor) <- 0
rownames(prsCor) <- prsLabels[rownames(prsCor)]
colnames(prsCor) <- prsLabels[colnames(prsCor)]
png("pgsCorrelation_20210601.png", width = 1000, height = 1000)
heatmap3(prsCor, balanceColor = T, margins = c(21,21), scale = "none")
dev.off()