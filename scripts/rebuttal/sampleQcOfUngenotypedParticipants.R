# Sample quality control

library(readr)


load_pheno_file <- function(phenoPath) {
  if(!file.exists(paste0(phenoPath, ".rds"))) {
    #Run once to conver pheno data to RDS format
    pheno <- read_delim(paste0(phenoPath, ".txt"), delim = "\t", quote = "", guess_max = 100000)
    dim(pheno)
    pheno2 <- as.data.frame(pheno)
    row.names(pheno2) <- pheno2[,1]
    
    colnames(pheno2)[1] <- "PROJECT_PSEUDO_ID"
    
    saveRDS(pheno2, paste0(phenoPath, ".rds"))
  } else {
    pheno2 <- readRDS(paste0(phenoPath, ".rds"))
  }
  
}

## Paths
workdir <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/analyses/20211125/sample_qc/samples_all"
phenoPath <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/combined_questionnaires_v10_06-05-2021_genome_fitered/questionaire_df_subset_participants_with_genome_data_06-05-2021"
phenoPathFull <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/combined_questionnaires_v10_06-05-2021/covid_export_questionaire_1-17_inclusive_correction_columns_06-05-2021"
qOverviewFile <- "quest_overview_nl_new_quest17_codes_updated_14-days-include-complete-qof.txt"
longitudinalSamplesFile <- "longitudinalSamples_NoGenotypeData.txt"
validationSamplesFile <- "validationSamples_NoGenotypeData.txt"
inclusionPrVlFile <- "inclusionPerVl_NoGenotypeData.txt"

setwd(workdir)

pheno2 <- load_pheno_file(phenoPathFull)
phenoGeno <- load_pheno_file(phenoPath)

pheno2 <- pheno2[!pheno2$PROJECT_PSEUDO_ID %in% phenoGeno$PROJECT_PSEUDO_ID,]

dim(pheno2)

qOverview <- as.matrix(read.delim(qOverviewFile, stringsAsFactors = F, row.names = 1))
vls <- colnames(qOverview)[-c(20,21)]

missingVraagList <- sapply(vls, function(vl){
  
  cat("----", vl, "\n")
  
  vlq <- qOverview[,vl]
  vlq <- vlq[vlq!=""]
  vlq <- vlq[vlq%in%colnames(pheno2)]
  
  cat(length(vlq), "\n")
  
  
  vlDate <- qOverview["responsdatum covid-vragenlijst",vl]
  
  #participants of this vl
  vlp <- !is.na(pheno2[,vlDate])
  
  print(table(vlp))
  
  #calculate percentage only among participants this vl
  vlMissing <- apply(pheno2[vlp,vlq], 2, function(x){sum(is.na(x)*100/length(x))})
  
  print(str(vlMissing))
  
  return(vlMissing)
  
})
pheno2[ "00067887-50eb-4d61-9f5b-820de4c18c26",1:10]

missingVraagList[["X14.0"]]

range(sapply(missingVraagList, length))

pdf("missingPerVl_NoGenotypeData.pdf")
sapply(vls, function(vl){
  hist(missingVraagList[[vl]], main = vl, breaks = 100)
})
dev.off()

qForSampleQc <- do.call("c", sapply(missingVraagList, function(x){names(x)[x<=5]}))

missing <- sapply(vls, function(vl){
  
  vlq <- qOverview[,vl]
  vlq <- vlq[vlq!=""]
  vlq <- vlq[vlq%in%qForSampleQc]
  
  vlMissing <- apply(pheno2[,vlq], 1, function(x){sum(is.na(x)*100/length(x))})
  
  return(vlMissing)
  
})

str(missing)

hist(missing, breaks = 100, ylim = c(0,1000))
dev.off()

ppvl <- apply(missing, 2, function(x){
  sum(x <= 5)
})
barplot(ppvl)
dev.off()

qpp <- apply(missing, 1, function(x){
  sum(x <= 5)
})
table(qpp)



colnames(missing)[1:12]

participatedFirstHalf <- apply(missing[,c(1:12)], 1, function(x){
  any(x <= 5)
})

participatedSecondHalf <- apply(missing[,c(13:19)], 1, function(x){
  any(x <= 5)
})

sum(participatedFirstHalf)
sum(participatedSecondHalf)


participatedBothHalf <- participatedFirstHalf & participatedSecondHalf
sum(participatedBothHalf)



qpp <- apply(missing[participatedBothHalf,], 1, function(x){
  sum(x <= 5)
})
table(qpp)


median(qpp)

inclusionPerVl <- missing <= 5
str(inclusionPerVl)

write.table(inclusionPerVl, file = inclusionPrVlFile, quote = F, sep = "\t", col.names = NA)


rpng(width = 1000, height = 1000)
barplot(table(qpp))
dev.off()




write.table(names(participatedBothHalf)[participatedBothHalf], file = longitudinalSamplesFile, quote = F, row.names = F, col.names = F);



validationSet <- apply(missing[,c("X4.0","X9.0", "X14.0", "X17.0")], 1, function(x){
  all(x <= 5)
})
table(validationSet)

write.table(names(validationSet)[validationSet], file = validationSamplesFile, quote = F, row.names = F, col.names = F);


pdf("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/rebuttal/pcPerWeek_NoGenotypeData.pdf")

sapply(vls, function(vl){


  
vlq <- qOverview[,vl]
vlq <- vlq[vlq!=""]
vlq <- vlq[vlq%in%qForSampleQc]


vlClean <- pheno2[missing[,vl] <= 5, vlq]

types <- sapply(vlClean, class)
vlClean <- as.matrix(vlClean[,types == "numeric" | types == "logical"])


vlCleanColMean <- apply(vlClean, 2, mean, na.rm = T)
for(c in 1:ncol(vlClean)){
  vlClean[is.na(vlClean[,c ]),c] <- vlCleanColMean[c]
}

vlClean <- vlClean[,apply(vlClean, 2, sd)>0]
dim(vlClean)
vlClean <- scale(vlClean)

pcaRes <- prcomp(t(vlClean), center = T, scale. = T)

#plot(pcaRes)
#dev.off()

png(paste0("/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/PRS_correlation/rebuttal/pcPerWeek_NoGenotypeData/", vl, ".png"))
plot(pcaRes$rotation[,c(1,2)], bg = adjustcolor("dodgerblue2", alpha.f = 0.1), pch = 21, col=adjustcolor("dodgerblue2", alpha.f = 0.1))
(sd1 <- sd(pcaRes$rotation[,1]))
abline(v=sd1*4, col = "firebrick")
abline(v=-sd1*4, col = "firebrick")
(sd2 <- sd(pcaRes$rotation[,2]))
abline(h=sd2*4, col = "firebrick")
abline(h=-sd2*4, col = "firebrick")
grDevices::dev.off()
})
