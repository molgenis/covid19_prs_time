# Plot 



library(heatmap3)
library(readr)
library(NMF)

## Paths
workdir <- "/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/pgs_correlations/"
inclusionPrVlFile <- "inclusionPerVl.txt"
preparedDataFile <- "longitudinal.RData"
qVsPrsBaselineZscoresFile <- "covid_combined_matrix_z_scores_first_values_25-05-2021.txt"
heatmapLabelsFiles <- "question_selection_heatmap_20210525.txt"

setwd(workdir)
load(preparedDataFile)

selectedQ <- read.table(longitudinalSelectionRecodingFile, sep = "\t",  stringsAsFactors = F, header = T, quote = "")
prsLabels <- as.matrix(read.delim(prsLabelFile, stringsAsFactors = F, row.names = 1))[,1]

inclusionPrVl <- read.delim(inclusionPrVlFile)
str(inclusionPrVl)



prs2 <- prs[prs$PROJECT_PSEUDO_ID %in% inclusionPrVl[,1],]

qVsPrsBaselineZscores <- read.delim(qVsPrsBaselineZscoresFile, row.names = 1)
str(qVsPrsBaselineZscores)

qVsPrsBaselineZscores[is.na(qVsPrsBaselineZscores)] <- 0
  
prs3 <- prs2[,colnames(qVsPrsBaselineZscores) ]

#300 is questions tested; 17 is PRS tested
bonThres <- abs(qnorm((0.05 / (300 * 17))/2))


anySig <- apply(qVsPrsBaselineZscores, 1, function(x){
  return(any(abs(x)>= bonThres))
})
sum(anySig)

apply(qVsPrsBaselineZscores, 2, function(x){
  return(any(abs(x)>= bonThres))
})

qVsPrsBaselineZscores <- qVsPrsBaselineZscores[anySig,]
dim(qVsPrsBaselineZscores)
row.names(qVsPrsBaselineZscores)[!row.names(qVsPrsBaselineZscores) %in% selectedQ$Question]
dim(qVsPrsBaselineZscores)

row.names(qVsPrsBaselineZscores) <- selectedQ$label_en[match(row.names(qVsPrsBaselineZscores), selectedQ$Question)]

all(colnames(qVsPrsBaselineZscores) %in% names(prsLabels))

colnames(qVsPrsBaselineZscores) <- prsLabels[match(colnames(qVsPrsBaselineZscores), names(prsLabels))]

colnames(prs3) <- prsLabels[match(colnames(prs3), names(prsLabels))]

dim(qVsPrsBaselineZscores)
save(bonThres, prs3, qVsPrsBaselineZscores, file = "heatmapData.Rdata")


library(pheatmap)



prsCor <- cor(prs3, use = "pa")
prsDist <- as.dist(1 - prsCor^2)
prsClust <- hclust(prsDist)
colPrsCor <-colorRampPalette(c("firebrick3", "white", "navy"))(200)
pheatmap(prsCor, scale = "none", color = colPrsCor, cluster_cols = prsClust, breaks = seq(-1,1,length.out = 201), cluster_rows = prsClust, show_rownames =T, show_colnames = F, cellwidth = 10, cellheight = 10, treeheight_row = 0, treeheight_col = 0, filename = "heatmapPrs.pdf")


#colorRampPalette(c("navy", "white", "firebrick3"))(200)
#colorRampPalette(c("navy", "white"))(50)
#colorRampPalette(c("white", "firebrick3"))(50)

colZscores <- c(colorRampPalette(c("firebrick3", "#F9E8E8"))(100), "white", colorRampPalette(c("#E4E4F2", "navy"))(100))
colBreaks <- c(min(qVsPrsBaselineZscores),seq(-10,-bonThres+0.01, length.out = 100), seq(bonThres-0.01, 10, length.out = 100),max(qVsPrsBaselineZscores))
length(colBreaks)
length(colZscores)



qClust <- hclust(as.dist(1 - cor(t(qVsPrsBaselineZscores), use = "pa")^2), method = "ward.D2")

qVsPrsBaselineZscoresBool <- (abs(qVsPrsBaselineZscores) >= bonThres)+0
qClust <- hclust(dist(qVsPrsBaselineZscoresBool, method = "manhattan"), method = "ward.D2")

qVsPrsBaselineZscoresSig <- qVsPrsBaselineZscores
qVsPrsBaselineZscoresSig[abs(qVsPrsBaselineZscoresSig) <= bonThres] <- 0
qVsPrsBaselineZscoresSig[qVsPrsBaselineZscoresSig <= 0] <- qVsPrsBaselineZscoresSig[qVsPrsBaselineZscoresSig <= 0] * -1
qClust <- hclust(dist(qVsPrsBaselineZscoresSig, method = "euclidean"), method = "ward.D2")





pheatmap(qVsPrsBaselineZscores, scale = "none", color = colZscores, cluster_cols = prsClust, breaks = colBreaks, cluster_rows = qClust, show_rownames =T, cellwidth = 10, cellheight = 10, filename = "heatmap.pdf", treeheight_row = 0, legend_breaks = c(-10, -round(bonThres, digits = 2), round(bonThres, digits = 2), 10,20,30,40))






for(class in unique(selectedQ$Class)){
  
  qVsPrsBaselineZscoresClass <- qVsPrsBaselineZscores[unique(selectedQ$label_en[selectedQ$label_en %in% row.names(qVsPrsBaselineZscores)  & selectedQ$Class == class]),]
  

  qVsPrsBaselineZscoresSig <- qVsPrsBaselineZscoresClass
  qVsPrsBaselineZscoresSig[abs(qVsPrsBaselineZscoresSig) <= bonThres] <- 0
  qVsPrsBaselineZscoresSig[qVsPrsBaselineZscoresSig <= 0] <- qVsPrsBaselineZscoresSig[qVsPrsBaselineZscoresSig <= 0] * -1
  qClustClass <- hclust(dist(qVsPrsBaselineZscoresSig, method = "euclidean"), method = "ward.D2")
  
  pheatmap(qVsPrsBaselineZscoresClass, scale = "none", color = colZscores, cluster_cols = prsClust, breaks = colBreaks, cluster_rows = qClustClass, show_rownames =T, cellwidth = 10, cellheight = 10, filename = paste0("heatmap_", class,".pdf"), treeheight_row = 0, legend_breaks = c(-10, -round(bonThres, digits = 2), round(bonThres, digits = 2), 10,20,30,40), main = class)  

  
}


