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


prsLabels <- as.matrix(read.delim(prsLabelFile, stringsAsFactors = F, row.names = 1))[,1]

inclusionPrVl <- read.delim(inclusionPrVlFile)
str(inclusionPrVl)

heatmapLabels <- read.delim(heatmapLabelsFiles, stringsAsFactors = F)


prs2 <- prs[prs$PROJECT_PSEUDO_ID %in% inclusionPrVl[,1],]



hmData <- heatmap3(cor(prs2[,-1]), scale = "none", balanceColor = T, keep.dendro = T, distfun = function(x) as.dist(1 - cor(t(x), use = "pa")^2))
dev.off()


str(hmData)

layout(matrix(c(1,2), nrow =2))

qVsPrsBaselineZscores <- read.delim(qVsPrsBaselineZscoresFile, row.names = 1)
str(qVsPrsBaselineZscores)

qVsPrsBaselineZscores[is.na(qVsPrsBaselineZscores)] <- 0

prs3 <- prs2[,colnames(qVsPrsBaselineZscores) ]

#300 is questions tested; 23 is PRS tested
bonThres <- abs(qnorm((0.05 / (300 * 17))/2))
#qWithHit <- apply(qVsPrsBaselineZscores, 1, function(x){any(abs(x) >= bonThres)})

anySig <- apply(qVsPrsBaselineZscores, 1, function(x){
  return(any(abs(x)>= bonThres))
})
sum(anySig)

qVsPrsBaselineZscores <- qVsPrsBaselineZscores[anySig,]

row.names(qVsPrsBaselineZscores) <- heatmapLabels$label_en[match(row.names(qVsPrsBaselineZscores), heatmapLabels$Question)]

all(colnames(qVsPrsBaselineZscores) %in% names(prsLabels))

colnames(qVsPrsBaselineZscores) <- prsLabels[match(colnames(qVsPrsBaselineZscores), names(prsLabels))]

colnames(prs3) <- prsLabels[match(colnames(prs3), names(prsLabels))]


save(bonThres, prs3, qVsPrsBaselineZscores, file = "heatmapData.Rdata")




library(pheatmap)



str(qVsPrsBaselineZscores)






prsCor <- cor(prs3, use = "pa")

prsDist <- as.dist(1 - prsCor^2)
prsClust <- hclust(prsDist)




col <-colorRampPalette(c("firebrick3", "white", "navy"))(200)


pheatmap(prsCor, scale = "none", color = col, cluster_cols = prsClust, breaks = seq(-1,1,length.out = 201), cluster_rows = prsClust, show_rownames =T, show_colnames = F, cellwidth = 10, cellheight = 10, treeheight_row = 0, treeheight_col = 0, filename = "heatmapPrs.pdf")


colorRampPalette(c("navy", "white", "firebrick3"))(200)
colorRampPalette(c("navy", "white"))(50)
colorRampPalette(c("white", "firebrick3"))(50)

col <- c(colorRampPalette(c("firebrick3", "#F9E8E8"))(100), "white", colorRampPalette(c("#E4E4F2", "navy"))(100))



colBreaks <- c(min(qVsPrsBaselineZscores),seq(-10,-bonThres+0.01, length.out = 100), seq(bonThres-0.01, 10, length.out = 100),max(qVsPrsBaselineZscores))
length(colBreaks)
length(col)



qClust <- hclust(as.dist(1 - cor(t(qVsPrsBaselineZscores), use = "pa")^2), method = "ward.D2")

qVsPrsBaselineZscoresBool <- (abs(qVsPrsBaselineZscores) >= bonThres)+0
qClust <- hclust(dist(qVsPrsBaselineZscoresBool, method = "manhattan"), method = "ward.D2")

qVsPrsBaselineZscoresSig <- qVsPrsBaselineZscores
qVsPrsBaselineZscoresSig[abs(qVsPrsBaselineZscoresSig) <= bonThres] <- 0
qVsPrsBaselineZscoresSig[qVsPrsBaselineZscoresSig <= 0] <- qVsPrsBaselineZscoresSig[qVsPrsBaselineZscoresSig <= 0] * -1
qClust <- hclust(dist(qVsPrsBaselineZscoresSig, method = "euclidean"), method = "ward.D2")



pheatmap(qVsPrsBaselineZscores, scale = "none", color = col, cluster_cols = prsClust, breaks = colBreaks, cluster_rows = qClust, show_rownames =F, cellwidth = 5, cellheight = 5, treeheight_row = 40, treeheight_col = 0, legend_breaks = c(-10, -round(bonThres, digits = 2), round(bonThres, digits = 2), 10,20,30,40))
rpng.off()





pheatmap(qVsPrsBaselineZscores, scale = "none", color = col, cluster_cols = prsClust, breaks = colBreaks, cluster_rows = qClust, show_rownames =T, cellwidth = 10, cellheight = 10, filename = "heatmap.pdf", treeheight_row = 0, legend_breaks = c(-10, -round(bonThres, digits = 2), round(bonThres, digits = 2), 10,20,30,40))
