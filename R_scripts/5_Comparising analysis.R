# Vergleichende Analysen der Proben
# Distanz-Analyse:
## Bei dieser Analyse wird dargestellt, wie sich die Proben ?hneln.
colors <- colorRampPalette(brewer.pal(4, "Blues"))(255) # Einstellung der Farbe mit Hilfe des packages "RColorBrewer"
poisd <- PoissonDistance(t(geneCounts_nonZeros)) # Hier werden nicht die geordneten Daten ben?tigt
poisd$dd

samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- paste(studyDesign$group)
colnames(samplePoisDistMatrix) <- paste(studyDesign$group)

pheatmap(samplePoisDistMatrix, 
         clustering_distance_rows = poisd$dd, 
         clustering_distance_cols = poisd$dd, 
         col = colors)

# Varianz-Analyse
## Hier wird die Varianz der Gene ?ber die Proben hinweg analysiert und grafisch dargestellt
vsd.deSeqRaw <- vst(object = deSeqRaw, blind = TRUE) # Datentransformation
topVarGenes <- head(order(rowVars(assay(vsd.deSeqRaw)), decreasing = TRUE), n = config$topVarianceGenes) # Filtern nach Genen mit der h?chsten Varianz 
mat <- assay(vsd.deSeqRaw)[topVarGenes,] # Erstellung einer Matrix mittels des Filters
mat <- mat - rowMeans(mat) # Schwankung um Mittelwert

## Erstellung der Heatmap
mat_breaks <- seq(min(mat), max(mat), length.out = 10)
anno <- as.data.frame(colData(vsd.deSeqRaw)[, c("group","replicate")]) # Hier k?nnen biliebige Spaltennamen eingetragen werden. M?ssen vorhanden sein!
variance_heatmap <- pheatmap(mat, 
                             annotation_col = anno, 
                             labels_row = ens.str.SYMBOL(mat),
                             labels_col = c(studyDesign$group, studyDesign$replicate), drop_levels = TRUE) # Die Reihenfolge kann ?ber die Faktorisierung modifiziert werden

## Falls Faktorisierung notwendeig
#colData(vsd.deSeqRaw)[,"replicate"] <- factor(colData(vsd.deSeqRaw)[,"replicate"])

# Principal component analysis
plotPCA(vsd.deSeqRaw, intgroup = c("group" , "replicate")) # Hiermit wird automatisch eine PCA Analyse durchgef?hrt und visualisiert

## Die Darstellung kann aber auch modifiziert werden. Hierzu wird das package ggplot2 verwendet
## Die Daten m?ssen dabei ausgegeben und in ein Objekt ?berf?hrt werden
pcaData <- plotPCA(vsd.deSeqRaw, intgroup = c("group" , "replicate"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

PCA_plot <- ggplot(pcaData, aes(PC1, PC2, color=group, shape=replicate)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

suppressWarnings((print(PCA_plot)))