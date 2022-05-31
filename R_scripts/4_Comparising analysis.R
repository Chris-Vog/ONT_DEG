# Vergleichende Analysen der Proben
# Distanz-Analyse:
## Bei dieser Analyse wird dargestellt, wie sich die Proben ähneln.
colors <- colorRampPalette(brewer.pal(4, "RdYlBu"))(255) # Einstellung der Farbe mit Hilfe des packages "RColorBrewer"
poisd <- PoissonDistance(t(geneCounts_nonZeros)) # Hier werden nicht die geordneten Daten benötigt
poisd$dd

samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- paste(studyDesign$group)
colnames(samplePoisDistMatrix) <- paste(studyDesign$group)

pheatmap(samplePoisDistMatrix, 
         clustering_distance_rows = poisd$dd, 
         clustering_distance_cols = poisd$dd, 
         col = colors)

# Varianz-Analyse
## Hier wird die Varianz der Gene über die Proben hinweg analysiert und grafisch dargestellt
vsd.deSeqRaw <- vst(object = deSeqRaw, blind = TRUE) # Datentransformation
topVarGenes <- head(order(rowVars(assay(vsd.deSeqRaw)), decreasing = TRUE), n = config$topVarianceGenes) # Filtern nach Genen mit der höchsten Varianz 
mat <- assay(vsd.deSeqRaw)[topVarGenes,] # Erstellung einer Matrix mittels des Filters
mat <- mat - rowMeans(mat) # Schwankung um Mittelwert

## Erstellung der Heatmap
mat_breaks <- seq(min(mat), max(mat), length.out = 10)
anno <- as.data.frame(colData(vsd.deSeqRaw)[, c("group","replicate")]) # Hier können biliebige Spaltennamen eingetragen werden. Müssen vorhanden sein!
variance_heatmap <- pheatmap(mat, 
                             annotation_col = anno, 
                             labels_row = ens.str.SYMBOL(mat),
                             labels_col = c(studyDesign$group, studyDesign$replicate), drop_levels = TRUE) # Die Reihenfolge kann über die Faktorisierung modifiziert werden

## Falls Faktorisierung notwendeig
# colData(vsd.deSeqRaw)[,"replicate"] <- factor(colData(vsd.deSeqRaw)[,"replicate"])

suppressWarnings(print(variance_heatmap))