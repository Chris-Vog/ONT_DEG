#Setting the working directory----
## Hier stellt ihr fest, in welchem Ordner ihr arbeitet
## Wenn nicht weiter spezifiziert, dann ist dies der Ausgangspunkt aller Dateipfade
setwd("/home/ag-rossi/Schreibtisch/CV043 - RNASeq MinION HaCaT/")
dir <- file.path("/home/ag-rossi/Schreibtisch/CV043 - RNASeq MinION HaCaT/")
list.files(dir)

#Loading required packages----
suppressMessages(library(digest))
suppressMessages(library(ShortRead))
suppressMessages(library(ggplot2))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(Rsubread))
suppressMessages(library(DESeq2))
suppressMessages(library(pcaMethods))
suppressMessages(library(caTools))                 
suppressMessages(library(writexl))                 
suppressMessages(library(yaml))
suppressMessages(library(session))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(RColorBrewer))
suppressMessages(library(PoiClaClu))
suppressMessages(library(pheatmap))
suppressMessages(library(vsn))
suppressMessages(library(karyoploteR))
suppressMessages(library(ensembldb))
suppressMessages(library(viridis))
suppressMessages(library(debrowser))
suppressMessages(library(edgeR))
suppressMessages(library(biomaRt))
suppressMessages(library(DOSE))
suppressMessages(library(pathview))
suppressMessages(library(clusterProfiler))
suppressMessages(library(tidyft))

#Loading custom functions and configurations----
## Selbstgeschriebene FUnktionen können in einer Extra-Datei gespeichert werden
## Diese werden hier geladen und stehen ab sofort zur Verfügung.
source("Static/R/common.R")
resultDir <- file.path("Analysis", "Results")

dir.create(resultDir, showWarnings = FALSE, recursive = TRUE) # Erstellung eines Ordners, in den die Ergebnisse gespeichert werden

## Laden der Konfiguration eurer Analyse, auf die im Laufe der Analyse zugegriffen wird
config <- yaml::yaml.load_file("config.yaml")
persistenceData <- file.path(resultDir, "CV043_RNASeq__HaCaT_ONT.Rdata")

#Creating the study design from the config-file----
## Das study design enthält die Informationen über eure Proben
studyDesign <- data_frame()
for (i in 1:length(config$Samples)) {
  studyDesign <- rbind(studyDesign,
                       data.frame(samples=names(config$Samples[[i]][[1]]),
                                  filename=unlist(config$Samples[[i]][[1]]),
                                  group=names(config$Samples[[i]])))
}

# Erweiterung des Study designs um die Anzahl der Replikate
studyDesign$replicate <-sapply(1:nrow(studyDesign), function(x)sum(studyDesign$group[1:x]==studyDesign$group[x]))

# Überprüfung, dass es sich um individuelle Dateien handelt
studyDesign$md5 <- lapply(as.character(studyDesign$filename), md5sum)

#Entfernung überflüssiger Spalten
studyDesign$samples <-NULL

# Darstellung des Study designs
knitr::kable(studyDesign, caption = "Study design for samples evaluated within report", booktabs=TRUE, table.envir='table*', linesep="") %>%
  kable_styling(latex_options=c("hold_position", "scale_down"))

#Quality control
#Creating a statistic-table of fastq-files----
processQCFastq <- function(rowname) {
  row <- which(row.names(studyDesign)==rowname)
  file <- as.character(studyDesign[row, "filename"])
  fastq <- readFastq(file, widthIds=FALSE)
  c(
    reads = formatC(length(fastq), big.mark=","), #number of total reads
    mbs = formatC(round(sum(width(fastq)) / 1000 / 1000, digits=1), big.mark=","), 
    min = min(width(fastq)), #minimum length
    max = max(width(fastq)), #maximum length
    mean = round(mean(width(fastq)), digits=1), #mean length
    median = round(median(width(fastq)), digits=0), #median of read length
    qval = round(mean(alphabetScore(fastq) / width(fastq)), digits=1), #mean quality score
    gc = round(mean(letterFrequency(sread(fastq), "GC")  / width(fastq)) * 100, digits=1), #Percent GC-content
    n50 = ncalc(width(fastq), n=0.5), #length of reads, which scaffold 50 % of assembled bases
    l50 = lcalc(width(fastq), n=0.5), #number of reads, which scaffold 50 % of assembled bases
    n90 = ncalc(width(fastq), n=0.9), #length of reads, which scaffold 90 % of assembled bases
    l90 = lcalc(width(fastq), n=0.9)  #number of reads, which scaffold 90 % of assembled bases
  )
}

data <- lapply(row.names(studyDesign), processQCFastq) #applies this function to a vector of samples
qcData <- data.frame(data) #transformation of the data into a data frame
colnames(qcData) <- row.names(studyDesign)

# Darstellung der Daten
knitr::kable(qcData, caption="Summary statistics for the cDNA libraries imported", booktabs=TRUE, table.envir='table*', linesep="")  %>%
  kable_styling(latex_options=c("hold_position", font_size=9))

#Extraction and visualization of read length----
extractLengths <- function(rowname) {
  row <- which(row.names(studyDesign)==rowname)
  file <- as.character(studyDesign[row, "filename"])
  fastq <- readFastq(file)
  t(as.matrix(width(fastq)))
}
lengthData <- mapply(extractLengths, row.names(studyDesign))
lengthDataMatrix <- data.frame(t(plyr::rbind.fill.matrix(lengthData)))
colnames(lengthDataMatrix) <-  row.names(studyDesign)

lengthMatrixMelt <- reshape2::melt(lengthDataMatrix, na.rm=TRUE, measure.vars=row.names(studyDesign))
lengthMatrixMelt <- cbind(lengthMatrixMelt, group=studyDesign[match(lengthMatrixMelt$variable,  rownames(studyDesign)), "group"])

plot <- ggplot(lengthMatrixMelt, aes(x=variable, y=value, fill=group)) + 
  geom_violin() + 
  scale_y_continuous(limits=c(0, as.numeric(quantile(lengthMatrixMelt$value, probs=c(0.975))))) + 
  xlab("study sample") +  
  ylab("Distribution of Read Lengths (bp)") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_brewer(palette="Paired") + 
  labs(title="Violin plot showing distribution of read lengths across samples")

suppressWarnings(print(plot))

#Extraction and visualization of quality scores----
extractQualities <- function(rowname) {
  row <- which(row.names(studyDesign)==rowname)
  file <- as.character(studyDesign[row, "filename"])
  fastq <- readFastq(file)
  t(as.matrix(alphabetScore(fastq) / width(fastq)))
}
qualityData <- mapply(extractQualities, row.names(studyDesign))
qualityDataMatrix <- data.frame(t(plyr::rbind.fill.matrix(qualityData)))
colnames(qualityDataMatrix) <-  row.names(studyDesign)

qualityMatrixMelt <- reshape2::melt(qualityDataMatrix, na.rm=TRUE, measure.vars=row.names(studyDesign))
qualityMatrixMelt <- cbind(qualityMatrixMelt, group=studyDesign[match(qualityMatrixMelt$variable,  rownames(studyDesign)), "group"])

plotQ <- ggplot(qualityMatrixMelt, aes(x=variable, y=value, fill=group)) + geom_violin() + scale_y_continuous(limits=c(min(qualityMatrixMelt$value), max(qualityMatrixMelt$value))) + xlab("study sample") +  ylab("Distribution of Read Qualities (QV)") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Paired") + labs(title="Violin plot showing distribution of read qualities across samples")

suppressWarnings(print(plotQ))

#Analysis of flagstat results----
flagstatTargets <- file.path("Analysis", "flagstat", 
                             paste(tools::file_path_sans_ext(basename(as.character(studyDesign$filename)), compression=TRUE),".txt",sep=""))

loadFlagstat <- function(file) {
  x <- read.table(file, header=FALSE, sep=" ", fill=NA)[c(1:5),c(1,3)]
  x[,1]
}

flagstatRes <- data.frame(matrix(unlist(lapply(flagstatTargets, loadFlagstat)), ncol=length(flagstatTargets)), stringsAsFactors = FALSE)
colnames(flagstatRes) <- rownames(studyDesign)
rownames(flagstatRes) <- c("read mappings", "Secondary", "Supplementary", "Duplicates", "Mapped")

flagstatRes[nrow(flagstatRes)+1,] <- as.numeric(gsub(",","",t(qcData)[, "reads"]))
rownames(flagstatRes)[6] <- "nreads"

getVal <- function(word) {
  sum(as.numeric(unlist(strsplit(word, "/"))))
}

zreads <- unlist(lapply(flagstatRes["read mappings", ], getVal)) -
  unlist(lapply(flagstatRes["Secondary", ], getVal)) -
  unlist(lapply(flagstatRes["Supplementary", ], getVal)) - 
  unlist(lapply(flagstatRes["Duplicates", ], getVal)) 

flagstatRes[nrow(flagstatRes)+1,] <- round(as.numeric(flagstatRes["Mapped", ]) / as.numeric(flagstatRes["nreads", ]) * 100, digits = 2)

flagstatRes <- flagstatRes[c(6,1,2,3,4,7),]
rownames(flagstatRes)[6] <- "%mapping"

knitr::kable(flagstatRes, caption="Summary statistics from the minimap2 long read spliced mapping.", booktabs=TRUE, table.envir='table*', linesep="")  %>%
  kable_styling(latex_options=c("hold_position", font_size=11)) %>%
  add_footnote(c("information from samtools flagstat"))

readCountTargets <- file.path("Analysis", "Minimap", 
                              paste(tools::file_path_sans_ext(basename(as.character(studyDesign$filename)), compression=TRUE), ".bam", sep=""))

ExternalAnnotation = file.path("ReferenceData", basename(config$genome_annotation))

#Assign mapped sequencing reads to specified genomic features----
geneCounts <- featureCounts(files=readCountTargets,
                            annot.ext=ExternalAnnotation,
                            isGTFAnnotationFile=TRUE,
                            GTF.featureType="exon",
                            GTF.attrType="gene_id",
                            isLongRead=TRUE,
                            largestOverlap=TRUE,
                            useMetaFeatures=TRUE,
                            nthreads = 16)$counts

colnames(geneCounts) <- rownames(studyDesign)

#Extraction of the top 10 genes according to the total number of reads----
knitr::kable(geneCounts[order(rowSums(geneCounts), decreasing=TRUE)[1:10],], caption="Table showing the 10 annotated gene features with the highest number of mapped reads", booktabs=TRUE, table.envir='table*', linesep="") %>%
  kable_styling(latex_options=c("hold_position", font_size=11)) %>%
  add_footnote(c("This is raw count data and no normalisation or transformation has been performed"))

#Removal of all genes with zero counts----
geneCounts_nonZeros <- geneCounts[which(rowSums(geneCounts) > 0),]

geneCounts_nonZeros <- as.data.frame(geneCounts_nonZeros)

ens.geneCounts_nonZeros <- rownames(geneCounts_nonZeros)
Symbols <- mapIds(org.Hs.eg.db, keys = ens.geneCounts_nonZeros, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
geneCounts_nonZeros$Symbols <- Symbols
dim(geneCounts_nonZeros)

geneCounts_nonZeros <- drop_na(geneCounts_nonZeros)
dim(geneCounts_nonZeros)

geneCounts_nonZeros <- geneCounts_nonZeros[!duplicated(geneCounts_nonZeros$Symbols),]

rownames(geneCounts_nonZeros) <- geneCounts_nonZeros$Symbols
geneCounts_nonZeros$Symbols <- NULL
dim(geneCounts_nonZeros)

#Change the names of the columns----
colnames(geneCounts_nonZeros) <- paste0(studyDesign$group, "_", studyDesign$replicate) #It is important to choose different column names

#Ordering of the data frame
geneCounts_ordered <- geneCounts_nonZeros
geneCounts_ordered$RowSums <- rowSums(geneCounts_ordered)

# As rownames will not be exported, we create an additional column containing the gene_id
geneCounts_ordered$gene_id <- rownames(geneCounts_ordered)

geneCounts_ordered <- geneCounts_ordered[!duplicated(geneCounts_ordered$gene_id),]
# Order the samples with the arrange function of dplyr----
geneCounts_ordered <- geneCounts_ordered %>%
  dplyr::arrange(desc(RowSums))

#Save the data set----
xlsExpressedGenes <- file.path("ExpressedGenes.xlsx")
write_xlsx(x = geneCounts_ordered, path = xlsExpressedGenes)

#Creation of an DESeq-Data set to store input values
deSeqRaw <- DESeqDataSetFromMatrix(countData = geneCounts, colData = studyDesign, design = ~group)

#Analysing the variance----
vsd.deSeqRaw <- vst(object = deSeqRaw, blind = TRUE)
topVarGenes <- head(order(rowVars(assay(vsd.deSeqRaw)), decreasing = TRUE), n = 50)
mat <- assay(vsd.deSeqRaw)[topVarGenes,]
mat <- mat - rowMeans(mat)

#Creation of a Poisson distance matrix----
colors <- colorRampPalette(brewer.pal(4, "RdYlBu"))(255)
poisd <- PoissonDistance(t(geneCounts_nonZeros))
poisd$dd

samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- paste(studyDesign$group)
colnames(samplePoisDistMatrix) <- paste(studyDesign$group)

pheatmap(samplePoisDistMatrix, clustering_distance_rows = poisd$dd, clustering_distance_cols = poisd$dd, col = colors)

#Variance heatmap (unnormalized data)----
#Creation of a heatmap depicting the genes with the highes variance
colData(vsd.deSeqRaw)[,"replicate"] <- factor(colData(vsd.deSeqRaw)[,"replicate"])
mat_breaks <- seq(min(mat), max(mat), length.out = 10)
anno <- as.data.frame(colData(vsd.deSeqRaw)[, c("group","replicate")])
variance_heatmap <- pheatmap(mat, 
         annotation_col = anno, 
         labels_row = mapIds(org.Hs.eg.db, keys = substr(rownames(mat),1,15), 
                             column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"),
         labels_col = c(studyDesign$group, studyDesign$replicate), drop_levels = TRUE)

suppressWarnings(print(variance_heatmap))

pdf(file = "Analysis/Results/Variance_heatmap.pdf", width = 5, height = 9)
suppressWarnings(print(variance_heatmap))
dev.off()

#Analyzing differentially expressed genes----
#DE-analysis based on the negative binomial distribution
dds <- DESeq(deSeqRaw)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#BaP vs DMSO----
ens.str.SYMBOL <- function(result.table) {
  ens.str.vector <- sub("\\.*","",rownames(result.table))
  mapIds(org.Hs.eg.db, keys = ens.str.vector, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
}

ens.str.ENTREZ <- function(result.table) {
  ens.str.vector <- sub("\\.*","",rownames(result.table))
  mapIds(org.Hs.eg.db, keys = ens.str.vector, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
}

res_BaP_DMSO <- results(dds, contrast = c("group", "BaP", "DMSO"), alpha = 0.1)

DESeq2::plotMA(res_DMSO_BaP, ylim = c(-5,5),
               colNonSig = "gray32", colSig = "red3",
               log = "x")+
  title("DMSO vs BaP")

summary(res_DMSO_BaP)
sink(file = "Analysis/DE_Analysis/BaP_DMSO_LFC_statistics.txt")
print(summary(res_DMSO_BaP))
sink()

res_BaP_DMSO$SYMBOL <- ens.str.SYMBOL(result.table = res_BaP_DMSO)
res_BaP_DMSO$ENTREZ <- ens.str.ENTREZ(result.table = res_BaP_DMSO)

res_BaP_DMSO <- res_BaP_DMSO[order(res_BaP_DMSO$padj, decreasing = FALSE),]
write.table(x = res_BaP_DMSO, file = "Analysis/DE_Analysis/BaP_DMSO.txt", sep = "\t")

#PCB vs DMSO----
res_PCB_DMSO <- results(dds, contrast = c("group","PCB","DMSO"), alpha =  0.1)
DESeq2::plotMA(res_DMSO_PCB, ylim=c(-5,5),
               colNonSig = "gray32", colSig = "red3",
               log = "x")+
  title("DMSO vs PCB")

summary(res_PCB_DMSO)
sink(file = "Analysis/DE_Analysis/PCB_DMSO_LFC_statistics.txt")
print(summary(res_PCB_DMSO))
sink()

res_PCB_DMSO$SYMBOL <- ens.str.SYMBOL(result.table = res_PCB_DMSO)
res_PCB_DMSO$ENTREZ <- ens.str.ENTREZ(result.table = res_PCB_DMSO)

res_PCB_DMSO <- res_PCB_DMSO[order(res_PCB_DMSO$padj, decreasing = FALSE),]
write.table(x = res_PCB_DMSO, file = "Analysis/DE_Analysis/PCB_DMSO.txt")

#BaP_PD153035 vs DMSO----
res_BaP_PD153035_DMSO <- results(dds, contrast = c("group","BaP_PD153035", "DMSO"), alpha = 0.1)
DESeq2::plotMA(res_BaP_PD153035_DMSO, ylim = c(-5,5),
               colNonSig = "gray32", colSig = "red3",
               log = "x")+
  title("DMSO vs BaP + PD153035")

summary(res_BaP_PD153035_DMSO)
sink(file = "Analysis/DE_Analysis/BaP_PD153035_DMSO_LFC_statistics.txt")
print(summary(res_BaP_PD153035_DMSO))
sink()

res_BaP_PD153035_DMSO$SYMBOL <- ens.str.SYMBOL(result.table = res_BaP_PD153035_DMSO)
res_BaP_PD153035_DMSO$ENTREZ <- ens.str.ENTREZ(result.table = res_BaP_PD153035_DMSO)

res_BaP_PD153035_DMSO <- res_BaP_PD153035_DMSO[order(res_BaP_PD153035_DMSO$padj, decreasing = FALSE),]
write.table(x = res_BaP_PD153035_DMSO, file = "Analysis/DE_Analysis/BaP_PD153035_DMSO.txt", sep = "\t")

#BaP vs BaP PD153035----
res_BaP_BaP_PD153035 <- results(dds, contrast = c("group","BaP","BaP_PD153035"), alpha = 0.1)

DESeq2::plotMA(res_BaP_BaP_PD153035, ylim =c(-5,5),
               colNonSig = "gray32", colSig = "red3",
               log = "x")+
  title("BaP vs BaP + PD153035")

summary(res_BaP_BaP_PD153035)
sink(file ="Analysis/DE_Analysis/BaP_PD153035_BaP_LFC_statistics.txt")
print(summary(res_BaP_BaP_PD153035))
sink()

res_BaP_BaP_PD153035$SYMBOL <- ens.str.SYMBOL(result.table = res_BaP_BaP_PD153035)
res_BaP_BaP_PD153035$ENTREZ <- ens.str.ENTREZ(result.table = res_BaP_BaP_PD153035)

res_BaP_BaP_PD153035[order(res_BaP_BaP_PD153035$padj, decreasing = FALSE),]
write.table(x = res_BaP_BaP_PD153035, file = "Analysis/DE_Analysis/BaP_BaP_PD153035.txt", sep = "\t")

# Venn diagram data vs DMSO----
PCB.total <- as.data.frame(rownames(res_PCB_DMSO[abs(res_PCB_DMSO$log2FoldChange) > 1.5 ,]))
PCB.up <- as.data.frame(rownames(res_PCB_DMSO[res_PCB_DMSO$log2FoldChange > 1.5 ,]))
PCB.down <- as.data.frame(rownames(res_PCB_DMSO[res_PCB_DMSO$log2FoldChange < -1.5 ,]))
dim(PCB.total)

BaP.total <- as.data.frame(rownames(res_BaP_DMSO[abs(res_BaP_DMSO$log2FoldChange) > 1.5, ]))
BaP.up <- as.data.frame(rownames(res_BaP_DMSO[res_BaP_DMSO$log2FoldChange > 1.5, ]))
BaP.down <- as.data.frame(rownames(res_BaP_DMSO[res_BaP_DMSO$log2FoldChange < -1.5, ]))
dim(BaP.total)

BaP_PD153035.total <- as.data.frame(rownames(res_BaP_PD153035_DMSO[abs(res_BaP_PD153035_DMSO$log2FoldChange) > 1.5, ]))
BaP_PD153035.up <- as.data.frame(rownames(res_BaP_PD153035_DMSO[res_BaP_PD153035_DMSO$log2FoldChange > 1.5 ,]))
BaP_PD153035.down <- as.data.frame(rownames(res_BaP_PD153035_DMSO[res_BaP_PD153035_DMSO$log2FoldChange < -1.5 ,]))
dim(BaP_PD153035.total)

sheets <- list("PCB.total" = PCB.total, 
               "PCB.up" = PCB.up,
               "PCB.down" = PCB.down,
               "BaP.total" = BaP.total,
               "BaP.up" = BaP.up,
               "BaP.down" = BaP.down,
               "BaP_PD153035.total" = BaP_PD153035.total,
               "BaP_PD153035.up" = BaP_PD153035.up,
               "BaP_PD153035.down" = BaP_PD153035.down)

write_xlsx(sheets, "Analysis/DE.complete.xlsx")

#Venn diagram of differential expressed genes
library(VennDiagram)

draw.triple.venn(area1 = dim(PCB.total)[1], #Genes PCB
                 area2 = dim(BaP.total)[1], #Genes BaP
                 area3 = dim(BaP_PD153035.total)[1], #Genes BaP_PD153035
                 n12 = length(base::intersect(PCB.total[[1]], BaP.total[[1]])),
                 n13 = length(base::intersect(PCB.total[[1]], BaP_PD153035.total[[1]])),
                 n23 = length(base::intersect(BaP.total[[1]], BaP_PD153035.total[[1]])),
                 n123 = length(base::intersect(intersect(PCB.total[[1]], BaP.total[[1]]), BaP_PD153035.total[[1]])),
                 category = c("PCB126 (910)", "B[a]P (1434)", "B[a]P + PD153035 (2062)"),
                 euler.d = TRUE,
                 fill = c("#F59B25", "#95E354", "#54A7E3"),
                 alpha = 0.5,
                 fontfamily = rep("Helvetica",7),
                 cat.fontfamily = rep("Helvetica",3),
                 cex = 0.9,
                 cat.cex = 0.9)

#Plotting Sample distance
vsd <- vst(dds, blind = FALSE)
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$group)
colnames(sampleDistMatrix) <- paste(vsd$group)
colors <- colorRampPalette(brewer.pal(4, "RdYlBu"))(255)

pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists)

#Variance heatmap (normalized data)----
topVarGenesDDS <- head(order(rowVars(assay(vsd)), decreasing = TRUE), n = 25)
mat_dds <- assay(vsd)[topVarGenesDDS,]

mat_dds <- mat_dds - rowMeans(mat_dds)
kable(mat_dds)

mat.ens.str <- sub("\\.*","",rownames(mat_dds))
colData(vsd)[,"replicate"] <- factor(colData(vsd)[,"replicate"])
anno <- as.data.frame(colData(vsd))[, c("group","replicate")]
pheatmap(mat_dds, annotation_col = anno, labels_row = mapIds(org.Hs.eg.db, keys = mat.ens.str, column = "SYMBOL",
                                                             keytype = "ENSEMBL", multiVals = "first"))

# Identification of genes, which are regulated by the non-genomic signaling pathway----
#PCB vs BaP - log2FoldChange > |2|
res_PCB_DMSO.unregulated <- res_PCB_DMSO[abs(res_PCB_DMSO$log2FoldChange) < 1,]

res_BaP_PCB <- results(dds, contrast = c("group","BaP","PCB"), independentFiltering = FALSE)
res_BaP_PCB <- res_BaP_PCB[order(res_BaP_PCB$log2FoldChange, decreasing = TRUE), ]
res_BaP_PCB.unregulated <- res_BaP_PCB[rownames(res_PCB_DMSO.unregulated),]

res_BaP_PCB.LFC2 <- res_BaP_PCB.unregulated[abs(res_BaP_PCB.unregulated$log2FoldChange) > 2,]
dim(res_BaP_PCB.LFC2)

#BaP vs BaP_PD153035 - Filtered by res_PCB_BaP.LFC2 and downregulated in BaP_PD153035
res_BaP_BaP_PD153035.LFC2 <- res_BaP_BaP_PD153035[rownames(res_BaP_PCB.LFC2),]
res_BaP_BaP_PD153035.LFC2 <- res_BaP_BaP_PD153035.LFC2[abs(res_BaP_BaP_PD153035.LFC2$log2FoldChange) > 1,]
dim(res_BaP_BaP_PD153035.LFC2)

# Filter PCB vs BaP for this results
res_BaP_PCB.refiltered <- res_BaP_PCB.LFC2[rownames(res_BaP_BaP_PD153035.LFC2),]
dim(res_BaP_PCB.refiltered)

#Depicting these genes in a variance matrix
mat_dds_filt <- assay(vsd)[rownames(res_BaP_PCB.refiltered),]
mat_dds_filt <- mat_dds_filt - rowMeans(mat_dds_filt)
mat_dds_filt <- as.data.frame(mat_dds_filt)
mat_dds_filt$SYMBOL <- ens.str.SYMBOL(result.table = mat_dds_filt)
NA.filter <- is.na(mat_dds_filt$SYMBOL)
mat_dds_filt <- mat_dds_filt[!NA.filter,]
mat_dds_filt$SYMBOL <- NULL

mat_dds_filt_topVar <- head(order(rowVars(as.matrix(mat_dds_filt)), decreasing = TRUE), n = 10)
mat_dds_filt <- mat_dds_filt[mat_dds_filt_topVar,]
breaksList <- seq(-1,1, by = 0.05)
mat.ens.str.filt <- sub("\\.*","",rownames(mat_dds_filt))
pheatmap(mat_dds_filt, annotation_col = anno, labels_row = mapIds(org.Hs.eg.db, keys = mat.ens.str.filt,
                                                                  column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList,
         labels_col = studyDesign$group)

# Principal Component Analysis
plotPCA(vsd, intgroup = c("group","replicate"))

#GO ontologies to attribute functions to genes
#Example BaP vs DMSO
library(DOSE)
library(pathview)
library(clusterProfiler)
library(tidyft)

#Prepare input
original_gene_list <- res_BaP_PCB.refiltered$log2FoldChange

names(original_gene_list) <- rownames(res_BaP_PCB.refiltered)

gene_list <- na.omit(original_gene_list)

gene_list <- sort(gene_list, decreasing = TRUE)

#Gene Set Enrichment
organism <- org.Hs.eg.db
keytypes(org.Hs.eg.db)

gse <- gseGO(geneList = gene_list,
             ont = "All",
             keyType = "ENSEMBL",
             minGSSize = 4,
             maxGSSize = 800,
             pvalueCutoff = 0.01,
             verbose = TRUE,
             OrgDb = organism,
             pAdjustMethod = "none",
             seed = TRUE)

gse.simple <- simplify(gseGO(geneList = gene_list,
               ont = "BP",
               keyType = "ENSEMBL",
               minGSSize = 4,
               maxGSSize = 800,
               pvalueCutoff = 0.01,
               verbose = TRUE,
               OrgDb = organism,
               pAdjustMethod = "none",
               seed = TRUE), cutoff = 0.9,by = "p.adjust", select_fun = min)

dotplot(gse, showCategory = 20, x = "GeneRatio", font.size = 12)+
  facet_grid(.~.sign)+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 105))
gse <- dropGO(gse, term = "GO:0016705")

gene.df <- bitr(names(gene_list), fromType = "ENSEMBL",
                toType = c("ENTREZID","SYMBOL"),
                OrgDb = org.Hs.eg.db)
gene.df <- na.omit(gene.df)
dim(gene.df)

ggo <- groupGO(gene = gene.df$ENTREZID,
               OrgDb = org.Hs.eg.db,
               ont = "CC",
               level = 3,
               readable = TRUE)
head(ggo)

ego <- enrichGO(gene = gene.df$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.1,
                readable = TRUE)
head(ego)

heatplot(gse, foldChange = gene_list)
#Enrichment map
emapplot(gse, showCategory = 20)

#Ridgeplot
ridgeplot(gse)+
  labs(x = "enrichment distribution")

#Gene set enrichment analysis using clusterProfiler and Pathview
res_BaP_PCB.refiltered$ENTREZ <- ens.str.ENTREZ(result.table = res_BaP_PCB.refiltered)
res_BaP_PCB.refiltered

res_BaP_PCB.refiltered <- res_BaP_PCB.refiltered[which(duplicated(res_BaP_PCB.refiltered$ENTREZ) == F),]

ENTREZ.NA.filter <- is.na(res_BaP_PCB.refiltered$ENTREZ)
res_BaP_PCB.refiltered <- res_BaP_PCB.refiltered[!ENTREZ.NA.filter,]

foldchanges <- res_BaP_PCB.refiltered$log2FoldChange
names(foldchanges) <- res_BaP_PCB.refiltered$ENTREZ
foldchanges <- sort(foldchanges, decreasing = TRUE)

gseaKEGG <- gseKEGG(geneList = foldchanges,
                    organism = "hsa",
                    nPerm = 1000,
                    minGSSize = 5,
                    pvalueCutoff = 0.05,
                    verbose = FALSE)


gseaKEGG_result <- gseaKEGG@result
View(gseaKEGG_result)

