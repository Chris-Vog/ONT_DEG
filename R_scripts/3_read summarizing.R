# Read summarizing----
readCountTargets <- file.path("Analysis", "Minimap2", 
                              paste(tools::file_path_sans_ext(basename(as.character(studyDesign$filename)), compression=TRUE), ".sorted.bam", sep=""))

ExternalAnnotation = file.path("ReferenceData", basename(config$genome_annotation))

#Assign mapped sequencing reads to specified genomic features
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

#Extraction of the top 10 genes according to the total number of reads
knitr::kable(geneCounts[order(rowSums(geneCounts), decreasing=TRUE)[1:10],], caption="Table showing the 10 annotated gene features with the highest number of mapped reads", booktabs=TRUE, table.envir='table*', linesep="") %>%
  kable_styling(latex_options=c("hold_position", font_size=11)) %>%
  add_footnote(c("This is raw count data and no normalisation or transformation has been performed"))

#Removal of all genes with zero counts
geneCounts_nonZeros <- geneCounts[which(rowSums(geneCounts) > 0),]

geneCounts_nonZeros <- as.data.frame(geneCounts_nonZeros)

ens.geneCounts_nonZeros <- rownames(geneCounts_nonZeros)
Symbols <- ens.str.SYMBOL(geneCounts_nonZeros)
geneCounts_nonZeros$Symbols <- Symbols
dim(geneCounts_nonZeros)

geneCounts_nonZeros <- drop_na(geneCounts_nonZeros)
dim(geneCounts_nonZeros)

geneCounts_nonZeros <- geneCounts_nonZeros[!duplicated(geneCounts_nonZeros$Symbols),]

rownames(geneCounts_nonZeros) <- geneCounts_nonZeros$Symbols
geneCounts_nonZeros$Symbols <- NULL
dim(geneCounts_nonZeros)

#Change the names of the columns
colnames(geneCounts_nonZeros) <- paste0(studyDesign$group, "_", studyDesign$replicate) #It is important to choose different column names

#Ordering of the data frame
geneCounts_ordered <- geneCounts_nonZeros
geneCounts_ordered$RowSums <- rowSums(geneCounts_ordered)

# As rownames will not be exported, we create an additional column containing the gene_id
geneCounts_ordered$gene_id <- rownames(geneCounts_ordered)

geneCounts_ordered <- geneCounts_ordered[!duplicated(geneCounts_ordered$gene_id),]
# Order the samples with the arrange function of dplyr
geneCounts_ordered <- geneCounts_ordered %>%
  dplyr::arrange(desc(RowSums))

#Save the data set
xlsExpressedGenes <- file.path("Analysis/Results/ExpressedGenes.xlsx")
write_xlsx(x = geneCounts_ordered, path = xlsExpressedGenes)
