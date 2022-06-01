# Zuordnung der Reads
## Beim Mapping wurden die Reads nur Bereichen im Genom zugeordnet.
## Im folgendem Schritt werden diese Bereiche Genen zugeordnet und diese gezählt.

# Pfad zu den .bam-Dateien, die die Informationen der gemappten Reads enthalten
readCountTargets <- file.path("Analysis", "Minimap", 
                              paste(tools::file_path_sans_ext(basename(as.character(studyDesign$filename)), compression=TRUE), ".bam", sep=""))

# Pfad zu der Gen-Annotation
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

colnames(geneCounts) <- rownames(studyDesign) # Änderung der Spaltennamen

# Entfernung der Gene ohne Reads
geneCounts_nonZeros <- geneCounts[which(rowSums(geneCounts) > 0) , ]

# Änderungen der Spaltennamen----
colnames(geneCounts_nonZeros) <- paste0(studyDesign$group, "_", studyDesign$replicate) #It is important to choose different column names

# Ordnen der Gene nach kumulative Read-Anzahl der Gene über alle Proben
geneCounts_nonZeros$RowSums <- rowSums(geneCounts_nonZeros)
geneCounts_ordered <- geneCounts_nonZeros %>%
  dplyr::arrange(desc(RowSums))

geneCounts_ordered$RowSums <- NULL

# Speichern des Datensatzes
#Save the data set----
xlsExpressedGenes <- file.path("ExpressedGenes.xlsx")
write_xlsx(x = geneCounts_ordered, path = xlsExpressedGenes)