# Qualitäts-Kontrolle ONT generierter Daten
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