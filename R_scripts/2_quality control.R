#Quality control----
#Creating a statistic-table of fastq-files
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
    n50 = ncalc(width(fastq), n=0.5), #number of reads, which scaffold 50 % of assembled bases
    l50 = lcalc(width(fastq), n=0.5), #length of reads, which scaffold 50 % of assembled bases
    n90 = ncalc(width(fastq), n=0.9), #number of reads, which scaffold 90 % of assembled bases
    l90 = lcalc(width(fastq), n=0.9)  #length of reads, which scaffold 90 % of assembled bases
  )
}

data <- lapply(row.names(studyDesign), processQCFastq) #applies this function to a vector of samples
qcData <- data.frame(data) #transformation of the data into a data frame
colnames(qcData) <- row.names(studyDesign)

# Darstellung der Daten
knitr::kable(qcData, caption="Summary statistics for the cDNA libraries imported", booktabs=TRUE, table.envir='table*', linesep="")  %>%
  kable_styling(latex_options=c("hold_position", font_size=9))

#Extraction and visualization of read length
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
lengthMatrixMelt <- cbind(lengthMatrixMelt, group=studyDesign[match(lengthMatrixMelt$variable,  rownames(studyDesign)), "group"],
                          replicate=studyDesign[match(lengthMatrixMelt$variable,  rownames(studyDesign)), "replicate"])
lengthMatrixMelt$replicate <- factor(lengthMatrixMelt$replicate)

plot <- ggplot(lengthMatrixMelt, aes(x=variable, y=value, fill=replicate)) + # fill = replicate kann durch fill = group ersetzt werden.
  geom_violin() + 
  scale_y_continuous(limits=c(0, as.numeric(quantile(lengthMatrixMelt$value, probs=c(0.975))))) + 
  xlab("Samples") +  
  ylab("Distribution Read length (bp)") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_brewer(palette="Paired") + 
  labs(title="Violin plot showing distribution of read lengths across samples")

suppressWarnings(print(plot))

#Extraction and visualization of quality scores
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
qualityMatrixMelt <- cbind(qualityMatrixMelt, group=studyDesign[match(qualityMatrixMelt$variable,  rownames(studyDesign)), "group"],
                           replicate=studyDesign[match(qualityMatrixMelt$variable,  rownames(studyDesign)), "replicate"])
qualityMatrixMelt$replicate <- factor(qualityMatrixMelt$replicate)

plotQ <- ggplot(qualityMatrixMelt, aes(x=variable, y=value, fill=replicate)) + # fill = replicate kann durch fill = group ersetzt werden.
  geom_violin() + 
  scale_y_continuous(limits=c(min(qualityMatrixMelt$value), max(qualityMatrixMelt$value))) + 
  xlab("study sample") +  ylab("Distribution of Read Qualities (QV)") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_brewer(palette="Paired") + 
  labs(title="Violin plot showing distribution of read qualities across samples")

suppressWarnings(print(plotQ))

#Analysis of flagstat results
flagstatTargets <- file.path("Analysis", "flagstat", 
                             paste(tools::file_path_sans_ext(basename(as.character(studyDesign$filename)), compression=TRUE),".txt",sep=""))

loadFlagstat <- function(file) {
  x <- read.table(file, header=FALSE, sep=" ", fill=NA)[c(1:5,7),c(1,3)]
  x[,1]
}

flagstatRes <- data.frame(matrix(unlist(lapply(flagstatTargets, loadFlagstat)), ncol=length(flagstatTargets)), stringsAsFactors = FALSE)
colnames(flagstatRes) <- rownames(studyDesign)
rownames(flagstatRes) <- c("total reads", "Primary" , "Secondary", "Supplementary", "Duplicates", "Mapped")

getVal <- function(word) {
  sum(as.numeric(unlist(strsplit(word, "/"))))
}

zreads <- unlist(lapply(flagstatRes["total reads", ], getVal)) -
  unlist(lapply(flagstatRes["Primary", ], getVal)) -
  unlist(lapply(flagstatRes["Secondary", ], getVal)) - 
  unlist(lapply(flagstatRes["Supplementary", ], getVal)) -
  unlist(lapply(flagstatRes["Duplicates", ], getVal))

flagstatRes[nrow(flagstatRes)+1,] <- round(as.numeric(flagstatRes["Mapped", ]) / as.numeric(flagstatRes["total reads", ]) * 100, digits = 2)
rownames(flagstatRes)[7] <- "%mapping"

knitr::kable(flagstatRes, caption="Summary statistics from the minimap2 long read spliced mapping.", booktabs=TRUE, table.envir='table*', linesep="")  %>%
  kable_styling(latex_options=c("hold_position", font_size=11)) %>%
  add_footnote(c("information from samtools flagstat"))