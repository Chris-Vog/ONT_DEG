# Check md5sum for uniqueness
md5sum <- function(filename) digest::digest(filename, algo="md5", file=TRUE)

# Load flagstat targets
loadFlagstat <- function(file) {
  x <- read.table(file, header = FALSE, sep= " ", fill =NA)[c(1:8),c(1:3)]
  x[,1]
}

ncalc <- function(len.vector, n) {
  # N50 - length such that scaffolds of this length or longer include half the bases of the assembly
  len.sorted <- rev(sort(len.vector))
  len.sorted[cumsum(len.sorted) >= sum(len.sorted)*n][1]
}


lcalc <- function(len.vector, n) {
  len.sorted <- rev(sort(len.vector))
  which(cumsum(len.sorted) >= sum(len.sorted)*n)[1]
}

#Gene name conversion
ah <- AnnotationHub()
ahDb <- query(ah, pattern = c("Homo sapiens","EnsDb", 104))
edb.v104 <- ahDb[[1]]

#Conversion form ENSEMBL ID to SYMBOL
ens.str.SYMBOL <- function(result.table) {
  ens.str.vector <- sub("\\..*","",rownames(result.table))
  mapIds(edb.v104, keys = ens.str.vector, column = "SYMBOL", keytype = "GENEID", multiVals = "first")
}

#Conversion from ENSEMBL ID to ENTREZ
ens.str.ENTREZ <- function(result.table) {
  ens.str.vector <- sub("\\..*","",rownames(result.table))
  mapIds(edb.v104, keys = ens.str.vector, column = "ENTREZID", keytype = "GENEID", multiVals = "first")
}
