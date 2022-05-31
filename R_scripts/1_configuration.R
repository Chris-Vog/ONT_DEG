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
