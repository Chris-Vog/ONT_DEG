# Differential gene expression (DGE) Analyse----
## Eine DGE-Analyse kann immer nur zwischen zwei Bedingungen durchgef?hrt werden.

# Erstellung eines DESeq2-Datensatzes um die importierten Daten zu speichern
## Dieses Datei-Objekt wird f?r die folgende DGE-Analyse ben?tigt und
## ist auch Grundlage f?r weitere funktionale Analysen.
## Hier wird zun?chst ?berpr?ft, ob der Datensatz eine N-Zahl > 1 aufweist. 
group_size <- 2
if (length(unique(studyDesign$group)) == length(studyDesign$group)){
  group_size <- 1
}

if (group_size > 1) { 
  deSeqRaw <- DESeqDataSetFromMatrix(countData=geneCounts, colData=studyDesign, design= ~ group)
} else {
  deSeqRaw <- DESeqDataSetFromMatrix(countData=geneCounts, colData=studyDesign, design= ~ 1)  
}

## countData = Matrix mit den Gene counts
## colData = studyDesign -> Informationen ?ber die Proben
## design = Eine Formel, die die Abh?ngigkeit der Counts zu den Proben beschreibt

#Differential gene expression - Analyse basierend auf der negativen Binomialverteilung
dds <- DESeq(deSeqRaw)
keep <- rowSums(counts(dds)) >= config$cutoff_geneCount # Dieser Wert ist abh?ngig von der Anzahl von reads und der Probenzahl
dds <- dds[keep,] # Hier wird der Datensatz durch den Filter `keep` indexiert

# Faktor level
## Voreingestellt wählt R ein Referenzlevel basierend auf der alphabetischen Reihenfolge aus.
## Daher muss dem Programm mitgeteilt werden, welche Proben die Kontrollgruppe repr?sentieren.
## Hierzu gibt es verschiedene M?glichkeiten, die beide gleich korrekt sind.
## 1. Releveln (empfohlen)
dds$group <- relevel(dds$group, ref = config$referenceGroup)
dds <- nbinomWaldTest(dds)
## 2. Faktorisierung
#dds$group <- factor(dds$group, levels = c("DMSO" , "PCB" , "BaP"))

resultsNames(dds) # Die Ergebnisse ändern sich, je nachdem welche Gruppe als Referenzgruppe verwendet wird.

res <- results(dds, name = "group_BaP_vs_DMSO")

## Bei mehr als einer Bedingung wird es etwas komplizierter in der Berechnung. 
## Die zu vergleichenden Gruppen können auch mit Hilfe des Attributs 'contrast' spezifiziert werden.
#res2 <- results(dds, contrast = c("group" , "BaP" , "PCB")) # Hierbei wird BaP über PCB verglichen. In diesem Fall wird der LFC von PCB auf 0 gesetzt.

# Beispiel-Analyse
res <- results(dds, contrast = c("group" , config$treatedGroup , config$referenceGroup)) # 1. Erster Characterstring = Spalte in studyDesign! 2. und 3. Characterstring = Gruppen aus dieser Spalte

summary(res) # Zusammenfassende Statistiken der Analyse

## Bisher sind die Gene lediglich mit ihrer ENSEMBL ID angegeben. Diese kann in die ENTREZ ID oder in das RefSeq Symbol ?berf?hrt werden.
## Die dazu n?tigen Funktionen sind im common.R script vorhanden.
res$SYMBOL <- ens.str.SYMBOL(result.table = res)
res$ENTREZ <- ens.str.ENTREZ(result.table = res)
res$ENSEMBLID <- rownames(res)

## Speichern der DEG-Analyse
res <- res[order(res$padj, decreasing = FALSE),]

res.filePath.table <- paste0("Analysis/Results/", config$treatedGroup , "_" , config$referenceGroup , ".txt")
write.table(x = res, file = res.filePath.table, sep = "\t")

res.filePath.xlsx <- paste0("Analysis/Results/", config$treatedGroup , "_" , config$referenceGroup , ".xlsx")
write_xlsx(x = as.data.frame(res) , path = res.filePath.xlsx)

#Visualisierung
## LFC Minderung f?r Visualisierungen
resLFC <- lfcShrink(dds, contrast = c("group" , config$treatedGroup , config$referenceGroup), type = "ashr") # Hierzu muss das package "ashr installiert sein.
### Alternativ kann auch "apeglm" verwendet werden. Dieser Estimator funktioniert aber nicht mit dem Argument 'contrast', sondern nur mit 'coef'.

## MA-Plot
## Hier werden die LFCs (y-Achse) gegen die Anzahl der counts (x-Achse) aufgetragen.
DESeq2::plotMA(resLFC, ylim = c(-2,2),
               colNonSig = "gray32", colSig = "red3",
               log = "x")+
  title(paste0(config$treatedGroup , " vs. " , config$referenceGroup))