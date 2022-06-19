# Differential gene expression (DGE) Analyse
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
  deSeqRaw <- DESeqDataSetFromMatrix(countData=geneCounts, colData=studyDesign, design= ~group)
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

# Exemplarische DEG-Analyse für BaP vs. DMSO
res_BaP_DMSO <- results(dds, contrast = c("group" , "BaP" , "DMSO")) # 1. Erster Characterstring = Spalte in studyDesign! 2. und 3. Characterstring = Gruppen aus dieser Spalte

summary(res_BaP_DMSO) # Zusammenfassende Statistiken der Analyse

## Bisher sind die Gene lediglich mit ihrer ENSEMBL ID angegeben. Diese kann in die ENTREZ ID oder in das RefSeq Symbol ?berf?hrt werden.
## Die dazu n?tigen Funktionen sind im common.R script vorhanden.
res_BaP_DMSO$SYMBOL <- ens.str.SYMBOL(result.table = res_BaP_DMSO)
res_BaP_DMSO$ENTREZ <- ens.str.ENTREZ(result.table = res_BaP_DMSO)

## Speichern der DEG-Analyse
res_BaP_DMSO <- res_BaP_DMSO[order(res_BaP_DMSO$padj, decreasing = FALSE),]
write.table(x = res_BaP_DMSO, file = "Analysis/BaP_DMSO.txt", sep = "\t")

#Visualisierung
## LFC Minderung f?r Visualisierungen
resLFC <- lfcShrink(dds, contrast = c("group" , "BaP" , "DMSO"), type = "ashr")# Hierzu muss das package "ashr installiert sein.
### Alternativ kann auch "apeglm" verwendet werden. Dieser Estimator funktioniert aber nicht mit dem Argument 'contrast', sondern nur mit 'coef'.

## MA-Plot
## Hier werden die LFCs (y-Achse) gegen die Anzahl der counts (x-Achse) aufgetragen.
DESeq2::plotMA(resLFC, ylim = c(-5,5),
               colNonSig = "gray32", colSig = "red3",
               log = "x")+
  title("B[a]P vs. DMSO")








