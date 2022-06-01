# Differential gene expression (DGE) Analyse
## Eine DGE-Analyse kann immer nur zwischen zwei Bedingungen durchgeführt werden.

# Erstellung eines DESeq2-Datensatzes um die importierten Daten zu speichern
## Dieses Datei-Objekt wird für die folgende DGE-Analyse benötigt und
## ist auch Grundlage für weitere funktionale Analysen.
## Hier wird zunächst überprüft, ob der Datensatz eine N-Zahl > 1 aufweist. 
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
## colData = studyDesign -> Informationen über die Proben
## design = Eine Formel, die die Abhängigkeit der Counts zu den Proben beschreibt

#Differential gene expression - Analyse basierend auf der negativen Binomialverteilung
dds <- DESeq(deSeqRaw)
keep <- rowSums(counts(dds)) >= config$cutoff_geneCount # Dieser Wert ist abhängig von der Anzahl von reads und der Probenzahl
dds <- dds[keep,] # Hier wird der Datensatz durch den Filter `keep` indexiert

# Faktor level
## Voreingestellt wählt R ein Referenzlevel basierend auf der alphabetischen Reihenfolge aus.
## Daher muss dem Programm mitgeteilt werden, welche Proben die Kontrollgruppe repräsentieren.
## Hierzu gibt es verschiedene Möglichkeiten, die beide gleich korrekt sind.
## 1. Releveln (empfohlen)
dds$group <- relevel(dds$group, ref = config$referenceGroup)

## 2. Faktorisierung
#dds$group <- factor(dds$group, levels = c("DMSO" , "PCB" , "BaP"))

resultsNames(dds) # Die Ergebnisse ändern sich, je nachdem welche Gruppe als Referenzgruppe verwendet wird.

res <- results(dds, name = "group_PCB_vs_DMSO")

## Bei mehr als einer Bedingung wird es etwas komplizierter in der Berechnung. 
## Die zu vergleichenden Gruppen können auch mit Hilfe des Attributs 'contrast' spezifiziert werden.
res2 <- results(dds, contrast = c("group" , "BaP" , "PCB")) # Hierbei wird BaP über PCB verglichen. In diesem Fall wird der LFC von PCB auf 0 gesetzt.

## Exemplarische DEG-Analyse für BaP vs. DMSO
res_BaP_DMSO <- results(dds, contrast = c("group" , "BaP" , "DMSO")) # 1. Erster Characterstring = Spalte in studyDesign! 2. und 3. Characterstring = Gruppen aus dieser Spalte

## MA-Plot
DESeq2::plotMA(res_BaP_DDMSO, ylim = c(-5,5),
               colNonSig = "gray32", colSig = "red3",
               log = "x")+
  title("B[a]P vs. DMSO")





