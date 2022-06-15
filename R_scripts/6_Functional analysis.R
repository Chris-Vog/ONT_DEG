# Funktionelle Analysen
## Analysen wie die Gene Set Enrichment (GSE) Analyse helfen die gewonnen Daten in einem Kontext zu betrachten
## Als Grundlage für diese Analyse dienen die L2FC Daten aus der DEG Analyse
## Wir bleiben bei dem BaP vs DMSO Beispiel
# Gene Set enrichment analysis----
res_BaP_DMSO.filtered <- res_BaP_DMSO %>%
  as.data.frame() %>%
  dplyr::filter(abs(log2FoldChange) > config$lfcThreshold)

##Erstellung eines Datensatzes bestehend aus Gen und L2FC
geneList <- res_BaP_DMSO.filtered$log2FoldChange
names(geneList) <- rownames(res_BaP_DMSO.filtered)
geneList <- sort(geneList, decreasing = TRUE)

# Gene Set Enrichment Analysis of Gene Ontology
## Die erstellte Genliste dient nun als Grundlage für die Funktion 'gseGO'.
gse <- gseGO(geneList = geneList,
             ont = "ALL", # Sub-Geneontology
             keyType = "ENSEMBL", # Datenbank 
             nPerm = 10000, # Number of Permutations. Je höher, desto akkurater wird der p-Wert berechnet. Höhere Nummern benötigen mehr Rechenleistung. 
             minGSSize = config$minGSSize,
             maxGSSize = config$maxGSSize,
             pvalueCutoff = 0.05, # min p-Wert
             verbose = TRUE, # Gibt einen Output in der Konsole aus
             OrgDb = organism, # Angabe der Datenbank
             pAdjustMethod = "none")

##Visualisierung
dotplot(gse, showCategory = 35, split =".sign", orderBy = "x")+
  facet_grid(.~.sign)

## Azeigen der gesamten Ergebnisse
## Da in der Visualisierung nur einige Ergebnisse dargestellt sind, 
## ist es sinnvoll sich alle Ergebnisse in einer Tabelle anzeigen zu lassen.
View(as.data.frame(gse))

## Bei den GO-Terms gibt es oft sehr hohe Überschneidungen.
## Diese können reduziert und in einem GO-Term zusammengefasst werden.
gse.simple<- simplify(x = gseGO(geneList = geneList,
                                ont = "ALL", 
                                keyType = "ENSEMBL",  
                                nPerm = 10000, 
                                minGSSize = config$minGSSize,
                                maxGSSize = config$maxGSSize,
                                pvalueCutoff = 0.05, 
                                verbose = TRUE, 
                                OrgDb = organism,
                                pAdjustMethod = "none"),
                             cutoff = 0.7, # Prozentuale Überschneidung der Gene
                             by = "p.adjust",
                             select_fun = min)

dotplot(gse.simple, showCategory = 35, split =".sign", orderBy = "x")+
  facet_grid(.~.sign)

# Speichern der Ergebnisse
## Damit die Ergebnisse gespeichert werden können, müssen diese in ein DataFrame umgewandelt werden
gse.simple.df <- as.data.frame(gse.simple)
pathway <- file.path("Analysis/Results/pathway.xlsx")
write_xlsx(x = gse.simple.df, path = pathway)