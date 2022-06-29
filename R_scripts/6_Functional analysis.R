# Funktionelle Analysen----
## Analysen wie die Gene Set Enrichment (GSE) Analyse helfen die gewonnen Daten in einem Kontext zu betrachten
## Als Grundlage f?r diese Analyse dienen die L2FC Daten aus der DEG Analyse
## Wir bleiben bei dem BaP vs DMSO Beispiel
# Gene Set enrichment analysis
res_BaP_DMSO.filtered <- res_BaP_DMSO %>%
  as.data.frame() %>%
  dplyr::filter(abs(log2FoldChange) > config$lfcThreshold)

##Erstellung eines Datensatzes bestehend aus Gen und L2FC
geneList <- res_BaP_DMSO.filtered$log2FoldChange
names(geneList) <- rownames(res_BaP_DMSO.filtered)
geneList <- sort(geneList, decreasing = TRUE)

organism <- org.Hs.eg.db
keytypes(organism)
# Gene Set Enrichment Analysis of Gene Ontology
## Die erstellte Genliste dient nun als Grundlage f?r die Funktion 'gseGO'.
gse <- gseGO(geneList = geneList,
             ont = "ALL", # Sub-Geneontology
             keyType = "ENSEMBL", # Datenbank 
             nPerm = 10000, # Number of Permutations. Je h?her, desto akkurater wird der p-Wert berechnet. H?here Nummern ben?tigen mehr Rechenleistung. 
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

## Bei den GO-Terms gibt es oft sehr hohe ?berschneidungen.
## Diese k?nnen reduziert und in einem GO-Term zusammengefasst werden.
gse.simple<- simplify(x = gseGO(geneList = geneList,
                                ont = "MF", # simplify-Methode funktioniert nicht mit "ALL"
                                keyType = "ENSEMBL",  
                                nPerm = 10000, 
                                minGSSize = config$minGSSize,
                                maxGSSize = config$maxGSSize,
                                pvalueCutoff = 0.05, 
                                verbose = TRUE, 
                                OrgDb = organism,
                                pAdjustMethod = "none"),
                      cutoff = 0.7, # Prozentuale ?berschneidung der Gene
                      by = "p.adjust",
                      select_fun = min)

dotplot(gse.simple, showCategory = 10, split =".sign", orderBy = "x")+
  facet_grid(.~.sign)

# Speichern der Ergebnisse
## Damit die Ergebnisse gespeichert werden k?nnen, m?ssen diese in ein DataFrame umgewandelt werden
gse.simple.df <- as.data.frame(gse.simple)
pathway <- file.path("Analysis/Results/pathway.xlsx")
write_xlsx(x = gse.simple.df, path = pathway)

# Umwandlung der ENSEMBL-ID in die entrez-ID----
geneList.entrez <- geneList
de_entrez <- geneList[abs(geneList.entrez) > 1]
de_entrez <- as.data.frame(de_entrez)

entr1 <- ens.str.ENTREZ(de_entrez)
entr1 <- as.data.frame(entr1)

de_entrez$entrez <- entr1$entr1
de_entrez <- na.omit(de_entrez)

de1 <- de_entrez$de_entrez
names(de1) <- de_entrez$entrez

# DGN Network
# Erstellung eines Balkendiagramms mit assoziierten Erkrankungen nach der Disease Gene Network (DGN) Datenbank
edo <- enrichDGN(names(de1))
barplot(edo, showCategory = 20)