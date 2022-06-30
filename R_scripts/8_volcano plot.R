# Volcano plot----
## Im Volcano plot sind die adjusted p-Values (y-Achse) gegen die log2 Fold Changes (x-Achse) aufgetragen
## Farblich markiert sind Datenpunkte (Gene), die:
## einen gewissen L2FC aufweisen
## einen definierten adjusted p-Value unterschreiten
## beide Bedingungen erfüllen

## Definition der Gene, die im Volcano plot beschriftet werden
volcano.genes <- head(res[, "SYMBOL"], n = config$genes_volcano)
EnhancedVolcano(res,
                lab = res$SYMBOL,
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = volcano.genes,
                xlab = bquote(~Log[2]~ 'fold change'),
                title = paste0(config$treatedGroup, " vs. " , config$referenceGroup),
                pCutoff = as.numeric(config$pCutoff_Volcano),
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 5.0,
                colAlpha = 0.5,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                colConnectors = "black")
