# Venn diagram---- 
## Hier wurden die Daten jeweils mit DMSO verglichen
res_PCB_DMSO <- readxl::read_xlsx("Analysis/Results/PCB_DMSO.xlsx")
res_BaP_DMSO <- readxl::read_xlsx("Analysis/Results/BaP_DMSO.xlsx")

PCB.total <- res_PCB_DMSO[abs(res_PCB_DMSO$log2FoldChange) > config$lfcThreshold , "ENSEMBLID"] # Muss individuell angepasst werden!
#PCB.up <- res_PCB_DMSO[res_PCB_DMSO$log2FoldChange > config$lfcThreshold , "ENSEMBLID"]
#PCB.down <- res_PCB_DMSO[res_PCB_DMSO$log2FoldChange < -config$lfcThreshold , "ENSEMBLID"]
#dim(PCB.total)

BaP.total <- res_BaP_DMSO[abs(res_BaP_DMSO$log2FoldChange) > config$lfcThreshold, "ENSEMBLID"] # Muss individuell angepasst werden!
#BaP.up <- res_BaP_DMSO[res_BaP_DMSO$log2FoldChange > config$lfcThreshold, "ENSEMBLID"]
#BaP.down <- res_BaP_DMSO[res_BaP_DMSO$log2FoldChange < -config$lfcThreshold, "ENSEMBLID"]
#dim(BaP.total)

#sheets <- list("PCB.total" = PCB.total, 
# "PCB.up" = PCB.up,
# "PCB.down" = PCB.down,
# "BaP.total" = BaP.total,
# "BaP.up" = BaP.up,
# "BaP.down" = BaP.down)

#write_xlsx(sheets, "Analysis/Results/DE.complete.xlsx")

#Venn diagram of differential expressed genes
dev.off()
draw.pairwise.venn(area1 = dim(PCB.total)[1], #Genes PCB
                   area2 = dim(BaP.total)[1], #Genes BaP
                   cross.area = length(base::intersect(PCB.total[[1]], BaP.total[[1]])),
                   category = c("PCB126 (714)", "B[a]P (1041)"),
                   euler.d = T,
                   fill = c("#F59B25", "#95E354"),
                   alpha = 0.5,
                   fontfamily = rep("Helvetica",3), # Zahl = Anzahl der Regionen
                   cat.fontfamily = rep("Helvetica",2), # Zahl = Anzahl der Bedingungen
                   cex = 0.9,
                   cat.cex = 0.9,
                   cat.dist = 0.1,
                   ext.text = FALSE)
