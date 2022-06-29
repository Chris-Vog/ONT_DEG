# Venn diagram---- 
## Hier wurden die Daten jeweils mit DMSO verglichen
PCB.total <- as.data.frame(rownames(res_PCB_DMSO[abs(res_PCB_DMSO$log2FoldChange) > 1.5 ,]))
PCB.up <- as.data.frame(rownames(res_PCB_DMSO[res_PCB_DMSO$log2FoldChange > 1.5 ,]))
PCB.down <- as.data.frame(rownames(res_PCB_DMSO[res_PCB_DMSO$log2FoldChange < -1.5 ,]))
dim(PCB.total)

BaP.total <- as.data.frame(rownames(res_BaP_DMSO[abs(res_BaP_DMSO$log2FoldChange) > 1.5, ]))
BaP.up <- as.data.frame(rownames(res_BaP_DMSO[res_BaP_DMSO$log2FoldChange > 1.5, ]))
BaP.down <- as.data.frame(rownames(res_BaP_DMSO[res_BaP_DMSO$log2FoldChange < -1.5, ]))
dim(BaP.total)

BaP_PD153035.total <- as.data.frame(rownames(res_BaP_PD153035_DMSO[abs(res_BaP_PD153035_DMSO$log2FoldChange) > 1.5, ]))
BaP_PD153035.up <- as.data.frame(rownames(res_BaP_PD153035_DMSO[res_BaP_PD153035_DMSO$log2FoldChange > 1.5 ,]))
BaP_PD153035.down <- as.data.frame(rownames(res_BaP_PD153035_DMSO[res_BaP_PD153035_DMSO$log2FoldChange < -1.5 ,]))
dim(BaP_PD153035.total)

sheets <- list("PCB.total" = PCB.total, 
               "PCB.up" = PCB.up,
               "PCB.down" = PCB.down,
               "BaP.total" = BaP.total,
               "BaP.up" = BaP.up,
               "BaP.down" = BaP.down,
               "BaP_PD153035.total" = BaP_PD153035.total,
               "BaP_PD153035.up" = BaP_PD153035.up,
               "BaP_PD153035.down" = BaP_PD153035.down)

write_xlsx(sheets, "Analysis/DE.complete.xlsx")

#Venn diagram of differential expressed genes
library(VennDiagram)

draw.triple.venn(area1 = dim(PCB.total)[1], #Genes PCB
                 area2 = dim(BaP.total)[1], #Genes BaP
                 area3 = dim(BaP_PD153035.total)[1], #Genes BaP_PD153035
                 n12 = length(base::intersect(PCB.total[[1]], BaP.total[[1]])),
                 n13 = length(base::intersect(PCB.total[[1]], BaP_PD153035.total[[1]])),
                 n23 = length(base::intersect(BaP.total[[1]], BaP_PD153035.total[[1]])),
                 n123 = length(base::intersect(intersect(PCB.total[[1]], BaP.total[[1]]), BaP_PD153035.total[[1]])),
                 category = c("PCB126 (910)", "B[a]P (1434)", "B[a]P + PD153035 (2062)"),
                 euler.d = TRUE,
                 fill = c("#F59B25", "#95E354", "#54A7E3"),
                 alpha = 0.5,
                 fontfamily = rep("Helvetica",7),
                 cat.fontfamily = rep("Helvetica",3),
                 cex = 0.9,
                 cat.cex = 0.9)
