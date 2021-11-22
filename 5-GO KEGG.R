rm(list = ls())
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)


load("../1-markers.rda")
table(liver.markers$cluster)

genes <- liver.markers[liver.markers$cluster%in%"0",]$gene
GO_C0 <- enrichGO(gene          =genes,
                    #universe     = row.names(dge.celltype),
                    OrgDb         = 'org.Hs.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
dotplot(GO_C0)
GO_C0_meta <- as.data.frame(GO_C0)
GO_C0_meta$data <- "GO"

#KEGG
genelist <- bitr(genes, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
genelist <- pull(genelist,ENTREZID)               
KEGG_C0 <- enrichKEGG(gene = genelist, organism = 'hsa')
dotplot(KEGG_C0, showCategory=10)
KEGG_C0_meta <- as.data.frame(KEGG_C0)
KEGG_C0_meta$ONTOLOGY <- "NA"
KEGG_C0_meta$data <- "KEGG"

save(GO_C0,KEGG_C0,file = "C0.rda")
load("C0.rda")

C0_meta <- rbind(GO_C0_meta,KEGG_C0_meta)
write.csv(C0_meta,file = "TABLES/C0.csv")
