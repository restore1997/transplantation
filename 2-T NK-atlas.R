rm(list = ls())

library(Seurat)
library(dplyr)
library(ggplot2)


##############
load("liver.rda")
liver
table(liver$orig.ident)

liver$status <- "NA"
liver$status[liver$data%in%c("L1","L2","L4","L5")] <- "transplantation"
liver$status[liver$data%in%c("N1","N2","N3")] <- "health"


liver <- RunPCA(liver, features = VariableFeatures(object = liver))
DimPlot(liver, reduction = "pca")
DimPlot(liver, reduction = "pca",group.by = "diagnosis")
DimPlot(liver, reduction = "pca",group.by = "status")



liver <- FindNeighbors(liver, dims = 1:10)
liver <- FindClusters(liver, resolution = 0.8)

table(liver$seurat_clusters)

liver <- RunTSNE(liver, dims = 1:10)
liver <- RunUMAP(liver, dims = 1:10)

DimPlot(liver, reduction = "tsne",label = T)
DimPlot(liver, reduction = "tsne",group.by = "status")
DimPlot(liver, reduction = "tsne",group.by = "diagnosis")


###
VlnPlot(liver,features = "CD3D",pt.size = 0)## T cell cluster1,3,4,5,6,7,8,9,10,11,12
VlnPlot(liver,features = "CD3G",pt.size = 0)## T cell cluster1,3,4,5,6,7,8,9,10,11,12

VlnPlot(liver,features = "FGFBP2",pt.size = 0)## NK cell cluster2,5,9
VlnPlot(liver,features = "KLRF1",pt.size = 0)## NK cell cluster0,2,4,5,9,12
VlnPlot(liver,features = "GNLY",pt.size = 0)## NK cell cluster2,5,9,11
VlnPlot(liver,features = "CD7",pt.size = 0)## NK cell cluster0,1,2,3,4,5,6,7,8,9,10,11,12


VlnPlot(liver,features = "CD3D",pt.size = 0)
VlnPlot(liver,features = "CD4",pt.size = 0)
VlnPlot(liver,features = "CD8A",pt.size = 0)
VlnPlot(liver,features = "KLRF1",pt.size = 0)

VlnPlot(liver,features = c("CD3D","CD4","CD8A","KLRF1"),pt.size = 0)


FeaturePlot(liver,features = "CD8A",reduction = "tsne")
FeaturePlot(liver,features = "CD4",reduction = "tsne")
FeaturePlot(liver,features = "KLRF1",reduction = "tsne")

liver$MM2 <- "NA"
liver$MM2[liver$seurat_clusters%in%"0"] <- "NK_1"
liver$MM2[liver$seurat_clusters%in%"1"] <- "CD4T_1"
liver$MM2[liver$seurat_clusters%in%"2"] <- "NK_2"
liver$MM2[liver$seurat_clusters%in%"3"] <- "CD8T_1"
liver$MM2[liver$seurat_clusters%in%"4"] <- "CD8T_2"
liver$MM2[liver$seurat_clusters%in%"5"] <- "CD8T_3"
liver$MM2[liver$seurat_clusters%in%"6"] <- "CD4T_2"
liver$MM2[liver$seurat_clusters%in%"7"] <- "CD8T_4"
liver$MM2[liver$seurat_clusters%in%"8"] <- "CD4T_3"
liver$MM2[liver$seurat_clusters%in%"9"] <- "CD8T_5"
liver$MM2[liver$seurat_clusters%in%"10"] <- "CD8T_6"
liver$MM2[liver$seurat_clusters%in%"11"] <- "CD8T_7"
liver$MM2[liver$seurat_clusters%in%"12"] <- "CD4T_4"
liver$MM2[liver$seurat_clusters%in%"13"] <- "other"
liver$MM2[liver$seurat_clusters%in%"14"] <- "other"

DimPlot(liver, group.by="MM2", repel=T, reduction='tsne',label = T)

liver$MM3 <- "NA"
liver$MM3[liver$seurat_clusters%in%"0"] <- "XCL1+ NK"    ####
liver$MM3[liver$seurat_clusters%in%"1"] <- "IL7R+ CD4T"
liver$MM3[liver$seurat_clusters%in%"2"] <- "FGFBP2+ NK"
liver$MM3[liver$seurat_clusters%in%"3"] <- "GZMK+ CD8T"
liver$MM3[liver$seurat_clusters%in%"4"] <- "CRTAM+ CD8T"
liver$MM3[liver$seurat_clusters%in%"5"] <- "GZMH+ CD8T"
liver$MM3[liver$seurat_clusters%in%"6"] <- "RPS4Y1+ CD4T"  ####
liver$MM3[liver$seurat_clusters%in%"7"] <- "HSPA1B+ CD8T"    ####
liver$MM3[liver$seurat_clusters%in%"8"] <- "CCR6+ CD4T"
liver$MM3[liver$seurat_clusters%in%"9"] <- "HBB+ CD8T"     ###
liver$MM3[liver$seurat_clusters%in%"10"] <- "CTLA4+ CD8T"
liver$MM3[liver$seurat_clusters%in%"11"] <- "MKI67+ CD8T"
liver$MM3[liver$seurat_clusters%in%"12"] <- "CCL3+ CD8T"
liver$MM3[liver$seurat_clusters%in%"13"] <- "other"
liver$MM3[liver$seurat_clusters%in%"14"] <- "other"
DimPlot(liver, group.by="MM3", repel=T, reduction='tsne',label = T)


liver$MM4 <- "NA"
liver$MM4[liver$seurat_clusters%in%c("0","2")] <- "NK"
liver$MM4[liver$seurat_clusters%in%c("1","6","8")] <- "CD4T"
liver$MM4[liver$seurat_clusters%in%c("3","4","5","7","9","11","12")] <- "CD8T"
DimPlot(liver, group.by="MM4", repel=T, reduction='tsne',label = T)



liver.markers <- FindAllMarkers(liver, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
liver.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top5 <- liver.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(liver, features = top5$gene) #+ NoLegend()



save(liver,file = "1-atlas.rda")

meta <- liver@meta.data
save(meta,file = "1-meta.rda")

load("1-atlas.rda")

VlnPlot(liver,features = c("IL4","EDF1","IL10","IL17A","IFNG","TGFB1","FOXP3","RORC"),pt.size = 0)
#IL-4，IL-5，IL-10，IL-17，IFN-γ，TGF-β，FOXP3，RORγt


liver.markers <- FindAllMarkers(liver, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
liver.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top5 <- liver.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(liver, features = top5$gene) #+ NoLegend()
save(liver.markers,file = "1-markers.rda")
load("1-markers.rda")
write.csv(liver.markers,file = "Table3_T_NK_markers.csv")

#RidgePlot(liver,features = "CD3D")
DotPlot(liver,features = unique(top5$gene))

load("1-markers.rda")
pdf("top5.pdf")
top5 <- liver.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(liver, features = top5$gene) #+ NoLegend()
dev.off()

table(liver.markers$cluster)
top20 <- liver.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

top20_C0 <- top20[top20$cluster%in%"0",]
VlnPlot(liver,features = top20_C0$gene,pt.size = 0)

top20_C1 <- top20[top20$cluster%in%"1",]
VlnPlot(liver,features = top20_C1$gene,pt.size = 0)

top20_C2 <- top20[top20$cluster%in%"2",]
VlnPlot(liver,features = top20_C2$gene,pt.size = 0)

top20_C3 <- top20[top20$cluster%in%"3",]
VlnPlot(liver,features = top20_C3$gene,pt.size = 0)

top20_C4 <- top20[top20$cluster%in%"4",]
VlnPlot(liver,features = top20_C4$gene,pt.size = 0)

top20_C5 <- top20[top20$cluster%in%"5",]
VlnPlot(liver,features = top20_C5$gene,pt.size = 0)

top20_C6 <- top20[top20$cluster%in%"6",]
VlnPlot(liver,features = top20_C6$gene,pt.size = 0)

top20_C7 <- top20[top20$cluster%in%"7",]
VlnPlot(liver,features = top20_C7$gene,pt.size = 0)

top20_C8 <- top20[top20$cluster%in%"8",]
VlnPlot(liver,features = top20_C8$gene,pt.size = 0)

top20_C9 <- top20[top20$cluster%in%"9",]
VlnPlot(liver,features = top20_C9$gene,pt.size = 0)

top20_C10 <- top20[top20$cluster%in%"10",]
VlnPlot(liver,features = top20_C10$gene,pt.size = 0)

top20_C11 <- top20[top20$cluster%in%"11",]
VlnPlot(liver,features = top20_C11$gene,pt.size = 0)

top20_C12 <- top20[top20$cluster%in%"12",]
VlnPlot(liver,features = top20_C12$gene,pt.size = 0)

##exhaust
VlnPlot(liver,features = "CTLA4",pt.size = 0)
VlnPlot(liver,features = "PDCD1",pt.size = 0)
VlnPlot(liver,features = "HAVCR2",pt.size = 0)




load("1-atlas.rda")
table(liver$status)
transplantation <- subset(liver,subset = status%in%"transplantation");save(transplantation,file = "transplantation.rda")

VlnPlot(liver,features = c("ITGAM","ITGAX"),pt.size = 0)
VlnPlot(liver,features = c("CD4","IL2RA","FOXP3"))
VlnPlot(liver,features = c("CD4","CD8A"),pt.size = 0)
FeaturePlot(liver,features = c("CD4","CD8A","IL2RA","FOXP3"),reduction = "tsne")

table(liver$seurat_clusters)
C10 <- subset(liver,subset = seurat_clusters%in%"10")
save(C10,file = "C10.rda")
