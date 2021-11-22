rm(list = ls())

library(Seurat)
library(dplyr)
library(ggplot2)


data1 <- Read10X(data.dir = "../analysis/data/B2/filtered_feature_bc_matrix/data/")
data1 <- CreateSeuratObject(counts = data1, project = "B2",
                   min.cells = 3, min.features = 200)

data1

data2 <- Read10X(data.dir = "../analysis/data/B3/filtered_feature_bc_matrix/data/")
data2 <- CreateSeuratObject(counts = data2, project = "B3",
                         min.cells = 3, min.features = 200)
data2

data3 <- Read10X(data.dir = "../analysis/data/B4/filtered_feature_bc_matrix/data/")
data3 <- CreateSeuratObject(counts = data3, project = "B4",
                            min.cells = 3, min.features = 200)
data3

data4 <- Read10X(data.dir = "../analysis/data/B5/filtered_feature_bc_matrix/data/")
data4 <- CreateSeuratObject(counts = data4, project = "B5",
                            min.cells = 3, min.features = 200)
data4


data5 <- Read10X(data.dir = "../analysis/data/L1/filtered_feature_bc_matrix/data/")
data5 <- CreateSeuratObject(counts = data5, project = "L1",
                            min.cells = 3, min.features = 200)
data5

data6 <- Read10X(data.dir = "../analysis/data/L2/filtered_feature_bc_matrix/data/")
data6 <- CreateSeuratObject(counts = data6, project = "L2",
                            min.cells = 3, min.features = 200)
data6


data7 <- Read10X(data.dir = "../analysis/data/L4/data/")
data7 <- CreateSeuratObject(counts = data7, project = "L4",
                            min.cells = 3, min.features = 200)
data7

data8 <- Read10X(data.dir = "../analysis/data/L5/data/")
data8 <- CreateSeuratObject(counts = data8, project = "L5",
                            min.cells = 3, min.features = 200)
data8

##health
data9 <- Read10X(data.dir = "../analysis/data/N1/data/")
data9 <- CreateSeuratObject(counts = data9, project = "N1",
                            min.cells = 3, min.features = 200)
data9


data10 <- Read10X(data.dir = "../analysis/data/N2/data/")
data10 <- CreateSeuratObject(counts = data10, project = "N2",
                             min.cells = 3, min.features = 200)
data10

data11 <- Read10X(data.dir = "../analysis/data/N3/data/")
data11 <- CreateSeuratObject(counts = data11, project = "N3",
                             min.cells = 3, min.features = 200)
data11

total <- merge(data1,y=c(data2,data3,data4,data5,data6,data7,data8,
                         data9,data10,data11),project = "total")

total
table(total$orig.ident)
total[["percent.mt"]] <- PercentageFeatureSet(total, pattern = "^MT-")
VlnPlot(total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)

plot1 <- FeatureScatter(total, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(total, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

total
total <- subset(total, subset = nFeature_RNA > 500  & nCount_RNA >2000 & percent.mt < 25)

total
table(total$orig.ident)

save(total,file = "raw_total.rda")


##############
load("raw_total.rda")
total
table(total$orig.ident)

###meta data
total$data <- total$orig.ident

table(total$data)
total$patient <- "NA"
total$patient[total$data%in%"B2"] <- "patient1"
total$patient[total$data%in%"B3"] <- "patient2"
total$patient[total$data%in%"B4"] <- "patient3"
total$patient[total$data%in%"B5"] <- "patient4"
total$patient[total$data%in%"L1"] <- "patient1"
total$patient[total$data%in%"L2"] <- "patient4"
total$patient[total$data%in%"L4"] <- "patient5"
total$patient[total$data%in%"L5"] <- "patient6"
total$patient[total$data%in%"N1"] <- "patient7"
total$patient[total$data%in%"N2"] <- "patient8"
total$patient[total$data%in%"N3"] <- "patient9"

total$pre_diagnosis <- "normal"
total$pre_diagnosis[total$orig.ident%in%c("B2","B3","B4","L1")] <- "tumor"
total$pre_diagnosis[total$orig.ident%in%c("B5","L2","L4","L5")] <- "biliary_atresia"

total$age_year <- "NA"
total$age_year[total$patient%in%"patient1"] <- "56"
total$age_year[total$patient%in%"patient2"] <- "66"
total$age_year[total$patient%in%"patient3"] <- "50"
total$age_year[total$patient%in%"patient4"] <- "3"
total$age_year[total$patient%in%"patient5"] <- "8"
total$age_year[total$patient%in%"patient6"] <- "6"


total$tissue <- "NA"
total$tissue[total$orig.ident%in%c("B2","B3","B4","B5")] <- "blood"
total$tissue[total$orig.ident%in%c("L1","L2","L4","L5","N1","N2","N3")] <- "liver"
table(total$tissue)

total$diagnosis <- "NA"
total$diagnosis[total$patient%in%c("patient1","patient3","patient6")] <- "rejection"
total$diagnosis[total$patient%in%c("patient2","patient4","patient5")] <- "non_rejection"
total$diagnosis[total$patient%in%c("patient7","patient8","patient9")] <- "health"

total$status <- "NA"
total$status[total$data%in%c("data1","data2","data3","data4","data5","data6","data7","data8")] <- "transplantation"
total$status[total$data%in%c("data9","data10","data11")] <- "health"



total <- NormalizeData(total, normalization.method = "LogNormalize", scale.factor = 10000)

total <- FindVariableFeatures(total, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(total), 10)

plot1 <- VariableFeaturePlot(total)
LabelPoints(plot = plot1, points = top10, repel = TRUE)


all.genes <- rownames(total)
total <- ScaleData(total, features = all.genes)

total <- RunPCA(total, features = VariableFeatures(object = total))
print(total[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(total, dims = 1:2, reduction = "pca")
DimPlot(total, reduction = "pca")

total <- JackStraw(total, num.replicate = 100)
total <- ScoreJackStraw(total, dims = 1:20)
JackStrawPlot(total, dims = 1:15)
ElbowPlot(total)

total <- FindNeighbors(total, dims = 1:10)
total <- FindClusters(total, resolution = 0.5)

table(total$seurat_clusters)

total <- RunTSNE(total, dims = 1:10)
total <- RunUMAP(total, dims = 1:10)


DimPlot(total, reduction = "tsne",pt.size = 1.6,label = T)
DimPlot(total, reduction = "tsne",pt.size = 1.6,group.by = "orig.ident")

DimPlot(total, reduction = "umap",label = T,pt.size = 1.6)
DimPlot(total, reduction = "umap",pt.size = 1.6,group.by = "orig.ident")



library(SingleR)

load("../package/SingleR/hpca.rda")
refdata <- hpca.se #使用hpca参考数据库鉴定

testdata <- GetAssayData(total, slot="data")
clusters <- total@meta.data$seurat_clusters

cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")

celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)

total@meta.data$celltype_hpca = "NA"
for(i in 1:nrow(celltype)){
  total@meta.data[which(total@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype_hpca'] <- celltype$celltype[i]}

DimPlot(total, group.by="celltype_hpca", repel=T, label=T, reduction='umap',pt.size = 1.6)
DimPlot(total, group.by="celltype_hpca", repel=T, label=T, reduction='tsne',pt.size = 1.6)

table(total$celltype_hpca)

save(total,file = "1-atlas.rda")
load("1-atlas.rda")

###
VlnPlot(total,features = "CD3D",pt.size = 0)## T cell cluster0,1,2,4,10,11,13,18
VlnPlot(total,features = "CD3E",pt.size = 0)## T cell cluster0,1,2,4,10,11,13,18
VlnPlot(total,features = "CD3G",pt.size = 0)## T cell cluster0,1,2,4,10,11,18
VlnPlot(total,features = "IL7R",pt.size = 0)## T cell cluster0,1,2,10,12,14,18

VlnPlot(total,features = "CD7",pt.size = 0)## NK cell cluster0,1,2,4,10,11,14,15,18
VlnPlot(total,features = "FGFBP2",pt.size = 0)## NK cell cluster1,10
VlnPlot(total,features = "KLRF1",pt.size = 0)## NK cell cluster0,4,10
VlnPlot(total,features = "GNLY",pt.size = 0)## NK cell cluster1,2,10,11,14

VlnPlot(total,features = "CD79A",pt.size = 0)## B cell cluster7,9
VlnPlot(total,features = "CD79B",pt.size = 0)## B cell cluster7,9,11,13
VlnPlot(total,features = "MS4A1",pt.size = 0)## B cell cluster7,9
VlnPlot(total,features = "MZB1",pt.size = 0)## B cell cluster7,9,17

VlnPlot(total,features = "CD68",pt.size = 0)## Myeloid cell cluster3,5,6,15,16,17,19
VlnPlot(total,features = "CD14",pt.size = 0)## Myeloid cell cluster3,5,6,13,15
VlnPlot(total,features = "CD163",pt.size = 0)## Myeloid cell cluster3,5,6,15
VlnPlot(total,features = "CSF1R",pt.size = 0)## Myeloid cell cluster3,5,6,15,16
VlnPlot(total,features = "S100A8",pt.size = 0)## Neutrophil cluster3,5,6,12,16,18

VlnPlot(total,features = "ALB",pt.size = 0)##Hepatocyte cluster10,13,14,15
VlnPlot(total,features = "TF",pt.size = 0)##Hepatocyte cluster
VlnPlot(total,features = "TTR",pt.size = 0)##Hepatocyte cluster10,14

VlnPlot(total,features = "PECAM1",pt.size = 0)##Endothelial cell cluster3,5,6,8,12,13,15,16,17,19
VlnPlot(total,features = "CDH5",pt.size = 0)##Endothelial cell cluster8,13,15
VlnPlot(total,features = "ICAM2",pt.size = 0)##Endothelial cell cluster1,3,4,5,6,8,9,11,12,13,15,16,19
VlnPlot(total,features = "KDR",pt.size = 0)##Endothelial cell cluster8,13
VlnPlot(total,features = "ERG",pt.size = 0)##Endothelial cell cluster8,13
VlnPlot(total,features = "VWF",pt.size = 0)##Endothelial cell cluster8
VlnPlot(total,features = "ENG",pt.size = 0)##Endothelial cell cluster8,13,15

VlnPlot(total,features = "COL1A2",pt.size = 0)##Fibroblast cluster12
VlnPlot(total,features = "FAP",pt.size = 0)##Fibroblast
VlnPlot(total,features = "PDPN",pt.size = 0)##Fibroblast
VlnPlot(total,features = "DCN",pt.size = 0)##Fibroblast cluster12
VlnPlot(total,features = "COL3A1",pt.size = 0)##Fibroblast cluster12
VlnPlot(total,features = "COL6A1",pt.size = 0)##Fibroblast cluster8,12
VlnPlot(total,features = "ACTA2",pt.size = 0)##Fibroblast cluster12

VlnPlot(total,features = "CD3D",pt.size = 0)## T cell cluster0,1,2,4,10,11,13,18
VlnPlot(total,features = "CD4",pt.size = 0)##CD4
VlnPlot(total,features = "PTPRC",pt.size = 0)##ILC CD45 cluster0,1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19
VlnPlot(total,features = "GATA3",pt.size = 0)##ILC cluster0,4,11
VlnPlot(total,features = "KLRF1",pt.size = 0)##ILC cluster0,4,10
VlnPlot(total,features = "KLRC1",pt.size = 0)##ILC cluster
VlnPlot(total,features = "GZMA",pt.size = 0)##ILC cluster0,1,2,4,10,11,12,13,14,15,18
VlnPlot(total,features = "GZMB",pt.size = 0)##ILC cluster0,1,2,10,11,17
VlnPlot(total,features = "NKG7",pt.size = 0)##ILC cluster0,1,2,4,5,9,10,11,12,13,14,15,16,17,18,19

##MDSC
VlnPlot(total,features = "PTPRC",pt.size = 0)#CD45
VlnPlot(total,features = "HLA-DRA",pt.size = 0)#HLA-DR
VlnPlot(total,features = "ITGAM",pt.size = 0)#CD11b CD18
VlnPlot(total,features = "CD33",pt.size = 0)#CD33


total$MM <- "NA"
total$MM[total$seurat_clusters%in%"0"] <- "T/NK_cell"
total$MM[total$seurat_clusters%in%"1"] <- "T/NK_cell"
total$MM[total$seurat_clusters%in%"2"] <- "T/NK_cell"
total$MM[total$seurat_clusters%in%"3"] <- "Myeloid_cell"
total$MM[total$seurat_clusters%in%"4"] <- "T/NK_cell"
total$MM[total$seurat_clusters%in%"5"] <- "Myeloid_cell"
total$MM[total$seurat_clusters%in%"6"] <- "Myeloid_cell"
total$MM[total$seurat_clusters%in%"7"] <- "B_cell"
total$MM[total$seurat_clusters%in%"8"] <- "Endothelial_cell"
total$MM[total$seurat_clusters%in%"9"] <- "B_cell"
total$MM[total$seurat_clusters%in%"10"] <- "T/NK_cell"
total$MM[total$seurat_clusters%in%"11"] <- "T/NK_cell"
total$MM[total$seurat_clusters%in%"12"] <- "Fibroblast"
total$MM[total$seurat_clusters%in%"13"] <- "Endothelial_cell"
total$MM[total$seurat_clusters%in%"14"] <- "mix"
total$MM[total$seurat_clusters%in%"15"] <- "Myeloid_cell"
total$MM[total$seurat_clusters%in%"16"] <- "Myeloid_cell"
total$MM[total$seurat_clusters%in%"17"] <- "mix"
total$MM[total$seurat_clusters%in%"18"] <- "T/NK_cell"
total$MM[total$seurat_clusters%in%"19"] <- "mix"
table(total$MM)

DimPlot(total, group.by="MM", repel=T, reduction='tsne',pt.size = 1.6)
DimPlot(total, group.by="MM", repel=T, reduction='umap',pt.size = 1.6)


total.markers <- FindAllMarkers(total, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
total.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top5 <- total.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(total, features = top5$gene) #+ NoLegend()
save(total.markers,file = "1-totalmarkers.rda")
load("1-total_select.markers.rda")
write.csv(total_select.markers,file = "Table2_total_select_markers.csv")

table(total$MM)
total_select <- subset(total,subset=MM!="mix")
DimPlot(total_select, group.by="MM", repel=T, reduction='tsne',pt.size = 1.6)
DimPlot(total_select, group.by="MM", repel=T, reduction='umap',pt.size = 1.6)


save(total_select,file = "1-atlas2.rda")
load("1-atlas2.rda")
DimPlot(total_select, repel=T, reduction='tsne',pt.size = 1.6)
DimPlot(total_select, group.by="MM", repel=T, reduction='tsne',pt.size = 1.6)
DiscretePalette(5, palette = NULL)
DimPlot(total_select, group.by="MM", repel=T, reduction='tsne',pt.size = 1.6,
        cols = c("#F0A0FF", "#0075DC", "#993F00", "#93669d" ,"#00733d"))
table(total_select$diagnosis)
DimPlot(total_select, group.by="diagnosis", repel=T, reduction='tsne',pt.size = 1.6)


###各样本的细胞比例
pred <- total_select@meta.data
pred <- pred[,c("orig.ident","MM")]

sort(table(pred$MM))
pred$value <- length(pred$MM[pred$MM%in%"Fibroblast"])
pred$value[pred$MM%in%"Endothelial_cell"] <- length(pred$MM[pred$MM%in%"Endothelial_cell"])
pred$value[pred$MM%in%"B_cell"] <- length(pred$MM[pred$MM%in%"B_cell"])
pred$value[pred$MM%in%"Myeloid_cell"] <- length(pred$MM[pred$MM%in%"Myeloid_cell"])
pred$value[pred$MM%in%"T/NK_cell"] <- length(pred$MM[pred$MM%in%"T/NK_cell"])
sort(table(pred$value))

##画图
library(ggplot2)
ggplot(pred, aes(fill=MM, y=value, x=orig.ident),order_by(pred)) + 
  geom_bar(position="fill", stat="identity")+theme_bw(base_size = 20)

table(total_select$MM)
if(T){
  B_cell.markers <- FindMarkers(total_select,ident.1 = "B_cell", min.pct = 0.25,group.by = "MM")
  B_cell.markers$cluster <- "B_cell"
  
  Endothelial_cell.markers <- FindMarkers(total_select,ident.1 = "Endothelial_cell", min.pct = 0.25,group.by = "MM")
  Endothelial_cell.markers$cluster <- "Endothelial_cell"
  
  Fibroblast.markers <- FindMarkers(total_select,ident.1 = "Fibroblast", min.pct = 0.25,group.by = "MM")
  Fibroblast.markers$cluster <- "Fibroblast"
  
  Myeloid_cell.markers <- FindMarkers(total_select,ident.1 = "Myeloid_cell", min.pct = 0.25,group.by = "MM")
  Myeloid_cell.markers$cluster <- "Myeloid_cell"
  
  T_NK_cell.markers <- FindMarkers(total_select,ident.1 = "T/NK_cell", min.pct = 0.25,group.by = "MM")
  T_NK_cell.markers$cluster <- "T/NK_cell"
  
  total_select.markers <- rbind(B_cell.markers,Endothelial_cell.markers,Fibroblast.markers,
                          Myeloid_cell.markers,T_NK_cell.markers)
  total_select.markers$gene <- rownames(total_select.markers)
  
  top5 <- total_select.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
  #DoHeatmap(total_select, features = top5$gene,group.by = "MM",slot = "counts")
  #DoHeatmap(total_select, features = top5$gene,group.by = "MM",slot = "scale.data")
  save(total_select.markers,file = "1-total_select.markers.rda")
  
}
load("1-total_select.markers.rda")
top5 <- total_select.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
top10 <- total_select.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

top10_genes <- top10$gene
pdf("top10.pdf",width = 8,height = 12)
top10_genes <- c("IGKC"    ,  "IGHM"   ,   "IGHG1"  ,   "CD79A"   ,  "IGLC2"   ,  "IGLC3"    , "IGHG3"  ,  
                 "MS4A1"   ,  "IGLC1"  ,   "IGHA1"  ,   "CCL14"   ,  "FCN3"   ,   "DNASE1L3" , "CRHBP"  ,  
                 "IFI27"   ,  "TM4SF1"  ,  "MGP"    ,   "SLC9A3R2" , "PLPP3"  ,   "HES1"    ,  "TAGLN"   , 
                 "MYL9"    ,  "ADIRF"  ,  "ACTA2"  ,   "IGFBP7"  , "MYH11"  ,   "TPM2"   ,   "CALD1"   ,
                 "C11orf96" ,"IGFBP5"  ,  "LYZ"   ,   "S100A9" ,  "CTSS"  ,   "S100A8" ,  "CST3"    ,
                 "FCN1"   ,   "VCAN"    ,  "IFI30" ,   "AIF1"   ,  "MNDA"   ,   "CCL5"   ,  "NKG7"   , 
                 "GNLY"  ,   "IL32"  ,   "KLRB1"  ,  "GZMA"   ,  "CST7"  ,   "IL7R"   ,  "CD3D"  ,  
                 "KLRD1" )
DoHeatmap(total_select, features = top10_genes,group.by = "MM")
total_select2 <- total_select
total_select2@active.ident <- total_select2$diagnosis
DoHeatmap(subset(total_select,downsample=100), features = top10_genes)
DoHeatmap(subset(total_select,downsample=100), features = top10_genes,group.by = "MM")

#DoHeatmap(total_select, features = top10$gene,group.by = "MM",slot = "counts")
dev.off()

top5_genes <- top5$gene
top5_genes <- c("IGKC"   ,  "IGHM"  ,   "IGHG1"   , "CD79A" ,   "IGLC2"   , "CCL14"   , "FCN3"    , "DNASE1L3",
                "CRHBP"   , "IFI27"  ,  "TAGLN"    ,"MYL9"   ,  "ADIRF" ,  "ACTA2"   , "IGFBP7" , "LYZ"    ,
                "S100A9" , "CTSS"   , "S100A8" , "CST3"   , "CCL5" ,   "NKG7"   , "GNLY"   , "IL32"   ,
                "KLRB1"  )
DotPlot(total_select,features = top5_genes,group.by = "MM")

DoHeatmap(total_select, features = top5$gene,group.by = "MM",slot = "counts")
DoHeatmap(total, features = top5$gene,group.by = "MM")

DoHeatmap(total_select, features = top5$gene,group.by = "MM",slot = "data")

dev.off()
DotPlot(total_select, features = top5$gene,group.by = "MM")


VlnPlot(total_select,features = "CD79A",group.by = "MM",pt.size = 0)
VlnPlot(total_select,features = "CDH5",group.by = "MM",pt.size = 0)
VlnPlot(total_select,features = "COL1A2",group.by = "MM",pt.size = 0)
VlnPlot(total_select,features = "CD68",group.by = "MM",pt.size = 0)
VlnPlot(total_select,features = "CD3D",group.by = "MM",pt.size = 0)
FeaturePlot(total_select,features = "CD79A",reduction = "tsne",cols = c("grey","blue"))
FeaturePlot(total_select,features = "CDH5",reduction = "tsne")
FeaturePlot(total_select,features = "COL1A2",reduction = "tsne")
FeaturePlot(total_select,features = "CD68",reduction = "tsne")
FeaturePlot(total_select,features = "CD3D",reduction = "tsne")
FeaturePlot(total_select,features = "KLRF1",reduction = "tsne")

T_NK_cell <- subset(total_select,subset = MM%in%"T/NK_cell");save(T_NK_cell,file = "T_NK_cell.rda")
Myeloid_cell <- subset(total_select,subset = MM%in%"Myeloid_cell");save(Myeloid_cell,file = "Myeloid_cell.rda")
B_cell <- subset(total_select,subset = MM%in%"B_cell");save(B_cell,file = "B_cell.rda")


load("1-atlas.rda")
VlnPlot(total,features = c("CD4","CD8A","PTGDR2","GATA3","KLRB1","IL13"),pt.size = 0)
