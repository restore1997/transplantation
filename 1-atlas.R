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

total1 <- merge(data1,y=c(data2,data3,data4),project = "total1")

total1
table(total1$orig.ident)
total1[["percent.mt"]] <- PercentageFeatureSet(total1, pattern = "^MT-")
VlnPlot(total1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)

plot1 <- FeatureScatter(total1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(total1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

total1
total1 <- subset(total1, subset = nFeature_RNA > 500  & nCount_RNA >2000 & percent.mt < 25)

total1
table(total1$orig.ident)
VlnPlot(total1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)


###add healthy bloods

load("../../单细胞/外周血单细胞/2020-A single-cell atlas of severe COVID-19/raw_total.rda")
total2 <- total

total <- merge(total1,total2)
table(total$orig.ident)
VlnPlot(total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)


table(total$orig.ident)
total <- subset(total,subset = orig.ident%in%c("B2","B3","B4","B5",
                                               "health2","health3","health4","health5"))

total
table(total$orig.ident)
VlnPlot(total, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)

save(total,file = "raw_total.rda")

###meta data
total$data <- total$orig.ident

total$patient <- "NA"
total$patient[total$data%in%"B2"] <- "patient1"
total$patient[total$data%in%"B3"] <- "patient2"
total$patient[total$data%in%"B4"] <- "patient3"
total$patient[total$data%in%"B5"] <- "patient4"
total$patient[total$data%in%"L1"] <- "patient1"
total$patient[total$data%in%"L2"] <- "patient4"
total$patient[total$data%in%"L4"] <- "patient5"
total$patient[total$data%in%"L5"] <- "patient6"

total$pre_diagnosis <- "NA"
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
total$tissue[total$orig.ident%in%c("L1","L2","L4","L5")] <- "liver"

total$diagnosis <- "NA"
total$diagnosis[total$patient%in%c("patient1","patient3","patient6")] <- "rejection"
total$diagnosis[total$patient%in%c("patient2","patient4","patient5")] <- "non_rejection"


total <- NormalizeData(total, normalization.method = "LogNormalize", scale.factor = 10000)

total <- FindVariableFeatures(total, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(total), 10)

plot1 <- VariableFeaturePlot(total)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

total <- ScaleData(total)

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

load("../../data/SingleR/hpca.rda")
refdata <- ref.data #使用hpca参考数据库鉴定

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
VlnPlot(total,features = "CD3D",pt.size = 0)## T cell cluster0,1,3,5,10,11,12,15,16
VlnPlot(total,features = "CD3E",pt.size = 0)## T cell cluster0,1,3,5,8,10,11,12,15,16,17
VlnPlot(total,features = "CD3G",pt.size = 0)## T cell cluster0,1,3,5,10,11,12,16
VlnPlot(total,features = "IL7R",pt.size = 0)## T cell cluster0,1,5,11,13,14,15,16,18

VlnPlot(total,features = "CD7",pt.size = 0)## NK cell cluster0,1,3,5,8,10,11,12,13,14,15,16,17,18
VlnPlot(total,features = "FGFBP2",pt.size = 0)## NK cell cluster3,5,8,11,16
VlnPlot(total,features = "KLRF1",pt.size = 0)## NK cell cluster3,8,11
VlnPlot(total,features = "GNLY",pt.size = 0)## NK cell cluster3,5,8,11,12,13,14,15,16,17,18

VlnPlot(total,features = "CD79A",pt.size = 0)## B cell cluster6,15,19,22
VlnPlot(total,features = "CD79B",pt.size = 0)## B cell cluster6,9,12,16,22
VlnPlot(total,features = "MS4A1",pt.size = 0)## B cell cluster6,15,22
VlnPlot(total,features = "MZB1",pt.size = 0)## B cell cluster6,15,17,19

VlnPlot(total,features = "CD68",pt.size = 0)## Myeloid cell cluster2,4,7,9,16,17,19,22,23
VlnPlot(total,features = "CD14",pt.size = 0)## Myeloid cell cluster2,4,7,9,16,22
VlnPlot(total,features = "CD163",pt.size = 0)## Myeloid cell cluster2,4,9,16,22
VlnPlot(total,features = "CSF1R",pt.size = 0)## Myeloid cell cluster2,4,7,9,16,22
VlnPlot(total,features = "S100A8",pt.size = 0)## Neutrophil cluster2,4,7,9,16,19,21,22

VlnPlot(total,features = "ALB",pt.size = 0)##Hepatocyte cluster11,13,14,18,20
VlnPlot(total,features = "TF",pt.size = 0)##Hepatocyte cluster20
VlnPlot(total,features = "TTR",pt.size = 0)##Hepatocyte cluster11,13,20

VlnPlot(total,features = "PECAM1",pt.size = 0)##Endothelial cell cluster2,4,7,9,14,16,17,19,22,23
VlnPlot(total,features = "CDH5",pt.size = 0)##Endothelial cell cluster14
VlnPlot(total,features = "ICAM2",pt.size = 0)##Endothelial cell cluster0,4,7,8,9,12,14,15,16,22,23
VlnPlot(total,features = "KDR",pt.size = 0)##Endothelial cell cluster14
VlnPlot(total,features = "ERG",pt.size = 0)##Endothelial cell cluster14
VlnPlot(total,features = "VWF",pt.size = 0)##Endothelial cell cluster14
VlnPlot(total,features = "ENG",pt.size = 0)##Endothelial cell cluster14

VlnPlot(total,features = "COL1A2",pt.size = 0)##Fibroblast cluster18
VlnPlot(total,features = "FAP",pt.size = 0)##Fibroblast
VlnPlot(total,features = "PDPN",pt.size = 0)##Fibroblast
VlnPlot(total,features = "DCN",pt.size = 0)##Fibroblast cluster18
VlnPlot(total,features = "COL3A1",pt.size = 0)##Fibroblast cluster18
VlnPlot(total,features = "COL6A1",pt.size = 0)##Fibroblast cluster18
VlnPlot(total,features = "ACTA2",pt.size = 0)##Fibroblast cluster18

##MDSC
VlnPlot(total,features = "PTPRC",pt.size = 0)#CD45
VlnPlot(total,features = "HLA-DRA",pt.size = 0)#HLA-DR
VlnPlot(total,features = "ITGAM",pt.size = 0)#CD11b CD18
VlnPlot(total,features = "CD33",pt.size = 0)#CD33


total$MM <- "NA"
total$MM[total$seurat_clusters%in%"0"] <- "T/NK_cell"
total$MM[total$seurat_clusters%in%"1"] <- "T/NK_cell"
total$MM[total$seurat_clusters%in%"2"] <- "Myeloid_cell"
total$MM[total$seurat_clusters%in%"3"] <- "T/NK_cell"
total$MM[total$seurat_clusters%in%"4"] <- "Myeloid_cell"
total$MM[total$seurat_clusters%in%"5"] <- "T/NK_cell"
total$MM[total$seurat_clusters%in%"6"] <- "B_cell"
total$MM[total$seurat_clusters%in%"7"] <- "Myeloid_cell"
total$MM[total$seurat_clusters%in%"8"] <- "T/NK_cell"
total$MM[total$seurat_clusters%in%"9"] <- "Myeloid_cell"
total$MM[total$seurat_clusters%in%"10"] <- "T/NK_cell"
total$MM[total$seurat_clusters%in%"11"] <- "T/NK_cell"
total$MM[total$seurat_clusters%in%"12"] <- "T/NK_cell"
total$MM[total$seurat_clusters%in%"13"] <- "T/NK_cell"
total$MM[total$seurat_clusters%in%"14"] <- "Endothelial_cell"
total$MM[total$seurat_clusters%in%"15"] <- "mix"
total$MM[total$seurat_clusters%in%"16"] <- "mix"
total$MM[total$seurat_clusters%in%"17"] <- "mix"
total$MM[total$seurat_clusters%in%"18"] <- "Fibroblast"
total$MM[total$seurat_clusters%in%"19"] <- "mix"
total$MM[total$seurat_clusters%in%"20"] <- "Hepatocyte"
total$MM[total$seurat_clusters%in%"21"] <- "Neutrophil"
total$MM[total$seurat_clusters%in%"22"] <- "mix"
total$MM[total$seurat_clusters%in%"23"] <- "mix"
table(total$MM)

DimPlot(total, group.by="MM", repel=T, reduction='tsne',pt.size = 1.6)
DimPlot(total, group.by="MM", repel=T, reduction='umap',pt.size = 1.6)


total.markers <- FindAllMarkers(total, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
total.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top5 <- total.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(total, features = top5$gene) #+ NoLegend()
save(total.markers,file = "1-totalmarkers.rda")
load("1-total_select.markers.rda")

table(total$MM)
total_select <- subset(total,subset=MM!="mix")
DimPlot(total_select, group.by="MM", repel=T, reduction='tsne',pt.size = 1.6)
DimPlot(total_select, group.by="MM", repel=T, reduction='umap',pt.size = 1.6)


save(total_select,file = "1-atlas2.rda")
load("1-atlas2.rda")


###各样本的细胞比例
pred <- total_select@meta.data
pred <- pred[,c("orig.ident","MM")]

sort(table(pred$MM))
pred$value <- length(pred$MM[pred$MM%in%"Neutrophil"])
pred$value[pred$MM%in%"Hepatocyte"] <- length(pred$MM[pred$MM%in%"Hepatocyte"])
pred$value[pred$MM%in%"Fibroblast"] <- length(pred$MM[pred$MM%in%"Fibroblast"])
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
if(F){
  B_cell.markers <- FindMarkers(total_select,ident.1 = "B_cell", min.pct = 0.25,group.by = "MM")
  B_cell.markers$cluster <- "B_cell"
  
  Endothelial_cell.markers <- FindMarkers(total_select,ident.1 = "Endothelial_cell", min.pct = 0.25,group.by = "MM")
  Endothelial_cell.markers$cluster <- "Endothelial_cell"
  
  Fibroblast.markers <- FindMarkers(total_select,ident.1 = "Fibroblast", min.pct = 0.25,group.by = "MM")
  Fibroblast.markers$cluster <- "Fibroblast"
  
  Hepatocyte.markers <- FindMarkers(total_select,ident.1 = "Hepatocyte", min.pct = 0.25,group.by = "MM")
  Hepatocyte.markers$cluster <- "Hepatocyte"
  
  Myeloid_cell.markers <- FindMarkers(total_select,ident.1 = "Myeloid_cell", min.pct = 0.25,group.by = "MM")
  Myeloid_cell.markers$cluster <- "Myeloid_cell"
  
  Neutrophil.markers <- FindMarkers(total_select,ident.1 = "Neutrophil", min.pct = 0.25,group.by = "MM")
  Neutrophil.markers$cluster <- "Neutrophil"
  
  T_NK_cell.markers <- FindMarkers(total_select,ident.1 = "T/NK_cell", min.pct = 0.25,group.by = "MM")
  T_NK_cell.markers$cluster <- "T/NK_cell"
  
  total_select.markers <- rbind(B_cell.markers,Endothelial_cell.markers,Fibroblast.markers,Hepatocyte.markers,
                          Myeloid_cell.markers,Neutrophil.markers,T_NK_cell.markers)
  total_select.markers$gene <- rownames(total_select.markers)
  
  top5 <- total_select.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
  DoHeatmap(total_select, features = top5$gene,group.by = "MM",slot = "counts")
  #DoHeatmap(total_select, features = top5$gene,group.by = "MM",slot = "scale.data")
  save(total_select.markers,file = "1-total_select.markers.rda")
  
}

VlnPlot(total_select,features = "CD79A",group.by = "MM",pt.size = 0)
VlnPlot(total_select,features = "CDH5",group.by = "MM",pt.size = 0)
VlnPlot(total_select,features = "COL1A2",group.by = "MM",pt.size = 0)
VlnPlot(total_select,features = "TF",group.by = "MM",pt.size = 0)
VlnPlot(total_select,features = "CD68",group.by = "MM",pt.size = 0)
VlnPlot(total_select,features = "S100A8",group.by = "MM",pt.size = 0)
VlnPlot(total_select,features = "CD3D",group.by = "MM",pt.size = 0)


T_NK_cell <- subset(total_select,subset = MM%in%"T/NK_cell");save(T_NK_cell,file = "T_NK_cell.rda")
Myeloid_cell <- subset(total_select,subset = MM%in%"Myeloid_cell");save(Myeloid_cell,file = "Myeloid_cell.rda")
