rm(list = ls())
library(Seurat)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

load("../T_NK_cell/liver/C10/1-atlas.rda")
table(liver$status)
T_NK_cell <- subset(liver,subset=status%in%"transplantation")
table(T_NK_cell$seurat_clusters)
T_NK_cell$MM2[T_NK_cell$seurat_clusters%in%"0"] <- "CD8T_0"
T_NK_cell$MM2[T_NK_cell$seurat_clusters%in%"1"] <- "CD8T_1"
T_NK_cell$MM2[T_NK_cell$seurat_clusters%in%"2"] <- "CD8T_2"

load("../Myeloid_cell/liver/1-atlas.rda")
table(liver$status)
Myeloid_cell <- subset(liver,subset=status%in%"transplantation")
table(Myeloid_cell$MM2)

T_NK_cell
Myeloid_cell

total <- merge(T_NK_cell,Myeloid_cell)
table(total$MM)
table(total$MM2)

rm(list = setdiff(ls(), "total"))

data.input <- GetAssayData(total, assay = "RNA", slot = "data") # normalized data matrix
meta <- total@meta.data

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "MM2")
rm(list = setdiff(ls(), "cellchat"))

#> Create a CellChat object from a data matrix
#> Set cell identities for the new CellChat object
#> The cell groups used for CellChat analysis are  APOE+ FIB FBN1+ FIB COL11A1+ FIB Inflam. FIB cDC1 cDC2 LC Inflam. DC TC Inflam. TC CD40LG+ TC NKT

cellchat <- setIdent(cellchat, ident.use = "MM2") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
table(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
CellChatDB_interaction <- CellChatDB$interaction
table(CellChatDB_interaction$annotation)

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multiprocess", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)


##Part II: Inference of cell-cell communication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)
sort(table(df.net$pathway_name))
sort(table(df.net$ligand))

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
sort(table(df.net$pathway_name))
save(cellchat,file = "1-cellchat.rda")



