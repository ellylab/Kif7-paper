options(stringsAsFactors = F)
# options(warn = -1, verbose = F)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(monocle)
library(Seurat)
library(cowplot)
library(RColorBrewer)
library(pheatmap)

source("https://raw.githubusercontent.com/leezx/Toolsets/master/R/Toolsets.R")
source("https://raw.githubusercontent.com/leezx/Toolsets/master/R/Plot.R")

# myColors <- brewer.pal(8,"Set2")
# myColors

# lineage
myColors5 <- c('#FDB462','#B3DE69','#FCCDE5',"#ff309a",'#D9D9D9')
# names(myColors5) <- c("GP", "BP", "NP", "late NP", "MF")
# names(myColors5) <- c("GP", "BP", "NPearly","NPlate","ENMFB")
names(myColors5) <- c("GP","BP","NPearly","NPlate","ENMFB")

# data1
print(load("data/organogenesis/ENS_cell_annotate.Rdata"))
organMatrix.anno <- ENS_cell_annotate
dim(organMatrix.anno)
cds <- readRDS("data/organogenesis/cds_cleaned.RDS")
organMatrix <- exprs(cds)[,ENS_cell_annotate$sample]
organMatrix <- as.matrix(organMatrix)
dim(organMatrix)

gene_short_name <- cds@featureData@data$gene_short_name
gene_short_name2 <- gene_short_name[!duplicated(gene_short_name)]
organMatrix <- organMatrix[!duplicated(gene_short_name),]
rownames(organMatrix) <- gene_short_name2
dim(organMatrix)

# data2
print(load("data/mouse_brain/Mouse_brain.ENS.Rdata"))
dim(mMatrix)
mMatrix.anno <- read.csv("data/mouse_brain/metaDf.csv", header = T, row.names = 1)
dim(mMatrix.anno)
colnames(mMatrix)[1:5]
mMatrix.anno$CellID[1:5]
colnames(mMatrix)[duplicated(colnames(mMatrix))]
mMatrix.anno$CellID[duplicated(mMatrix.anno$CellID)]
mMatrix.anno$CellID2 <- colnames(mMatrix)
rownames(mMatrix.anno) <- mMatrix.anno$CellID2
mMatrix.anno2 <- read.csv("data/mouse_brain/cellID.clusterName.csv", header = T, row.names = 1)
dim(mMatrix.anno2)
mMatrix.anno2 <- mMatrix.anno2[!duplicated(mMatrix.anno2$CellID),]
rownames(mMatrix.anno2) <- mMatrix.anno2$CellID
head(mMatrix.anno2)
mMatrix.anno$ClusterName <- mMatrix.anno2[mMatrix.anno$CellID,]$ClusterName
print(dim(mMatrix))
print(dim(mMatrix.anno))
mMatrix.anno$ClusterName <- as.character(mMatrix.anno$ClusterName)
# exclude irrelevant cells
mMatrix.anno <- subset(mMatrix.anno, Class %in% c("Neurons", "PeripheralGlia", "Vascular"))

# prepare data
print(load("data/organogenesis/Mouse_organogenesis.ENS.Rdata"))
dim(organMatrix)
dim(organMatrix.anno)
organMatrix.anno$sample <- as.character(organMatrix.anno$sample)
organMatrix.anno$Main_cell_type <- as.character(organMatrix.anno$Main_cell_type)
# creat standard format df
organMatrix.anno.st <- data.frame(cellName=organMatrix.anno$sample, origin.cluster=organMatrix.anno$Main_cell_type, 
                              stage=organMatrix.anno$development_stage, platform="sci-RNA-seq3")
organMatrix.anno.st$stage <- paste("E", organMatrix.anno.st$stage, sep="")
organMatrix <- organMatrix[,organMatrix.anno.st$cellName]

print(load("data/mouse_brain/Mouse_brain.ENS.Rdata"))
mMatrix <- mMatrix[,rownames(mMatrix.anno)]
dim(mMatrix)
dim(mMatrix.anno)
#sum(table(mMatrix.anno$ClusterName))
mMatrix.anno <- mMatrix.anno[!is.na(mMatrix.anno$ClusterName),]
mMatrix.anno$ClusterName <- as.character(mMatrix.anno$ClusterName)
dim(mMatrix.anno)
mMatrix <- mMatrix[,rownames(mMatrix.anno)]
dim(mMatrix)
mMatrix.anno.st <- data.frame(cellName=rownames(mMatrix.anno), origin.cluster=mMatrix.anno$ClusterName, stage=mMatrix.anno$Age, 
                              platform="10x Genomics")
mMatrix.anno.st[grep("p19", mMatrix.anno.st$stage),]$stage <- "P19"
mMatrix.anno.st[grep("p23", mMatrix.anno.st$stage),]$stage <- "P23"

# Elly lab data
Ctrl_YFP_ENCC_report <- "/Users/zxli/Dropbox/Projects/EllyLab/mouse/singleCell/control/ENCC/raw/Ctrl_YFP_ENCC_report/outs/filtered_gene_bc_matrices/mm10"
Ctrl_E165_YFP_ENCC_report <- "/Users/zxli/Dropbox/Projects/EllyLab/mouse/singleCell/control/ENCC/raw/Control-E165-YFP-ENCC/outs/filtered_gene_bc_matrices/mm10/"
Kif7_YFP_ENCC_report <- "/Users/zxli/Dropbox/Projects/EllyLab/mouse/singleCell/case/Kif7_ENCC/raw/Kif7l-YFP-ENCC/outs/filtered_gene_bc_matrices/mm10"
GBS_GFPneg_ENCC_report <- "/Users/zxli/Dropbox/Projects/EllyLab/mouse/singleCell/case/GBS_ENCC/raw/GBS_GFPneg_ENCC_report/outs/filtered_gene_bc_matrices/mm10"
GBS_GFPpos_ENCC_report <- "/Users/zxli/Dropbox/Projects/EllyLab/mouse/singleCell/case/GBS_ENCC/raw/GBS_GFPpos_ENCC_report/outs/filtered_gene_bc_matrices/mm10"
Vcl_YFP_ENCC_report <- "/Users/zxli/Dropbox/Projects/EllyLab/mouse/singleCell/case/Vcl_ENCC/raw/Vcl_YFP_ENCC_report/outs/filtered_gene_bc_matrices/mm10"

ctrl_rawdata <- Read10X(Ctrl_YFP_ENCC_report)
colnames(ctrl_rawdata) <- paste("ctrl", colnames(ctrl_rawdata), sep="_")
print(dim(ctrl_rawdata))
ctrl_rawdata[1:3,1:3]

ctr2_rawdata <- Read10X(Ctrl_E165_YFP_ENCC_report)
colnames(ctr2_rawdata) <- paste("ctr2", colnames(ctr2_rawdata), sep="_")
print(dim(ctr2_rawdata))
ctr2_rawdata[1:3,1:3]

kif7_rawdata <- Read10X(Kif7_YFP_ENCC_report)
colnames(kif7_rawdata) <- paste("kif7", colnames(kif7_rawdata), sep="_")
print(dim(kif7_rawdata))
kif7_rawdata[1:3,1:3]

vcl_rawdata <- Read10X(Vcl_YFP_ENCC_report)
colnames(vcl_rawdata) <- paste("vcl", colnames(vcl_rawdata), sep="_")
print(dim(vcl_rawdata))
vcl_rawdata[1:3,1:3]

GBS_GFPpos_rawdata <- Read10X(GBS_GFPpos_ENCC_report)
colnames(GBS_GFPpos_rawdata) <- paste("pos", colnames(GBS_GFPpos_rawdata), sep="_")
print(dim(GBS_GFPpos_rawdata))
GBS_GFPpos_rawdata <- GBS_GFPpos_rawdata[,GBS_GFPpos_rawdata["Phox2b",]>0]
GBS_GFPpos_rawdata[1:3,1:3]
print(dim(GBS_GFPpos_rawdata))

GBS_GFPneg_rawdata <- Read10X(GBS_GFPneg_ENCC_report)
colnames(GBS_GFPneg_rawdata) <- paste("neg", colnames(GBS_GFPneg_rawdata), sep="_")
print(dim(GBS_GFPneg_rawdata))
GBS_GFPneg_rawdata <- GBS_GFPneg_rawdata[,GBS_GFPneg_rawdata["Phox2b",]>0]
GBS_GFPneg_rawdata[1:3,1:3]
print(dim(GBS_GFPneg_rawdata))

E135Matrix.anno.st <- data.frame(cellName=colnames(ctrl_rawdata), origin.cluster="control", 
                              stage="E13.5", platform="10x Genomics")

head(E135Matrix.anno.st)

E135Matrix <- ctrl_rawdata[,E135Matrix.anno.st$cellName]

# print(load("ellylab/E165_ENCC_2.Rdata"))

# seuset@raw.data[1:5,1:5]

E165Matrix.anno.st <- data.frame(cellName=colnames(ctr2_rawdata), origin.cluster="control", 
                              stage="E16.5", platform="10x Genomics")

head(E165Matrix.anno.st)

E165Matrix <- ctr2_rawdata[,E165Matrix.anno.st$cellName]

# E165Matrix <- as.matrix(seuset@raw.data)

# head(all_tsne)

kif7Matrix.anno.st <- data.frame(cellName=colnames(kif7_rawdata), origin.cluster="Kif7 cKO", 
                              stage="E13.5", platform="10x Genomics")

head(kif7Matrix.anno.st)

rownames(kif7Matrix.anno.st) <- kif7Matrix.anno.st$cellName

set.seed(49)
kif7Matrix.anno.st <- kif7Matrix.anno.st[sample(kif7Matrix.anno.st$cellName, 5000),]

kif7Matrix <- kif7_rawdata[,kif7Matrix.anno.st$cellName]

dim(kif7Matrix)

# merge
metadata <- rbind(mMatrix.anno.st, organMatrix.anno.st, E135Matrix.anno.st, E165Matrix.anno.st, 
                  kif7Matrix.anno.st, vclMatrix.anno.st, negMatrix.anno.st, posMatrix.anno.st)
dim(metadata)

all.count <- cbind(mMatrix[genes,], organMatrix[genes,], E135Matrix[genes,], E165Matrix[genes,], 
                   kif7Matrix[genes,], vclMatrix[genes,], negMatrix[genes,], posMatrix[genes,])

metadata <- metadata[!duplicated(metadata$cellName),]
dim(metadata)



# Integrate rawData using seurat
print(load("all.raw.integration.Rdata"))

table(metadata$stage)

table(metadata$platform)

# this cells are wired, remove them
metadata <- metadata[!(metadata$stage=="E13.5" & metadata$platform=="sci-RNA-seq3"),]

all.count <- all.count[,rownames(metadata)]

ens.integration <- CreateSeuratObject(counts = all.count, meta.data = metadata)
ens.list <- SplitObject(object = ens.integration, split.by = "platform")
# ens.list <- SplitObject(object = ens.integration, split.by = "stage")

names(ens.list)

options(warn = -1)

for (i in 1:length(ens.list)) {
    ens.list[[i]] <- SCTransform(ens.list[[i]], verbose = FALSE)
}

# save(ens.list, file = "all.ens.list.Rdata")

print(load("ens.list.Rdata"))

options(future.globals.maxSize= 5000*1024^2)

ens.features <- SelectIntegrationFeatures(object.list = ens.list, nfeatures = 10000)

print(load("cor.genes.Rdata"))

ens.features <- ens.features[ens.features %in% cor.genes]
print(length(ens.features))

# save(ens.features, file="ens.features.Rdata")

ens.list <- PrepSCTIntegration(object.list = ens.list, anchor.features = ens.features, verbose = FALSE)

ens.anchors <- FindIntegrationAnchors(object.list = ens.list, normalization.method = "SCT", 
    anchor.features = ens.features, verbose = FALSE)
ens.integrated <- IntegrateData(anchorset = ens.anchors, normalization.method = "SCT", verbose = FALSE)

# switch to integrated assay. The variable features of this assay are
# automatically set during IntegrateData
DefaultAssay(object = ens.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
ens.integrated <- ScaleData(object = ens.integrated, verbose = FALSE)
ens.integrated <- RunPCA(object = ens.integrated, npcs = 30, verbose = FALSE)
ens.integrated <- RunUMAP(object = ens.integrated, reduction = "pca", dims = 1:30)

ens.integrated <- RunTSNE(object = ens.integrated, reduction = "pca", dims = 1:30)

ens.integrated <- FindNeighbors(ens.integrated, reduction = "pca", dims = 1:20)

# 0.2 -> 10 clusters
# 0.3 -> 11 clusters
# 0.4 -> 14 clusters
ens.integrated <- FindClusters(ens.integrated, resolution = 0.15)
table(ens.integrated$seurat_clusters)

table(ens.integrated$origin.cluster)

# options(repr.plot.width=6, repr.plot.height=4)
# p1 <- DimPlot(object = ens.integrated, reduction = "umap", group.by = "platform")
# p1

# options(repr.plot.width=6, repr.plot.height=4)
# p1 <- DimPlot(object = ens.integrated, reduction = "umap", group.by = "stage")
# p1

# options(repr.plot.width=5, repr.plot.height=4)
# p2 <- DimPlot(object = ens.integrated, reduction = "umap", group.by = "origin.cluster", 
#     label = TRUE, repel = TRUE) + NoLegend()
# p2

# plot_grid(p1, p2)

# options(repr.plot.width=5, repr.plot.height=4)
# p3 <- DimPlot(object = ens.integrated, reduction = "umap", group.by = "seurat_clusters", 
#     label = TRUE, repel = TRUE) + NoLegend()
# p3

# options(repr.plot.width=5, repr.plot.height=4)
# p3 <- DimPlot(object = ens.integrated, reduction = "pca", group.by = "seurat_clusters", 
#     label = TRUE, repel = TRUE) + NoLegend()
# p3

# options(repr.plot.width=5, repr.plot.height=4)
# p3 <- DimPlot(object = ens.integrated, reduction = "pca", group.by = "origin.cluster", 
#     label = TRUE, repel = TRUE) + NoLegend()
# p3

table(ens.integrated$platform)

table(ens.integrated$stage)

# options(repr.plot.width=5, repr.plot.height=4)
# DimPlot(object = ens.integrated, reduction = "tsne", group.by = "origin.cluster", 
#     label = TRUE, repel = TRUE) + NoLegend()

# options(repr.plot.width=5, repr.plot.height=4)
# DimPlot(object = ens.integrated, reduction = "tsne", group.by = "seurat_clusters", 
#     label = TRUE, repel = TRUE) + NoLegend()

table(ens.integrated$origin.cluster)

unique(ens.integrated$origin.cluster)

# ens.integrated$origin.cluster2 <- factor(ens.integrated$origin.cluster, 
#                                         levels = rev(c('Sensory neurons', 'Schwann cell precursor', '',
#                                             'ENT9', 'ENT8', 'ENT7', 'ENT6', 'ENT4', 'ENT5', 'ENTG7',
#                                                    'ENTG5', 'ENTG2', 'ENT3', 'ENT2', 'ENTG4', 'ENT1', 'ENTG6',
#                                                    'ENTG3', 'ENTG1'))
#                                        )

# options(repr.plot.width=5, repr.plot.height=4)
# DimPlot(object = ens.integrated, reduction = "pca", group.by = "origin.cluster2", 
#     label = TRUE, repel = TRUE) + NoLegend()
