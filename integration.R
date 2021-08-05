

# integration

ens.integration <- CreateSeuratObject(counts = all.count, meta.data = metadata)
# ens.list <- SplitObject(object = ens.integration, split.by = "platform")
ens.list <- SplitObject(object = ens.integration, split.by = "stage")

for (i in 1:length(x = ens.list)) {
    ens.list[[i]] <- NormalizeData(object = ens.list[[i]], verbose = FALSE)
    ens.list[[i]] <- FindVariableFeatures(object = ens.list[[i]], 
        selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

ens.anchors <- FindIntegrationAnchors(object.list = ens.list, dims = 1:30, max.features = 1000)
ens.integrated <- IntegrateData(anchorset = ens.anchors, dims = 1:30)

options(future.globals.maxSize= 5000*1024^2)
ens.features <- SelectIntegrationFeatures(object.list = ens.list, nfeatures = 1500)
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

# 0.5, okay
ens.integrated <- FindClusters(ens.integrated, resolution = 0.5)
table(ens.integrated$seurat_clusters)

# options(repr.plot.width=6, repr.plot.height=4)
# p1 <- DimPlot(object = ens.integrated, reduction = "tsne", group.by = "seurat_clusters", cols = brewer.pal(10,"Set3")) # seurat_clusters, platform
# p1

ens.integrated$cluster <- plyr::mapvalues(ens.integrated$seurat_clusters, from = 0:9, 
                            to=c("c7","c1","c4","c2","c5","c5","c6","c3","c8","c9"))

ens.integrated$cluster <- factor(ens.integrated$cluster, levels = paste("c",9:1,sep=""))
