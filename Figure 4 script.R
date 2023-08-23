##Figure 4 Scripts

library(Seurat)
library(dplyr)
library(ggplot2)
list <- readRDS("~/Desktop/Preprocessing-finalized-nf-1-23.RDS")
for (i in 1:length(list)) {
  z <- quantile(list[[i]]@meta.data[["nFeature_RNA"]], probs = c(0.05, 0.1, 0.15, 0.9, 0.95))
  y <- quantile(list[[i]]@meta.data[["nCount_RNA"]], probs = c(0.05, 0.1, 0.15, 0.9, 0.95))
  list[[i]][["percent.mt"]] <- PercentageFeatureSet(list[[i]], pattern = "^mt-")
  m <- quantile(list[[i]]@meta.data[["percent.mt"]], probs = c(0.05, 0.1, 0.15, 0.9, 0.95))
  list[[i]] <- subset(list[[i]], subset = percent.mt < 5)
  list[[i]] <-  subset(list[[i]], subset = nFeature_RNA <= z[[5]])
  list[[i]] <-  subset(list[[i]], subset = nCount_RNA <= y[[5]])
  list[[i]] <-  subset(list[[i]], subset = percent.mt <= m[[5]])
  list[[i]] <- SCTransform(list[[i]], vst.flavor = "v2", vars.to.regress = "percent.mt")
  list[[i]] <- RunPCA(list[[i]])
  list[[i]] <- RunUMAP(list[[i]], dims = 1:30)
  list[[i]] <- FindNeighbors(list[[i]], dims = 1:30)
  list[[i]] <- FindClusters(list[[i]], resolution = 0.5)
}


list_mov <- c(list[3], list[4], list[5])
features <- SelectIntegrationFeatures(object.list = list_mov, nfeatures = 3000)
list_mov <- PrepSCTIntegration(object.list = list_mov, anchor.features = features)

immune.anchors <- FindIntegrationAnchors(object.list = list_mov, normalization.method = "SCT",
                                         anchor.features = features)


int_mov <- IntegrateData(immune.anchors, normalization.method = "SCT")
DefaultAssay(int_mov) <- "integrated"
int_mov <- RunPCA(int_mov, npcs = 60, verbose = TRUE)
int_mov <- FindNeighbors(int_mov, reduction = "pca", dims = 1:35)
int_mov <- FindClusters(int_mov, resolution = 1)
int_mov <- RunUMAP(int_mov, reduction = "pca", dims = 1:35)
DimPlot(int_mov, split.by = "orig.ident", label = TRUE)




Idents(int_mov) <- "doublet"

int_mov <- subset(int_mov, idents = c("Singlet"))
DefaultAssay(int_mov) <- "integrated"
int_mov <- RunPCA(int_mov, npcs = 60, verbose = TRUE)
int_mov <- FindNeighbors(int_mov, reduction = "pca", dims = 1:35)
int_mov <- FindClusters(int_mov, resolution = 0.5)
int_mov <- RunUMAP(int_mov, reduction = "pca",dims = 1:35)
DimPlot(int_mov, split.by = "orig.ident", label = TRUE)
DefaultAssay(int_mov) <- "SCT"
DimPlot(int_mov, pt.size = 1) + NoLegend() + NoAxes()
ggsave("mov-umap-before-subset.pdf", width = 10, height = 10, dpi = 1000)
DefaultAssay(int_mov) <- "SCT"
FeaturePlot(int_mov, features = c("Myh2"), pt.size = 1) + NoAxes() + NoLegend()
ggsave("mov-umap-before-subset-myh2-feature.pdf", width = 10, height = 10, dpi = 1000)
FeaturePlot(int_mov, features = c("Ttn"), pt.size = 1) + NoAxes() + NoLegend()
ggsave("mov-umap-before-subset-Ttn-feature.pdf", width = 10, height = 10, dpi = 1000)
FeaturePlot(int_mov, features = c("Myh4"), pt.size = 1) + NoAxes() + NoLegend()
ggsave("mov-umap-before-subset-myh4-feature.pdf", width = 10, height = 10, dpi = 1000)
FeaturePlot(int_mov, features = c("Myh1"), pt.size = 1) + NoAxes() + NoLegend()
ggsave("mov-umap-before-subset-myh1-feature.pdf", width = 10, height = 10, dpi = 1000)
FeaturePlot(int_mov, features = c("Ckm"), pt.size = 1) + NoAxes() + NoLegend()
ggsave("mov-umap-before-subset-ckm-feature.pdf", width = 10, height = 10, dpi = 1000)
FeaturePlot(int_mov, features = c("Ebf1"), pt.size = 1) + NoLegend() + NoAxes()
ggsave("mov-umap-before-subset-ebf1-feature.pdf", width = 10, height = 10, dpi = 1000)

write.csv(table(Idents(int_mov)), "int-mov-raw-cluster-counts-all-cells.csv")
Idents(int_mov) <- "orig.ident"
x <- subset(int_mov, idents = c("Mov1"))
Idents(x) <- "seurat_clusters"
write.csv(table(Idents(x)), "int-mov-raw-cluster-counts-mov1.csv")
x <- subset(int_mov, idents = c("Mov2"))
Idents(x) <- "seurat_clusters"
write.csv(table(Idents(x)), "int-mov-raw-cluster-counts-mov2.csv")
x <- subset(int_mov, idents = c("Sed"))
Idents(x) <- "seurat_clusters"
write.csv(table(Idents(x)), "int-mov-raw-cluster-counts-sham.csv")
Idents(int_mov) <- "seurat_clusters"
DimPlot(int_mov, label = TRUE, pt.size = 1) + NoLegend() + NoAxes()
ggsave("int-mov-umap-before-subset-labelled.pdf", width = 10, height = 10, dpi = 1000)

myogenic_mov <- subset(int_mov, idents = c("8", "11", "10", "7"), invert = TRUE)
x <- which(myogenic_mov@meta.data[["orig.ident"]] == "Mov1")
myogenic_mov@meta.data[["orig.ident"]][x] <- "Mov"
x <- which(myogenic_mov@meta.data[["orig.ident"]] == "Mov2")
myogenic_mov@meta.data[["orig.ident"]][x] <- "Mov"


DefaultAssay(myogenic_mov) <- "integrated"
myogenic_mov <- RunPCA(myogenic_mov, npcs = 60, verbose = TRUE)
myogenic_mov <- FindNeighbors(myogenic_mov, reduction = "pca", dims = 1:20)
myogenic_mov <- FindClusters(myogenic_mov, resolution = 0.5)
myogenic_mov <- RunUMAP(myogenic_mov, reduction = "pca",dims = 1:20)
DimPlot(myogenic_mov, label = TRUE)
myogenic_mov <- RenameIdents(myogenic_mov, "3" = "Type IIb Myonuclei", "2" = "Type IIb Myonuclei", "4" = "Type IIb Myonuclei", "8" = "Type IIb Myonuclei", 
                             "0" = "Type IIx Myonuclei", "6" = "Type IIx Myonuclei", "5" = "Atf3+ Myonuclei", "1" = "Type IIa Myonuclei", "11" = "Type I Myonuclei", "9" = "Type IIa Myonuclei", "10" = "Unclassified B", "7" = "Unclassified A")


DimPlot(myogenic_mov, label = FALSE, pt.size = 1) + NoLegend() + NoAxes()
ggsave("Figures-May3/mov-int-umap-myogenic.pdf", width = 10, height = 10, dpi = 1000)

DimPlot(myogenic_mov, label = TRUE, pt.size = 1) + NoLegend() + NoAxes()
ggsave("Figures-May3/mov-int-umap-myogenic-labelled.pdf", width = 10, height = 10, dpi = 1000)





DimPlot(myogenic_mov, label = TRUE, split.by = "orig.ident") + NoLegend() + NoAxes()
ggsave("Figures-May3/mov-int-umap-myogenic-labelled-split-by.pdf", width = 20, height = 10, dpi = 1000)
DimPlot(myogenic_mov, label = TRUE, group.by = "orig.ident") + NoLegend() + NoAxes()
ggsave("Figures-May3/mov-int-umap-myogenic-labelled-group-by.pdf", width = 10, height = 10, dpi = 1000)


myogenic_mov$celltype <- Idents(myogenic_mov)
Idents(myogenic_mov) <- "orig.ident"
x <- subset(myogenic_mov, idents = c("Mov"))
Idents(x) <- "celltype"
DimPlot(x, label = FALSE, pt.size = 1) + NoAxes() + NoLegend()
ggsave("Figures-May3/mov-samples-only-split-integrated-myogenic-umap.pdf", width = 10, height = 10, dpi = 1000)
x <- subset(myogenic_mov, idents = c("Sed"))
Idents(x) <- "celltype"
DimPlot(x, label = FALSE, pt.size = 1) + NoAxes() + NoLegend()
ggsave("Figures-May3/sham-sample-only-split-integrated-myogenic-umap.pdf", width = 10, height = 10, dpi = 1000)

Idents(myogenic_mov) <- "celltype"
DimPlot(myogenic_mov,  pt.size = 1, label = FALSE) + NoAxes() + NoLegend()
ggsave("Figures-May3/myogenic-mov-integrated.pdf", width = 10, height = 10, dpi = 1000)


myogenic_mov <- PrepSCTFindMarkers(myogenic_mov)
DefaultAssay(myogenic_mov) <- "SCT"
myogenic_mov_markers <- FindAllMarkers(myogenic_mov, only.pos = TRUE)
write.csv(myogenic_mov_markers, "Figures-May3/all-myogenic-mov-markers-integrated-sct.csv")
DefaultAssay(myogenic_mov) <- "RNA"
myogenic_mov_markers_rna <- FindAllMarkers(myogenic_mov, only.pos = TRUE)
write.csv(myogenic_mov_markers_rna, "Figures-May3/all-myogenic-mov-markers-integrated-rna.csv")
DefaultAssay(p10_fivemonth_myogenic) <- "SCT"
Idents(p10_fivemonth_myogenic) <- "celltype"
FeaturePlot(p10_fivemonth_myogenic, features = c("Ebf1"), pt.size = 1) + NoLegend() + NoAxes()
ggsave("Figures-May3/p10-fivemonth-myogenic-ebf1-feature.pdf", width = 10, height = 10, dpi = 1000)
FeaturePlot(p10_fivemonth_myogenic, features = c("Col1a1"), pt.size = 1) + NoLegend() + NoAxes()
ggsave("Figures-May3/p10-fivemonth-myogenic-col1a1-feature.pdf", width = 10, height = 10, dpi = 1000)



x<- FindAllMarkers(p10_fivemonth_myogenic, only.pos = TRUE)
write.csv(x, "Figures-May3/all-myogenic-p10-fivemonth-markers-integrated-sct.csv")
DefaultAssay(p10_fivemonth_myogenic) <- "RNA"

Idents(p10_fivemonth_myogenic) <- "celltype"
x <- FindAllMarkers(p10_fivemonth_myogenic, only.pos = TRUE)
write.csv(x, "Figures-May3/all-myogenic-p10-fivemonth-markers-integrated-rna.csv")

DefaultAssay(myogenic_mov) <- "RNA"
x<- FindMarkers(myogenic_mov, ident.1 = c("Unclassified A"), ident.2 = c("Unclassified B"))
write.csv(x, "unclassified-A-vs-B-RNA.csv")
DefaultAssay(myogenic_mov) <- "SCT"
x<- FindMarkers(myogenic_mov, ident.1 = c("Unclassified A"), ident.2 = c("Unclassified B"))
write.csv(x, "unclassified-A-vs-B-SCT.csv")

myogenic_mov$celltype.orig <- paste(Idents(myogenic_mov), myogenic_mov$orig.ident, sep = "_")
Idents(myogenic_mov) <- "celltype.orig"
DefaultAssay(myogenic_mov) <- "SCT"
