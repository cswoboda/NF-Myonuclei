library(Seurat)
library(dplyr)
library(ggplot2)
library(reticulate)
reticulate::py_discover_config(required_module = "phate")
reticulate::import("phate")
library(phateR)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(stringr)
##loaded data is processed with both cellbender and solo
list <- readRDS("~/Desktop/Preprocessing-finalized-nf-1-23.RDS")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes <- tolower(s.genes)
g2m.genes <- tolower(g2m.genes)
s.genes <- str_to_title(s.genes)
g2m.genes <- str_to_title(g2m.genes)
###data processing
for (i in 1:length(list)) {
  z <- quantile(list[[i]]@meta.data[["nFeature_RNA"]], probs = c(0.05, 0.1, 0.15, 0.9, 0.95))
  y <- quantile(list[[i]]@meta.data[["nCount_RNA"]], probs = c(0.05, 0.1, 0.15, 0.9, 0.95))
  list[[i]][["percent.mt"]] <- PercentageFeatureSet(list[[i]], pattern = "^mt-")
  m <- quantile(list[[i]]@meta.data[["percent.mt"]], probs = c(0.05, 0.1, 0.15, 0.9, 0.95))
  list[[i]] <- subset(list[[i]], subset = percent.mt < 5)
  list[[i]] <-  subset(list[[i]], subset = nFeature_RNA <= z[[5]])
  list[[i]] <-  subset(list[[i]], subset = nCount_RNA <= y[[5]])
  list[[i]] <-  subset(list[[i]], subset = percent.mt <= m[[5]])
  list[[i]] <- NormalizeData(list[[i]])
  list[[i]] <- CellCycleScoring(list[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  list[[i]] <- SCTransform(list[[i]], vst.flavor = "v2", vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
  list[[i]] <- RunPCA(list[[i]])
  list[[i]] <- RunUMAP(list[[i]], dims = 1:30)
  list[[i]] <- FindNeighbors(list[[i]], dims = 1:30)
  list[[i]] <- FindClusters(list[[i]], resolution = 0.5)
}
###pull out just P10 object from list
p10 <- list[[1]]

DimPlot(p10, label = TRUE)
Idents(p10) <- "doublet"
p10 <- subset(p10, idents = c("Singlet"))
DimPlot(p10, label = TRUE)
p10 <- RunUMAP(p10, dims = 1:30)
p10 <- FindNeighbors(p10, dims = 1:30)
p10 <- FindClusters(p10, resolution = 1)
Idents(p10) <- "seurat_clusters"
###generate markers for cluster identification
all.markers.p10.clean <- FindAllMarkers(p10, only.pos = TRUE)
DimPlot(p10, label = TRUE)
###subset out myogenic populations
p10_myogenic <- subset(p10, idents = c("11", "14", "5", "2", "3", "0", "1", "12", "16", "7"))
###reperform dimensionality reduction
p10_myogenic<- RunUMAP(p10_myogenic, dims = 1:30)
p10_myogenic <- FindNeighbors(p10_myogenic, dims = 1:30)
p10_myogenic <- FindClusters(p10_myogenic, resolution = 1)
###generate all markers for myogenic populations
all.markers.p10.myogenic <- FindAllMarkers(p10_myogenic, only.pos = TRUE)
DimPlot(p10_myogenic, pt.size = 1) + NoAxes() + NoLegend()
dir.create("Figures-May3")
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-umap-unlabelled-raw-clusters.jpeg", dpi = 1000, height = 10, width = 10)
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-umap-unlabelled-raw-clusters.pdf", dpi = 1000, height = 10, width = 10)

DimPlot(p10_myogenic, pt.size = 1, label = TRUE) + NoAxes() + NoLegend()
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-umap-labelled-raw-clusters.jpeg", dpi = 1000, height = 10, width = 10)
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-umap-labelled-raw-clusters.pdf", dpi = 1000, height = 10, width = 10)
##Sanity check for fiber type
FeaturePlot(p10_myogenic, features = c("Myh4"), pt.size = 1) + NoAxes()
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-feature-myh4.jpeg", dpi = 1000, height = 10, width = 10)
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-feature-myh4.pdf", dpi = 1000, height = 10, width = 10)


FeaturePlot(p10_myogenic, features = c("Myh1"), pt.size = 1) + NoAxes()
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-feature-myh1.jpeg", dpi = 1000, height = 10, width = 10)
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-feature-myh1.pdf", dpi = 1000, height = 10, width = 10)

FeaturePlot(p10_myogenic, features = c("Myh2"), pt.size = 1) + NoAxes()
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-feature-myh2.jpeg", dpi = 1000, height = 10, width = 10)
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-feature-myh2.pdf", dpi = 1000, height = 10, width = 10)

FeaturePlot(p10_myogenic, features = c("Myh7"), pt.size = 1) + NoAxes()
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-feature-myh7.jpeg", dpi = 1000, height = 10, width = 10)
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-feature-myh7.pdf", dpi = 1000, height = 10, width = 10)

FeaturePlot(p10_myogenic, features = c("Pax7"), pt.size = 1) + NoAxes()
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-feature-pax7.jpeg", dpi = 1000, height = 10, width = 10)
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-feature-pax7.pdf", dpi = 1000, height = 10, width = 10)

FeaturePlot(p10_myogenic, features = c("Myog"), pt.size = 1) + NoAxes()
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-feature-myog.jpeg", dpi = 1000, height = 10, width = 10)
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-feature-myog.pdf", dpi = 1000, height = 10, width = 10)
Idents(p10_myogenic) <- "seurat_clusters"
###Name clusters
p10_myogenic <- RenameIdents(p10_myogenic, "0" = "Type IIb Myonuclei", "1" = "Type IIb Myonuclei", "5" = "Type IIb Myonuclei",
                             "7" = "Type IIb Myonuclei", "9" = "Type IIb Myonuclei", "6" = "Type IIb Myonuclei", "3" = "Type IIb Myonuclei", "4" = "Type IIx Myonuclei", 
                             "2" = "Type IIx Myonuclei", "10" = "Type IIa Myonuclei", "13" = "Type I Myonuclei", "12" = "Dev 2", "11" = "Dev 1", "14" = "MUSCs", "8" = "MUSCs")
###Figure 1A
DimPlot(p10_myogenic, pt.size = 1, label = TRUE) + NoAxes() + NoLegend()
###Figure 1B
genes <- c("Pax7", "Myog", "Ttn", "Myh4", "Myh1", "Myh2", "Myh7")
DotPlot(p10_myogenic, features = genes) 

###Phate Analysis
reticulate::py_discover_config(required_module = "phate")
library(phateR)
DefaultAssay(p10_myogenic) <- "RNA"
p10_myogenic <- NormalizeData(p10_myogenic)
seurat_data <- as.data.frame(p10_myogenic@assays$RNA@data)
phate_data_input <- t(seurat_data)
phate_output <- phate(phate_data_input, t = 30) 
ggplot(phate_output, aes(x=PHATE1, y=PHATE2)) +
  geom_point()
p10_myogenic[["PHATE"]] <- CreateDimReducObject(embeddings = phate_output$embedding, key = "PHATE_", assay = DefaultAssay(p10_myogenic))

##Figure 1C
DimPlot(p10_myogenic , reduction = "PHATE", label = TRUE) + ggtitle(label = "PHATE") + NoLegend()
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-phate-no-legend.pdf", dpi = 1000, height = 10, width = 15)
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-phate-no-legend.jpeg", dpi = 1000, height = 10, width = 15)

cells <- Cells(subset(p10_myogenic, idents = c("Dev 1")))
DimPlot(p10_myogenic , reduction = "PHATE", cells.highlight = cells) + ggtitle(label = "PHATE") + NoLegend()
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-phate-no-legend-unclassified1-highlight.pdf", dpi = 1000, height = 10, width = 15)
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-phate-no-legend-unclassified1-highlight.jpeg", dpi = 1000, height = 10, width = 15)

cells <- Cells(subset(p10_myogenic, idents = c("Dev 2")))
DimPlot(p10_myogenic , reduction = "PHATE", cells.highlight = cells) + ggtitle(label = "PHATE") + NoLegend()
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-phate-no-legend-unclassified2-highlight.pdf", dpi = 1000, height = 10, width = 15)
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-phate-no-legend-unclassified2-highlight.jpeg", dpi = 1000, height = 10, width = 15)

##Figure 1D
DefaultAssay(p10_myogenic) <- "SCT"
p10_myogenic_markers_clusters_labelled <- FindAllMarkers(p10_myogenic, only.pos = TRUE)
p10_myogenic_markers_clusters_labelled %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10

DoHeatmap(p10_myogenic, features = top10$gene) + NoLegend()

FeaturePlot(p10_myogenic, features = c("Ebf1"), pt.size = 1) + NoLegend() + NoAxes()
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-ebf1-feature.pdf", width = 10, height = 10, dpi = 1000)
FeaturePlot(p10_myogenic, features = c("Col1a1"), pt.size = 1) + NoLegend() + NoAxes()
ggsave("~/Desktop/Research/Bulk-RNA-Seq/Figures-May3/p10-myogenic-col1a1-feature.pdf", width = 10, height = 10, dpi = 1000)

###Integration for Figure 1E
list_p10_fivemonth <- c(list[1], list[2])
features <- SelectIntegrationFeatures(object.list = list_p10_fivemonth, nfeatures = 3000)
list_p10_fivemonth <- PrepSCTIntegration(object.list = list_p10_fivemonth, anchor.features = features)

immune.anchors <- FindIntegrationAnchors(object.list = list_p10_fivemonth, normalization.method = "SCT",
                                         anchor.features = features)


int_p10_fivemonth <- IntegrateData(immune.anchors, normalization.method = "SCT")
DefaultAssay(int_p10_fivemonth) <- "integrated"
int_p10_fivemonth <- RunPCA(int_p10_fivemonth, npcs = 60, verbose = TRUE)
int_p10_fivemonth <- FindNeighbors(int_p10_fivemonth, reduction = "pca", dims = 1:35)
int_p10_fivemonth <- FindClusters(int_p10_fivemonth, resolution = 1)
int_p10_fivemonth <- RunUMAP(int_p10_fivemonth, reduction = "pca", dims = 1:35)
DimPlot(int_p10_fivemonth, split.by = "orig.ident", label = TRUE)
int_p10_fivemonth$int_cluster <- Idents(int_p10_fivemonth)
Idents(int_p10_fivemonth) <- "doublet"
int_p10_fivemonth <- subset(int_p10_fivemonth, idents = c("Singlet"))
int_p10_fivemonth <- FindNeighbors(int_p10_fivemonth, reduction = "pca", dims = 1:50)
int_p10_fivemonth <- FindClusters(int_p10_fivemonth, resolution = 1)
int_p10_fivemonth <- RunUMAP(int_p10_fivemonth, reduction = "pca", dims = 1:50)
DimPlot(int_p10_fivemonth, pt.size = 1) + NoLegend() + NoAxes()

int_p10_fivemonth <- RenameIdents(int_p10_fivemonth, "3" = "Type IIb Myonuclei", "11" = "Type IIb Myonuclei", "7" = "Type IIb Myonuclei", "2" = "Type IIb Myonuclei", "9" = "Type IIb Myonuclei", "0" = "Type IIb Myonuclei", "5" = "Type IIx Myonuclei", "1" = "Type IIx Myonuclei", "4" = "Type IIb Myonuclei", "10" = "Type IIa Myonuclei", "20" = "Type I Myonuclei", 
                                  "12" = "Myotendinous Junction", "18" = "Developing Myonuclei", "16" = "Neuromuscular Junction", "14" = "MUSCs")
DefaultAssay(int_p10_fivemonth) <- "integrated"
p10_fivemonth_myogenic <- subset(int_p10_fivemonth, idents = c("Type IIb Myonuclei", "Type IIx Myonuclei", "Type IIa Myonuclei", "Type I Myonuclei", "MUSCs", "Developing Myonuclei"))
p10_fivemonth_myogenic <- RunUMAP(p10_fivemonth_myogenic, reduction = "pca", dims = 1:60)
##Figure 1E
DimPlot(p10_fivemonth_myogenic, split.by = "orig.ident")  + NoAxes()
p10_fivemonth_myogenic$celltype <-  Idents(p10_fivemonth_myogenic) 
p10_fivemonth_myogenic$celltype.orig <- paste(Idents(p10_fivemonth_myogenic), p10_fivemonth_myogenic$orig.ident, sep = "_")


Idents(p10_fivemonth_myogenic) <- "celltype"
p10_fivemonth_myogenic <- RenameIdents(p10_fivemonth_myogenic, "Type IIb Myonuclei" = "Body Myonuclei", "Type IIx Myonuclei" = "Body Myonuclei", "Type IIa Myonuclei" = "Body Myonuclei", "Type I Myonuclei" = "Body Myonuclei")
p10_fivemonth_myogenic$simple.orig <- paste(Idents(p10_fivemonth_myogenic), p10_fivemonth_myogenic$orig.ident, sep = "_")
Idents(p10_fivemonth_myogenic) <- "simple.orig"

dotplot_obj <- subset(p10_fivemonth_myogenic, idents = c("Body Myonuclei_P10", "Body Myonuclei_Five Months", "Developing Myonuclei_P10", "MUSCs_P10"))
DefaultAssay(p10_fivemonth_myogenic) <- "SCT"
dotplot_obj <- PrepSCTFindMarkers(dotplot_obj)
Idents(dotplot_obj) <- "simple.orig"
genes <- c('Neat1', 'Gm42418', 'Gm26917', 'Clcn5', 'Gm47283', 'Zswim6' ,'Megf10','Tcf7l1','Spats2l','Sema6a','Gab2','Tshz3','Slc9a9', 'Pcdh7', 'Elmo1', 'Sptbn1', 'Myzap', 'Dpp6', 'Utrn', 'Arhgap18', 'Ankrd12', 'Cobl', 'Casz1', 'Nkain2', 'Gm45159', '5430431A17Rik', 'Camta1', 'Dlg2', 'Igf2', 'H19')
###Figure 1F
DotPlot(dotplot_obj, features = unique(genes)) + coord_flip()




