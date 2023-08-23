
##Load in DESEQ2 object
ddsTxi <- readRDS("all_transcript_counts.RDS")
#load in sample metadata
samples <- read.csv("sample-metadata.csv")
idx <- which(samples$condition %in% c("DEV_POS", "DEV_NEG"))
ddsTxi <- ddsTxi[,idx]
ddsTxi$condition <- droplevels(ddsTxi$condition)
ddsTxi$mouse<- droplevels(ddsTxi$mouse)
ddsTxi$batch<- droplevels(ddsTxi$batch)
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "DEV_NEG")
##include mouse in design to control for individual mouse effects
design(ddsTxi) <- formula(~mouse + condition)
ddsTxi <- DESeq(ddsTxi)
res <- results(ddsTxi, name = "condition_DEV_POS_vs_DEV_NEG")
resultsNames(ddsTxi)
res_df <- as.data.frame(res)
vsd <- vst(ddsTxi, blind=FALSE)
plotPCA(vsd, intgroup = "condition")
ggsave("Figures-May3/DEV-PCA-by-Condition.jpeg", width = 10, height = 10, dpi = 1000)

res_df$gene_symbol <- gsub("\\..*","",rownames(res_df))
x <- which(res_df$log2FoldChange >= 0.5)
up_regulated_dev<- res_df[x, ]
write.csv(up_regulated_dev, "upregulated_dev_pos_vs_dev_neg_paired_final.csv")
x <- which(res_df$log2FoldChange <= -0.5)
down_regulated_dev<- res_df[x, ]
write.csv(down_regulated_dev, "downregulated_dev_pos_vs_dev_neg_paired_final.csv")
library(ggplot2)
library(ggrepel)
pval_threshold <- 1e-2
logfc_threshold1 <- 1
logfc_threshold2 <- -1
logfc_threshold <- 1
deseq.threshold1 <- which(res_df$log2FoldChange >= logfc_threshold1 & 
                            res_df$padj < pval_threshold)
deseq.threshold2 <- which(res_df$log2FoldChange <= logfc_threshold2 & 
                            res_df$padj < pval_threshold)
tmp <- rep("black", nrow(res_df))
tmp[deseq.threshold1] <- "red" 
tmp[deseq.threshold2] <- "blue" 
res_df$color <- tmp
tmp <- rep("NS", nrow(res_df))
tmp[deseq.threshold1] <- "Upregulated" 
tmp[deseq.threshold2] <- "Downregulated" 
res_df$label <- tmp
tmp <- as.factor(tmp)
g = ggplot(data=res_df, 
           aes(x=log2FoldChange, y=-log10(padj), 
               color=label)) +
  geom_point(alpha=0.4, size=1.75) +
  geom_vline(xintercept = logfc_threshold) +
  geom_vline(xintercept = -logfc_threshold) +
  geom_hline(yintercept = -log10(pval_threshold)) +
  xlab("log2 fold change") + ylab("-log10 FDR") +
  ggtitle("Differentially Expressed Genes DEV+ vs DEV-")
myColors <- c("blue", "black", "red")
res_df$label <- as.factor(res_df$label)
names(myColors) <- levels(res_df$label)
colScale <- scale_colour_manual(name = "label", values = myColors)
## Figure 3C
g + colScale
ggsave("PCA-dev-pos-vs-dev-neg-paired.pdf", width = 10, height = 10, dpi = 1000)


ddsTxi <- readRDS("all_transcript_counts.RDS")
idx <- which(samples$condition %in% c("MYMK_KO"))
ddsTxi <- ddsTxi[,-idx]
ddsTxi$condition <- droplevels(ddsTxi$condition)
ddsTxi$mouse<- droplevels(ddsTxi$mouse)
ddsTxi$batch<- droplevels(ddsTxi$batch)
ddsTxi <- DESeq(ddsTxi)
vsd <- vst(ddsTxi, blind=FALSE)
##Figure 3F
plotPCA(vsd, intgroup = "condition")


ddsTxi <- readRDS("all_transcript_counts.RDS")
idx <- which(samples$sample %in% c("Mov_Neg1", "Mov_Neg2", "Mov_Neg3", "Sham1", "Sham2", "Sham3", "Sham4"))
ddsTxi <- ddsTxi[,idx]
ddsTxi$condition <- droplevels(ddsTxi$condition)
ddsTxi$mouse<- droplevels(ddsTxi$mouse)
ddsTxi$batch<- droplevels(ddsTxi$batch)
ddsTxi$group <- droplevels(ddsTxi$group)
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "SHAM")
design(ddsTxi) <- formula(~condition)
ddsTxi <- DESeq(ddsTxi)
vsd <- vst(ddsTxi, blind=FALSE)
plotPCA(vsd)

resultsNames(ddsTxi)
res <- results(ddsTxi, name = "condition_MOV_NEG_vs_SHAM")
res_df <- as.data.frame(res)
res_df$gene_symbol <- gsub("\\..*","",rownames(res_df))
x <- which(res_df$log2FoldChange >= 0.5)
up_regulated_dev<- res_df[x, ]
write.csv(up_regulated_dev, "upregulated_mov_neg_13_vs_sham_final.csv")
x <- which(res_df$log2FoldChange <= -0.5)
down_regulated_dev<- res_df[x, ]
write.csv(down_regulated_dev, "downregulated_mov_neg_13_vs_sham_final.csv")

library(ggplot2)
library(ggrepel)
pval_threshold <- 1e-2
logfc_threshold1 <- 1
logfc_threshold2 <- -1
logfc_threshold <- 1


deseq.threshold1 <- which(res_df$log2FoldChange >= logfc_threshold1 & 
                            res_df$padj < pval_threshold)
deseq.threshold2 <- which(res_df$log2FoldChange <= logfc_threshold2 & 
                            res_df$padj < pval_threshold)




test <- rep("black", nrow(res_df))
test[deseq.threshold1] <- "red" 
test[deseq.threshold2] <- "blue" 

res_df$color <- test

test <- rep("black", nrow(res_df))
test[deseq.threshold1] <- "red" 
test[deseq.threshold2] <- "blue" 

res_df$color <- test
test <- rep("NS", nrow(res_df))
test[deseq.threshold1] <- "Upregulated" 
test[deseq.threshold2] <- "Downregulated" 

res_df$label <- test
test <- as.factor(test)
g = ggplot(data=res_df, 
           aes(x=log2FoldChange, y=-log10(padj), 
               # Colour based on the threshold defined before
               color=label)) +
  # Define the look of the points
  geom_point(alpha=0.4, size=1.75) +
  # Hide the legend
  # Apply another theme
  # Add the lines separating the DEGs
  geom_vline(xintercept = logfc_threshold) +
  geom_vline(xintercept = -logfc_threshold) +
  geom_hline(yintercept = -log10(pval_threshold)) +
  xlab("log2 fold change") + ylab("-log10 FDR") +
  ggtitle("Differentially Expressed Genes MOV NEG vs SHAM")
myColors <- c("blue", "black", "red")
res_df$label <- as.factor(res_df$label)
names(myColors) <- levels(res_df$label)
colScale <- scale_colour_manual(name = "label", values = myColors)
#Figure 3G (MOV NEG)
g + colScale

ddsTxi <- readRDS("all_transcript_counts.RDS")
idx <- which(samples$sample %in% c("Mov_Pos1", "Mov_Pos2", "Mov_Pos3", "Sham1", "Sham2", "Sham3", "Sham4"))
ddsTxi <- ddsTxi[,idx]
ddsTxi$condition <- droplevels(ddsTxi$condition)
ddsTxi$mouse<- droplevels(ddsTxi$mouse)
ddsTxi$batch<- droplevels(ddsTxi$batch)
ddsTxi$group <- droplevels(ddsTxi$group)
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "SHAM")
design(ddsTxi) <- formula(~condition)
ddsTxi <- DESeq(ddsTxi)
vsd <- vst(ddsTxi, blind=FALSE)
plotPCA(vsd)
resultsNames(ddsTxi)
res <- results(ddsTxi, name = "condition_MOV_POS_vs_SHAM")
res_df <- as.data.frame(res)
res_df$gene_symbol <- gsub("\\..*","",rownames(res_df))
x <- which(res_df$log2FoldChange >= 0.5)
up_regulated_dev<- res_df[x, ]
write.csv(up_regulated_dev, "upregulated_mov_pos_13_vs_sham_final.csv")
x <- which(res_df$log2FoldChange <= -0.5)
down_regulated_dev<- res_df[x, ]
write.csv(down_regulated_dev, "downregulated_mov_pos_13_vs_sham_final.csv")

pval_threshold <- 1e-2
logfc_threshold1 <- 1
logfc_threshold2 <- -1
logfc_threshold <- 1

deseq.threshold1 <- which(res_df$log2FoldChange >= logfc_threshold1 & 
                            res_df$padj < pval_threshold)
deseq.threshold2 <- which(res_df$log2FoldChange <= logfc_threshold2 & 
                            res_df$padj < pval_threshold)

test <- rep("black", nrow(res_df))
test[deseq.threshold1] <- "red" 
test[deseq.threshold2] <- "blue" 

res_df$color <- test

test <- rep("black", nrow(res_df))
test[deseq.threshold1] <- "red" 
test[deseq.threshold2] <- "blue" 

res_df$color <- test
test <- rep("NS", nrow(res_df))
test[deseq.threshold1] <- "Upregulated" 
test[deseq.threshold2] <- "Downregulated" 

res_df$label <- test
test <- as.factor(test)
g = ggplot(data=res_df, 
           aes(x=log2FoldChange, y=-log10(padj), 
               color=label)) +
  geom_point(alpha=0.4, size=1.75) +
  geom_vline(xintercept = logfc_threshold) +
  geom_vline(xintercept = -logfc_threshold) +
  geom_hline(yintercept = -log10(pval_threshold)) +
  xlab("log2 fold change") + ylab("-log10 FDR") +
  ggtitle("Differentially Expressed Genes MOV POS vs SHAM")
myColors <- c("blue", "black", "red")
res_df$label <- as.factor(res_df$label)
names(myColors) <- levels(res_df$label)
colScale <- scale_colour_manual(name = "label", values = myColors)
#Figure 3G (MOV Pos)
g + colScale

##GFP+ VS GFP- MOV DEGs
ddsTxi <- readRDS("all_transcript_counts.RDS")
idx <- which(samples$sample %in% c("Mov_Pos1", "Mov_Pos2", "Mov_Pos3", "Mov_Neg1", "Mov_Neg2", "Mov_Neg3"))
ddsTxi <- ddsTxi[,idx]
ddsTxi$condition <- droplevels(ddsTxi$condition)
ddsTxi$mouse<- droplevels(ddsTxi$mouse)
ddsTxi$batch<- droplevels(ddsTxi$batch)
ddsTxi$group <- droplevels(ddsTxi$group)
ddsTxi$condition <- relevel(ddsTxi$condition, ref = "MOV_NEG")
design(ddsTxi) <- formula(~mouse + condition)
ddsTxi <- DESeq(ddsTxi)
vsd <- vst(ddsTxi, blind=FALSE)
plotPCA(vsd)

resultsNames(ddsTxi)
res <- results(ddsTxi, name = "condition_MOV_POS_vs_MOV_NEG")
res_df <- as.data.frame(res)
res_df$gene_symbol <- gsub("\\..*","",rownames(res_df))
x <- which(res_df$log2FoldChange >= 0.5)
up_regulated_dev<- res_df[x, ]
write.csv(up_regulated_dev, "Figures-May3/upregulated_mov_pos13_vs_mov_neg13_paired_final.csv")
x <- which(res_df$log2FoldChange <= -0.5)
down_regulated_dev<- res_df[x, ]
write.csv(down_regulated_dev, "Figures-May3/downregulated_mov_pos13_vs_mov_neg13_paired_final.csv")

library(ggplot2)
library(ggrepel)
pval_threshold <- 1e-2
logfc_threshold1 <- 1
logfc_threshold2 <- -1
logfc_threshold <- 1
deseq.threshold1 <- which(res_df$log2FoldChange >= logfc_threshold1 & 
                            res_df$padj < pval_threshold)
deseq.threshold2 <- which(res_df$log2FoldChange <= logfc_threshold2 & 
                            res_df$padj < pval_threshold)
test <- rep("black", nrow(res_df))
test[deseq.threshold1] <- "red" 
test[deseq.threshold2] <- "blue" 
res_df$color <- test
test <- rep("black", nrow(res_df))
test[deseq.threshold1] <- "red" 
test[deseq.threshold2] <- "blue" 
res_df$color <- test
test <- rep("NS", nrow(res_df))
test[deseq.threshold1] <- "Upregulated" 
test[deseq.threshold2] <- "Downregulated" 
res_df$label <- test
test <- as.factor(test)
g = ggplot(data=res_df, 
           aes(x=log2FoldChange, y=-log10(padj), 
               color=label)) +
  geom_point(alpha=0.4, size=1.75) +
  geom_vline(xintercept = logfc_threshold) +
  geom_vline(xintercept = -logfc_threshold) +
  geom_hline(yintercept = -log10(pval_threshold)) +
  xlab("log2 fold change") + ylab("-log10 FDR") +
  ggtitle("Differentially Expressed Genes MOV13 POS and Neg 13")
myColors <- c("blue", "black", "red")
res_df$label <- as.factor(res_df$label)
names(myColors) <- levels(res_df$label)
colScale <- scale_colour_manual(name = "label", values = myColors)
#Figure 3 H
g + colScale


