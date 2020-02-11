library(devtools)
library(Seurat)
library(Matrix)
library(xlsx)
library(umap)
library(reshape2)
library(dplyr)
library(ggplot2)
library(cowplot)

## Setup the Seurat Object-WT
wt6h.data <- Read10X(data.dir = "../WT_hCG6h/outs/filtered_feature_bc_matrix/")
wt6h.atleastone <- apply(wt6h.data, 2, function(x) sum(x>0))
hist(wt6h.atleastone, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")
wt6h.tmp <- apply(wt6h.data, 1, function(x) sum(x>0))
table(wt6h.tmp>=3)
wt6h.keep <- wt6h.tmp>=3
wt6h.tmp <- wt6h.data[wt6h.keep,]
wt6h.atleastone <- apply(wt6h.tmp, 2, function(x) sum(x>0))
summary(wt6h.atleastone)
rownames(x = wt6h.data) <- gsub(pattern = '_', replacement = '-', x = rownames(x = wt6h.data))
WT6 <- CreateSeuratObject(counts = wt6h.data, project = "6h_WT", min.cells = 5, min.features = 200)
WT6$geno <- "WT"
WT6 <- subset(x = WT6, subset = nFeature_RNA > 500)
WT6 <- NormalizeData(object = WT6, verbose = FALSE)
WT6 <- FindVariableFeatures(object = WT6, selection.method = "vst", nfeatures = 2000)

## Setup the Seurat Object-KO
ko6h.data <- Read10X(data.dir = "../KO_hCG6h/outs/filtered_feature_bc_matrix/")
ko6h.atleastone <- apply(ko6h.data, 2, function(x) sum(x>0))
hist(ko6h.atleastone, breaks = 100,
     main = "Distribution of detected genes",
     xlab = "Genes with at least one tag")
ko6h.tmp <- apply(ko6h.data, 1, function(x) sum(x>0))
table(ko6h.tmp>=3)
ko6h.keep <- ko6h.tmp>=3
ko6h.tmp <- ko6h.data[ko6h.keep,]
ko6h.atleastone <- apply(ko6h.tmp, 2, function(x) sum(x>0))
summary(ko6h.atleastone)
rownames(x = ko6h.data) <- gsub(pattern = '_', replacement = '-', x = rownames(x = ko6h.data))
KO6 <- CreateSeuratObject(counts = ko6h.data, project = "6h_KO", min.cells = 5, min.features = 200)
KO6$geno <- "KO"
KO6 <- subset(x = KO6, subset = nFeature_RNA > 500)
KO6 <- NormalizeData(object = KO6, verbose = FALSE)
KO6 <- FindVariableFeatures(object = KO6, selection.method = "vst", nfeatures = 2000)

## Find Anchors
anchors <- FindIntegrationAnchors(object.list = list(WT6, KO6), dims = 1:20)
combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = combined) <- "integrated"

## visualization and clustering
combined <- ScaleData(object = combined, verbose = FALSE)
combined <- RunPCA(object = combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(object = combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(object = combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)
DefaultAssay(object = combined) <- "RNA"

## Plotting
p1 <- DimPlot(object = combined, reduction = "umap", group.by = "geno")
p2 <- DimPlot(object = combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(object = combined, reduction = "umap", split.by = "geno", label = TRUE, label.size = 5)
VlnPlot(object = combined, features = c("Pgr"),  
        pt.size = 0, combine = FALSE)

## Plotting - Localization of PGE2-receptor and cytokine genes
combined.Ptger2 <- subset(x = combined, subset = Ptger2 > 0.2)
combined.Il6 <- subset(x = combined, subset = Il6 > 0.2)
combined.Il1b <- subset(x = combined, subset = Il1b > 0.2)
combined.Tnf <- subset(x = combined, subset = Tnf > 0.2)
combined.Ccl2 <- subset(x = combined, subset = Ccl2 > 0.2)
DimPlot(object = combined.Ptger2, reduction = "umap", split.by = "geno", label.size = 5)
DimPlot(object = combined, reduction = "umap", label.size = 5)
DimPlot(object = combined.Ptger2, reduction = "umap", label.size = 5)
DimPlot(object = combined.Il6, reduction = "umap", label.size = 5)
DimPlot(object = combined.Il1b, reduction = "umap", label.size = 5)
DimPlot(object = combined.Tnf, reduction = "umap", label.size = 5)
DimPlot(object = combined.Ccl2, reduction = "umap", label.size = 5)
FeaturePlot(object = int.Ptger2, features = c("Il11", "Ptger2"), 
            cols = c("grey", "red"))

## Sub-clustering (granulosa cells)
cells.GCs <- subset(x = combined, idents = c("0", "2", "3", "4", "5", "6", "11"))
cells.GCs
cells.GCs$geno <- factor(x = cells.GCs$geno, levels = c('WT', 'KO'))
cells.GCs <- ScaleData(object = cells.GCs, verbose = FALSE)
cells.GCs <- RunPCA(object = cells.GCs, npcs = 30, verbose = FALSE)
cells.GCs <- RunUMAP(object = cells.GCs, reduction = "pca", dims = 1:25)
cells.GCs <- FindNeighbors(object = cells.GCs, reduction = "pca", dims = 1:25)
cells.GCs <- FindClusters(cells.GCs, resolution = 0.7)

## Plotting
p1 <- DimPlot(object = cells.GCs, reduction = "umap", group.by = "geno", label.size = 5)
p2 <- DimPlot(object = cells.GCs, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1)
plot_grid(p1, p2)
DimPlot(object = cells.GCs, reduction = "umap", split.by = "geno", label.size = 5)

# cumulus cell marker - Ube2c
FeaturePlot(object = cells.GCs, features = c("Ube2c"), split.by = "geno", 
            cols = c("grey", "red"), pt.size = 0.5)
# ovulatory mural cell marker - Pgr
FeaturePlot(object = cells.GCs, features = c("Pgr"), split.by = "geno", 
            cols = c("grey", "red"), pt.size = 0.5)

## Sub-clustering - mural granulosa cells
mGCs <- subset(x = cells.GCs, idents = c("1", "2", "3", "4", "5", "6", "8", "9", "10"))
mGCs
mGCs$geno <- factor(x = mGCs$geno, levels = c('WT', 'KO'))
DimPlot(object = mGCs, reduction = "umap", split.by = "geno", label.size = 5)

## Differential gene expression - mural granulosa cells
VlnPlot(object = mGCs, features = c("Pgr"),  
        pt.size = 0, combine = FALSE)
Idents(object = mGCs) <- "geno"
avg.mGCs <- log1p(x = AverageExpression(object = mGCs, verbose = FALSE)$RNA)
avg.mGCs$gene <- rownames(x = avg.mGCs)
Idents(object = mGCs) <- "celltype.geno"
mGCs$celltype <- Idents(object = mGCs)
mGCs$celltype.geno <- paste(Idents(object = mGCs), mGCs$geno, sep = "_")
head(mGCs$celltype.geno)
Idents(object = mGCs) <- "celltype.geno"
KO.response.mGCs <- FindMarkers(object = mGCs, ident.1 = c("celltype.geno_KO"), ident.2 = c("celltype.geno_WT"), verbose = FALSE)
head(x = KO.response.mGCs, n = 15)

## Sub-clustering (ovulatory mural; Pgr-positive mural granulosa cells)
ovul.mGCs <- subset(x = cells.GCs, idents = c("1", "2", "5", "10"))
ovul.mGCs
p1 <- DimPlot(object = ovul.mGCs, reduction = "umap", group.by = "geno", label.size = 5)
p2 <- DimPlot(object = ovul.mGCs, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1)
plot_grid(p1, p2)
DimPlot(object = ovul.mGCs, reduction = "umap", split.by = "geno", label = TRUE, label.size = 5)

## Differential gene expression in ovulatory mural GC
Idents(object = ovul.mGCs) <- "geno"
avg.ovul.mGCs <- log1p(x = AverageExpression(object = ovul.mGCs, verbose = FALSE)$RNA)
avg.ovul.mGCs$gene <- rownames(x = avg.ovul.mGCs)
Idents(object = ovul.mGCs) <- "celltype.geno"
ovul.mGCs$celltype <- Idents(object = ovul.mGCs)
ovul.mGCs$celltype.geno <- paste(Idents(object = ovul.mGCs), ovul.mGCs$geno, sep = "_")
head(ovul.mGCs$celltype.geno)
Idents(object = ovul.mGCs) <- "celltype.geno"
KO.response.ovul.mGCs <- FindMarkers(object = ovul.mGCs, ident.1 = c("celltype.geno_KO"), ident.2 = c("celltype.geno_WT"), verbose = FALSE)
head(x = KO.response.ovul.mGCs, n = 15)

## Plotting
genes.to.label = c("Sprr2g", "Mt2", "Gm8730", "Slc25a33", "Gm26917", "Susd3", "Mt1", "Etl4", "S100a6", "Socs2",
                   "F3", "Areg", "Tnfaip6", "Ptgs2", "Cited4", "Ereg", "Inhba", "Pde10a", "Ddit4l", "Btc")
table(KO.response.ovul.mGCs$p_val_adj<0.05)
# only mark genes with FDR p-value < 0.05
sig_genes <- KO.response.ovul.mGCs[KO.response.ovul.mGCs$p_val_adj<0.05,]  
#make common coloum "gene" for joining
sig_genes$gene <- rownames(sig_genes) 
#need to merge test results and epression table together so everything is in the same order
avg.ovul.mGCs <- avg.ovul.mGCs %>% left_join(sig_genes) 
avg.ovul.mGCs$status <- "NoSig"
avg.ovul.mGCs$status[avg.ovul.mGCs$avg_logFC>(0) & avg.ovul.mGCs$p_val_adj<0.05] <- "Up"   
avg.ovul.mGCs$status[avg.ovul.mGCs$avg_logFC<(0) & avg.ovul.mGCs$p_val_adj<0.05] <- "Down"
avg.ovul.mGCs$status <- factor(avg.ovul.mGCs$status, levels = c("NoSig", "Up", "Down"))
#re-arrange the df so "NoSig" points got drew first
avg.ovul.mGCs <- avg.ovul.mGCs %>% arrange(status) 
#re-order the levels to make the legend show in ggplot
avg.ovul.mGCs$status <- factor(avg.ovul.mGCs$status, levels = c( "Up","NoSig", "Down")) 
#for label points to work
rownames(avg.ovul.mGCs) <- avg.ovul.mGCs$gene 
pp1 <- ggplot(avg.ovul.mGCs, aes(WT, KO, color = status)) + geom_point() + ggtitle("ovulatory mural")+
  scale_color_manual(values=c("red",  "grey", "blue")) 
pp1 <- LabelPoints(plot = pp1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)
x11();plot_grid(pp1)

## Sub-clustering each cell type
# ovulatory mural-1
FeaturePlot(object = combined, features = c("Pgr"), cols = c("grey", "red"), label = TRUE, label.size = 5)
ov_mural_1 <- subset(x = combined, idents = "4")
Idents(object = ov_mural_1) <- "geno"
avg.ov_mural_1 <- log1p(x = AverageExpression(object = ov_mural_1, verbose = FALSE)$RNA)
avg.ov_mural_1$gene <- rownames(x = avg.ov_mural_1)
# ovulatory mural-2
ov_mural_2 <- subset(x = combined, idents = "0")
Idents(object = ov_mural_2) <- "geno"
avg.ov_mural_2 <- log1p(x = AverageExpression(object = ov_mural_2, verbose = FALSE)$RNA)
avg.ov_mural_2 <- rownames(x = avg.ov_mural_2)
# ovulatory cumulus
FeaturePlot(object = combined, features = c("Ube2c"), cols = c("grey", "red"), label = TRUE, label.size = 5)
ov_cumulus <- subset(x = combined, idents = "3")
Idents(object = ov_cumulus) <- "geno"
avg.ov_cumulus <- log1p(x = AverageExpression(object = ov_cumulus, verbose = FALSE)$RNA)
avg.ov_cumulus$gene <- rownames(x = avg.ov_cumulus)
# immune cells
FeaturePlot(object = combined, features = c("Ptprc"), cols = c("grey", "red"), label = TRUE, label.size = 5)
immune <- subset(x = combined, idents = "12")
Idents(object = immune) <- "geno"
avg.immune <- log1p(x = AverageExpression(object = immune, verbose = FALSE)$RNA)
avg.immune$gene <- rownames(x = avg.immune)
# theca cells
FeaturePlot(object = combined, features = c("Cyp17a1"), cols = c("grey", "red"), label = TRUE, label.size = 5)
theca <- subset(x = combined, idents = "10")
Idents(object = theca) <- "geno"
avg.theca <- log1p(x = AverageExpression(object = theca, verbose = FALSE)$RNA)
avg.theca$gene <- rownames(x = avg.theca)
# endothelial cells
FeaturePlot(object = combined, features = c("Icam2"), cols = c("grey", "red"), label = TRUE, label.size = 5)
endo <- subset(x = combined, idents = "7")
Idents(object = endo) <- "geno"
avg.endo <- log1p(x = AverageExpression(object = endo, verbose = FALSE)$RNA)
avg.endo$gene <- rownames(x = avg.endo)
# early granulosa cells
FeaturePlot(object = combined, features = c("Foxl2"), cols = c("grey", "red"), label = TRUE, label.size = 5)
earlygc <- subset(x = combined, idents = c("2","5", "6", "11"))
Idents(object = earlygc) <- "geno"
avg.earlygc <- log1p(x = AverageExpression(object = earlygc, verbose = FALSE)$RNA)
avg.earlygc$gene <- rownames(x = avg.earlygc)
# interstitial cells
FeaturePlot(object = combined, features = c("Col1a2"), cols = c("grey", "red"), label = TRUE, label.size = 5)
int <- subset(x = combined, idents = c("1", "8", "9", "13", "15", "17"))
Idents(object = int) <- "geno"
avg.int <- log1p(x = AverageExpression(object = int, verbose = FALSE)$RNA)
avg.int$gene <- rownames(x = avg.int)

## Differential gene expression in 'ov_mural_1'
Idents(object = ov_mural_1) <- "geno"
avg.ov_mural_1 <- log1p(x = AverageExpression(object = ov_mural_1, verbose = FALSE)$RNA)
avg.ov_mural_1$gene <- rownames(x = avg.ov_mural_1)
Idents(object = ov_mural_1) <- "celltype.geno"
ov_mural_1$celltype <- Idents(object = ov_mural_1)
ov_mural_1$celltype.geno <- paste(Idents(object = ov_mural_1), ov_mural_1$geno, sep = "_")
head(ov_mural_1$celltype.geno)
Idents(object = ov_mural_1) <- "celltype.geno"
KO.response.ov_mural_1 <- FindMarkers(object = ov_mural_1, ident.1 = c("celltype.geno_KO"), ident.2 = c("celltype.geno_WT"), verbose = FALSE)
head(x = KO.response.ov_mural_1, n = 40)
genes.to.label.4 = c("S100a6", "Gsta4", "Tmx1", "Susd3", "Sfrp4", "Tsc22d3", "Tmsb4x", "Mt2", "Hsd17b11", "S100a10", "Adamts1" 
                     , "Areg", "F3", "Tnfaip6", "Ptgs2", "Nap1l5", "Inhba", "Cdkn1a", "Ereg", "Btg1", "Gja1")
table(KO.response.ov_mural_1$p_val_adj<0.05)
sig_genes <- KO.response.ov_mural_1[KO.response.ov_mural_1$p_val_adj<0.05,]  
sig_genes$gene <- rownames(sig_genes) 
avg.ov_mural_1 <- avg.ov_mural_1 %>% left_join(sig_genes) 
avg.ov_mural_1$status <- "NoSig"
avg.ov_mural_1$status[avg.ov_mural_1$avg_logFC>(0) & avg.ov_mural_1$p_val_adj<0.05] <- "Up"   
avg.ov_mural_1$status[avg.ov_mural_1$avg_logFC<(0) & avg.ov_mural_1$p_val_adj<0.05] <- "Down"
avg.ov_mural_1$status <- factor(avg.ov_mural_1$status, levels = c("NoSig", "Up", "Down"))
avg.ov_mural_1 <- avg.ov_mural_1 %>% arrange(status) 
avg.ov_mural_1$status <- factor(avg.ov_mural_1$status, levels = c( "Up","NoSig", "Down")) 
rownames(avg.ov_mural_1) <- avg.ov_mural_1$gene 
pp1 <- ggplot(avg.ov_mural_1, aes(WT, KO, color = status)) + geom_point() + ggtitle("Ov_mural-1")+
  scale_color_manual(values=c("red",  "grey", "blue")) 
pp1 <- LabelPoints(plot = pp1, points = genes.to.label.4, repel = TRUE, xnudge = 0, ynudge = 0)
x11();plot_grid(pp1)

## Differential gene expression in 'ov_mural_2'
Idents(object = ov_mural_2) <- "geno"
avg.ov_mural_2 <- log1p(x = AverageExpression(object = ov_mural_2, verbose = FALSE)$RNA)
avg.ov_mural_2$gene <- rownames(x = avg.ov_mural_2)
Idents(object = ov_mural_2) <- "celltype.geno"
ov_mural_2$celltype <- Idents(object = ov_mural_2)
ov_mural_2$celltype.geno <- paste(Idents(object = ov_mural_2), ov_mural_2$geno, sep = "_")
head(ov_mural_2$celltype.geno)
Idents(object = ov_mural_2) <- "celltype.geno"
KO.response.ov_mural_2 <- FindMarkers(object = ov_mural_2, ident.1 = c("celltype.geno_KO"), ident.2 = c("celltype.geno_WT"), verbose = FALSE)
head(x = KO.response.ov_mural_2, n = 50)
genes.to.label.0 = c("Sprr2g", "Mt2", "Mt1", "Slc25a33", "S100a6", "Etl4", "Susd3", "Socs2", "Hsd17b11", "Lyve1", "Adamts1" 
                     , "F3", "Cited4", "Areg", "Tnfaip6", "Ptgs2", "Inhba", "Ereg", "Rgcc", "Ier3", "Ddit4l", "Btc", "Nfkbia")
table(KO.response.ov_mural_2$p_val_adj<0.05)
sig_genes <- KO.response.ov_mural_2[KO.response.ov_mural_2$p_val_adj<0.05,]  
sig_genes$gene <- rownames(sig_genes) 
avg.ov_mural_2 <- avg.ov_mural_2 %>% left_join(sig_genes) 
avg.ov_mural_2$status <- "NoSig"
avg.ov_mural_2$status[avg.ov_mural_2$avg_logFC>(0) & avg.ov_mural_2$p_val_adj<0.05] <- "Up"   
avg.ov_mural_2$status[avg.ov_mural_2$avg_logFC<(0) & avg.ov_mural_2$p_val_adj<0.05] <- "Down"
avg.ov_mural_2$status <- factor(avg.ov_mural_2$status, levels = c("NoSig", "Up", "Down"))
avg.ov_mural_2 <- avg.ov_mural_2 %>% arrange(status) 
avg.ov_mural_2$status <- factor(avg.ov_mural_2$status, levels = c( "Up","NoSig", "Down")) 
rownames(avg.ov_mural_2) <- avg.ov_mural_2$gene 
pp1 <- ggplot(avg.ov_mural_2, aes(WT, KO, color = status)) + geom_point() + ggtitle("Ov_mural_2")+
  scale_color_manual(values=c("red",  "grey", "blue")) 
pp1 <- LabelPoints(plot = pp1, points = genes.to.label.0, repel = TRUE, xnudge = 0, ynudge = 0)
x11();plot_grid(pp1)

## Differential gene expression in 'ov_cumulus'
Idents(object = ov_cumulus) <- "geno"
avg.ov_cumulus <- log1p(x = AverageExpression(object = ov_cumulus, verbose = FALSE)$RNA)
avg.ov_cumulus$gene <- rownames(x = avg.ov_cumulus)
Idents(object = ov_cumulus) <- "celltype.geno"
ov_cumulus$celltype <- Idents(object = ov_cumulus)
ov_cumulus$celltype.geno <- paste(Idents(object = ov_cumulus), ov_cumulus$geno, sep = "_")
head(ov_cumulus$celltype.geno)
Idents(object = ov_cumulus) <- "celltype.geno"
KO.response.ov_cumulus <- FindMarkers(object = ov_cumulus, ident.1 = c("celltype.geno_KO"), ident.2 = c("celltype.geno_WT"), verbose = FALSE)
head(x = KO.response.ov_cumulus, n = 15)
genes.to.label.3 = c("Mt2", "Socs2", "Gsta4", "S100a10", "Sprr2g", "Tmsb4x", "Mt1", "Tmx1", "Susd3", "S100a6", "Hsd17b11" 
                     , "F3", "Tnfaip6", "Areg", "Ptgs2", "Nap1l5", "Ereg", "Btc", "Odc1", "Inhba", "Cited4")
table(KO.response.ov_cumulus$p_val_adj<0.05)
sig_genes <- KO.response.ov_cumulus[KO.response.ov_cumulus$p_val_adj<0.05,]  
sig_genes$gene <- rownames(sig_genes) 
avg.ov_cumulus <- avg.ov_cumulus %>% left_join(sig_genes) 
avg.ov_cumulus$status <- "NoSig"
avg.ov_cumulus$status[avg.ov_cumulus$avg_logFC>(0) & avg.ov_cumulus$p_val_adj<0.05] <- "Up"   
avg.ov_cumulus$status[avg.ov_cumulus$avg_logFC<(0) & avg.ov_cumulus$p_val_adj<0.05] <- "Down"
avg.ov_cumulus$status <- factor(avg.ov_cumulus$status, levels = c("NoSig", "Up", "Down"))
avg.ov_cumulus <- avg.ov_cumulus %>% arrange(status) 
avg.ov_cumulus$status <- factor(avg.ov_cumulus$status, levels = c( "Up","NoSig", "Down")) 
rownames(avg.ov_cumulus) <- avg.ov_cumulus$gene 
pp1 <- ggplot(avg.ov_cumulus, aes(WT, KO, color = status)) + geom_point() + ggtitle("Ov_cumulus")+
  scale_color_manual(values=c("red",  "grey", "blue")) 
pp1 <- LabelPoints(plot = pp1, points = genes.to.label.3, repel = TRUE, xnudge = 0, ynudge = 0)
x11();plot_grid(pp1)

## Differential gene expression in 'immune'
Idents(object = immune) <- "geno"
avg.immune <- log1p(x = AverageExpression(object = immune, verbose = FALSE)$RNA)
avg.immune$gene <- rownames(x = avg.immune)
Idents(object = immune) <- "celltype.geno"
immune$celltype <- Idents(object = immune)
immune$celltype.geno <- paste(Idents(object = immune), immune$geno, sep = "_")
head(immune$celltype.geno)
Idents(object = immune) <- "celltype.geno"
KO.response.immune <- FindMarkers(object = immune, ident.1 = c("celltype.geno_KO"), ident.2 = c("celltype.geno_WT"), verbose = FALSE)
head(x = KO.response.immune, n = 15)
genes.to.label.12 = c("Ctsl", "Chil3", "Igkc", "Ccl3", "Slpi", "Ccl2", "Ccl5", "Ccl4", "Cd74", "Ctla2a", "Ccl2")
table(KO.response.immune$p_val_adj<0.05)
sig_genes <- KO.response.immune[KO.response.immune$p_val_adj<0.05,]  
sig_genes$gene <- rownames(sig_genes) 
avg.immune <- avg.immune %>% left_join(sig_genes) 
avg.immune$status <- "NoSig"
avg.immune$status[avg.immune$avg_logFC>(0) & avg.immune$p_val_adj<0.05] <- "Up"   
avg.immune$status[avg.immune$avg_logFC<(0) & avg.immune$p_val_adj<0.05] <- "Down"
avg.immune$status <- factor(avg.immune$status, levels = c("NoSig", "Up", "Down"))
avg.immune <- avg.immune %>% arrange(status) 
avg.immune$status <- factor(avg.immune$status, levels = c("NoSig", "Down")) 
rownames(avg.immune) <- avg.immune$gene 
pp1 <- ggplot(avg.immune, aes(WT, KO, color = status)) + geom_point() + ggtitle("Immune")+
  scale_color_manual(values=c("grey", "blue")) 
pp1 <- LabelPoints(plot = pp1, points = genes.to.label.12, repel = TRUE, xnudge = 0, ynudge = 0)
x11();plot_grid(pp1)

## Differential gene expression in 'theca'
Idents(object = theca) <- "geno"
avg.theca <- log1p(x = AverageExpression(object = theca, verbose = FALSE)$RNA)
avg.theca$gene <- rownames(x = avg.theca)
Idents(object = theca) <- "celltype.geno"
theca$celltype <- Idents(object = theca)
theca$celltype.geno <- paste(Idents(object = theca), theca$geno, sep = "_")
head(theca$celltype.geno)
Idents(object = theca) <- "celltype.geno"
KO.response.theca <- FindMarkers(object = theca, ident.1 = c("celltype.geno_KO"), ident.2 = c("celltype.geno_WT"), verbose = FALSE)
head(x = KO.response.theca, n = 15)
genes.to.label.10 = c("Gm8730", "Mt1", "Mt2", "Glul", "Sphk1", "Zbtb16", "B4galt5", "Hsd17b11", "Hipk2", "Wbp1l", "Susd3", "Smarca1", "Cln5", "Zcchc24", "Pla2g15")
table(KO.response.theca$p_val_adj<0.05)
sig_genes <- KO.response.theca[KO.response.theca$p_val_adj<0.05,]  
sig_genes$gene <- rownames(sig_genes) 
avg.theca <- avg.theca %>% left_join(sig_genes) 
avg.theca$status <- "NoSig"
avg.theca$status[avg.theca$avg_logFC>(0) & avg.theca$p_val_adj<0.05] <- "Up"   
avg.theca$status[avg.theca$avg_logFC<(0) & avg.theca$p_val_adj<0.05] <- "Down"
avg.theca$status <- factor(avg.theca$status, levels = c("NoSig", "Up", "Down"))
avg.theca <- avg.theca %>% arrange(status) 
avg.theca$status <- factor(avg.theca$status, levels = c( "Up","NoSig", "Down")) 
rownames(avg.theca) <- avg.theca$gene 
pp1 <- ggplot(avg.theca, aes(WT, KO, color = status)) + geom_point() + ggtitle("Theca")+
  scale_color_manual(values=c("red",  "grey", "blue")) 
pp1 <- LabelPoints(plot = pp1, points = genes.to.label.10, repel = TRUE, xnudge = 0, ynudge = 0)
x11();plot_grid(pp1)

## Differential gene expression in 'endo'
Idents(object = endo) <- "geno"
avg.endo <- log1p(x = AverageExpression(object = endo, verbose = FALSE)$RNA)
avg.endo$gene <- rownames(x = avg.endo)
Idents(object = endo) <- "celltype.geno"
endo$celltype <- Idents(object = endo)
endo$celltype.geno <- paste(Idents(object = endo), endo$geno, sep = "_")
head(endo$celltype.geno)
Idents(object = endo) <- "celltype.geno"
KO.response.endo <- FindMarkers(object = endo, ident.1 = c("celltype.geno_KO"), ident.2 = c("celltype.geno_WT"), verbose = FALSE)
head(x = KO.response.endo, n = 15)
genes.to.label.7 = c("Gm8730", "Pdlim5", "Tnfrsf12a", "Tbc1d1", "Actr3", "Gnl3", "Lrrfip1", "Lmna", "Actg1", "Rnd1", "Sp100", "Hif1a", "Cav1")
table(KO.response.endo$p_val_adj<0.05)
sig_genes <- KO.response.endo[KO.response.endo$p_val_adj<0.05,]  
sig_genes$gene <- rownames(sig_genes) 
avg.endo <- avg.endo %>% left_join(sig_genes) 
avg.endo$status <- "NoSig"
avg.endo$status[avg.endo$avg_logFC>(0) & avg.endo$p_val_adj<0.05] <- "Up"   
avg.endo$status[avg.endo$avg_logFC<(0) & avg.endo$p_val_adj<0.05] <- "Down"
avg.endo$status <- factor(avg.endo$status, levels = c("NoSig", "Up", "Down"))
avg.endo <- avg.endo %>% arrange(status) 
avg.endo$status <- factor(avg.endo$status, levels = c( "Up","NoSig", "Down")) 
rownames(avg.endo) <- avg.endo$gene 
pp1 <- ggplot(avg.endo, aes(WT, KO, color = status)) + geom_point() + ggtitle("Endo")+
  scale_color_manual(values=c("red",  "grey", "blue")) 
pp1 <- LabelPoints(plot = pp1, points = genes.to.label.7, repel = TRUE, xnudge = 0, ynudge = 0)
x11();plot_grid(pp1)

## Differential gene expression in 'earlygc'
Idents(object = earlygc) <- "geno"
avg.earlygc <- log1p(x = AverageExpression(object = earlygc, verbose = FALSE)$RNA)
avg.earlygc$gene <- rownames(x = avg.earlygc)
Idents(object = earlygc) <- "celltype.geno"
earlygc$celltype <- Idents(object = earlygc)
earlygc$celltype.geno <- paste(Idents(object = earlygc), earlygc$geno, sep = "_")
head(earlygc$celltype.geno)
Idents(object = earlygc) <- "celltype.geno"
KO.response.earlygc <- FindMarkers(object = earlygc, ident.1 = c("celltype.geno_KO"), ident.2 = c("celltype.geno_WT"), verbose = FALSE)
head(x = KO.response.earlygc, n = 15)
genes.to.label.egc = c("Gm8730", "Grb14", "Inhbb", "Prss23", "Pik3ip1", "Cdc42ep3", "Kcnk2", "Pik3r1", "Bmp2", "Tmx1", "Ccnd2", "Sprr2g", "Bcat1", "Tdrd5", "Asb4")
table(KO.response.earlygc$p_val_adj<0.05)
sig_genes <- KO.response.earlygc[KO.response.earlygc$p_val_adj<0.05,]  
sig_genes$gene <- rownames(sig_genes) 
avg.earlygc <- avg.earlygc %>% left_join(sig_genes) 
avg.earlygc$status <- "NoSig"
avg.earlygc$status[avg.earlygc$avg_logFC>(0) & avg.earlygc$p_val_adj<0.05] <- "Up"   
avg.earlygc$status[avg.earlygc$avg_logFC<(0) & avg.earlygc$p_val_adj<0.05] <- "Down"
avg.earlygc$status <- factor(avg.earlygc$status, levels = c("NoSig", "Up", "Down"))
avg.earlygc <- avg.earlygc %>% arrange(status) 
avg.earlygc$status <- factor(avg.earlygc$status, levels = c( "Up","NoSig", "Down")) 
rownames(avg.earlygc) <- avg.earlygc$gene 
pp1 <- ggplot(avg.earlygc, aes(WT, KO, color = status)) + geom_point() + ggtitle("Early GC")+
  scale_color_manual(values=c("red",  "grey", "blue")) 
pp1 <- LabelPoints(plot = pp1, points = genes.to.label.egc, repel = TRUE, xnudge = 0, ynudge = 0)
x11();plot_grid(pp1)

## Differential gene expression in 'int'
Idents(object = int) <- "geno"
avg.int <- log1p(x = AverageExpression(object = int, verbose = FALSE)$RNA)
avg.int$gene <- rownames(x = avg.int)
Idents(object = int) <- "celltype.geno"
int$celltype <- Idents(object = int)
int$celltype.geno <- paste(Idents(object = int), int$geno, sep = "_")
head(int$celltype.geno)
Idents(object = int) <- "celltype.geno"
KO.response.int <- FindMarkers(object = int, ident.1 = c("celltype.geno_KO"), ident.2 = c("celltype.geno_WT"), verbose = FALSE)
head(x = KO.response.int, n = 15)
genes.to.label.int = c("Gm8730", "Ptpn2", "Slc2a1", "Pmp22", "Tnfaip6", "Rdh10", "Pde10a", "Crem", "Angpt4", "Cebpb", "Gm12840", "Scd2", "Btg3", "Pid1", "Tns1")
table(KO.response.int$p_val_adj<0.05)
sig_genes <- KO.response.int[KO.response.int$p_val_adj<0.05,]  
sig_genes$gene <- rownames(sig_genes) 
avg.int <- avg.int %>% left_join(sig_genes) 
avg.int$status <- "NoSig"
avg.int$status[avg.int$avg_logFC>(0) & avg.int$p_val_adj<0.05] <- "Up"   
avg.int$status[avg.int$avg_logFC<(0) & avg.int$p_val_adj<0.05] <- "Down"
avg.int$status <- factor(avg.int$status, levels = c("NoSig", "Up", "Down"))
avg.int <- avg.int %>% arrange(status) 
avg.int$status <- factor(avg.int$status, levels = c( "Up","NoSig", "Down")) 
rownames(avg.int) <- avg.int$gene 
pp1 <- ggplot(avg.int, aes(WT, KO, color = status)) + geom_point() + ggtitle("Interstitial")+
  scale_color_manual(values=c("red",  "grey", "blue")) 
pp1 <- LabelPoints(plot = pp1, points = genes.to.label.int, repel = TRUE, xnudge = 0, ynudge = 0)
x11();plot_grid(pp1)

save.image(file = "D:/Single-cell_ovary_hCG/Cell_Rep/ovary_Seurat.RData")

## Plotting for manuscript revision
DimPlot(object = combined, reduction = "umap", label = TRUE, label.size = 10)
combined$geno <- factor(x = combined$geno, levels = c('WT', 'KO'))
FeaturePlot(object = combined, features = c("Ptgs2"), split.by = "geno", max.cutoff = 3, 
            cols = c("grey", "red"))
FeaturePlot(object = combined, features = c("Ptger2"), split.by = "geno", max.cutoff = 3, 
            cols = c("grey", "red"))
FeaturePlot(object = combined, features = c("Il1b"), split.by = "geno", max.cutoff = 3, 
            cols = c("grey", "red"))
FeaturePlot(object = combined, features = c("Il6"), split.by = "geno", max.cutoff = 3, 
            cols = c("grey", "red"))
FeaturePlot(object = combined, features = c("Il11"), split.by = "geno", max.cutoff = 3, 
            cols = c("grey", "red"))
FeaturePlot(object = combined, features = c("Cxcl12"), split.by = "geno", max.cutoff = 3, 
            cols = c("grey", "red"))
FeaturePlot(object = combined, features = c("Abcc4"), split.by = "geno", max.cutoff = 3, 
            cols = c("grey", "red"))
FeaturePlot(object = combined, features = c("Il11"), cols = c("grey", "red"))
FeaturePlot(object = combined, features = c("Ptger2"), cols = c("grey", "red"))
FeaturePlot(object = combined, features = c("Abcc4"), cols = c("grey", "red"))
FeaturePlot(object = combined, features = c("Il1b"), cols = c("grey", "red"))
FeaturePlot(object = combined, features = c("Notch3"), cols = c("grey", "red"))
FeaturePlot(object = combined, features = c("Itga2b"), cols = c("grey", "red"))
FeaturePlot(object = combined, features = c("Selp"), cols = c("grey", "red"))

FeaturePlot(object = ovul.mGCs, features = c("Nfkb1", "Rela"), split.by = "geno", max.cutoff = 3, 
            cols = c("grey", "red"))
FeaturePlot(object = cells.GCs, features = c("Pgr"), pt.size = 1,
            cols = c("grey", "red"))
FeaturePlot(object = ovul.mGCs, features = c("Pgr"), pt.size = 1,
            cols = c("grey", "red"))
FeaturePlot(object = ovul.mGCs, features = c("Nfkb1", "Rela"), dims = c(1, 2), cells = NULL,
            cols = c("lightgrey", "blue"), pt.size = 1, order = FALSE,
            min.cutoff = 0.5, max.cutoff = 3, reduction = NULL,
            split.by = "geno", shape.by = NULL, blend = FALSE,
            blend.threshold = 0.5, label = FALSE, label.size = 4,
            ncol = NULL, combine = TRUE, coord.fixed = FALSE, by.col = TRUE)
FeaturePlot(object = ovul.mGCs, features = c("Chuk", "Ikbkb"), dims = c(1, 2), cells = NULL,
            cols = c("lightgrey", "blue"), pt.size = 1, order = FALSE,
            min.cutoff = 0.5, max.cutoff = 3, reduction = NULL,
            split.by = "geno", shape.by = NULL, blend = FALSE,
            blend.threshold = 0.5, label = FALSE, label.size = 4,
            ncol = NULL, combine = TRUE, coord.fixed = FALSE, by.col = TRUE)
FeatureScatter(object = ovul.mGCs, feature1 = "Nfkbia", feature2 = "Pgr", split.by = "geno")

table(ovul.mGCs$celltype.geno)
ovul.mGCs.Pgr <- subset(x = ovul.mGCs, subset = Pgr > 1)

ovul.mGCs.WT <- subset(x = ovul.mGCs, ident = c("celltype.geno_WT"))
ovul.mGCs.WT
FeatureScatter(object = ovul.mGCs.WT, feature1 = "Nfkbia", feature2 = "Pgr")
FeatureScatter(object = ovul.mGCs.WT, feature1 = "Nfkbia", feature2 = "Ptgs2")
ovul.mGCs.WT.Pgr <- subset(x = ovul.mGCs.WT, subset = Pgr > 1)
ovul.mGCs.WT.Pgr 
FeatureScatter(object = ovul.mGCs.WT.Pgr, feature1 = "Nfkbia", feature2 = "Pgr", pt.size = 1)

ovul.mGCs.KO <- subset(x = ovul.mGCs, ident = c("celltype.geno_KO"))
ovul.mGCs.KO
FeatureScatter(object = ovul.mGCs.KO, feature1 = "Nfkbia", feature2 = "Pgr")
FeatureScatter(object = ovul.mGCs.KO, feature1 = "Nfkbia", feature2 = "Ptgs2", pt.size = 1)
ovul.mGCs.KO.Pgr <- subset(x = ovul.mGCs.KO, subset = Pgr > 1)
ovul.mGCs.KO.Pgr 
FeatureScatter(object = ovul.mGCs.KO.Pgr, feature1 = "Nfkbia", feature2 = "Pgr", pt.size = 1)

VlnPlot(object = ovul.mGCs, features = c("Pgr"),  
        group.by = "geno", pt.size = 1, combine = FALSE)

VlnPlot(object = ovul.mGCs, features = c("Nfkb1"),  
        group.by = "geno", pt.size = 1, combine = FALSE)
VlnPlot(object = ovul.mGCs, features = c("Rela"),  
        group.by = "geno", pt.size = 1, combine = FALSE)
VlnPlot(object = ovul.mGCs, features = c("Nfkbia"),  
        group.by = "geno", pt.size = 1, combine = FALSE)
VlnPlot(object = ovul.mGCs, features = c("Chuk"),  
        group.by = "geno", pt.size = 1, combine = FALSE)
VlnPlot(object = ovul.mGCs, features = c("Ikbkb"),  
        group.by = "geno", pt.size = 1, combine = FALSE)
VlnPlot(object = ovul.mGCs, features = c("Ikbkg"),  
        group.by = "geno", pt.size = 1, combine = FALSE)

ovul.mGCs.WT <- subset(x = ovul.mGCs, ident = c("celltype.geno_WT"))
ovul.mGCs.WT
FeatureScatter(object = ovul.mGCs.WT, feature1 = "Nfkbia", feature2 = "Pgr")
FeatureScatter(object = ovul.mGCs.WT, feature1 = "Nfkbia", feature2 = "Ptgs2")
ovul.mGCs.WT.Pgr <- subset(x = ovul.mGCs.WT, subset = Pgr > 1)
ovul.mGCs.WT.Pgr 
FeatureScatter(object = ovul.mGCs.WT.Pgr, feature1 = "Nfkbia", feature2 = "Pgr", pt.size = 1)

ovul.mGCs.KO <- subset(x = ovul.mGCs, ident = c("celltype.geno_KO"))
ovul.mGCs.KO
FeatureScatter(object = ovul.mGCs.KO, feature1 = "Nfkbia", feature2 = "Pgr")
FeatureScatter(object = ovul.mGCs.KO, feature1 = "Nfkbia", feature2 = "Ptgs2", pt.size = 1)
ovul.mGCs.KO.Pgr <- subset(x = ovul.mGCs.KO, subset = Pgr > 1)
ovul.mGCs.KO.Pgr 
FeatureScatter(object = ovul.mGCs.KO.Pgr, feature1 = "Nfkbia", feature2 = "Pgr", pt.size = 1)

VlnPlot(object = ovul.mGCs, features = c("Pgr"),  
        group.by = "geno", pt.size = 1, combine = FALSE)

VlnPlot(object = ovul.mGCs, features = c("Nfkb1"),  
        group.by = "geno", pt.size = 1, combine = FALSE)
VlnPlot(object = ovul.mGCs, features = c("Rela"),  
        group.by = "geno", pt.size = 1, combine = FALSE)
VlnPlot(object = ovul.mGCs, features = c("Nfkbia"),  
        group.by = "geno", pt.size = 1, combine = FALSE)
VlnPlot(object = ovul.mGCs, features = c("Chuk"),  
        group.by = "geno", pt.size = 1, combine = FALSE)
VlnPlot(object = ovul.mGCs, features = c("Ikbkb"),  
        group.by = "geno", pt.size = 1, combine = FALSE)
VlnPlot(object = ovul.mGCs, features = c("Ikbkg"),  
        group.by = "geno", pt.size = 1, combine = FALSE)

immune$geno <- factor(x = immune$geno, levels = c('WT', 'KO'))
FeaturePlot(object = immune, features = c("Il1b", "Tnf"), split.by = "geno", 
            cols = c("grey", "red"))
VlnPlot(object = immune, features = c("Il6"),  
        group.by = "geno", pt.size = 1, combine = FALSE)
VlnPlot(object = immune, features = c("Tnf"),  
        group.by = "geno", pt.size = 1, combine = FALSE)

