library(monocle)
library(devtools)
library(Seurat)
library(cowplot)
library(ggplot2)
library(umap)
library(dplyr)
library(xlsx)
library(Matrix)
library(DDRTree)
library(pheatmap)
library(reshape2)
library(clues)
library(stringr)
library(argparse)
library(gridExtra)
library(scales)
library(viridis)
library(data.table)
library(readr)
library(reticulate)
library(VGAM)
library(ggrastr)
library(L1Graph)
library(tidyverse)

## Setup the Seurat Object
# Setup the Seurat Object-WT
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
wt <- CreateSeuratObject(counts = wt6h.data, project = "6h_WT", min.cells = 5, min.features = 200)
wt$stim <- "WT"
wt <- subset(x = wt, subset = nFeature_RNA > 500)
wt <- NormalizeData(object = wt, verbose = FALSE)
wt <- FindVariableFeatures(object = wt, selection.method = "vst", nfeatures = 2000)
# Setup the Seurat Object-KO
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
ko <- CreateSeuratObject(counts = ko6h.data, project = "6h_KO", min.cells = 5, min.features = 200)
ko$stim <- "KO"
ko <- subset(x = ko, subset = nFeature_RNA > 500)
ko <- NormalizeData(object = ko, verbose = FALSE)
ko <- FindVariableFeatures(object = ko, selection.method = "vst", nfeatures = 2000)
# Find Anchors
anchors <- FindIntegrationAnchors(object.list = list(wt, ko), dims = 1:20)
combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(object = combined) <- "integrated"
# visualization and clustering
combined <- ScaleData(object = combined, verbose = FALSE)
combined <- RunPCA(object = combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(object = combined, reduction = "pca", dims = 1:20)
combined <- FindNeighbors(object = combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)
DefaultAssay(object = combined) <- "RNA"
# Plotting
p1 <- DimPlot(object = combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(object = combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
DimPlot(object = combined, reduction = "umap", split.by = "stim")
# Sub-clustering (granulosa cells)
cells.GCs <- subset(x = combined, idents = c("0", "2", "3", "4", "5", "6", "11"))
cells.GCs
cells.GCs$stim <- factor(x = cells.GCs$stim, levels = c('WT', 'KO'))
# Standard workflow for visualization and clustering
cells.GCs <- ScaleData(object = cells.GCs, verbose = FALSE)
cells.GCs <- RunPCA(object = cells.GCs, npcs = 30, verbose = FALSE)
cells.GCs <- RunUMAP(object = cells.GCs, reduction = "pca", dims = 1:25)
cells.GCs <- FindNeighbors(object = cells.GCs, reduction = "pca", dims = 1:25)
cells.GCs <- FindClusters(cells.GCs, resolution = 0.7)
# Plotting
p1 <- DimPlot(object = cells.GCs, reduction = "umap", group.by = "stim", label.size = 5)
p2 <- DimPlot(object = cells.GCs, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1)
plot_grid(p1, p2)
DimPlot(object = cells.GCs, reduction = "umap", split.by = "stim", label.size = 5)
DimPlot(object = cells.GCs, reduction = "umap", label.size = 5)
FeaturePlot(object = cells.GCs, features = c("Ube2c"), split.by = "stim", 
            cols = c("grey", "red"), blend.threshold = 0.5, pt.size = 0.5)
FeaturePlot(object = cells.GCs, features = c("Ube2c"),  
            cols = c("grey", "red"), blend.threshold = 0.5, pt.size = 0.5)
# mural granulosa cells
mGCs <- subset(x = cells.GCs, idents = c("0", "2", "3", "4", "5", "6", "8", "9", "10"))
mGCs
mGCs$stim <- factor(x = mGCs$stim, levels = c('WT', 'KO'))
DimPlot(object = mGCs, reduction = "umap", split.by = "stim", label.size = 5)
FeaturePlot(object = mGCs, features = c("Nfkbia"), split.by = "stim", 
            cols = c("grey", "red"), blend.threshold = 0.5, pt.size = 1)
FeaturePlot(object = mGCs, features = c("Nfkbia"),  
            cols = c("grey", "red"), pt.size = 1)

## Monocle - WT
class(mGCs)
slotNames(mGCs)
head(mGCs@meta.data)

mGCWT <- SubsetData(mGCs, subset.name = "stim", accept.value = c("WT"), assay = 'RNA')
DimPlot(object = mGCWT, reduction = "umap", split.by = "stim", label.size = 5)

mGCKO <- SubsetData(mGCs, subset.name = "stim", accept.value = c("KO"), assay = 'RNA')
DimPlot(object = mGCKO, reduction = "umap", split.by = "stim", label.size = 5)

newimport <- function(otherCDS, import_all = FALSE) {
  if(class(otherCDS)[1] == 'Seurat') {
    requireNamespace("Seurat")
    data <- otherCDS@assays$RNA@counts
    
    if(class(data) == "data.frame") {
      data <- as(as.matrix(data), "sparseMatrix")
    }
    
    pd <- tryCatch( {
      pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
      pd
    }, 
    #warning = function(w) { },
    error = function(e) { 
      pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
      pd <- new("AnnotatedDataFrame", data = pData)
      
      message("This Seurat object doesn't provide any meta data");
      pd
    })
    
    # remove filtered cells from Seurat
    if(length(setdiff(colnames(data), rownames(pd))) > 0) {
      data <- data[, rownames(pd)]  
    }
    
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    lowerDetectionLimit <- 0
    
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    
    valid_data <- data[, row.names(pd)]
    
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)
    
    if(import_all) {
      if("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
        
        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS
        
      } else {
        # mist_list <- list(ident = ident) 
        mist_list <- otherCDS
      }
    } else {
      mist_list <- list()
    }
    
    if(1==1) {
      var.genes <- setOrderingFilter(monocle_cds, otherCDS@assays$RNA@var.features)
      
    }
    monocle_cds@auxClusteringData$seurat <- mist_list
    
  } else if (class(otherCDS)[1] == 'SCESet') {
    requireNamespace("scater")
    
    message('Converting the exprs data in log scale back to original scale ...')    
    data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset
    
    fd <- otherCDS@featureData
    pd <- otherCDS@phenoData
    experimentData = otherCDS@experimentData
    if("is.expr" %in% slotNames(otherCDS))
      lowerDetectionLimit <- otherCDS@is.expr
    else 
      lowerDetectionLimit <- 1
    
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    
    if(import_all) {
      # mist_list <- list(iotherCDS@sc3,
      #                   otherCDS@reducedDimension)
      mist_list <- otherCDS 
      
    } else {
      mist_list <- list()
    }
    
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)
    # monocle_cds@auxClusteringData$sc3 <- otherCDS@sc3
    # monocle_cds@auxOrderingData$scran <- mist_list
    
    monocle_cds@auxOrderingData$scran <- mist_list
    
  } else {
    stop('the object type you want to export to is not supported yet')
  }
  
  return(monocle_cds)
}

cdsWT <- newimport(mGCWT, import_all = TRUE)
class(cdsWT)
dim(exprs(cdsWT))
exprs(cdsWT)[1:5, 1:5]
dim(pData(cdsWT))
head(pData(cdsWT))
dim(fData(cdsWT))
head(fData(cdsWT))
head(cdsWT)

my_featWT <- fData(cdsWT)

my_cdsWT <- newCellDataSet(exprs(cdsWT),
                           phenoData = new("AnnotatedDataFrame", data = pData(cdsWT)),
                           featureData = new("AnnotatedDataFrame", data = my_featWT),
                           lowerDetectionLimit = 0.5,
                           expressionFamily = negbinomial.size())
my_cdsWT
my_cdsWT <- estimateSizeFactors(my_cdsWT)
my_cdsWT <- estimateDispersions(my_cdsWT)

my_cdsWT <- detectGenes(my_cdsWT, min_expr = 0.1)
head(fData(my_cdsWT))
summary(fData(my_cdsWT)$num_cells_expressed)
sum((exprs(my_cdsWT['Pgr',])))
sum((exprs(my_cdsWT['Ptgs2',])))
head(pData(my_cdsWT))
sum((exprs(my_cdsWT)[,"AAAGAACGTCGAGTGA_1"])>0)
summary(pData(my_cdsWT)$num_genes_expressed)

xpWT <- pData(my_cdsWT)$num_genes_expressed
xp_1WT <- (xpWT - mean(xpWT)) / sd(xpWT)
summary(xp_1WT)

dfpWT <- data.frame(xpWT = xp_1WT)
ggplot(dfpWT, aes(xpWT)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')
pData(my_cdsWT)$UMI <- Matrix::colSums(exprs(my_cdsWT))
ggplot(pData(my_cdsWT), aes(num_genes_expressed, UMI)) + geom_point()
disp_tableWT <- dispersionTable(my_cdsWT)
head(disp_tableWT)
table(disp_tableWT$mean_expression>=0.1)
unsup_clustering_genesWT <- subset(disp_tableWT, mean_expression >= 0.1)
my_cdsWT <- setOrderingFilter(my_cdsWT, unsup_clustering_genesWT$gene_id)
plot_ordering_genes(my_cdsWT)
plot_pc_variance_explained(my_cdsWT, return_all = FALSE)

my_cdsWT <- reduceDimension(my_cdsWT, max_components = 2, num_dim = 25,
                            reduction_method = 'tSNE', verbose = TRUE)
my_cdsWT <- clusterCells(my_cdsWT, num_clusters = 10)
plot_cell_clusters(my_cdsWT)
plot_cell_clusters(my_cdsWT, color_by = "stim")
head(pData(my_cdsWT))

my_vectorWT <- rep('no', nrow(pData(my_cdsWT)))
my_vectorWT[pData(my_cdsWT)$Cluster == 1] <- rep('yes', sum(pData(my_cdsWT)$Cluster == 1))
pData(my_cdsWT)$test <- my_vectorWT
head(pData(my_cdsWT))
length(unsup_clustering_genesWT$gene_id)
de_cluster_oneWT <- differentialGeneTest(my_cdsWT[unsup_clustering_genesWT$gene_id,],
                                         fullModelFormulaStr = '~test',
                                         cores = 8)
dim(de_cluster_oneWT)
head(pData(my_cdsWT))
de_cluster_oneWT %>% arrange(qval) %>% head()
plot_genes_jitter(my_cdsWT['Lyve1',], grouping = "Cluster")

## Constructing Single Cell Trajectories
expressed_genesWT <- row.names(subset(fData(my_cdsWT), num_cells_expressed >= 10))

clustering_DEG_genesWT <- differentialGeneTest(my_cdsWT,
                                               fullModelFormulaStr = '~Cluster',
                                               cores = 8)
dim(clustering_DEG_genesWT)
clustering_DEG_genesWT %>% arrange(qval) %>% head()

my_ordering_genesWT <- row.names(clustering_DEG_genesWT)[order(clustering_DEG_genesWT$qval)][1:1000]
my_cdsWT <- setOrderingFilter(my_cdsWT, ordering_genes = my_ordering_genesWT)
my_cdsWT <- reduceDimension(my_cdsWT, method = 'DDRTree', norm_method = 'log')

my_cdsWT <- orderCells(my_cdsWT)
plot_cell_trajectory(my_cdsWT, color_by = "Cluster")
plot_cell_trajectory(my_cdsWT)

# pseudotime is now a column in the phenotypic data as well as the cell state
head(pData(my_cdsWT))
my_pseudotime_deWT <- differentialGeneTest(my_cdsWT, fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                           cores = 4)
my_pseudotime_deWT %>% arrange(qval) %>% head()

plot_cell_trajectory(my_cdsWT, markers="Pgr", color_by = "Pseudotime")
plot_cell_trajectory(my_cdsWT, markers="Esr2", color_by = "Cluster", cell_size = 0.1)
plot_cell_trajectory(my_cdsWT, markers="Lyve1", color_by = "Cluster", cell_size = 0.1)
plot_cell_trajectory(my_cdsWT, markers="Ptgs2", color_by = "Cluster")
plot_cell_trajectory(my_cdsWT, markers="Foxl2", color_by = "Cluster", cell_size = 0.1)
plot_cell_trajectory(my_cdsWT, markers="Esr2", color_by = "Cluster")
head(my_pseudotime_deWT)

# save the top 6 genes
my_pseudotime_deWT %>% arrange(qval) %>% head() %>% select(gene_short_name) -> my_pseudotime_geneWT
my_pseudotime_geneWT <- my_pseudotime_geneWT$gene_short_name

plot_genes_in_pseudotime(my_cdsWT[my_pseudotime_geneWT,])

# cluster the top 50 genes
my_pseudotime_deWT %>% arrange(qval) %>% head(50) %>% select(gene_short_name) -> gene_to_clusterWT
gene_to_clusterWT <- gene_to_clusterWT$gene_short_name

my_pseudotime_clusterWT <- plot_pseudotime_heatmap(my_cdsWT[gene_to_clusterWT,],
                                                   num_clusters = 9,
                                                   cores = 8,
                                                   show_rownames = TRUE,
                                                   return_heatmap = TRUE)
newdataWT <- data.frame(Pseudotime = seq(min(pData(my_cdsWT)$Pseudotime), 
                                         max(pData(my_cdsWT)$Pseudotime), length.out = 100))

my_clusterWT <- cutree(my_pseudotime_clusterWT$tree_row, 9)
my_clusterWT 

# genes in cluster 1
my_pseudotime_deWT[names(my_clusterWT[my_clusterWT == 1]),"gene_short_name"]
# genes in cluster 2
my_pseudotime_deWT[names(my_clusterWT[my_clusterWT == 2]),"gene_short_name"]
# genes in cluster 3
my_pseudotime_deWT[names(my_clusterWT[my_clusterWT == 3]),"gene_short_name"]
# genes in cluster 4
my_pseudotime_deWT[names(my_clusterWT[my_clusterWT == 4]),"gene_short_name"]
# genes in cluster 5
my_pseudotime_deWT[names(my_clusterWT[my_clusterWT == 5]),"gene_short_name"]
# genes in cluster 6
my_pseudotime_deWT[names(my_clusterWT[my_clusterWT == 6]),"gene_short_name"]
# genes in cluster 7
my_pseudotime_deWT[names(my_clusterWT[my_clusterWT == 7]),"gene_short_name"]
# genes in cluster 8
my_pseudotime_deWT[names(my_clusterWT[my_clusterWT == 8]),"gene_short_name"]
# genes in cluster 9
my_pseudotime_deWT[names(my_clusterWT[my_clusterWT == 9]),"gene_short_name"]

## Analyzing Branches in Single-Cell Trajectories
plot_cell_trajectory(my_cdsWT, x = 1, y = 2, color_by = "Pseudotime", 
                     markers="Pgr", cell_name_size = 1, cell_link_size = 1, 
                     show_tree=TRUE, show_backbone=TRUE, show_branch_points = TRUE, theta = 90)

# Set the root state (using Pgr)
plot_cell_trajectory(my_cdsWT, x = 1, y = 2, color_by = "State", markers="Pgr", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, 
                     show_backbone=TRUE, show_branch_points = FALSE, theta = 90)
my_cdsWT <- orderCells(my_cdsWT, root_state = 3)

pp1 <- plot_cell_trajectory(my_cdsWT, x = 1, y = 2, color_by = "Pseudotime", 
                            markers="Pgr",
                            show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp2 <- plot_cell_trajectory(my_cdsWT, x = 1, y = 2, color_by = "Pseudotime", 
                            markers="Foxo1", cell_name_size = 1, cell_link_size = 1, 
                            show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp3 <- plot_cell_trajectory(my_cdsWT, x = 1, y = 2, color_by = "Pseudotime", 
                            markers="Amh", cell_name_size = 1, cell_link_size = 1, 
                            show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp4 <- plot_cell_trajectory(my_cdsWT, x = 1, y = 2, color_by = "Pseudotime", 
                            markers="Lyve1", cell_name_size = 1, cell_link_size = 1, 
                            show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp5 <- plot_cell_trajectory(my_cdsWT, x = 1, y = 2, color_by = "Pseudotime", 
                            markers="Snap25", cell_name_size = 1, cell_link_size = 1, 
                            show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp6 <- plot_cell_trajectory(my_cdsWT, x = 1, y = 2, color_by = "Pseudotime", 
                            markers="Gsta4", cell_name_size = 1, cell_link_size = 1, 
                            show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp7 <- plot_cell_trajectory(my_cdsWT, x = 1, y = 2, color_by = "Pseudotime", 
                            markers="Cited4", cell_name_size = 1, cell_link_size = 1, 
                            show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp8 <- plot_cell_trajectory(my_cdsWT, x = 1, y = 2, color_by = "Pseudotime", 
                            markers="S100a10", cell_name_size = 1, cell_link_size = 1, 
                            show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp9 <- plot_cell_trajectory(my_cdsWT, x = 1, y = 2, color_by = "Pseudotime", 
                            markers="Star", cell_name_size = 1, cell_link_size = 1, 
                            show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp10 <- plot_cell_trajectory(my_cdsWT, x = 1, y = 2, color_by = "Pseudotime", 
                             markers="Cyp11a1", cell_name_size = 1, cell_link_size = 1, 
                             show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp11 <- plot_cell_trajectory(my_cdsWT, x = 1, y = 2, color_by = "Pseudotime", 
                             markers="Ptgs2", cell_name_size = 1, cell_link_size = 1, 
                             show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp12 <- plot_cell_trajectory(my_cdsWT, x = 1, y = 2, color_by = "Pseudotime", 
                             markers="Nfkbia", cell_name_size = 1, cell_link_size = 1,
                             show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')

plot_cell_trajectory(my_cdsWT, x = 1, y = 2, color_by = "Pseudotime", 
                     cell_name_size = 1, cell_link_size = 1, 
                     show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 90)
plot_grid(pp1, pp2, pp3, pp4)
plot_grid(pp1, pp2, pp3, pp4, pp5, pp6)
plot_grid(pp7, pp8, pp9, pp10, pp11, pp12)
plot_grid(pp1, pp4, pp11, pp12)
plot_grid(pp1, pp2, pp3, pp4, pp5, pp6, pp7, pp8, pp9, pp10, pp11, pp12)

plot_cell_clusters(my_cdsWT, color_by="Cluster")
plot_cell_clusters(my_cdsWT, color_by="Pseudotime")

plot_cell_clusters(my_cdsWT, markers="Pgr")
plot_cell_clusters(my_cdsWT, markers="Ptgs2")
plot_cell_clusters(my_cdsWT, markers="Nfkbia")

GM_stateWT <- function(my_cdsWT){
  if (length(unique(pData(my_cdsWT)$State)) > 1){
    T0_counts <- table(pData(my_cdsWT)$State, pData(my_cdsWT)$Hours)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
plot_cell_trajectory(my_cdsWT, x = 1, y = 2, color_by = "Pseudotime", use_color_gradient = FALSE, show_branch_points = FALSE, theta = 90)
plot_cell_trajectory(my_cdsWT, x = 1, y = 2, color_by = "State", use_color_gradient = FALSE, show_branch_points = FALSE, theta = 90) +
  facet_wrap(~State, nrow = 1)

blast_genesWT <- row.names(subset(fData(my_cdsWT),
                                  gene_short_name %in% c("Pgr", "Ptgs2", "Nfkbia", "Lyve1", "Snap25", "Amh")))
plot_genes_jitter(my_cdsWT[blast_genesWT,],
                  grouping = "State",
                  min_expr = 0.1)
plot_genes_jitter(my_cdsWT[blast_genesWT,],
                  grouping = "Pseudotime",
                  min_expr = 0.1)

## Monocle - KO
cdsKO <- newimport(mGCKO, import_all = TRUE)
class(cdsKO)
dim(exprs(cdsKO))
exprs(cdsKO)[1:5, 1:5]
dim(pData(cdsKO))
head(pData(cdsKO))
dim(fData(cdsKO))
head(fData(cdsKO))
head(cdsKO)
my_featKO <- fData(cdsKO)
my_cdsKO <- newCellDataSet(exprs(cdsKO),
                           phenoData = new("AnnotatedDataFrame", data = pData(cdsKO)),
                           featureData = new("AnnotatedDataFrame", data = my_featKO),
                           lowerDetectionLimit = 0.5,
                           expressionFamily = negbinomial.size())
my_cdsKO
my_cdsKO <- estimateSizeFactors(my_cdsKO)
my_cdsKO <- estimateDispersions(my_cdsKO)
my_cdsKO <- detectGenes(my_cdsKO, min_expr = 0.1)
head(fData(my_cdsKO))
summary(fData(my_cdsKO)$num_cells_expressed)
sum((exprs(my_cdsKO['Pgr',])))
sum((exprs(my_cdsKO['Ptgs2',])))
head(pData(my_cdsKO))
sum((exprs(my_cdsKO)[,"AAAGAACGTCGAGTGA_1"])>0)
summary(pData(my_cdsKO)$num_genes_expressed)
xpKO <- pData(my_cdsKO)$num_genes_expressed
xp_1KO <- (xpKO - mean(xpKO)) / sd(xpKO)
summary(xp_1KO)
dfpKO <- data.frame(xpKO = xp_1KO)
ggplot(dfpKO, aes(xpKO)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')
pData(my_cdsKO)$UMI <- Matrix::colSums(exprs(my_cdsKO))
ggplot(pData(my_cdsKO), aes(num_genes_expressed, UMI)) + geom_point()
disp_tableKO <- dispersionTable(my_cdsKO)
head(disp_tableKO)
table(disp_tableKO$mean_expression>=0.1)
unsup_clustering_genesKO <- subset(disp_tableKO, mean_expression >= 0.1)
my_cdsKO <- setOrderingFilter(my_cdsKO, unsup_clustering_genesKO$gene_id)
plot_ordering_genes(my_cdsKO)
plot_pc_variance_explained(my_cdsKO, return_all = FALSE)

my_cdsKO <- reduceDimension(my_cdsKO, max_components = 2, num_dim = 25,
                            reduction_method = 'tSNE', verbose = TRUE)
my_cdsKO <- clusterCells(my_cdsKO, num_clusters = 10)
plot_cell_clusters(my_cdsKO)
plot_cell_clusters(my_cdsKO, color_by = "stim")
head(pData(my_cdsKO))

my_vectorKO <- rep('no', nrow(pData(my_cdsKO)))
my_vectorKO[pData(my_cdsKO)$Cluster == 1] <- rep('yes', sum(pData(my_cdsKO)$Cluster == 1))
pData(my_cdsKO)$test <- my_vectorKO
head(pData(my_cdsKO))
length(unsup_clustering_genesKO$gene_id)
de_cluster_oneKO <- differentialGeneTest(my_cdsKO[unsup_clustering_genesKO$gene_id,],
                                         fullModelFormulaStr = '~test',
                                         cores = 8)
dim(de_cluster_oneKO)
head(pData(my_cdsKO))
de_cluster_oneKO %>% arrange(qval) %>% head()
plot_genes_jitter(my_cdsKO['Ptgs2',], grouping = "Cluster")

## Constructing Single Cell Trajectories
expressed_genesKO <- row.names(subset(fData(my_cdsKO), num_cells_expressed >= 10))

clustering_DEG_genesKO <- differentialGeneTest(my_cdsKO,
                                               fullModelFormulaStr = '~Cluster',
                                               cores = 8)
dim(clustering_DEG_genesKO)
clustering_DEG_genesKO %>% arrange(qval) %>% head()

my_ordering_genesKO <- row.names(clustering_DEG_genesKO)[order(clustering_DEG_genesKO$qval)][1:1000]
my_cdsKO <- setOrderingFilter(my_cdsKO, ordering_genes = my_ordering_genesKO)
my_cdsKO <- reduceDimension(my_cdsKO, method = 'DDRTree', norm_method = 'log')

my_cdsKO <- orderCells(my_cdsKO)
plot_cell_trajectory(my_cdsKO, color_by = "Cluster")
plot_cell_trajectory(my_cdsKO)

head(pData(my_cdsKO))
my_pseudotime_deKO <- differentialGeneTest(my_cdsKO, fullModelFormulaStr = "~sm.ns(Pseudotime)", 
                                           cores = 4)
my_pseudotime_deKO %>% arrange(qval) %>% head()
plot_cell_clusters(my_cdsKO, markers="Pgr")
plot_cell_clusters(my_cdsKO, markers="Nfkbia")
plot_cell_clusters(my_cdsKO, markers="Ptgs2")
plot_cell_clusters(my_cdsKO, color_by = "Pseudotime")
plot_cell_trajectory(my_cdsKO, markers="Pgr", color_by = "Pseudotime")
plot_cell_trajectory(my_cdsKO, markers="Esr2", color_by = "Cluster", cell_size = 0.1)
plot_cell_trajectory(my_cdsKO, markers="Lyve1", color_by = "Cluster", cell_size = 0.1)
plot_cell_trajectory(my_cdsKO, markers="Ptgs2", color_by = "Cluster")
plot_cell_trajectory(my_cdsKO, markers="Foxl2", color_by = "Cluster", cell_size = 0.1)
plot_cell_trajectory(my_cdsKO, markers="Esr2", color_by = "Cluster")
head(my_pseudotime_deKO)

# top 6 genes
my_pseudotime_deKO %>% arrange(qval) %>% head() %>% select(gene_short_name) -> my_pseudotime_geneKO
my_pseudotime_geneKO <- my_pseudotime_geneKO$gene_short_name
plot_genes_in_pseudotime(my_cdsKO[my_pseudotime_geneKO,])

# cluster the top 50 genes 
my_pseudotime_deKO %>% arrange(qval) %>% head(50) %>% select(gene_short_name) -> gene_to_clusterKO
gene_to_clusterKO <- gene_to_clusterKO$gene_short_name
my_pseudotime_clusterKO <- plot_pseudotime_heatmap(my_cdsKO[gene_to_clusterKO,],
                                                   num_clusters = 9,
                                                   cores = 8,
                                                   show_rownames = TRUE,
                                                   return_heatmap = TRUE)

newdataKO <- data.frame(Pseudotime = seq(min(pData(my_cdsKO)$Pseudotime), 
                                         max(pData(my_cdsKO)$Pseudotime), length.out = 100))
my_clusterKO <- cutree(my_pseudotime_clusterKO$tree_row, 9)
my_clusterKO 
# genes in cluster 1
my_pseudotime_deKO[names(my_clusterKO[my_clusterKO == 1]),"gene_short_name"]
# genes in cluster 2
my_pseudotime_deKO[names(my_clusterKO[my_clusterKO == 2]),"gene_short_name"]
# genes in cluster 3
my_pseudotime_deKO[names(my_clusterKO[my_clusterKO == 3]),"gene_short_name"]
# genes in cluster 4
my_pseudotime_deKO[names(my_clusterKO[my_clusterKO == 4]),"gene_short_name"]
# genes in cluster 5
my_pseudotime_deKO[names(my_clusterKO[my_clusterKO == 5]),"gene_short_name"]
# genes in cluster 6
my_pseudotime_deKO[names(my_clusterKO[my_clusterKO == 6]),"gene_short_name"]
# genes in cluster 7
my_pseudotime_deKO[names(my_clusterKO[my_clusterKO == 7]),"gene_short_name"]
# genes in cluster 8
my_pseudotime_deKO[names(my_clusterKO[my_clusterKO == 8]),"gene_short_name"]
# genes in cluster 9
my_pseudotime_deKO[names(my_clusterKO[my_clusterKO == 9]),"gene_short_name"]

## Analyzing Branches in Single-Cell Trajectories
plot_cell_trajectory(my_cdsKO, x = 1, y = 2, color_by = "Pseudotime", 
                     markers="Pgr", cell_name_size = 1, cell_link_size = 1, 
                     show_tree=TRUE, show_backbone=TRUE, show_branch_points = TRUE, theta = 45)

# Set the root state (using Pgr)
plot_cell_trajectory(my_cdsKO, x = 1, y = 2, color_by = "State", markers="Pgr", 
                     cell_name_size = 1, cell_link_size = 1, show_tree=TRUE, 
                     show_backbone=TRUE, show_branch_points = FALSE, theta = 0)
my_cdsKO <- orderCells(my_cdsKO, root_state = 2)

pp1 <- plot_cell_trajectory(my_cdsKO, x = 1, y = 2, color_by = "Pseudotime", 
                            markers="Pgr",
                            show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 45) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp2 <- plot_cell_trajectory(my_cdsKO, x = 1, y = 2, color_by = "Pseudotime", 
                            markers="Foxo1", cell_name_size = 1, cell_link_size = 1, 
                            show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 45) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp3 <- plot_cell_trajectory(my_cdsKO, x = 1, y = 2, color_by = "Pseudotime", 
                            markers="Amh", cell_name_size = 1, cell_link_size = 1, 
                            show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 45) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp4 <- plot_cell_trajectory(my_cdsKO, x = 1, y = 2, color_by = "Pseudotime", 
                            markers="Lyve1", cell_name_size = 1, cell_link_size = 1, 
                            show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 45) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp5 <- plot_cell_trajectory(my_cdsKO, x = 1, y = 2, color_by = "Pseudotime", 
                            markers="Snap25", cell_name_size = 1, cell_link_size = 1, 
                            show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 45) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp6 <- plot_cell_trajectory(my_cdsKO, x = 1, y = 2, color_by = "Pseudotime", 
                            markers="Gsta4", cell_name_size = 1, cell_link_size = 1, 
                            show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 45) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp7 <- plot_cell_trajectory(my_cdsKO, x = 1, y = 2, color_by = "Pseudotime", 
                            markers="Cited4", cell_name_size = 1, cell_link_size = 1, 
                            show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 45) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp8 <- plot_cell_trajectory(my_cdsKO, x = 1, y = 2, color_by = "Pseudotime", 
                            markers="S100a10", cell_name_size = 1, cell_link_size = 1, 
                            show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 45) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp9 <- plot_cell_trajectory(my_cdsKO, x = 1, y = 2, color_by = "Pseudotime", 
                            markers="Star", cell_name_size = 1, cell_link_size = 1, 
                            show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 45) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp10 <- plot_cell_trajectory(my_cdsKO, x = 1, y = 2, color_by = "Pseudotime", 
                             markers="Cyp11a1", cell_name_size = 1, cell_link_size = 1, 
                             show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 45) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp <- plot_cell_trajectory(my_cdsKO, x = 1, y = 2, color_by = "Pseudotime", 
                           markers="Esr1", cell_name_size = 1, cell_link_size = 1, 
                           show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 45)
pp <- plot_cell_trajectory(my_cdsKO, x = 1, y = 2, color_by = "Pseudotime", 
                           markers="Esr2", cell_name_size = 1, cell_link_size = 1, 
                           show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 45)
pp11 <- plot_cell_trajectory(my_cdsKO, x = 1, y = 2, color_by = "Pseudotime", 
                             markers="Ptgs2", cell_name_size = 1, cell_link_size = 1, 
                             show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 45) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')
pp12 <- plot_cell_trajectory(my_cdsKO, x = 1, y = 2, color_by = "Pseudotime", 
                             markers="Nfkbia", cell_name_size = 1, cell_link_size = 1,
                             show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 45) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position='right')

plot_cell_trajectory(my_cdsKO, x = 1, y = 2, color_by = "Pseudotime", 
                     cell_name_size = 1, cell_link_size = 1, 
                     show_tree=TRUE, show_backbone=TRUE, show_branch_points = FALSE, theta = 45)
plot_grid(pp1, pp2, pp3, pp4)
plot_grid(pp1, pp2, pp3, pp4, pp5, pp6)
plot_grid(pp7, pp8, pp9, pp10, pp11, pp12)
plot_grid(pp1, pp4, pp11, pp12)
plot_grid(pp1, pp2, pp3, pp4, pp5, pp6, pp7, pp8, pp9, pp10, pp11, pp12)

plot_cell_clusters(my_cdsKO, color_by="Cluster")
plot_cell_clusters(my_cdsKO, color_by="Pseudotime")

plot_cell_clusters(my_cdsKO, markers="Pgr")
plot_cell_clusters(my_cdsKO, markers="Ptgs2")
plot_cell_clusters(my_cdsKO, markers="Nfkbia")

plot_cell_trajectory(my_cdsKO, markers = c("Pgr", "Ptgs2", "Nfkbia"), use_color_gradient = TRUE, theta = 45)

GM_stateKO <- function(my_cdsKO){
  if (length(unique(pData(my_cdsKO)$State)) > 1){
    T0_counts <- table(pData(my_cdsKO)$State, pData(my_cdsKO)$Hours)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
plot_cell_trajectory(my_cdsKO, x = 1, y = 2, color_by = "Pseudotime", use_color_gradient = FALSE, show_branch_points = FALSE, theta = 45)
plot_cell_trajectory(my_cdsKO, x = 1, y = 2, color_by = "State", use_color_gradient = FALSE, show_branch_points = FALSE, theta = 45) +
  facet_wrap(~State, nrow = 1)

blast_genes <- row.names(subset(fData(my_cdsKO),
                                gene_short_name %in% c("Pgr", "Ptgs2", "Nfkbia", "Lyve1", "Snap25", "Amh")))
plot_genes_jitter(my_cdsKO[blast_genes,],
                  grouping = "State",
                  min_expr = 0.1)
plot_genes_jitter(my_cdsKO[blast_genes,],
                  grouping = "Pseudotime",
                  min_expr = 0.1)

## plotting
# plotting WT
expressed_genesWT <-  row.names(subset(fData(my_cdsWT),
                                       num_cells_expressed >= 10))
filteredWT <- my_cdsWT[expressed_genesWT,]
my_genesWT <- row.names(subset(fData(filteredWT),
                               gene_short_name %in% c("Pgr", "Ptgs2", "Snap25", "Amh", "Lyve1")))
cds_subsetWT <- filteredWT[my_genesWT,]

branchTest(my_cdsWT, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch",
           reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)",
           branch_states = NULL, branch_point = 1, relative_expr = TRUE,
           cores = 8, branch_labels = NULL, verbose = FALSE)

plot.1 <- plot_genes_branched_pseudotime(cds_subsetWT, branch_point = 1, 
                                         color_by = "Pseudotime", 
                                         ncol = 1, cell_size = 1.5,
                                         panel_order = c("Pgr", "Ptgs2", "Snap25", "Amh", "Lyve1")
)

# plotting KO
expressed_genesKO <-  row.names(subset(fData(my_cdsKO),
                                       num_cells_expressed >= 10))
filteredKO <- my_cdsKO[expressed_genesKO,]
my_genesKO <- row.names(subset(fData(filteredKO),
                               gene_short_name %in% c("Pgr", "Ptgs2", "Snap25", "Amh", "Lyve1")))
cds_subsetKO <- filteredKO[my_genesKO,]
branchTest(my_cdsKO, fullModelFormulaStr = "~sm.ns(Pseudotime, df = 3)*Branch",
           reducedModelFormulaStr = "~sm.ns(Pseudotime, df = 3)",
           branch_states = NULL, branch_point = 1, relative_expr = TRUE,
           cores = 8, branch_labels = NULL, verbose = FALSE)

plot.2 <- plot_genes_branched_pseudotime(cds_subsetKO, branch_point = 1, color_by = "Pseudotime", 
                                         ncol = 1, 
                                         cell_size = 1.5,
                                         panel_order = c("Pgr", "Ptgs2", "Snap25", "Amh", "Lyve1"))


# From https://apps.fishandwhistle.net/archives/1344


scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}

facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}

CustomFacetWrap <- ggproto(
  "CustomFacetWrap", FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)


facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}

# manually reset scales
# "Pgr", "Ptgs2", "Snap25", "Amh", "Lyve1"

x11(); plot.1 + facet_wrap_custom(~feature_label, scales = "free", ncol = 1, scale_overrides = list(
  scale_override(1, scale_y_log10(breaks = c(1,3,10,30), limits= c(0.5,30))), #set Pgr
  scale_override(2, scale_y_log10(breaks = c(1,10,100), limits= c(0.5,100))), #set for Ptgs2
  scale_override(3, scale_y_log10(breaks = c(1,3,10,30), limits= c(0.5,30))), #Snap25
  scale_override(4, scale_y_log10(breaks = c(1,10,100), limits= c(0.5,100))),  #"Amh"
  scale_override(5, scale_y_log10(breaks = c(1,10,100), limits= c(0.5,100))) #Lyve1
)) +ggtitle("WT")
plot.1 + facet_wrap_custom(~feature_label, scales = "free", ncol = 1, scale_overrides = list(
  scale_override(1, scale_y_log10(breaks = c(1,3,10,30), limits= c(0.5,30))), #set Pgr
  scale_override(2, scale_y_log10(breaks = c(1,10,100), limits= c(0.5,100))), #set for Ptgs2
  scale_override(3, scale_y_log10(breaks = c(1,3,10,30), limits= c(0.5,30))), #Snap25
  scale_override(4, scale_y_log10(breaks = c(1,10,100), limits= c(0.5,100))),  #"Amh"
  scale_override(5, scale_y_log10(breaks = c(1,10,100), limits= c(0.5,100))) #Lyve1
)) +ggtitle("WT")

x11(); plot.2 + facet_wrap_custom(~feature_label, scales = "free", ncol = 1, scale_overrides = list(
  scale_override(1, scale_y_log10(breaks = c(1,3,10,30), limits= c(0.5,30))), #set Pgr
  scale_override(2, scale_y_log10(breaks = c(1,10,100), limits= c(0.5,100))), #set for Ptgs2
  scale_override(3, scale_y_log10(breaks = c(1,3,10,30), limits= c(0.5,30))), #Snap25
  scale_override(4, scale_y_log10(breaks = c(1,10,100), limits= c(0.5,100))),  #"Amh"
  scale_override(5, scale_y_log10(breaks = c(1,10,100), limits= c(0.5,100))) #Lyve1
))  +ggtitle("Esr2-PgrKO")
plot.2 + facet_wrap_custom(~feature_label, scales = "free", ncol = 1, scale_overrides = list(
  scale_override(1, scale_y_log10(breaks = c(1,3,10,30), limits= c(0.5,30))), #set Pgr
  scale_override(2, scale_y_log10(breaks = c(1,10,100), limits= c(0.5,100))), #set for Ptgs2
  scale_override(3, scale_y_log10(breaks = c(1,3,10,30), limits= c(0.5,30))), #Snap25
  scale_override(4, scale_y_log10(breaks = c(1,10,100), limits= c(0.5,100))),  #"Amh"
  scale_override(5, scale_y_log10(breaks = c(1,10,100), limits= c(0.5,100))) #Lyve1
))  +ggtitle("Esr2-PgrKO")