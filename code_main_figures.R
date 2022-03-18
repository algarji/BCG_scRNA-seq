library(dplyr)
library(tidyverse)
library(patchwork)
library(Seurat)
library(dittoSeq)
library(RColorBrewer)
library(EnhancedVolcano)
library(pheatmap)
library(dorothea)
library(tibble)
library(tidyr)
library(viper)
library(ggraph)
library(sceasy)
library(monocle)
library(escape)
memory.size()
memory.limit(size=56000)

ibcg <- readRDS("D:/sc/blood/IBCGSelectedCells_Regressed.rds")
oncotice <- readRDS("D:/sc/blood/OncoTSelectedCells_Regressed.rds")

n <- 36
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

pbmc.merged <- merge(ibcg, y = oncotice, add.cell.ids = c("ibcg", "oncotice"), project = "blood")
pbmc.merged <- NormalizeData(pbmc.merged)
pbmc.merged <- FindVariableFeatures(pbmc.merged, selection.method = "mean.var.plot")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
pbmc.merged <- CellCycleScoring(pbmc.merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
pbmc.merged <- ScaleData(pbmc.merged, vars.to.regress = c("S.Score", "G2M.Score", "HTO_maxID"), features = rownames(pbmc.merged), block.size = 100) 
pbmc.merged <- RunPCA(pbmc.merged, npcs = 30, verbose = FALSE)
pbmc.merged <- FindNeighbors(pbmc.merged, dims = 1:30)
pbmc.merged <- FindClusters(pbmc.merged, resolution = 0.5)
pbmc.merged <- RunUMAP(pbmc.merged, dims = 1:30)

pbmc.markers <- FindAllMarkers(pbmc.merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc.merged, features = top10$gene, group.colors = col_vector, size = 2, group.bar.height = 0.01) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))

cluster0.markers <- FindMarkers(pbmc.merged, ident.1 = 0, min.pct = 0.25)
cluster0 <- filter(cluster0.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster1.markers <- FindMarkers(pbmc.merged, ident.1 = 1, min.pct = 0.25)
cluster1 <- filter(cluster1.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster2.markers <- FindMarkers(pbmc.merged, ident.1 = 2, min.pct = 0.25)
cluster2 <- filter(cluster2.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster3.markers <- FindMarkers(pbmc.merged, ident.1 = 3, min.pct = 0.25)
cluster3 <- filter(cluster3.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster4.markers <- FindMarkers(pbmc.merged, ident.1 = 4, min.pct = 0.25)
cluster4 <- filter(cluster4.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster5.markers <- FindMarkers(pbmc.merged, ident.1 = 5, min.pct = 0.25)
cluster5 <- filter(cluster5.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster6.markers <- FindMarkers(pbmc.merged, ident.1 = 6, min.pct = 0.25)
cluster6 <- filter(cluster6.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster7.markers <- FindMarkers(pbmc.merged, ident.1 = 7, min.pct = 0.25)
cluster7 <- filter(cluster7.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster8.markers <- FindMarkers(pbmc.merged, ident.1 = 8, min.pct = 0.25)
cluster8 <- filter(cluster8.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster9.markers <- FindMarkers(pbmc.merged, ident.1 = 9, min.pct = 0.25)
cluster9 <- filter(cluster9.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster10.markers <- FindMarkers(pbmc.merged, ident.1 = 10, min.pct = 0.25)
cluster10 <- filter(cluster10.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster11.markers <- FindMarkers(pbmc.merged, ident.1 = 11, min.pct = 0.25)
cluster11 <- filter(cluster11.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster12.markers <- FindMarkers(pbmc.merged, ident.1 = 12, min.pct = 0.25)
cluster12 <- filter(cluster12.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster13.markers <- FindMarkers(pbmc.merged, ident.1 = 13, min.pct = 0.25)
cluster13 <- filter(cluster13.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)

azimuth <- DietSeurat(pbmc.merged)
saveRDS(azimuth, file = "../azimuth.rds")	# run azimuth and download tsv predictions
azimuth_pred <- read.table('azimuth_pred.tsv', header=TRUE, sep='\t')
identical(azimuth_pred$cell, rownames(pbmc.merged@meta.data))
pbmc.merged@meta.data$azimuth = azimuth_pred$predicted.celltype.l2
DimPlot(pbmc.merged, group.by='azimuth', label=TRUE, label.size=2, repel=TRUE, cols=col_vector)
matrix_azimuth <- data.frame(row.names=azimuth_pred$cell, pred = azimuth_pred$predicted.celltype.l2.score)
matrix_azimuth <- t(matrix_azimuth)
pbmc.merged[["azimuth_pred"]] <- CreateAssayObject(counts = matrix_azimuth)
FeaturePlot(pbmc.merged, features = "azimuthpred_pred", cols = rev(brewer.pal(n = 11, name = "RdBu")))
pbmc.merged[["azimuth_pred"]] <- NULL

new.cluster.ids <- c("gd T cell", "CD4_basal", "CD4_expanded", "NK cell", "TEM", "CD8_basal", "B cell", "gd_cycling", "gd T cell", "Plasma cell", "TEM", "CD4_expanded", "Apoptotic cell", "Monocyte")
names(new.cluster.ids) <- levels(pbmc.merged)
pbmc.merged <- RenameIdents(pbmc.merged, new.cluster.ids)
DimPlot(pbmc.merged, split.by='condition', cols=col_vector)
annotation <- Idents(pbmc.merged)
annotation <- unname(annotation)
pbmc.merged@meta.data$annotation = annotation
dittoBarPlot(pbmc.merged, 'annotation', group.by='condition')

pbmc.markers <- FindAllMarkers(pbmc.merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc.merged, features = top10$gene) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))

sceasy::convertFormat(pbmc.merged, from="seurat", to="anndata", outFile='scanpy_input.h5ad')
# python
import scanpy as sc
adata =sc.read_h5ad("scanpy_input.h5ad")
marker_genes=['CD3E','CD4', 'CD8A', 'TRDV2', 'TNFRSF4', 'CCR7', 'TYROBP', 'GZMA', 'MS4A1', 'MKI67', 'IGHA1', 'TRAV1-2', 'percent.mt', 'SERPING1']
sc.pl.dotplot(adata, marker_genes, groupby='annotation', standard_scale='var', dendrogram=True, cmap='Purples', save='.pdf')

Idents(pbmc.merged) <- 'condition'
pseudo <- AverageExpression(object = pbmc.merged, group.by='ident', add.ident='SampleName', features=pbmc.markers.cond$gene, assays = 'RNA', slot='scale.data')
pseudo <- pseudo$RNA
breaksList=seq(-1,0.5, by=0.05)
breaksList2=seq(0.6,0.7, by=0.05)
breaksList3=seq(0.8,0.9, by=0.05)
breaksList4=seq(0.9,1, by=0.05)
bk=c(breaksList, breaksList2, breaksList3, breaksList4)
paletteLength <- 40
myColor <- colorRampPalette(c("darkblue", "white", "red", "darkred"))(paletteLength)
pheatmap(pseudo, show_rownames = FALSE, border_color = 'black', cluster_cols = FALSE, cluster_rows=FALSE, color=myColor, cellwidth = 15, fontsize_row = 5)

gd <- subset(pbmc.merged, idents = c("gd T cell"))
gd <- DietSeurat(gd)
gd <- ScaleData(gd, vars.to.regress = c("S.Score", "G2M.Score", "HTO_maxID"), features = rownames(gd), block.size = 100)
gd <- RunPCA(gd, npcs = 20, verbose = FALSE)
gd <- FindNeighbors(gd, dims = 1:20)
gd <- FindClusters(gd, resolution = 0.5)
gd <- RunUMAP(gd, dims = 1:20)

poscells <- WhichCells(gd, expression = TRDV2 > 0 & TRGV9 > 0)
gd$gd_logical <- ifelse(colnames(gd) %in% poscells, "Pos", "Neg")
DimPlot(gd, group.by='gd_logical', split.by='condition')
gd <- subset(gd, idents=c('Pos'))
gd.markers <- FindAllMarkers(gd, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gd.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10_gd
DoHeatmap(gd, features = top10_gd$gene, group.colors = col_vector, size = 2, group.bar.height = 0.01) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))

cluster0.markers <- FindMarkers(gd, ident.1 = 0, min.pct = 0.25)
cluster0 <- filter(cluster0.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster1.markers <- FindMarkers(gd, ident.1 = 1, min.pct = 0.25)
cluster1 <- filter(cluster1.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster2.markers <- FindMarkers(gd, ident.1 = 2, min.pct = 0.25)
cluster2 <- filter(cluster2.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster3.markers <- FindMarkers(gd, ident.1 = 3, min.pct = 0.25)
cluster3 <- filter(cluster3.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster4.markers <- FindMarkers(gd, ident.1 = 4, min.pct = 0.25)
cluster4 <- filter(cluster4.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)
cluster5.markers <- FindMarkers(gd, ident.1 = 5, min.pct = 0.25)
cluster5 <- filter(cluster5.markers, avg_log2FC > 0.58 & p_val_adj < 0.05)

new.cluster.ids <- c("gd_deffective", "gd_TRDV2", "gd_TRDV1-like", "gd_activated", "gd_proliferative", "gd_eff")
names(new.cluster.ids) <- levels(gd)
gd <- RenameIdents(gd, new.cluster.ids)
DimPlot(gd, label = TRUE, label.box = TRUE, cols = "Set2", label.size = 2)

cell.list <- WhichCells(gd, idents = c("gd_deffective", "gd_TRDV2", "gd_TRDV1-like", "gd_activated", "gd_proliferative", "gd_eff"), downsample = 100)
gd_subset <- gd[, cell.list]
data <- as(as.matrix(gd_subset@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = gd_subset@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
HSMM <- newCellDataSet(data, phenoData = pd, featureData = fd, lowerDetectionLimit = 1, expressionFamily = uninormal())
var_genes <- gd_subset[["RNA"]]@var.features
ordering_genes <- var_genes
HSMM <- setOrderingFilter(HSMM, ordering_genes)
HSMM <- reduceDimension(HSMM,norm_method="none", reduction_method="DDRTree", max_components=2, scaling=TRUE, verbose=TRUE, pseudo_expr=0)
HSMM <- orderCells(HSMM,reverse = F)
plot_cell_trajectory(HSMM, color_by = "Pseudotime", theta = -15, show_branch_points = FALSE, show_tree = TRUE, cell_size = 1) + theme(legend.position = "right")
plot_cell_trajectory(HSMM, color_by = "gd_final", theta = -15, show_branch_points = FALSE, show_tree = TRUE, cell_size = 1) + theme(legend.position = "right") + scale_colour_brewer(palette = "Set2")

VlnPlot(gd, features=c('NCR3', 'KLRB1', 'KLRC1', 'KLRC4', 'KLRD1', 'KLRG1', 'KLRK1', 'FCGR3A'), stack = TRUE, flip = TRUE)

clusters_gd <- FindMarkers(gd, ident.1 = "gd_effector", ident.2 = "gd_TRDV2", min.pct = 0.25)
volcano_gd <- data.frame(row.names=row.names(clusters_gd), log2FoldChange=clusters_gd$avg_log2FC, pvalue=clusters_gd$p_val_adj)
EnhancedVolcano(volcano_gd, lab = rownames(volcano_gd), x = 'log2FoldChange', y = 'pvalue', title = 'gd_TRDV2_vs_gd_effector', pCutoff = 10e-4, FCcutoff = 0.58, pointSize = 3.0, labSize = 2.0, col=c('black', 'black', 'black', 'red3'), colAlpha = 1)

dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
regulon <- dorothea_regulon_human %>% dplyr::filter(confidence %in% c("A","B","C"))
pbmc.merged <- run_viper(pbmc.merged, regulon, options = list(method = "scale", minsize = 4, verbose = FALSE))
DefaultAssay(object = pbmc.merged) <- "dorothea"
pbmc.merged <- ScaleData(pbmc.merged, assay = 'dorothea')
sceasy::convertFormat(pbmc.merged, from="seurat", to="anndata", outFile='scanpy_input.h5ad')
# python
import scanpy as sc
adata =sc.read_h5ad("scanpy_input.h5ad")
marker_genes=['IRF4','BATF', 'TCF7', 'STAT5B', 'TBX21', 'EOMES', 'STAT3', 'SOX13']
sc.pl.matrixplot(adata, marker_genes, 'annotation', dendrogram=True, cmap='seismic', standard_scale='var', colorbar_title='column scaled\nexpression', save='.pdf')

DefaultAssay(object = pbmc.merged) <- "RNA"
sceasy::convertFormat(pbmc.merged, from="seurat", to="anndata", outFile='scanpy_input.h5ad')
#python
import scanpy as sc
adata =sc.read_h5ad("scanpy_input.h5ad")
sc.pl.umap(adata, color=['FASLG', 'GZMA', 'TNFSF10', 'PRF1'], size=20, color_map= 'seismic', save='.pdf')

gene.sets <- getGeneSets(library = "C5", gene.sets = 'GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY')
ES <- enrichIt(obj = pbmc.merged, gene.sets = gene.sets)
pbmc.merged <- AddMetaData(pbmc.merged, ES)
ES2 <- data.frame(pbmc.merged[[]], Idents(pbmc.merged))
colnames(ES2)[ncol(ES2)] <- "cluster"
ridgeEnrichment(ES2, gene.set = "GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY", group = "cluster", add.rug = TRUE)

nk <- subset(pbmc.merged, idents = c("NK cell"))
nk <- DietSeurat(nk)
nk <- ScaleData(nk, vars.to.regress = c("S.Score", "G2M.Score", "HTO_maxID"), features = rownames(nk), block.size = 100)
nk <- RunPCA(nk, npcs = 20, verbose = FALSE)
nk <- FindNeighbors(nk, dims = 1:20)
nk <- FindClusters(nk, resolution = 0.5)
nk <- RunUMAP(nk, dims = 1:20)
DimPlot(nk, split.by='condition', cols='Dark2')
VlnPlot(nk, features=c('NCAM1', 'FCGR3A', 'XCL1', 'XCL2', 'FCER1G', 'KLRB1', 'KLRC1', 'KLRC2', 'KLRD1', 'KLRG1', 'KLRK1', 'KIR2DL4', 'KIR2DL3', 'KIR3DL1', 'NCR1', 'NCR3'), stack = TRUE, flip = TRUE)
nk.markers <- FindAllMarkers(nk, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
nk.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10_nk
DoHeatmap(nk, features = top10_nk$gene, group.colors = Dark1, size = 2, group.bar.height = 0.01) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))
