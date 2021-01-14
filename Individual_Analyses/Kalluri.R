#Kalluri et al 2019 dataset - Analysis for Zeca et al
#Whole mice aortas
#####################################################################

#####################################################################
#Load libraries
######################################################################

reticulate::use_python('/usr/bin/python3')
library(loomR)
library(reticulate)
library(Matrix)
library(BUSpaRse)
library(Seurat)
library(SeuratWrappers)
library(BSgenome.Mmusculus.UCSC.mm10)
library(AnnotationHub)
library(zeallot) # For %<-% that unpacks lists in the Python manner
library(DropletUtils)
library(tidyverse)
library(uwot) # For umap
library(GGally) # For ggpairs
library(velocyto.R)
library(SingleR)
library(scales)
library(plotly)

setwd("~/Documents/Bioinfo/Ebru/Kalluri")
#####################################################################
#Read in the data
######################################################################

meta <- read.table('META_DATA_Chow_12PCs_outfile.txt', sep = '\t', header = T, row.names = 1)
raw <- read.table('Seurat_Chow_12PCs_outfile.mtx') #.mtx file extension but data in table format - how did this even get past publication?
raw <- as.matrix(raw)

genes <- raw[,1]
genes <- genes[-1]

cells <- raw[1,]
cells <- cells[-1]

raw <- raw[-1,-1]
colnames(raw) <- cells
rownames(raw) <- genes

#Create Seurat Object

dat <- CreateSeuratObject(counts = raw, project = 'Aorta', meta.data = meta, names.field = 1, names.delim = "-")


#####################################################################
#QC
######################################################################

#Mitochondrial Genes
mito.genes <- grep(pattern = "^mt-", x = rownames(dat@assays$RNA@counts), value = TRUE)
percent.mito <- Matrix::colSums(dat@assays$RNA@counts[mito.genes, ])/Matrix::colSums(dat@assays$RNA@counts)
dat <- AddMetaData(object = dat, metadata = percent.mito, col.name = "percent.mito")

#General Metrics
counts_per_cell <- Matrix::colSums(dat@assays$RNA@counts)
counts_per_gene <- Matrix::rowSums(dat)
genes_per_cell <- Matrix::colSums(dat@assays$RNA@counts >0)

#Metrics plotting
hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')
plot(counts_per_cell, genes_per_cell, log='xy', col='wheat', title('counts vs genes per cell'))
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell')
plot(sort(counts_per_cell), xlab='cell', log='y', main='counts per cell')

VlnPlot(object = dat, features = "nCount_RNA")
VlnPlot(object = dat, features = "nFeature_RNA")
VlnPlot(object = dat, features = "percent.mito")

dat <- subset(dat, subset = nCount_RNA > 2000 & nCount_RNA < 3500 & 
                            percent.mito < 0.04)

#####################################################################
#Processing 
######################################################################

dat <- SCTransform(dat,  return.only.var.genes = F, variable.features.n = 3500)

dat <- RunPCA(dat, npcs = 100)
ElbowPlot(dat, ndims = 100)
PCAPlot(dat)

dat <- FindNeighbors(dat, dims = 1:18)
dat <- FindClusters(dat, algorithm = 2)
PCAPlot(dat)

dat <- RunUMAP(dat, dims = 1:18, min.dist = 0.5)
UMAPPlot(dat, pt.size = 1, group.by = 'Sub.Cluster')

#####################################################################
#dbMAP
######################################################################
genes <- dat@assays$SCT@var.features
data <- dat@assays$SCT@data[genes,]
barcodes <- colnames(data)

write.table(genes, file = 'inputs/genes.txt', sep = ' ', col.names = F, row.names = F)
write.table(barcodes, file = 'inputs/barcodes.txt', sep = ' ', col.names = F, row.names = F)
writeMM(obj = t(data), file = 'inputs/matrix.mtx')

py_run_file('~/Documents/Bioinfo/Zeca/Diff_Maps.py') #Run the diffusion script

######################################################################
#Export python outputs to R 
######################################################################
dmap <- py$dmap
evals <- py$evals
plot(evals) #Chose 20 dims

mvals <- py$mvals
mms <- py$ms_data
mms <- as.matrix(mms)
plot(mvals) #Chose 20 dims

dat@reductions$dmap <- dat@reductions$pca
rownames(dmap) <- colnames(dat)
colnames(dmap) <- colnames(dat@reductions$pca[,1:100])
dat@reductions$dmap@cell.embeddings <- dmap

dat@reductions$mms <- dat@reductions$pca
rownames(mms) <- colnames(dat)
colnames(mms) <- colnames(dat@reductions$pca[,1:99])
dat@reductions$mms@cell.embeddings <- mms

######################################################################
#Run diffusion embeddings
######################################################################

dat <- RunUMAP(object = dat, dims = 1:20, min.dist = 1, reduction = 'dmap', reduction.key = 'dbMAP_', reduction.name = 'dbmap', learning.rate = 2)
DimPlot(dat, reduction = 'dbmap', pt.size = 1, group.by = 'Sub.Cluster')

dat <- RunUMAP(object = dat, dims = 1:20, min.dist = 1, reduction = 'mms', reduction.key = 'dbMAPa_', reduction.name = 'dbmapa', learning.rate = 2)
DimPlot(dat, reduction = 'dbmapa', pt.size = 1, group.by = 'Sub.Cluster')
DimPlot(dat, reduction = 'dbmapa', pt.size = 1, group.by = 'seurat_clusters')

dat <- FindNeighbors(dat, reduction = 'pca', dims = 1:20)
dat <- FindClusters(dat, algorithm = 2, resolution = 0.5)
DimPlot(dat, reduction = 'dbmapa', pt.size = 1, group.by = 'seurat_clusters')

saveRDS(dat, 'Kalluri.Rds')

######################################################################
#Heatmaps
######################################################################
library(dplyr)
dat <- RunALRA(dat)
ALRAChooseKPlot(dat)
dat <- ScaleData(dat)
dat.markers <- FindAllMarkers(dat, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25, assay = 'alra', )
dat.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top10 <- dat.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DoHeatmap(dat, features = top10$gene, assay = 'alra') + NoLegend()

######################################################################
#Export to cerebro
######################################################################

library(cerebroApp)
dat <- addPercentMtRibo(dat, organism = 'mm', gene_nomenclature = 'name')
dat <- getMostExpressedGenes(dat, column_sample = 'Sub.Cluster', column_cluster = 'seurat_clusters')
dat <- getMarkerGenes(dat, organism = 'mm', column_sample = 'Sub.Cluster', column_cluster = 'seurat_clusters')

dat <- getEnrichedPathways(dat, column_sample = 'Sub.Cluster', column_cluster = 'seurat_clusters',
                           databases = c('Transcription_Factor_PPIs', 'Chromossome_Location', 'Mouse_Gene_Atlas', 
                                         'MGI_Mammalian_Phenotype_Level_3', 'MGI_Mammalian_Phenotyoe_Level_4_2019', 'MGI_Mammalian_Phenotype_2017',
                                         'OMIM_Disease', 'OMIM_Expanded', 'MsigDB_Computational', 'UK_Biobank_GWAS_v1',
                                         'KEGG_2019_Mouse', 'PheWeb_2019', 'WikiPathways_2019_Mouse', 'GWAS_Catalog_2019',
                                         'GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018',
                                         'Panther_2016', 'BioCarta_2016', 'Reactome_2016', 'Kegg_2016'))



cerebro <- exportFromSeurat(dat, file = 'Kalluri.crb', experiment_name = 'Murine Aorta', organism = 'mm',
                            column_sample = 'Sub.Cluster', column_cluster = 'seurat_clusters',
                            column_nUMI = 'nCount_RNA', column_nGene = 'nFeature_RNA')

library(cerebroApp)
launchCerebro()


# Convert to h5Seurat, h5ad
library(SeuratDisk)


counts <- kal@assays$RNA@counts
meta <- kal@meta.data
dat <- CreateSeuratObject(counts, meta.data = meta)

as.h5Seurat(dat, 'Kalluri.h5Seurat', overwrite = T)
Convert('Kalluri.h5Seurat', dest = 'h5ad', overwrite = T, assay = 'RNA')


######################################################################
#Export to Arboreto - SCENIC pipeline
######################################################################

exprMat <- t(as.matrix(dat@assays$SCT@data))
write.table(exprMat, file = 'ExpMat.tsv', sep = '\t', col.names = T, row.names = F)

















