######################################################################
#Analysis of Pan et al mouse data (will use as mouse reference)
######################################################################

#Load packages 
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(reticulate)
library(cerebroApp)
library(plotly)
setwd("~/Documents/Bioinfo/Vascular/Pan/")

#Load data
# SMC-derived (Myh11-CreERT2+ = ZsGreen+ = cells that came from SMCs)
ZsPos_ApoeKO_8w <- as.matrix(read.delim('data/Mouse/ZsGreen1+_Apoe_KO_8w_WD_matrix.txt.gz', sep = '\t', header = T, row.names = 1))
ZsPos_ApoeKO_8w <- CreateSeuratObject(ZsPos_ApoeKO_8w, project = 'Myh11-CreERT2-ZsGreen+/Apoe KO/8w WD')
#
ZsPos_ApoeKO_16w <- as.matrix(read.delim('data/Mouse/ZsGreen1+_Apoe_KO_16w_WD_matrix.txt.gz', sep = '\t', header = T, row.names = 1))
ZsPos_ApoeKO_16w <- CreateSeuratObject(ZsPos_ApoeKO_16w, project = 'Myh11-CreERT2-ZsGreen+/Apoe KO/16w WD')
#
ZsPos_ApoeKO_22w <- as.matrix(read.delim('data/Mouse/ZsGreen1+_Apoe_KO_22w_WD_matrix.txt.gz', sep = '\t', header = T, row.names = 1))
ZsPos_ApoeKO_22w <- CreateSeuratObject(ZsPos_ApoeKO_22w, project = 'Myh11-CreERT2-ZsGreen+/Apoe KO/22w WD')
#
ZsPos_LdlrKO_0w <- as.matrix(read.delim('data/Mouse/ZsGreen1+_Ldlr_KO_0w_WD_matrix.txt.gz', sep = '\t', header = T, row.names = 1))
ZsPos_LdlrKO_0w <- CreateSeuratObject(ZsPos_LdlrKO_0w, project = 'Myh11-CreERT2-ZsGreen+/Ldlr KO/0w WD')
#
ZsPos_LdlrKO_8w <- as.matrix(read.delim('data/Mouse/ZsGreen1+_Ldlr_KO_8w_WD_matrix.txt.gz', sep = '\t', header = T, row.names = 1))
ZsPos_LdlrKO_8w <- CreateSeuratObject(ZsPos_LdlrKO_8w, project = 'Myh11-CreERT2-ZsGreen+/Ldlr KO/8w WD')
#
ZsPos_LdlrKO_16w <- as.matrix(read.delim('data/Mouse/ZsGreen1+_Ldlr_KO_16w_WD_matrix.txt.gz', sep = '\t', header = T, row.names = 1))
ZsPos_LdlrKO_16w <- CreateSeuratObject(ZsPos_LdlrKO_16w, project = 'Myh11-CreERT2-ZsGreen+/Ldlr KO/16w WD')
#
ZsPos_LdlrKO_26w <- as.matrix(read.delim('data/Mouse/ZsGreen1+_Ldlr_KO_26w_WD_matrix.txt.gz', sep = '\t', header = T, row.names = 1))
ZsPos_LdlrKO_26w <- CreateSeuratObject(ZsPos_LdlrKO_26w, project = 'Myh11-CreERT2-ZsGreen+/Ldlr KO/26w WD')

# NOT-SMC-derived (Myh11-CreERT2- = ZsGreen- = cells that DID NOT came from SMCs)
ZsNeg_ApoeKO_8w <- as.matrix(read.delim('data/Mouse/ZsGreen1-_Apoe_KO_8w_WD_matrix.txt.gz', sep = '\t', header = T, row.names = 1))
ZsNeg_ApoeKO_8w <- CreateSeuratObject(ZsNeg_ApoeKO_8w, project = 'Myh11-CreERT2-ZsGreen-/Apoe KO/8w WD')
#
ZsNeg_ApoeKO_16w <- as.matrix(read.delim('data/Mouse/ZsGreen1-_Apoe_KO_16w_WD_matrix.txt.gz', sep = '\t', header = T, row.names = 1))
ZsNeg_ApoeKO_16w <- CreateSeuratObject(ZsNeg_ApoeKO_16w, project = 'Myh11-CreERT2-ZsGreen-/Apoe KO/16w WD')
#
ZsNeg_ApoeKO_22w <- as.matrix(read.delim('data/Mouse/ZsGreen1-_Apoe_KO_22w_WD_matrix.txt.gz', sep = '\t', header = T, row.names = 1))
ZsNeg_ApoeKO_22w <- CreateSeuratObject(ZsNeg_ApoeKO_22w, project = 'Myh11-CreERT2-ZsGreen-/Apoe KO/22w WD')
#
ZsNeg_LdlrKO_0w <- as.matrix(read.delim('data/Mouse/ZsGreen1-_Ldlr_KO_0w_WD_matrix.txt.gz', sep = '\t', header = T, row.names = 1))
ZsNeg_LdlrKO_0w <- CreateSeuratObject(ZsNeg_LdlrKO_0w, project = 'Myh11-CreERT2-ZsGreen-/Ldlr KO/0w WD')
#
ZsNeg_LdlrKO_8w <- as.matrix(read.delim('data/Mouse/ZsGreen1-_Ldlr_KO_8w_WD_matrix.txt.gz', sep = '\t', header = T, row.names = 1))
ZsNeg_LdlrKO_8w <- CreateSeuratObject(ZsNeg_LdlrKO_8w, project = 'Myh11-CreERT2-ZsGreen-/Ldlr KO/8w WD')
#
ZsNeg_LdlrKO_16w <- as.matrix(read.delim('data/Mouse/ZsGreen1-_Ldlr_KO_16w_WD_matrix.txt.gz', sep = '\t', header = T, row.names = 1))
ZsNeg_LdlrKO_16w <- CreateSeuratObject(ZsNeg_LdlrKO_16w, project = 'Myh11-CreERT2-ZsGreen-/Ldlr KO/16w WD')
#
ZsNeg_LdlrKO_26w <- as.matrix(read.delim('data/Mouse/ZsGreen1-_Ldlr_KO_26w_WD_matrix.txt.gz', sep = '\t', header = T, row.names = 1))
ZsNeg_LdlrKO_26w <- CreateSeuratObject(ZsNeg_LdlrKO_26w, project = 'Myh11-CreERT2-ZsGreen-/Ldlr KO/26w WD')
#

# Myh11-CreERT2-ZsGreen+ / Apoe KO
ZsPos_Apoe <- list(ZsPos_ApoeKO_8w, ZsPos_ApoeKO_16w, ZsPos_ApoeKO_22w)
#
# Myh11-CreERT2-ZsGreen- / Apoe KO
ZsNeg_Apoe <- list(ZsNeg_ApoeKO_8w, ZsNeg_ApoeKO_16w, ZsNeg_ApoeKO_22w)
#
# Myh11-CreERT2-ZsGreen+ / Ldlr KO
ZsPos_Ldlr <- list(ZsPos_LdlrKO_0w, ZsPos_LdlrKO_8w,  ZsPos_LdlrKO_16w, ZsPos_LdlrKO_26w)
#
# Myh11-CreERT2-ZsGreen- / Ldlr KO
ZsNeg_Ldlr <- list(ZsNeg_LdlrKO_0w, ZsNeg_LdlrKO_8w,  ZsNeg_LdlrKO_16w, ZsNeg_LdlrKO_26w)
#

######################################################################
#QC and filtering 
######################################################################
dats <- list(ZsPos_Apoe, ZsNeg_Apoe, ZsPos_Ldlr, ZsNeg_Ldlr)

group <- list()
seurat <- list()

# This is the filtered public-released data. Let's check how filtered it is.
for(i in 1:length(dats)){
  group <- dats[[i]]
    for(i in 1:length(group)){
      seurat <- group[[i]]
      
      seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-", assay = 'RNA')
      counts_per_cell <- Matrix::colSums(seurat)
      counts_per_gene <- Matrix::rowSums(seurat)
      genes_per_cell <- Matrix::colSums(seurat@assays$RNA@counts >0)
      
      hist(log10(counts_per_cell+1),main='Counts per cell',col='chartreuse4')
      hist(log10(genes_per_cell+1), main='Genes per cell', col='chartreuse4')
      
      plot(counts_per_cell, genes_per_cell, type = 'p', log='xy', pch=19, col="firebrick3", xlab="Counts per cell (log)", ylab="Genes per cell (log)", title(paste('Counts vs Genes per Cell', as.character(Project(seurat), sep = ' - '))))
      plot(sort(genes_per_cell), pch=19, col="firebrick3", xlab="Cell rank", ylab="Genes per cell (log)", 
           log='y', main= paste('Gene detection rank', as.character(Project(seurat), sep = ' - ')))
      plot(sort(counts_per_cell), pch=19, col="firebrick3", xlab="Cell rank", ylab="Genes per cell (log)",
           log='y', main= paste('Count detection rank', as.character(Project(seurat), sep = ' - ')))
      
      group[[i]] <- seurat
    }  
}

# Actually looks quite good! Great job, Dr. Pan and colleagues !!!
# We won't filter this any further
# Merge data sets

ZsPos_Apoe <-  dats[[1]]
ZsPos_Apoe <- merge(ZsPos_Apoe[[1]], ZsPos_Apoe[-1]) # 9508 cells vs 14542 genes
ZsPos_Apoe[["percent.mt"]] <- PercentageFeatureSet(ZsPos_Apoe, pattern = "^mt-", assay = 'RNA')

ZsNeg_Apoe <-  dats[[2]]
ZsNeg_Apoe <-  merge(ZsNeg_Apoe[[1]], ZsNeg_Apoe[-1]) # 5572 cells vs 14182 genes
ZsNeg_Apoe[["percent.mt"]] <- PercentageFeatureSet(ZsNeg_Apoe, pattern = "^mt-", assay = 'RNA')

ZsPos_Ldlr <-  dats[[3]]
ZsPos_Ldlr <-  merge(ZsPos_Ldlr[[1]], ZsPos_Ldlr[-1]) # 15920 cells vs 14855 genes 
ZsPos_Ldlr[["percent.mt"]] <- PercentageFeatureSet(ZsPos_Ldlr, pattern = "^mt-", assay = 'RNA')

ZsNeg_Ldlr <-  dats[[4]]
ZsNeg_Ldlr <-  merge(ZsNeg_Ldlr[[1]], ZsNeg_Ldlr[-1]) # 11372 cells vs 14623 genes
ZsNeg_Ldlr[["percent.mt"]] <- PercentageFeatureSet(ZsNeg_Ldlr, pattern = "^mt-", assay = 'RNA')

VlnPlot(object = ZsPos_Apoe, features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), group.by = 'orig.ident')
VlnPlot(object = ZsNeg_Apoe, features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), group.by = 'orig.ident')
VlnPlot(object = ZsPos_Ldlr, features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), group.by = 'orig.ident')
VlnPlot(object = ZsNeg_Ldlr, features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), group.by = 'orig.ident')

dat <- merge(ZsPos_Apoe, list(ZsNeg_Apoe, ZsPos_Ldlr, ZsNeg_Ldlr))
# Total of 42,372 cells and 15,852 detected genes

VlnPlot(object = dat, features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), group.by = 'orig.ident')

######################################################################
#Add some meta-data
######################################################################

# Diet time-course WD meta-data
diet_l_0w = rep('No WD', times = (ncol(ZsPos_LdlrKO_0w) + ncol(ZsNeg_LdlrKO_0w)))
names(diet_l_0w) <- c(colnames(ZsPos_LdlrKO_0w), colnames(ZsNeg_LdlrKO_0w))
#
diet_l_8w = rep('8w WD', times = (ncol(ZsPos_ApoeKO_8w) + ncol(ZsNeg_ApoeKO_8w)
                                  + ncol(ZsPos_LdlrKO_8w) + ncol(ZsNeg_LdlrKO_8w)))
names(diet_l_8w) <- c(colnames(ZsPos_ApoeKO_8w), colnames(ZsNeg_ApoeKO_8w),
                      colnames(ZsPos_LdlrKO_8w), colnames(ZsNeg_LdlrKO_8w))
#
diet_l_16w = rep('16w WD', times = (ncol(ZsPos_ApoeKO_16w) + ncol(ZsNeg_ApoeKO_16w) +
                                ncol(ZsPos_LdlrKO_16w) + ncol(ZsNeg_LdlrKO_16w)))
names(diet_l_16w) <- c(colnames(ZsPos_ApoeKO_16w), colnames(ZsNeg_ApoeKO_16w), 
                 colnames(ZsPos_LdlrKO_16w), colnames(ZsNeg_LdlrKO_16w))
#
diet_l_22w = rep('22w WD', times = (ncol(ZsPos_ApoeKO_22w) + ncol(ZsNeg_ApoeKO_22w)))
names(diet_l_22w) <- c(colnames(ZsPos_ApoeKO_22w), colnames(ZsNeg_ApoeKO_22w))
#
diet_l_26w = rep('26w WD', times = (ncol(ZsPos_LdlrKO_26w) + ncol(ZsNeg_LdlrKO_26w)))
names(diet_l_26w) <- c(colnames(ZsPos_LdlrKO_26w), colnames(ZsNeg_LdlrKO_26w))
#
diet_label <- c(diet_l_0w, diet_l_8w, diet_l_16w, diet_l_22w, diet_l_26w)
dat <- AddMetaData(dat, diet_label, col.name = 'Diet')

# Genotype meta-data
gen_l_apoe <- rep('Apoe KO', times = ncol(ZsPos_ApoeKO_8w) + ncol(ZsNeg_ApoeKO_8w) 
                  + ncol(ZsPos_ApoeKO_16w) + ncol(ZsNeg_ApoeKO_16w)
                  + ncol(ZsPos_ApoeKO_22w) + ncol(ZsNeg_ApoeKO_22w)
                  )
names(gen_l_apoe) <- c(colnames(ZsPos_ApoeKO_8w), colnames(ZsNeg_ApoeKO_8w),
                       colnames(ZsPos_ApoeKO_16w), colnames(ZsNeg_ApoeKO_16w),
                       colnames(ZsPos_ApoeKO_22w), colnames(ZsNeg_ApoeKO_22w)
                       )
#
gen_l_ldlr <- rep('Ldlr KO', times = ncol(ZsPos_LdlrKO_0w) + ncol(ZsNeg_LdlrKO_0w) 
                  + ncol(ZsPos_LdlrKO_8w) + ncol(ZsNeg_LdlrKO_8w)
                  + ncol(ZsPos_LdlrKO_16w) + ncol(ZsNeg_LdlrKO_16w)
                  + ncol(ZsPos_LdlrKO_26w) + ncol(ZsNeg_LdlrKO_26w)
                  )
names(gen_l_ldlr) <- c(colnames(ZsPos_LdlrKO_0w), colnames(ZsNeg_LdlrKO_0w),
                       colnames(ZsPos_LdlrKO_8w), colnames(ZsNeg_LdlrKO_8w),
                       colnames(ZsPos_LdlrKO_16w), colnames(ZsNeg_LdlrKO_16w),
                       colnames(ZsPos_LdlrKO_26w), colnames(ZsNeg_LdlrKO_26w)
                       )
#
gen_label <- c(gen_l_apoe, gen_l_ldlr)
dat <- AddMetaData(dat, gen_label, col.name = 'Genotype')

# Sorting
# Myh11-CreERT2-ZsGreen+
sort_l_ZsPos <- rep('Myh11-CreERT2-ZsGreen Positive', times = ncol(ZsPos_LdlrKO_0w) 
                    + ncol(ZsPos_LdlrKO_8w) + ncol(ZsPos_LdlrKO_16w) 
                    + ncol(ZsPos_LdlrKO_26w) + ncol(ZsPos_ApoeKO_8w) 
                    + ncol(ZsPos_ApoeKO_16w) + ncol(ZsPos_ApoeKO_22w)
                    )
names(sort_l_ZsPos) <- c(colnames(ZsPos_LdlrKO_0w), colnames(ZsPos_LdlrKO_8w),
                       colnames(ZsPos_LdlrKO_16w), colnames(ZsPos_LdlrKO_26w),
                       colnames(ZsPos_ApoeKO_8w), colnames(ZsPos_ApoeKO_16w),
                       colnames(ZsPos_ApoeKO_22w)
                       )
#
sort_l_ZsNeg <- rep('Myh11-CreERT2-ZsGreen Negative', times = ncol(ZsNeg_LdlrKO_0w) 
                    + ncol(ZsNeg_LdlrKO_8w) + ncol(ZsNeg_LdlrKO_16w) 
                    + ncol(ZsNeg_LdlrKO_26w) + ncol(ZsNeg_ApoeKO_8w) 
                    + ncol(ZsNeg_ApoeKO_16w) + ncol(ZsNeg_ApoeKO_22w)
)
names(sort_l_ZsPos) <- c(colnames(ZsNeg_LdlrKO_0w), colnames(ZsNeg_LdlrKO_8w),
                         colnames(ZsNeg_LdlrKO_16w), colnames(ZsNeg_ApoeKO_8w),
                         colnames(ZsPos_ApoeKO_8w), colnames(ZsPos_ApoeKO_16w),
                         colnames(ZsPos_ApoeKO_22w)
)
#
sorting_label <- c(sort_l_ZsPos, sort_l_ZsNeg)
dat <- AddMetaData(dat, sorting_label, col.name = 'Sorting')

VlnPlot(object = dat, features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), group.by = 'Sorting')
VlnPlot(object = dat, features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), group.by = 'Genotype')
VlnPlot(object = dat, features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), group.by = 'Diet')


######################################################################
#Default workflow, integrate with CCA anchoring and cluster with graph Louvain
######################################################################
genes <- rownames(dat)
dat$Sample <- dat$orig.ident
sets <- SplitObject(dat, split.by = 'Sample')

for(i in 1:length(sets)){
  sets[[i]] <- NormalizeData(sets[[i]])
  sets[[i]] <- FindVariableFeatures(sets[[i]], selection.method = 'disp', nfeatures = 3000, num.bin = 100)
  sets[[i]] <- ScaleData(sets[[i]], features = rownames(sets[[i]]))
}

feat <- SelectIntegrationFeatures(object.list = sets, nfeatures = 3000)
anchors <- FindIntegrationAnchors(object.list = sets, anchor.features = feat, 
                                  normalization.method = 'LogNormalize', dims = 1:50, 
                                  max.features = 800, nn.method = 'annoy', reduction = 'cca')
dat <- IntegrateData(anchorset = anchors, features.to.integrate = genes, dims = 1:100, k.weight = 30)
dat <- ScaleData(dat)

dat <- RunPCA(dat, npcs = 100)
ElbowPlot(dat, ndims = 100)

dat <- FindNeighbors(object = dat, features = VariableFeatures(dat))
dat <- FindClusters(object = dat, features = VariableFeatures(dat), algorithm = 2,
                    resolution = 1.5)
dat$gene_clusters <- dat$seurat_clusters

dat <- FindNeighbors(object = dat, dims = 1:50)
dat <- FindClusters(object = dat, algorithm = 2, resolution = 1.5)
dat$pca_clusters <- dat$seurat_clusters


dat <- RunUMAP(object = dat, dims = 1:50)
UMAPPlot(dat, group.by = 'Sample') 
UMAPPlot(dat, group.by = 'Diet') 
UMAPPlot(dat, group.by = 'Genotype') 
UMAPPlot(dat, group.by = 'Sorting') 

UMAPPlot(dat, group.by = 'pca_clusters', pt.size = 0.5) 
UMAPPlot(dat, group.by = 'gene_clusters', pt.size = 0.5) 

FeaturePlot(dat, features = 'Foxp3', order = T, min.cutoff = 0, max.cutoff = 1, pt.size = 1)
FeaturePlot(dat, features = 'Ikzf2', order = T)
FeaturePlot(dat, features = 'Ccr2', order = T, min.cutoff = 2.5)
FeaturePlot(dat, features = 'Cd209a', order = T, min.cutoff = 0, max.cutoff = 1)
FeaturePlot(dat, features = 'Fabp5', order = T, min.cutoff = 0, max.cutoff = 4)

VlnPlot(dat, features = c('Ikzf2', 'Cd209a', 'Ccr2'), group.by = 'Group') 

###############################################################################
# Run dbMAP
###############################################################################
library(reticulate)
np <- reticulate::import("numpy")
pd <- reticulate::import("pandas")
sp <- reticulate::import("scipy")
dbmap <- reticulate::import('dbmap')

data <- t(dat@assays$integrated@data[VariableFeatures(dat),])
a <- r_to_py(data)
b <- a$tocsr()
diff <- dbmap$diffusion$Diffusor(n_components = as.integer(200), n_neighbors = as.integer(30),
                                 transitions = as.logical(F),
                                 norm = as.logical(F), ann_dist = as.character('cosine'),
                                 n_jobs = as.integer(10), kernel_use = as.character('simple')
)
diff = diff$fit(b)
sc = as.matrix(diff$transform(b))
res = diff$return_dict()

evals <- py_to_r(res$EigenValues)
plot(evals) #Visualize meaningful diffusion components.

rownames(sc) <- colnames(dat)
new_names <- list()
for(i in 1:length(colnames(sc))){
  new_names[i] <- paste('SC_' , as.integer(colnames(sc[i])) + 1, sep = '')
}
colnames(sc) <- as.vector(new_names)
names(evals) <- as.vector(new_names)

dat@reductions$sc <- dat@reductions$pca 
dat@reductions$sc@cell.embeddings <- as.matrix(sc)
dat@reductions$sc@feature.loadings <- matrix(data = c(0,0))
dat@reductions$sc@key <- 'SC_'

##################################################################################
# Post-diffusion processing (structure clustering, adaptive UMAP)
##################################################################################

dat <- FindNeighbors(dat, reduction = 'sc', dims = 1:ncol(dat@reductions$sc@cell.embeddings), annoy.metric = 'cosine', graph.name = 'dbgraph')
dat <- FindClusters(dat, resolution = 1, graph.name = 'dbgraph', algorithm = 2)

dat <- RunUMAP(dat, reduction = 'sc', dims = 1:ncol(dat@reductions$sc@cell.embeddings), n.neighbors = 30, init = 'spectral',
               min.dist = 0.6, spread = 1.5, learning.rate = 1.5, n.epochs = 200, reduction.key = 'dbMAP_', reduction.name = 'dbmap')

dat <- RunUMAP(dat, reduction = 'sc', n.components = 3, dims = 1:ncol(dat@reductions$sc@cell.embeddings),
               min.dist = 0.45, spread = 0.8, n.epochs = 200, reduction.key = 'dbMAP_3D_', reduction.name = 'dbmap3d')

DimPlot(dat, reduction = 'dbmap', group.by = 'seurat_clusters', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'gene_clusters', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'Group', pt.size = 1)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
dat <- CellCycleScoring(dat, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

saveRDS(dat, 'Pan_Mouse.Rds')

# Convert to h5Seurat, h5ad
library(SeuratDisk)
as.h5Seurat(dat, 'Pan_Mouse.h5Seurat', overwrite = T)
Convert('Pan_Mouse.h5Seurat', dest = 'h5ad', overwrite = T, assay = 'RNA')


######################################################################
#Plot and explore for annotation
########################################################### ###########

DimPlot(dat, reduction = 'dbmap', group.by = 'Group', pt.size = 0.5)
DimPlot(dat, reduction = 'dbmap', group.by = 'pca_clusters', pt.size = 0.5)
DimPlot(dat, reduction = 'dbmap', group.by = 'seurat_clusters', pt.size = 0.5)

# General markers
FeaturePlot(dat, reduction = 'dbmap', features = 'Cd7', pt.size = 0.5, min.cutoff = 0, max.cutoff = 1, order = T) 


# Functional macrophage markers
FeaturePlot(dat, reduction = 'dbmap', features = 'Mki67', pt.size = 0.5, min.cutoff = 0, max.cutoff = 1, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'Trem2', pt.size = 0.5, min.cutoff = 0, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'Ccl8', pt.size = 0.5, min.cutoff = 0, order = T) 

plot.data <- FetchData(object = dat, vars = c("dbMAP3D_1", "dbMAP3D_2", "dbMAP3D_3", 
                                              "pca_clusters",  'seurat_clusters', 'Group'
))
plot.data$label <- dat$seurat_clusters
plot_ly(data = plot.data, 
        x = ~dbMAP3D_1, y = ~dbMAP3D_2, z = ~dbMAP3D_3, 
        color = ~seurat_clusters, 
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 3, width=0.5),
        text=~seurat_clusters,
        hoverinfo="text", plot_bgcolor = 'black') 


library(cerebroApp)
dat <- addPercentMtRibo(dat, organism = 'mm', gene_nomenclature = 'name', assay = 'RNA')
dat <- getMarkerGenes(dat, organism = 'mm', groups = c('seurat_clusters', 'Diet', 'Genotype', 'Sorting', 'Sample'), min_pct = 0.3, assay = 'integrated', test = 'wilcox')


dat <- getEnrichedPathways(dat, 
                           databases = c('Transcription_Factor_PPIs', 'Chromossome_Location', 'Mouse_Gene_Atlas', 
                                         'MGI_Mammalian_Phenotype_Level_3', 'MGI_Mammalian_Phenotyoe_Level_4_2019', 'MGI_Mammalian_Phenotype_2017',
                                         'OMIM_Disease', 'OMIM_Expanded', 'MsigDB_Computational', 'UK_Biobank_GWAS_v1',
                                         'KEGG_2019_Mouse', 'PheWeb_2019', 'WikiPathways_2019_Mouse', 'GWAS_Catalog_2019',
                                         'GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018',
                                         'Panther_2016', 'BioCarta_2016', 'Reactome_2016', 'Kegg_2016'))
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('seurat_clusters', 'Diet', 'Genotype', 'Sorting', 'Sample'),
                                        GMT_file =  'h.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('seurat_clusters', 'Diet', 'Genotype', 'Sorting', 'Sample'),
                                        GMT_file =  'c2.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('seurat_clusters', 'Diet', 'Genotype', 'Sorting', 'Sample'),
                                        GMT_file =  'c3.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('seurat_clusters', 'Diet', 'Genotype', 'Sorting', 'Sample'),
                                        GMT_file =  'c5.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('seurat_clusters', 'Diet', 'Genotype', 'Sorting', 'Sample'),
                                        GMT_file =  'c3.tft.v7.1.symbols.gmt')

cerebro <- exportFromSeurat(dat, file = 'Pan_Mouse.crb', experiment_name = 'Pan et al - mouse data', organism = 'mm',
                            cell_cycle = 'Phase',  groups = c('seurat_clusters', 'Diet', 'Genotype', 'Sorting', 'Sample'),
                            nUMI = 'nCount_RNA', nGene = 'nFeature_RNA')
launchCerebro()


####################################################################
# Transfer labels
####################################################################
win <- readRDS('Winkels/Winkels.Rds')
transfer.anchors <- FindTransferAnchors(reference = win, query = dat, reduction = 'cca',
                                        dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = win$cluster, weight.reduction = 'cca',
                            dims = 1:30)
dat <- AddMetaData(dat, metadata = predictions)
dat$pred.win <- dat$predicted.id

coc <- readRDS('Cochain/Cochain.Rds')
transfer.anchors <- FindTransferAnchors(reference = coc, query = dat, 
                                        dims = 1:30, reduction = 'cca')
predictions <- TransferData(anchorset = transfer.anchors, refdata = coc$db_clusters, weight.reduction = 'cca',
                            dims = 1:30)
dat <- AddMetaData(dat, metadata = predictions)
dat$pred.coc <- dat$predicted.id



Idents(dat) <- dat$db_clusters
new.clusters.ids <- c('B-Cells (1)',
                      'T CD8/Lef1',
                      'Macrophages',
                      'T CD8/Gzmk',
                      'Monocyte',
                      'T pro-helper',
                      'B-Cells (2)',
                      'T intermediate',
                      'T CD8/Ccl5',
                      'B-Cells (3)',
                      'Ly-6C+ mono',
                      'T Progenitor',
                      'T Dusp10',
                      'Th2',
                      'Th17',
                      'NK',
                      'NA/IFIT1',
                      'NA/MYL10')
names(new.clusters.ids) <- levels(dat$db_clusters)
dat <- RenameIdents(dat, new.clusters.ids)
dat$db_clusters <- Idents(dat)

DimPlot(dat, reduction = 'dbmap', pt.size = 0.5)

