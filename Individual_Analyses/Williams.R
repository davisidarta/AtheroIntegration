######################################################################
#Analysis of Williams et al data
######################################################################

#Load packages 
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(reticulate)
library(cerebroApp)
setwd("~/Documents/Bioinfo/Vascular/Williams/")

#Load data
c <- Read10X('data/B6', gene.column = 2)
h <- Read10X('data/D21HFD', gene.column = 2)

c <- CreateSeuratObject(c)
h <- CreateSeuratObject(h)

######################################################################
#QC and filtering 
######################################################################

# Control
c[["percent.mt"]] <- PercentageFeatureSet(c, pattern = "^mt-", assay = 'RNA')

counts_per_cell <- Matrix::colSums(c)
counts_per_gene <- Matrix::rowSums(c)
genes_per_cell <- Matrix::colSums(c@assays$RNA@counts >0)

hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')

plot(counts_per_cell, genes_per_cell, log='xy', col='wheat', title('counts vs genes per cell'))
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
plot(sort(counts_per_cell), xlab='cell', log='y', main='counts per cell (ordered)')

VlnPlot(object = c, features = "nCount_RNA")
VlnPlot(object = c, features = "nFeature_RNA")
VlnPlot(object = c, features = "percent.mt")

c <- subset(c, subset = nFeature_RNA > 800 & nFeature_RNA < 4500 &
              nCount_RNA < 15000 & percent.mt < 10)

# High-fat diet (atherosclerosis)
h[["percent.mt"]] <- PercentageFeatureSet(h, pattern = "^mt-", assay = 'RNA')

counts_per_cell <- Matrix::colSums(h)
counts_per_gene <- Matrix::rowSums(h)
genes_per_cell <- Matrix::colSums(h@assays$RNA@counts >0)

hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')

plot(counts_per_cell, genes_per_cell, log='xy', col='wheat', title('counts vs genes per cell'))
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
plot(sort(counts_per_cell), xlab='cell', log='y', main='counts per cell (ordered)')

VlnPlot(object = h, features = "nCount_RNA")
VlnPlot(object = h, features = "nFeature_RNA")
VlnPlot(object = h, features = "percent.mt")

h <- subset(h, subset = nFeature_RNA > 800 & nFeature_RNA < 6000 &
              nCount_RNA < 25000 & percent.mt < 15)

#Add some metadata
l_c <- rep('Control', times = ncol(c))
l_h <- rep('HFD', times = ncol(h))
label <- c(l_c, l_h)

dat <- merge(c, h)
dat <- AddMetaData(dat, label, col.name = 'Group')

######################################################################
#Default workflow, integrate with CCA anchoring and cluster with graph Louvain
######################################################################
sets <- SplitObject(dat, split.by = 'Group')

for(i in 1:length(sets)){
  sets[[i]] <- NormalizeData(sets[[i]])
  sets[[i]] <- FindVariableFeatures(sets[[i]], selection.method = 'disp', nfeatures = 5000)
  sets[[i]] <- ScaleData(sets[[i]], features = rownames(sets[[i]]))
}

dat <- merge(sets[[1]], sets[[2]])
genes <- rownames(dat)

feat <- SelectIntegrationFeatures(object.list = sets, nfeatures = 1000)
anchors <- FindIntegrationAnchors(object.list = sets, anchor.features = feat, dims = 1:30, nn.method = 'annoy')
dat <- IntegrateData(anchorset = anchors, features.to.integrate = genes)
dat <- ScaleData(dat)

dat <- RunPCA(dat, npcs = 100)
ElbowPlot(dat, ndims = 100)

dat <- FindNeighbors(object = dat, dims = 1:30)
dat <- FindClusters(object = dat, algorithm = 2, resolution = 1.5)
dat$pca_clusters <- dat$seurat_clusters

dat <- FindNeighbors(object = dat, features = VariableFeatures(dat))
dat <- FindClusters(object = dat, features = VariableFeatures(dat), algorithm = 2,
                    resolution = 1.5)
dat$gene_clusters <- dat$seurat_clusters

dat <- RunUMAP(object = dat, dims = 1:50)
UMAPPlot(dat, group.by = 'Group') 
UMAPPlot(dat, group.by = 'pca_clusters', pt.size = 0.5) 

FeaturePlot(dat, features = 'S100a9', order = T, min.cutoff = 0, pt.size = 1)
FeaturePlot(dat, features = 'Ikzf2', order = T)
FeaturePlot(dat, features = 'Ccr2', order = T, min.cutoff = 2.5)
FeaturePlot(dat, features = 'Cd209a', order = T, min.cutoff = 2.5)
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
diff <- dbmap$diffusion$Diffusor(n_components = as.integer(60), n_neighbors = as.integer(15),
                                 transitions = as.logical(F),
                                 norm = as.logical(F), ann_dist = as.character('cosine'),
                                 n_jobs = as.integer(8), kernel_use = as.character('simple_adaptive')
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
# Post-diffusion processing (structure clustering, MAP)
##################################################################################

dat <- FindNeighbors(dat, reduction = 'sc', dims = 1:ncol(dat@reductions$sc@cell.embeddings), annoy.metric = 'cosine', graph.name = 'dbgraph')
dat <- FindClusters(dat, resolution = 1, graph.name = 'dbgraph', algorithm = 4)

dat <- RunUMAP(dat, reduction = 'sc', dims = 1:ncol(dat@reductions$sc@cell.embeddings), n.neighbors = 10, init = 'spectral',
               min.dist = 0.45, spread = 1.5, learning.rate = 1.5, n.epochs = 800, reduction.key = 'dbMAP_', reduction.name = 'dbmap')

dat <- RunUMAP(dat, reduction = 'sc', n.components = 3, dims = 1:ncol(dat@reductions$sc@cell.embeddings),
               min.dist = 0.45, spread = 0.8, n.epochs = 200, reduction.key = 'dbMAP_3D', reduction.name = 'dbmap3d')

DimPlot(dat, reduction = 'dbmap', group.by = 'seurat_clusters', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'gene_clusters', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'Group', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'Phase', pt.size = 1)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
dat <- CellCycleScoring(dat, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)


# Convert to h5Seurat, h5ad
library(SeuratDisk)
as.h5Seurat(dat, 'Williams.h5Seurat', overwrite = T)
Convert('Williams.h5Seurat', dest = 'h5ad', overwrite = T, assay = 'RNA')


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
FeaturePlot(dat, reduction = 'dbmap', features = 'Trem2', pt.size = 0.5, min.cutoff = 0, max.cutoff = 2, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'Cx3cr1', pt.size = 0.5, min.cutoff = 0, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'Cd209a', pt.size = 0.5, min.cutoff = 0, order = T) 

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
dat <- getMarkerGenes(dat, organism = 'mm', groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'pred.win', 'pred.coc', 'Group', 'Phase'), min_pct = 0.3, assay = 'integrated', test = 'wilcox')


dat <- getEnrichedPathways(dat, 
                           databases = c('Transcription_Factor_PPIs', 'Chromossome_Location', 'Mouse_Gene_Atlas', 
                                         'MGI_Mammalian_Phenotype_Level_3', 'MGI_Mammalian_Phenotyoe_Level_4_2019', 'MGI_Mammalian_Phenotype_2017',
                                         'OMIM_Disease', 'OMIM_Expanded', 'MsigDB_Computational', 'UK_Biobank_GWAS_v1',
                                         'KEGG_2019_Mouse', 'PheWeb_2019', 'WikiPathways_2019_Mouse', 'GWAS_Catalog_2019',
                                         'GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018',
                                         'Panther_2016', 'BioCarta_2016', 'Reactome_2016', 'Kegg_2016'))
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'Group', 'Phase'),
                                        GMT_file =  'h.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'Group', 'Phase'),
                                        GMT_file =  'c2.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'Group', 'Phase'),
                                        GMT_file =  'c3.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'Group', 'Phase'),
                                        GMT_file =  'c5.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'Group', 'Phase'),
                                        GMT_file =  'c3.tft.v7.1.symbols.gmt')

cerebro <- exportFromSeurat(dat, file = 'Kim.crb', experiment_name = 'Kim et al', organism = 'mm',
                            cell_cycle = 'Phase', groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'pred.win', 'pred.coc', 'Group', 'Phase'),
                            nUMI = 'nCount_RNA', nGene = 'nFeature_RNA')
launchCerebro()


####################################################################
# Annotation
####################################################################
win <- readRDS('Winkels/Winkels.Rds')
transfer.anchors <- FindTransferAnchors(reference = win, query = dat, reduction = 'cca',
                                        dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = win$cluster, weight.reduction = 'cca',
                            dims = 1:30)
dat <- AddMetaData(dat, metadata = predictions)
dat$pred.win <- dat$predicted.id
DimPlot(dat, reduction = 'dbmap', group.by = 'pred.win', pt.size = 0.5)

coc <- readRDS('Cochain/Cochain.Rds')
transfer.anchors <- FindTransferAnchors(reference = coc, query = dat, 
                                        dims = 1:30, reduction = 'cca')
predictions <- TransferData(anchorset = transfer.anchors, refdata = coc$db_clusters, weight.reduction = 'cca',
                            dims = 1:30)
dat <- AddMetaData(dat, metadata = predictions)
dat$pred.coc <- dat$predicted.id
DimPlot(dat, reduction = 'dbmap', group.by = 'pred.coc', pt.size = 0.5)

####################################################################
# Re-name Groups
####################################################################
Idents(dat) <- dat$Group
new.clusters.ids <- c('Chow_Williams', 'HFD_Williams')
names(new.clusters.ids) <- levels(dat$Group)
dat <- RenameIdents(dat, new.clusters.ids)
dat$Group <- Idents(dat)


######################################################################
#Lineage inference with Slingshot
######################################################################
library(slingshot)
#Subset macro
mac <- subset(dat, idents=c("Macrophages", "Monocyte", 'Ly-6C+ mono'))

mac <- RunUMAP(mac, reduction = 'db', dims = 1:(ncol(dat@reductions$db@cell.embeddings)), n.neighbors = 20,
               min.dist = 0.9, reduction.key = 'dbMAP_', reduction.name = 'dbmap')

DimPlot(mac, reduction = 'dbmap', pt.size = 2)
m <- getLineages(mac@reductions$dbmap@cell.embeddings, clusterLabels = mac$db_clusters, 
                 start.clus = 'Monocyte')
m <- getCurves(m, shrink = FALSE, extend = 'n', reassign = FALSE)

#Plotting colors
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = mac)))
names(x = ident.colors) <- levels(x = mac)
cell.colors <- ident.colors[Idents(object = mac)]
names(x = cell.colors) <- colnames(x = mac)

plot(m@reducedDim, col = cell.colors, pch = 16, cex = 1)
lines(m, lwd = 2, type = 'lineages', col = 'black')

plot(m@reducedDim, col = cell.colors, pch = 16, cex = 1)
lines(m, lwd = 2, col = 'black')



#Subset T cells
t_cells <- subset(dat, idents=c("T CD8/Lef1", "T CD8/Gzmk", 'T pro-helper',
                                "T Progenitor", "T intermediate", 'T CD8/Ccl5',
                                "T Dusp10", "Th2", "Th17", "NK"))

t_cells <- RunUMAP(t_cells, reduction = 'db', dims = 1:(ncol(dat@reductions$db@cell.embeddings)), n.neighbors = 30,
                   min.dist = 0.6, reduction.key = 'dbMAP_', reduction.name = 'dbmap')

DimPlot(t_cells, reduction = 'dbmap', pt.size = 1.5)

t <- slingshot(t_cells@reductions$dbmap@cell.embeddings, clusterLabels = t_cells$db_clusters, 
               start.clus = 'T Progenitor')

# Fit GAM
require(gam)
t <- sim$slingPseudotime_1

# for time, only look at the 100 most variable genes
Y <- log1p(assays(sim)$norm)
var100 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
Y <- Y[var100,]

# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  suppressWarnings({
    tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
  })
  p <- summary(tmp)[3][[1]][2,3]
  p
})

#Plotting colors
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = t_cells)))
names(x = ident.colors) <- levels(x = t_cells)
cell.colors <- ident.colors[Idents(object = t_cells)]
names(x = cell.colors) <- colnames(x = t_cells)

plot(t@reducedDim, col = cell.colors, pch = 16, cex = 1)
lines(t, lwd = 2, type = 'lineages', col = 'black')

plot(t@reducedDim, col = cell.colors, pch = 16, cex = 1)
lines(t, lwd = 2, col = 'black')
