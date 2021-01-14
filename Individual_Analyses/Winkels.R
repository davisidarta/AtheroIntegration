######################################################################
#Analysis of Winkels et al data
######################################################################

#Load packages 
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(reticulate)
library(cerebroApp)
setwd("~/Documents/Bioinfo/Vascular/Winkels")

#Load data
chow <- Read10X('chow', gene.column = 1)
hfd <- Read10X('hfd', gene.column = 1)

# Handle the goddamned all capital gene names (WTF, Winkels? Not cool.)
library('biomaRt')
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes_chow <- rownames(chow)
chow_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                          "external_gene_name", "description"),values=genes_chow,mart= mart)
new_chow_genes <- chow_list$external_gene_name
chow <- chow[chow_list$ensembl_gene_id,]
rownames(chow) <- new_chow_genes

genes_hfd <- rownames(hfd)
hfd_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                             "external_gene_name", "description"),values=genes_hfd,mart= mart)
new_hfd_genes <- hfd_list$external_gene_name
hfd <- hfd[hfd_list$ensembl_gene_id,]
rownames(hfd) <- new_hfd_genes

chow <- CreateSeuratObject(chow)
hfd <- CreateSeuratObject(hfd)
######################################################################
#QC and filtering 
######################################################################

#Chow
chow[["percent.mt"]] <- PercentageFeatureSet(chow, pattern = "^mt-", assay = 'RNA')

counts_per_cell <- Matrix::colSums(chow)
counts_per_gene <- Matrix::rowSums(chow)
genes_per_cell <- Matrix::colSums(chow@assays$RNA@counts >0)

hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')

plot(counts_per_cell, genes_per_cell, log='xy', col='wheat', title('counts vs genes per cell'))
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
plot(sort(counts_per_cell), xlab='cell', log='y', main='counts per cell (ordered)')

VlnPlot(object = chow, features = "nCount_RNA")
VlnPlot(object = chow, features = "nFeature_RNA")
VlnPlot(object = chow, features = "percent.mt")

chow <- subset(chow, subset = nFeature_RNA > 600 & nFeature_RNA < 3000 &
                nCount_RNA < 10000 & percent.mt < 10)

#HFD
hfd[["percent.mt"]] <- PercentageFeatureSet(hfd, pattern = "^mt-", assay = 'RNA')

counts_per_cell <- Matrix::colSums(hfd)
counts_per_gene <- Matrix::rowSums(hfd)
genes_per_cell <- Matrix::colSums(hfd@assays$RNA@counts >0)

hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')

plot(counts_per_cell, genes_per_cell, log='xy', col='wheat', title('counts vs genes per cell'))
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
plot(sort(counts_per_cell), xlab='cell', log='y', main='counts per cell (ordered)')

VlnPlot(object = hfd, features = "nCount_RNA")
VlnPlot(object = hfd, features = "nFeature_RNA")
VlnPlot(object = hfd, features = "percent.mt")

hfd <- subset(hfd, subset = nFeature_RNA > 700 & nFeature_RNA < 3000 &
                 nCount_RNA < 10000 & percent.mt < 10)


#Add some metadata
h <- rep('HFD', times = ncol(hfd))
ch <- rep('Chow', times = ncol(chow))
label <- c(h, ch)

dat <- merge(hfd, chow)
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

feat <- SelectIntegrationFeatures(object.list = sets, nfeatures = 3000)
anchors <- FindIntegrationAnchors(object.list = sets, anchor.features = feat)
dat <- IntegrateData(anchorset = anchors, features.to.integrate = genes)
dat <- ScaleData(dat, features = rownames(dat))

dat <- FindVariableFeatures(dat, selection.method = 'disp', nfeatures = 4000)
VariableFeaturePlot(dat)

dat <- RunPCA(dat, npcs = 100)
ElbowPlot(dat, ndims = 100)

dat <- FindNeighbors(object = dat, dims = 1:30)
dat <- FindClusters(object = dat, algorithm = 2)
dat$pca_clusters <- dat$seurat_clusters


dat <- RunUMAP(object = dat, dims = 1:30)
UMAPPlot(dat, group.by = 'Group') 
UMAPPlot(dat, group.by = 'pca_clusters', pt.size = 0.5) 

FeaturePlot(dat, features = 'Mki67', order = T)
FeaturePlot(dat, features = 'IFITM1', order = T)
FeaturePlot(dat, features = 'CCR7', order = T, min.cutoff = 2.5)


###############################################################################
# Run dbMAP
###############################################################################
library(reticulate)
np <- reticulate::import("numpy")
pd <- reticulate::import("pandas")
sp <- reticulate::import("scipy")
dm <- reticulate::import('dbmap')

genes <- dat@assays$integrated@var.features
data <- t(dat@assays$integrated@data[genes,])
a <- r_to_py(data)
b <- a$tocsr()

diff = dm$diffusion$Diffusor(n_components = as.integer(300),
                             n_neighbors = as.integer(15),
                             ann_dist = as.character('cosine'),
                             n_jobs = as.integer(10),
                             kernel_use = 'simple',
                             transitions = 'False',
                             norm = 'False')
diff = diff$fit(b)
mms = diff$transform(b)
res = diff$return_dict()

sc = py_to_r(res$StructureComponents)
ev = py_to_r(res$EigenValues)

################
rownames(sc) <- colnames(dat)
new_names <- list()
for(i in 1:length(sc)){
  new_names[i] <- paste('SC_' , as.integer(colnames(sc[i])) + 1, sep = '')
}
colnames(sc) <- as.vector(new_names)
names(ev) <- as.vector(new_names)
################

plot(ev)

dat@reductions$sc <- dat@reductions$pca 
dat@reductions$sc@cell.embeddings <- as.matrix(sc)
dat@reductions$sc@feature.loadings <- matrix(data = c(0,0))
dat@reductions$sc@key <- 'SC_'

######################################################################
# diffusion-based Manifold Approximation and Projection
######################################################################
dat <- RunUMAP(dat, reduction = 'sc', dims = 1:ncol(dat@reductions$sc@cell.embeddings), n.neighbors = 30, n.epochs = 1000, 
               min.dist = 0.6, spread = 1.2, learning.rate = 0.8, reduction.key = 'dbMAP_', reduction.name = 'dbmap')

dat <- RunUMAP(dat, reduction = 'sc', n.components = 3, dims = 1:ncol(dat@reductions$sc@cell.embeddings), n.neighbors = 15, n.epochs = 1000,
               min.dist = 0.8, spread = 1.5, learning.rate = 0.8, reduction.key = 'dbMAP3D_', reduction.name = 'dbmap3d')

# Constroi o grafo baseado nos componentes estruturais
dat <- FindNeighbors(dat, reduction = 'sc', dims = 1:(ncol(dat@reductions$sc@cell.embeddings)), 
                     graph.name = 'sc_graph', k.param = 10)
# Clusteriza usando o algoritmo de leiden nesse grafo
dat <- FindClusters(dat, method = 2, resolution = 0.3, graph.name = 'sc_graph')
dat$sc_clusters <- dat$seurat_clusters # Vamos guardar como sc_clusters
######################################################################
#Plot and save
######################################################################

DimPlot(dat, reduction = 'dbmap', group.by = 'Group', pt.size = 0.5)
DimPlot(dat, reduction = 'dbmap', group.by = 'pca_clusters', pt.size = 0.5)
DimPlot(dat, reduction = 'dbmap', group.by = 'db_clusters', pt.size = 0.5)
DimPlot(dat, reduction = 'dbmap', group.by = 'seurat_clusters', pt.size = 0.5)

FeaturePlot(dat, reduction = 'dbmap', features = 'MKI67', pt.size = 0.5, min.cutoff = 0, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'IFITM1', pt.size = 0.5, min.cutoff = 0, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'CCR7', pt.size = 0.5, min.cutoff = 2.5, order = T) 

  plot.data <- FetchData(object = dat, vars = c("dbMAP3D_1", "dbMAP3D_2", "dbMAP3D_3", 
                                                "pca_clusters", 'db_clusters', 'Group', 'seurat_clusters', 
                                                'MKI67', 'IFITM1', 'CCR7'
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
  


######################################################################
#Export to Cerebro
######################################################################
dat@assays$integrated@counts <- dat@assays$integrated@data
dat <- addPercentMtRibo(dat, organism = 'hg', gene_nomenclature = 'name', assay = 'integrated')
dat <- getMarkerGenes(dat, organism = 'hg', column_sample = 'Group', column_cluster = 'db_clusters', assay = 'integrated', min_pct = 0.7)
dat <- getEnrichedPathways(dat, column_sample = 'Group', column_cluster = 'db_clusters',
                           databases = c('Transcription_Factor_PPIs', 'Chromossome_Location', 'Mouse_Gene_Atlas', 
                                         'MGI_Mammalian_Phenotype_Level_3', 'MGI_Mammalian_Phenotyoe_Level_4_2019', 'MGI_Mammalian_Phenotype_2017',
                                         'OMIM_Disease', 'OMIM_Expanded', 'MsigDB_Computational', 'UK_Biobank_GWAS_v1',
                                         'KEGG_2019_Mouse', 'PheWeb_2019', 'WikiPathways_2019_Mouse', 'GWAS_Catalog_2019',
                                         'GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018',
                                         'Panther_2016', 'BioCarta_2016', 'Reactome_2016', 'Kegg_2016'))
cerebro <- exportFromSeurat(dat, file = 'Winkels.crb', experiment_name = 'Winkels et al', organism = 'hg',
                            column_sample = 'Group', column_cluster = 'db_clusters',
                            column_nUMI = 'nCount_RNA', column_nGene = 'nFeature_RNA', assay = 'integrated')
launchCerebro(maxFileSize = 500000)

#Save as RDS
saveRDS(dat, 'Winkels.Rds')

####################################################################
# Re-name clusters
####################################################################
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

####################################################################
# Re-name Groups
####################################################################
Idents(dat) <- dat$Group
new.clusters.ids <- c('Chow_Winkels', 'HFD_Winkels')
names(new.clusters.ids) <- levels(dat$Group)
dat <- RenameIdents(dat, new.clusters.ids)
dat$Group <- Idents(dat)


# Convert to h5Seurat, h5ad
library(SeuratDisk)
as.h5Seurat(dat, 'Winkels.h5Seurat', overwrite = T)
Convert('Winkels.h5Seurat', dest = 'h5ad', overwrite = T, assay = 'RNA')


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
