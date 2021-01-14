######################################################################
#Analysis of Lin et al data
######################################################################

#Load packages 
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(reticulate)
library(cerebroApp)
setwd("~/Documents/Bioinfo/Vascular/Lin")

#Load data
prog <- Read10X('progression')
reg <- Read10X('regression')

prog <- CreateSeuratObject(prog)
reg <- CreateSeuratObject(reg)

######################################################################
#QC and filtering 
######################################################################

# Progression group
prog[["percent.mt"]] <- PercentageFeatureSet(prog, pattern = "^mt-", assay = 'RNA')

counts_per_cell <- Matrix::colSums(prog)
counts_per_gene <- Matrix::rowSums(prog)
genes_per_cell <- Matrix::colSums(prog@assays$RNA@counts >0)

hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')

plot(counts_per_cell, genes_per_cell, log='xy', col='wheat', title('counts vs genes per cell'))
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
plot(sort(counts_per_cell), xlab='cell', log='y', main='counts per cell (ordered)')

VlnPlot(object = prog, features = "nCount_RNA")
VlnPlot(object = prog, features = "nFeature_RNA")
VlnPlot(object = prog, features = "percent.mt")

prog <- subset(prog, subset = nFeature_RNA > 600 & nFeature_RNA < 2500 &
                 nCount_RNA < 8000 & percent.mt < 10)

# Regression group
reg[["percent.mt"]] <- PercentageFeatureSet(reg, pattern = "^mt-", assay = 'RNA')

counts_per_cell <- Matrix::colSums(reg)
counts_per_gene <- Matrix::rowSums(reg)
genes_per_cell <- Matrix::colSums(reg@assays$RNA@counts >0)

hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')

plot(counts_per_cell, genes_per_cell, log='xy', col='wheat', title('counts vs genes per cell'))
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
plot(sort(counts_per_cell), xlab='cell', log='y', main='counts per cell (ordered)')

VlnPlot(object = reg, features = "nCount_RNA")
VlnPlot(object = reg, features = "nFeature_RNA")
VlnPlot(object = reg, features = "percent.mt")

reg <- subset(reg, subset = nFeature_RNA > 600 & nFeature_RNA < 2500 &
                  nCount_RNA < 8000 & percent.mt < 10)

#Add some metadata
lp <- rep('Progression', times = ncol(prog))
lr <- rep('Regression', times = ncol(reg))
label <- c(lp, lr)

dat <- merge(prog, reg)
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
UMAPPlot(dat, group.by = 'Group', pt.size = 0.5) 

FeaturePlot(dat, features = 'Mki67', order = T)
FeaturePlot(dat, features = 'Clec10a', order = T) #CD301
FeaturePlot(dat, features = 'Pdcd1lg2', order = T, min.cutoff = 0, max.cutoff = 0.01) #PDL2

######################################################################
#Embedd with dbMAP
######################################################################
dbmap <- reticulate::import('dbmap')
pd <- reticulate::import('pandas')
umap <- reticulate::import('umap')
genes <- dat@assays$integrated@var.features
data <- t(dat@assays$integrated@data[genes,])
data <- as.sparse(data)

data <- r_to_py(data)
data <- data$tocoo()

diff <- dbmap$diffusion$diffuse(data, n_components = as.integer(100), knn = as.integer(20))
evals <- diff$EigenValues
print(diff$Suggested_eigs)
plot(evals) #Select meaningful diffusion components. Used 60 (automated).
res <- dbmap$multiscale$Multiscale(diff)
db <- as.matrix(res)

#Add to Seurat
dat@reductions$db <- dat@reductions$pca
rownames(db) <- colnames(dat)
dat@reductions$db@cell.embeddings <- db

dat <- FindNeighbors(dat, reduction = 'db', dims = 1:(ncol(dat@reductions$db@cell.embeddings)), graph.name = 'dbgraph')
dat <- FindClusters(dat, resolution = 0.8, graph.name = 'dbgraph')
dat$db_clusters <- dat$seurat_clusters

######################################################################
#dbMAP  Adjustment
######################################################################
dat <- RunUMAP(dat, reduction = 'db', dims = 1:(ncol(dat@reductions$db@cell.embeddings)), min.dist = 0.45, reduction.key = 'dbMAP_', reduction.name = 'dbmap')

dat <- RunUMAP(dat, reduction = 'db', n.components = 3, dims = 1:ncol(dat@reductions$db@cell.embeddings), min.dist = 0.3, spread = 2, learning.rate = 2, reduction.key = 'dbMAP3D_', reduction.name = 'dbmap3d', init = 'spectral')

dat <- FindNeighbors(dat, reduction = 'dbmap3d', dims = 1:3, graph.name = 'db3d')
dat <- FindClusters(dat, resolution = 0.8, graph.name = 'db3d')

######################################################################
#Plot and save
######################################################################

DimPlot(dat, reduction = 'dbmap', group.by = 'Group', pt.size = 0.5)
DimPlot(dat, reduction = 'dbmap', group.by = 'pca_clusters', pt.size = 0.5)

DimPlot(dat, reduction = 'dbmap', group.by = 'db_clusters', pt.size = 0.5,
        label = F, repel = T, label.size = 5)

FeaturePlot(dat, features = 'Mki67', order = T, reduction = 'dbmap')
FeaturePlot(dat, features = 'Clec10a', order = T, reduction = 'dbmap') #CD301
FeaturePlot(dat, features = 'Pdcd1lg2', order = T, min.cutoff = 0, max.cutoff = 0.01 , reduction = 'dbmap') #PDL2

FeaturePlot(dat, features = 'Cd9', order = T, min.cutoff = 0, reduction = 'dbmap') #PDL2
FeaturePlot(dat, features = 'Pdcd1lg2', order = T, min.cutoff = 0, max.cutoff = 0.01 , reduction = 'dbmap') #PDL2
plot.data$label <- dat$seurat_clusters
plot_ly(data = plot.data, 
        x = ~dbMAP3D_1, y = ~dbMAP3D_2, z = ~dbMAP3D_3, 
        color = ~seurat_clusters, 
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 3, width=0.5),
        text=~seurat_clusters,
        hoverinfo="text", plot_bgcolor = 'black') 

# Convert to h5Seurat, h5ad
library(SeuratDisk)
as.h5Seurat(dat, 'Lin.h5Seurat', overwrite = T)
Convert('Lin.h5Seurat', dest = 'h5ad', overwrite = T, assay = 'integrated')

######################################################################
#Export to Cerebro
######################################################################
dat@assays$integrated@counts <- dat@assays$integrated@data
dat <- addPercentMtRibo(dat, organism = 'mm', gene_nomenclature = 'name', assay = 'integrated')
dat <- getMarkerGenes(dat, organism = 'mm', column_sample = 'Group', column_cluster = 'db_clusters', assay = 'integrated', min_pct = 0.7)
dat <- getEnrichedPathways(dat, column_sample = 'Group', column_cluster = 'db_clusters',
                           databases = c('Transcription_Factor_PPIs', 'Chromossome_Location', 'Mouse_Gene_Atlas', 
                                         'MGI_Mammalian_Phenotype_Level_3', 'MGI_Mammalian_Phenotyoe_Level_4_2019', 'MGI_Mammalian_Phenotype_2017',
                                         'OMIM_Disease', 'OMIM_Expanded', 'MsigDB_Computational', 'UK_Biobank_GWAS_v1',
                                         'KEGG_2019_Mouse', 'PheWeb_2019', 'WikiPathways_2019_Mouse', 'GWAS_Catalog_2019',
                                         'GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018',
                                         'Panther_2016', 'BioCarta_2016', 'Reactome_2016', 'Kegg_2016'))
cerebro <- exportFromSeurat(dat, file = 'Lin.crb', experiment_name = 'Lin et al', organism = 'mm',
                            column_sample = 'Group', column_cluster = 'db_clusters',
                            column_nUMI = 'nCount_RNA', column_nGene = 'nFeature_RNA', assay = 'integrated')

launchCerebro(maxFileSize = 500000)

#Save as RDS
saveRDS(dat, 'Lin.Rds')


####################################################################
# Re-name clusters
####################################################################
Idents(dat) <- dat$db_clusters
new.clusters.ids <- c('NMES', #0
                      'NA 1', #1
                      'Folr2-high 1', #2
                      'NA 2', #3
                      'Intermediate Stem-like', #4
                      'Cd74/Mhc2-high', #5
                      'Folr2-high 2', #6
                      'Res-high', #7
                      'HSP', #8
                      'Trem2-high 1', #9
                      'INF-sign high', #10
                      'Chemokine high', #11
                      'Trem2-high 2', #12
                      'NA 3', #13
                      'Stem', #14
                      'DNAse 1/3') #15
names(new.clusters.ids) <- levels(dat$db_clusters)
dat <- RenameIdents(dat, new.clusters.ids)
dat$db_clusters <- Idents(dat)

####################################################################
# Select only clusters of interest
####################################################################
dat1 <- dat

dat <- subset(dat, idents = 'Chemokine high')

RidgePlot(dat, features = c('Ldha', 'Ldhb', 'L2hgdh', 'Tet2'), group.by = 'Group')

VlnPlot(dat1,features = c('Ldha', 'Ldhb', 'L2hgdh', 'Tet2'), group.by = 'db_clusters', 
        ncol=2, same.y.lims = T)

DotPlot(dat1,features = c('Ldha', 'Ldhb', 'L2hgdh', 'Tet2'), group.by = 'db_clusters', 
        col.min = 0, dot.scale = 15)


#### Differential gene expression: progression vs. regression
dat <- readRDS("~/Documents/Bioinfo/Ebru/Lin/Lin.Rds")
DefaultAssay(dat) <- 'integrated'
Idents(dat) <- 'Group'
prog.markers <- FindMarkers(dat, ident.1 = 'Progression', ident.2 = 'Regression', logfc.threshold = 0.1, test.use = 'LR')
prog.markers$log2FC <- prog.markers$avg_logFC
write.csv(prog.markers, file = 'Prog_VS_Reg_markers.csv')
prog.markers <- prog.markers[-2]
write.table(prog.markers, file = 'Prog_VS_Reg_markers.tsv', sep = '\t', row.names = F)
prog.markers = prog.markers[order(prog.markers$log2FC),]

reg.markers <- FindMarkers(dat, ident.1 = 'Regression', ident.2 = 'Progression', logfc.threshold = 0.1)
reg.markers$log2FC <- reg.markers$avg_logFC
write.csv(reg.markers, file = 'Reg_VS_Prog_markers.csv')

DefaultAssay(dat) <- 'RNA'
prog <- subset(dat, idents = 'Progression')
reg <- subset(dat, idents = 'Regression')

phgdh_p <- FetchData(prog, vars = 'Phgdh')
phgdh_r <- FetchData(reg, vars = 'Phgdh')
p_value_phgdh <- wilcox.test(phgdh_p$Phgdh, phgdh_r$Phgdh)
p_value_phgdh

psat1_p <- FetchData(prog, vars = 'Psat1')
psat1_r <- FetchData(reg, vars = 'Psat1')
p_value_psat1 <- wilcox.test(psat1_p$Psat1, psat1_r$Psat1)
p_value_psat1

Psph_p <- FetchData(prog, vars = 'Psph')
Psph_r <- FetchData(reg, vars = 'Psph')
p_value_Psph <- wilcox.test(Psph_p$Psph, Psph_r$Psph)
p_value_Psph

Shmt1_p <- FetchData(prog, vars = 'Shmt1')
Shmt1_r <- FetchData(reg, vars = 'Shmt1')
p_value_Shmt1 <- wilcox.test(Shmt1_p$Shmt1, Shmt1_r$Shmt1)
p_value_Shmt1

Shmt2_p <- FetchData(prog, vars = 'Shmt2')
Shmt2_r <- FetchData(reg, vars = 'Shmt2')
p_value_Shmt2 <- wilcox.test(Shmt2_p$Shmt2, Shmt2_r$Shmt2)
p_value_Shmt2
  
FeaturePlot(dat, features = c('Phgdh', 'Psat1', 'Psph'), reduction = 'dbmap', order = T, min.cutoff = 0, max.cutoff = 1)
VlnPlot(dat, features = c('Phgdh', 'Psat1', 'Psph'), group.by = 'Group')

FeaturePlot(dat, features = c('Shmt1', 'Shmt2'), reduction = 'dbmap', order = T, min.cutoff = 0, max.cutoff = 1)
VlnPlot(dat, features = c('Shmt1', 'Shmt2'), group.by = 'Group')

FeaturePlot(dat, features = c('Phgdh', 'Psat1', 'Psph','Shmt1', 'Shmt2'), reduction = 'dbmap', order = T, min.cutoff = 0, max.cutoff = 1)
VlnPlot(dat, features = c('Phgdh', 'Psat1', 'Psph','Shmt1', 'Shmt2'), group.by = 'Group')

