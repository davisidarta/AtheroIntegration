######################################################################
#Analysis of Cochain et al data
######################################################################

#Load packages 
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(reticulate)
library(cerebroApp)
setwd("~/Documents/Bioinfo/Vascular/Cochain")

#Load data
chow <- Read10X('Cochain_2018_healthy_aorta')
hfd11 <- Read10X('Cochain_2018_11weeks_HFD')
hfd20 <- Read10X('Cochain_2018_20weeks_HFD')

chow <- CreateSeuratObject(chow)
hfd11 <- CreateSeuratObject(hfd11)
hfd20 <- CreateSeuratObject(hfd20)

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

#HFD 11w
hfd11[["percent.mt"]] <- PercentageFeatureSet(hfd11, pattern = "^mt-", assay = 'RNA')

counts_per_cell <- Matrix::colSums(hfd11)
counts_per_gene <- Matrix::rowSums(hfd11)
genes_per_cell <- Matrix::colSums(hfd11@assays$RNA@counts >0)

hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')

plot(counts_per_cell, genes_per_cell, log='xy', col='wheat', title('counts vs genes per cell'))
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
plot(sort(counts_per_cell), xlab='cell', log='y', main='counts per cell (ordered)')

VlnPlot(object = hfd11, features = "nCount_RNA")
VlnPlot(object = hfd11, features = "nFeature_RNA")
VlnPlot(object = hfd11, features = "percent.mt")

hfd11 <- subset(hfd11, subset = nFeature_RNA > 700 & nFeature_RNA < 3500 &
                nCount_RNA < 15000 & percent.mt < 10)

#HFD 20w
hfd20[["percent.mt"]] <- PercentageFeatureSet(hfd20, pattern = "^mt-", assay = 'RNA')

counts_per_cell <- Matrix::colSums(hfd20)
counts_per_gene <- Matrix::rowSums(hfd20)
genes_per_cell <- Matrix::colSums(hfd20@assays$RNA@counts >0)

hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')

plot(counts_per_cell, genes_per_cell, log='xy', col='wheat', title('counts vs genes per cell'))
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
plot(sort(counts_per_cell), xlab='cell', log='y', main='counts per cell (ordered)')

VlnPlot(object = hfd20, features = "nCount_RNA")
VlnPlot(object = hfd20, features = "nFeature_RNA")
VlnPlot(object = hfd20, features = "percent.mt")

hfd20 <- subset(hfd20, subset = nFeature_RNA > 400 & nFeature_RNA < 1800 &
                  nCount_RNA < 5000 & percent.mt < 10)

#Add some metadata
ch <- rep('Chow', times = ncol(chow))
h11 <- rep('HFD 11w', times = ncol(hfd11))
h20 <- rep('HFD 20w', times = ncol(hfd20))
label <- c(ch, h11, h20)

dat <- merge(chow, list(hfd11, hfd20))
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

dat <- merge(sets[[1]], list(sets[[2]], sets[[3]]))
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
FeaturePlot(dat, features = 'Ifitm1', order = T)
FeaturePlot(dat, features = 'Ccr7', order = T, min.cutoff = 2.5)

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

diff <- dbmap$diffusion$diffuse(data, n_components = as.integer(100), knn = as.integer(15))
evals <- diff$EigenValues
print(diff$Suggested_eigs)
plot(evals) #Select meaningful diffusion components. Used 64 (automated).
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
dat <- RunUMAP(dat, reduction = 'db', dims = 1:(ncol(dat@reductions$db@cell.embeddings)), min.dist = 0.8, reduction.key = 'dbMAP_', reduction.name = 'dbmap')

dat <- RunUMAP(dat, reduction = 'db', n.components = 3, dims = 1:ncol(dat@reductions$db@cell.embeddings), min.dist = 1, spread = 1.5, learning.rate = 2, reduction.key = 'dbMAP3D_', reduction.name = 'dbmap3d', init = 'spectral')

dat <- FindNeighbors(dat, reduction = 'dbmap3d', dims = 1:3, graph.name = 'db3d')
dat <- FindClusters(dat, resolution = 0.8, graph.name = 'db3d')

######################################################################
#Plot
######################################################################

DimPlot(dat, reduction = 'dbmap', group.by = 'Group', pt.size = 0.5)
DimPlot(dat, reduction = 'dbmap', group.by = 'pca_clusters', pt.size = 0.5)
DimPlot(dat, reduction = 'dbmap', group.by = 'db_clusters', pt.size = 1, label = T, repel = T, label.size = 5)

FeaturePlot(dat, features = 'Mki67', order = T, reduction = 'dbmap')
FeaturePlot(dat, features = 'Ifitm1', order = T, reduction = 'dbmap')
FeaturePlot(dat, features = 'Ccr7', order = T, min.cutoff = 1.5, reduction = 'dbmap')

plot.data <- FetchData(object = dat, vars = c("dbMAP3D_1", "dbMAP3D_2", "dbMAP3D_3", 
                                              "pca_clusters", 'db_clusters', 'Group', 'seurat_clusters', 
                                              'Mki67', 'Ifitm1', 'Ccr7', 'Pcna'
))
plot.data$label <- dat$db_clusters
plot_ly(data = plot.data, 
        x = ~dbMAP3D_1, y = ~dbMAP3D_2, z = ~dbMAP3D_3, 
        color = ~db_clusters, 
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 3, width=0.5),
        text=~db_clusters,
        hoverinfo="text") 

# Convert to h5Seurat, h5ad
library(SeuratDisk)
as.h5Seurat(dat, 'Cochain.h5Seurat', overwrite = T)
Convert('Cochain.h5Seurat', dest = 'h5ad', overwrite = T, assay = 'RNA')


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
cerebro <- exportFromSeurat(dat, file = 'Cochain.crb', experiment_name = 'Cochain et al', organism = 'mm',
                            column_sample = 'Group', column_cluster = 'db_clusters',
                            column_nUMI = 'nCount_RNA', column_nGene = 'nFeature_RNA', assay = 'integrated')
launchCerebro(maxFileSize = 500000)

#Save as RDS
saveRDS(dat, 'Cochain.Rds')


####################################################################
# Re-name clusters
####################################################################
Idents(dat) <- dat$db_clusters
new.clusters.ids <- c('Inflammatory Macro',
                      'T CD8/NK',
                      'Macro Res-like',
                      'Mono',
                      'T CXCR6-high',
                      'Mono-DC',
                      'T Progenitors',
                      'B Cells',
                      'T CD8',
                      'Macro Trem2-high',
                      'Mast',
                      'Macro FOLR2-high',
                      'Macro Resistin-high',
                      'Granulocytes',
                      'NA',
                      'DC',
                      'NK',
                      'DC activated')
names(new.clusters.ids) <- levels(dat$db_clusters)
dat <- RenameIdents(dat, new.clusters.ids)
dat$db_clusters <- Idents(dat)

####################################################################
# Re-name Groups
####################################################################
Idents(dat) <- dat$Group
new.clusters.ids <- c('Chow_Cochain', 'HFD_11w_Cochain', 'HFD_20w_Cochain')
names(new.clusters.ids) <- levels(dat$Group)
dat <- RenameIdents(dat, new.clusters.ids)
dat$Group <- Idents(dat)


######################################################################
#Lineage inference with Slingshot
######################################################################
library(slingshot)
#Subset macro
mac <- subset(dat, idents=c("Inflammatory Macro", "Macro Res-like", 'Mono', 'Macro Trem2-high',
                             'Mono-DC', 'Macro FOLR2-high', 'Macro Resistin-high',
                             'Granulocytes', 'NA', 'DC'))

mac <- RunUMAP(mac, reduction = 'db', dims = 1:(ncol(dat@reductions$db@cell.embeddings)), n.neighbors = 15,
                   min.dist = 0.6, spread = 1.2, reduction.key = 'dbMAP_', reduction.name = 'dbmap')

DimPlot(mac, reduction = 'dbmap', pt.size = 2, label = T, repel = T, label.size = 6)



# Glycolytic genes
FeaturePlot(mac, features = c('Aldoa', 'Aldoc', 'Bid', 'Cd4', 'Eno1', 
                              'Eno2', 'Eno3', 'Gpi1',
                              'Hk1', 'Pfkfb1', 'Pfkfb2', 'Pfkm', 'Pgk1',
                              'Pgm1', 'Pkm', 'Tpi1' ),
            reduction = 'dbmap', pt.size = 1, order = T)

DotPlot(mac, features= c('Aldoa', 'Aldoc', 'Bid', 'Cd4', 'Eno1', 
                         'Eno2', 'Eno3', 'Gpi1',
                         'Hk1', 'Pfkfb1', 'Pfkfb2', 'Pfkm', 'Pgk1',
                         'Pgm1', 'Pkm', 'Tpi1' ), col.min = 0)
# OXPHOS genes
DotPlot(mac, features= c('Atp5b', 'Atp5c1', 'Atp5d', 'Atp5e', 'Atp5g1', 'Atp5g2',
                         'Atp5g3', 'Atp5k', 'Atp5j2', 'Atp5l', 'Atp5h', 'Atp5j', 
                         'Atp5o', 'Bcl2l1', 'Bcs1l', 'Cbarp', 'Cnot7', 'Cox10', 
                         'Cox11', 'Cox15', 'Cox4i1', 'Cox5a', 'Cox5b', 'Cox6a1',
                         'Cox6a2', 'Cox6b1', 'Cox6c', 'Cox7a1', 'Cox7a2', 'Cox7a2l',
                         'Cox7b', 'Cox7c', 'Cox8a', 'Cyb5a', 'Cyc1', 'Cycs', 'Gatb',
                         'Mbd3', 'Ndufa1', 'Ndufa10', 'Ndufa13', 'Ndufa2', 'Ndufa3',
                         'Ndufa4', 'Ndufa4l2', 'Ndufa5', 'Ndufa6', 'Ndufa7', 'Ndufa8',
                         'Ndufa9', 'Ndufab1', 'Ndufb2', 'Ndufb3', 'Ndufb4', 'Ndufb5',
                         'Ndufb6', 'Ndufb7', 'Ndufb8', 'Ndufc1', 'Ndufc2', 'Ndufs1',
                         'Ndufs2', 'Ndufs3', 'Ndufs4', 'Ndufs5', 'Ndufs6', 'Ndufs7',
                         'Ndufs8', 'Sdha', 'Sdhb', 'Sdhc', 'Sdhd', 'Surf1', 'Surf2',
                         'Tbcb', 'Tmco6', 'Uqcr10', 'Uqcr11', 'Uqcrb', 'Uqcrc1', 
                         'Uqcrc2', 'Uqcrfs1', 'Uqcrh', 'Uqcrq'), col.min = 0, dot.scale = 3) + RotatedAxis()


# TCA genes
DotPlot(mac, features= c('Aco2', 'Cs', 'Fh1', 'Idh2', 'Idh3a', 'Idh3b', 'Idh3g',
                         'Mdh2', 'Ogdh', 'Sdha', 'Sdhb', 'Sdhc', 'Sdhd', 'Sucla2', 'Suclg1' ), col.min = 0, dot.scale = 6) + RotatedAxis()




m <- getLineages(mac@reductions$dbmap@cell.embeddings, clusterLabels = mac$db_clusters, 
                   start.clus = 'Macro FOLR2-high')
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
t_cells <- subset(dat, idents=c("T CD8/NK", "T CXCR6-high", 'T Progenitors',
                            'T CD8', 'NK'))

t_cells <- RunUMAP(t_cells, reduction = 'db', dims = 1:(ncol(dat@reductions$db@cell.embeddings)), n.neighbors = 30,
                   min.dist = 1.2, reduction.key = 'dbMAP_', reduction.name = 'dbmap')

DimPlot(t_cells, reduction = 'dbmap', pt.size = 2, label = T, repel = T, label.size = 6)
FeaturePlot(t_cells, features = 'Mki67', reduction = 'dbmap', pt.size = 2, label = T, repel = T, order = T)

t <- getLineages(t_cells@reductions$dbmap@cell.embeddings, clusterLabels = t_cells$db_clusters, 
                 start.clus = 'T Progenitors')
t <- getCurves(t, shrink = FALSE, extend = 'n', reassign = FALSE)

#Plotting colors
ident.colors <- (scales::hue_pal())(n = length(x = levels(x = t_cells)))
names(x = ident.colors) <- levels(x = t_cells)
cell.colors <- ident.colors[Idents(object = t_cells)]
names(x = cell.colors) <- colnames(x = t_cells)

plot(t@reducedDim, col = cell.colors, pch = 16, cex = 1)
lines(t, lwd = 2, type = 'lineages', col = 'black')

plot(t@reducedDim, col = cell.colors, pch = 16, cex = 1)
lines(t, lwd = 2, col = 'black')





