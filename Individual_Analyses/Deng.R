######################################################################
#Analysis of Deng et al data
######################################################################

#Load packages 
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(reticulate)
library(cerebroApp)
setwd("~/Documents/Bioinfo/Ebru/Deng/")

#Load data
wt_adv <- Read10X_h5('data/GSM4186882_WTADV_filtered_gene_bc_matrices_h5.h5')
wt_me <- Read10X_h5('data/GSM4186884_WTME_filtered_gene_bc_matrices_h5.h5')
apoeKO_adv <- Read10X_h5('data/GSM4186883_APOEADV_filtered_gene_bc_matrices_h5.h5')
apoeKO_me <- Read10X_h5('data/GSM4186885_APOEME_filtered_gene_bc_matrices_h5.h5')

wt_adv <- CreateSeuratObject(wt_adv,  project = 'WT Adventitia')
wt_me <- CreateSeuratObject(wt_me, project = 'WT Media')
apoeKO_adv <- CreateSeuratObject(apoeKO_adv, project = 'Apoe-KO Adventitia ')
apoeKO_me <- CreateSeuratObject(apoeKO_me, project = 'Apoe-KO Media')

######################################################################
#QC and filtering 
######################################################################
dats <- list(wt_adv, wt_me, apoeKO_adv, apoeKO_me)

# This is the filtered public-released data. Let's check how filtered it is.
for(i in 1:length(dats)){
  seurat <- dats[[i]]
  
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
  
  dats[[i]] <- seurat
}  

VlnPlot(object = dats[[1]], features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), pt.size = 0.01)
VlnPlot(object = dats[[2]], features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), pt.size = 0.01)
VlnPlot(object = dats[[3]], features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), pt.size = 0.01)
VlnPlot(object = dats[[4]], features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), pt.size = 0.01)

# WT Adventitia
dats[[1]] <- subset(dats[[1]], subset = nFeature_RNA > 500 & nFeature_RNA < 4000 &
                      nCount_RNA < 15000 & percent.mt < 10)
# WT Media
dats[[2]] <- subset(dats[[2]], subset = nFeature_RNA > 500 & nFeature_RNA < 5000 &
                      nCount_RNA < 15000 & percent.mt < 20)
# ApoeKO Aventitia
dats[[3]] <- subset(dats[[3]], subset = nFeature_RNA > 500 & nFeature_RNA < 3500 &
                      nCount_RNA < 15000 & percent.mt < 10)
# ApoeKO Media
dats[[4]] <- subset(dats[[4]], subset = nFeature_RNA > 500 & nFeature_RNA < 4000 &
                      nCount_RNA < 15000 & percent.mt < 20)

#Add some metadata
l_wt_ad <- rep('WT Adventitia', times = ncol(dats[[1]]))
l_wt_me <- rep('WT Media', times = ncol(dats[[2]]))
l_apoeko_ad <- rep('ApoeKO Adventitia', times = ncol(dats[[3]]))
l_apoeko_me <- rep('ApoeKO Media', times = ncol(dats[[4]]))
label <- c(l_wt_ad, l_wt_me, l_apoeko_ad, l_apoeko_me)

dat <- merge(dats[[1]], dats[-1])
dat <- AddMetaData(dat, label, col.name = 'Sample')

######################################################################
#Default workflow, integrate with CCA anchoring and cluster with graph Louvain
######################################################################
sets <- SplitObject(dat, split.by = 'Sample')
genes <- rownames(dat)

for(i in 1:length(sets)){
  sets[[i]] <- NormalizeData(sets[[i]])
  sets[[i]] <- FindVariableFeatures(sets[[i]], selection.method = 'disp', nfeatures = 3000, num.bin = 100)
  sets[[i]] <- ScaleData(sets[[i]], features = rownames(sets[[i]]))
}

feat <- SelectIntegrationFeatures(object.list = sets, nfeatures = 3000)
anchors <- FindIntegrationAnchors(object.list = sets, anchor.features = feat, 
                                  normalization.method = 'LogNormalize', dims = 1:50, 
                                  max.features = 800, nn.method = 'annoy', reduction = 'cca')
dat <- IntegrateData(anchorset = anchors, features.to.integrate = genes)
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

dat <- RunUMAP(object = dat, dims = 1:30)
UMAPPlot(dat, group.by = 'Sample', pt.size = 0.5) 

UMAPPlot(dat, group.by = 'pca_clusters', pt.size = 0.5) 
UMAPPlot(dat, group.by = 'gene_clusters', pt.size = 0.5) 

FeaturePlot(dat, features = 'Mki67', order = T, min.cutoff = 0, max.cutoff = 1, pt.size = 1)
FeaturePlot(dat, features = 'Ikzf2', order = T)
FeaturePlot(dat, features = 'Ccr2', order = T, min.cutoff = 2.5)
FeaturePlot(dat, features = 'Cd209a', order = T, min.cutoff = 0, max.cutoff = 2, pt.size = 1)
FeaturePlot(dat, features = 'Fabp5', order = T, min.cutoff = 0, max.cutoff = 4)

VlnPlot(dat, features = c('Ikzf2', 'Cd209a', 'Ccr2'), group.by = 'Sample', assay = 'RNA') 

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
diff <- dbmap$diffusion$Diffusor(n_components = as.integer(120), n_neighbors = as.integer(15),
                                 transitions = as.logical(F),
                                 norm = as.logical(F), ann_dist = as.character('cosine'),
                                 n_jobs = as.integer(10), kernel_use = as.character('decay_adaptive')
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
dat <- FindClusters(dat, resolution = 1, graph.name = 'dbgraph', algorithm = 2)

dat <- RunUMAP(dat, reduction = 'sc', dims = 1:ncol(dat@reductions$sc@cell.embeddings), n.neighbors = 15, init = 'spectral',
               min.dist = 0.85, spread = 1.5, learning.rate = 2, n.epochs = 3000, reduction.key = 'dbMAP_', reduction.name = 'dbmap')

dat <- RunUMAP(dat, reduction = 'sc', n.components = 3, dims = 1:ncol(dat@reductions$sc@cell.embeddings),
               min.dist = 0.85, spread = 1.5, learning.rate = 2, n.epochs = 3000, reduction.key = 'dbMAP_3D_', reduction.name = 'dbmap3d')

DimPlot(dat, reduction = 'dbmap', group.by = 'seurat_clusters', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'Sample', pt.size = 1)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
dat <- CellCycleScoring(dat, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)


plot.data <- FetchData(object = dat, vars = c("dbMAP3D_1", "dbMAP3D_2", "dbMAP3D_3", 
                                              "pca_clusters", 'gene_clusters', 'seurat_clusters', 'Sample', 'pred.kal'
))
plot.data$label <- dat$pred.kal
plot_ly(data = plot.data, 
        x = ~dbMAP3D_1, y = ~dbMAP3D_2, z = ~dbMAP3D_3, 
        color = ~pred.kal, 
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 3, width=0.5),
        text=~pred.kal,
        hoverinfo="text", plot_bgcolor = 'black') 


####################################################################
# Transfer labels from Kalluri data
####################################################################
kal <- readRDS("~/Documents/Bioinfo/Ebru/Kalluri/Kalluri.Rds")
transfer.anchors <- FindTransferAnchors(reference = kal, query = dat, reduction = 'cca',
                                        dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = kal$Sub.Cluster, weight.reduction = 'cca',
                            dims = 1:30)
dat <- AddMetaData(dat, metadata = predictions)
dat$pred.kal <- dat$predicted.id

coc <- readRDS("~/Documents/Bioinfo/Ebru/Cochain/Cochain.Rds")
transfer.anchors <- FindTransferAnchors(reference = coc, query = dat, reduction = 'cca',
                                        dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = coc$db_clusters, weight.reduction = 'cca',
                            dims = 1:30)
dat <- AddMetaData(dat, metadata = predictions)
dat$pred.coc <- dat$predicted.id

DimPlot(dat, reduction = 'dbmap', group.by = 'pred.kal',, label = T, repel = T, label.size = 5, pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'pred.coc',, label = T, repel = T, label.size = 5, pt.size = 1)

FeaturePlot(dat, reduction = 'dbmap', features = 'Cd209a', pt.size = 0.5, min.cutoff = 0, max.cutoff = 1, order = T) 

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
dat <- CellCycleScoring(dat, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

# Save RDS
saveRDS(dat, 'Deng.Rds')

# Convert to h5Seurat, h5ad
library(SeuratDisk)
as.h5Seurat(dat, 'Deng.h5Seurat', overwrite = T)
Convert('Deng.h5Seurat', dest = 'h5ad', overwrite = T, assay = 'RNA')


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
dat <- getMarkerGenes(dat, organism = 'mm', groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'pred.kal', 'pred.coc', 'Sample', 'Phase'), min_pct = 0.3, assay = 'integrated', test = 'wilcox')


dat <- getEnrichedPathways(dat, 
                           databases = c('Transcription_Factor_PPIs', 'Chromossome_Location', 'Mouse_Gene_Atlas', 
                                         'MGI_Mammalian_Phenotype_Level_3', 'MGI_Mammalian_Phenotyoe_Level_4_2019', 'MGI_Mammalian_Phenotype_2017',
                                         'OMIM_Disease', 'OMIM_Expanded', 'MsigDB_Computational', 'UK_Biobank_GWAS_v1',
                                         'KEGG_2019_Mouse', 'PheWeb_2019', 'WikiPathways_2019_Mouse', 'GWAS_Catalog_2019',
                                         'GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018',
                                         'Panther_2016', 'BioCarta_2016', 'Reactome_2016', 'Kegg_2016'))
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'pred.kal', 'pred.coc', 'Sample', 'Phase'),
                                        GMT_file =  'h.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'pred.kal', 'pred.coc', 'Sample', 'Phase'),
                                        GMT_file =  'c2.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'pred.kal', 'pred.coc', 'Sample', 'Phase'),
                                        GMT_file =  'c3.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'pred.kal', 'pred.coc', 'Sample', 'Phase'),
                                        GMT_file =  'c5.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'pred.kal', 'pred.coc', 'Sample', 'Phase'),
                                        GMT_file =  'c3.tft.v7.1.symbols.gmt')

cerebro <- exportFromSeurat(dat, file = 'Deng.crb', experiment_name = 'Deng et al', organism = 'mm',
                            cell_cycle = 'Phase', groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'pred.kal', 'pred.coc', 'Sample', 'Phase'),
                            nUMI = 'nCount_RNA', nGene = 'nFeature_RNA')
launchCerebro()

