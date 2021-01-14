######################################################################
#Analysis of Alencar et al mouse data
######################################################################

#Load packages 
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(reticulate)
library(cerebroApp)
setwd("~/Documents/Bioinfo/Vascular/Alencar/")

#Load data
# Unsorted
wt_unsorted1 <- Read10X('data/GSM4555583_KLF4_WT_unsorted01_filtered_gene_bc_matrices')
wt_unsorted1 <- CreateSeuratObject(wt_unsorted1, project = 'Unsorted WT 1')
#
wt_unsorted2 <- Read10X('data/GSM4555584_KLF4_WT_unsorted02_filtered_gene_bc_matrices')
wt_unsorted2 <- CreateSeuratObject(wt_unsorted2, project = 'Unsorted WT 2')
#
ko_unsorted1 <- Read10X('data/GSM4555587_KLF4KO_unsorted01_filtered_gene_bc_matrices')
ko_unsorted1 <- CreateSeuratObject(ko_unsorted1, project = 'Unsorted Klf4 KO 1')
#
ko_unsorted2 <- Read10X('data/GSM4555588_KLF4KO_unsorted02_filtered_gene_bc_matrices')
ko_unsorted2 <- CreateSeuratObject(ko_unsorted2, project = 'Unsorted Klf4 KO 2')
#
wt_unsorted3 <- Read10X('data/GSM4555591_SMC_WT_eYFP_NEG_unsorted01_filtered_gene_bc_matrices')
wt_unsorted3 <- CreateSeuratObject(wt_unsorted3, project = 'Unsorted WT 3 (eYFP-)')
#
endo_unsorted <- Read10X('data/GSM4555601_endothelial_lineage_unsorted_control_filtered_gene_bc_matrices')
endo_unsorted <- CreateSeuratObject(endo_unsorted, project = 'Cdh5+ Traced Unsorted')
#
bca_unsorted1 <- Read10X('data/GSM4555598_bca_unsorted_01_filtered_gene_bc_matrices')
bca_unsorted1 <- CreateSeuratObject(bca_unsorted1, project = 'BCA Unsorted 1')
#
bca_unsorted2 <- Read10X('data/GSM4555592_SMC_WT_eYFP_POSITIVE_filtered_gene_bc_matrices')
bca_unsorted2 <- CreateSeuratObject(bca_unsorted2, project = 'BCA Unsorted 2')
#
aorta_unsorted <- Read10X('data/GSM4555600_aorta_03_filtered_gene_bc_matrices')
aorta_unsorted <- CreateSeuratObject(aorta_unsorted, project = 'Aorta Unsorted')
#
endo_unsorted <- Read10X('data/GSM4555601_endothelial_lineage_unsorted_control_filtered_gene_bc_matrices')
endo_unsorted <- CreateSeuratObject(endo_unsorted, project = 'Cdh5+ Traced Unsorted')


# Sorted
wt_smc1 <- Read10X('data/GSM4555585_KLF4WTeYFPpos01_filtered_gene_bc_matrices')
wt_smc1 <- CreateSeuratObject(wt_smc1, project = 'Sorted eYFP+ WT SMCs 1')
#
wt_smc2 <- Read10X('data/GSM4555586_KLF4WTeYFPpos02_filtered_gene_bc_matrices')
wt_smc2 <- CreateSeuratObject(wt_smc2, project = 'Sorted eYFP+ WT SMCs 2')
#
ko_smc1 <- Read10X('data/GSM4555589_KLF4KOeYFPpos01_filtered_gene_bc_matrices')
ko_smc1 <- CreateSeuratObject(ko_smc1, project = 'Sorted eYFP+ Klf4 KO SMCs 1')
#
ko_smc2 <- Read10X('data/GSM4555590_KLF4KOeYFPpos02_filtered_gene_bc_matrices')
ko_smc2 <- CreateSeuratObject(ko_smc2, project = 'Sorted eYFP+ Klf4 KO SMCs 2')
#
wt_smc3 <- Read10X('data/GSM4555592_SMC_WT_eYFP_POSITIVE_filtered_gene_bc_matrices')
wt_smc3 <- CreateSeuratObject(wt_smc3, project = 'Sorted eYFP+ WT SMCs 3')
#
# Dual reporter model
# GFP+ (Myh11-DreERT2 Lgals3-Cre Rosa-tdTomato-eGFP apoE-/- ---> Express eGFP upon Lgals3 expression)
gfp_smc1 <- Read10X('data/GSM4555593_green_01_filtered_gene_bc_matrices')
gfp_smc1 <- CreateSeuratObject(gfp_smc1, project = 'Sorted eGFP+ Transitioned SMCs 1')
#
gfp_smc2 <- Read10X('data/GSM4555594_green_02_filtered_gene_bc_matrices')
gfp_smc2 <- CreateSeuratObject(gfp_smc2, project = 'Sorted eGFP+ Transitioned SMCs 2')
#
gfp_smc3 <- Read10X('data/GSM4555595_green_03_filtered_gene_bc_matrices')
gfp_smc3 <- CreateSeuratObject(gfp_smc3, project = 'Sorted eGFP+ Transitioned SMCs 3')
#
td_smc1 <- Read10X('data/GSM4555596_tdtomato_01_filtered_gene_bc_matrices')
td_smc1 <- CreateSeuratObject(td_smc1, project = 'Sorted tdTomato+ SMCs 1')
#
td_smc2 <- Read10X('data/GSM4555597_tdtomato_03_filtered_gene_bc_matrices')
td_smc2 <- CreateSeuratObject(td_smc2, project = 'Sorted tdTomato+ SMCs 2')
#
# Endothelial lineage
endo_sorted <- Read10X('data/GSM4555601_endothelial_lineage_unsorted_control_filtered_gene_bc_matrices')
endo_sorted <- CreateSeuratObject(endo_sorted, project = 'Cdh5+ Traced Sorted')


### Unsorted vs Sorted QCs
unsorted <- list(bca_unsorted1, bca_unsorted2, wt_unsorted1, wt_unsorted2, wt_unsorted3, endo_unsorted, ko_unsorted1, ko_unsorted2)
sorted <- list(wt_smc1, wt_smc2, wt_smc3, ko_smc1, ko_smc2, td_smc1, td_smc2, gfp_smc1, gfp_smc2, gfp_smc3) # without `endo_sorted` !!!

######################################################################
#QC and filtering 
######################################################################

for(i in 1:length(unsorted)){
  seurat <- unsorted[[i]]
  
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-", assay = 'RNA')
  counts_per_cell <- Matrix::colSums(seurat)
  counts_per_gene <- Matrix::rowSums(seurat)
  genes_per_cell <- Matrix::colSums(seurat@assays$RNA@counts >0)
  
  plot(counts_per_cell, genes_per_cell, type = 'p', log='xy', pch=19, col="firebrick3", xlab="Counts per cell (log)", ylab="Genes per cell (log)", main = paste('Counts vs. Genes per cell (log)', as.character(Project(seurat), sep = ' - ')))
  plot(sort(genes_per_cell), pch=19, col="firebrick3", xlab="Cell rank", ylab="Genes per cell (log)", 
       log='y', main= paste('Gene detection rank', as.character(Project(seurat), sep = ' - ')))
  plot(sort(counts_per_cell), pch=19, col="firebrick3", xlab="Cell rank", ylab="Genes per cell (log)",
       log='y', main= paste('Count detection rank', as.character(Project(seurat), sep = ' - ')))
  
  unsorted[[i]] <- seurat
}  

for(i in 1:length(sorted)){
  seurat <- sorted[[i]]
  
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-", assay = 'RNA')
  counts_per_cell <- Matrix::colSums(seurat)
  counts_per_gene <- Matrix::rowSums(seurat)
  genes_per_cell <- Matrix::colSums(seurat@assays$RNA@counts >0)
  
  plot(counts_per_cell, genes_per_cell, type = 'p', log='xy', pch=19, col="firebrick3", xlab="Counts per cell (log)", ylab="Genes per cell (log)", main = paste('Counts vs. Genes per cell (log)', as.character(Project(seurat), sep = ' - ')))
  plot(sort(genes_per_cell), pch=19, col="firebrick3", xlab="Cell rank", ylab="Genes per cell (log)", 
       log='y', main= paste('Gene detection rank', as.character(Project(seurat), sep = ' - ')))
  plot(sort(counts_per_cell), pch=19, col="firebrick3", xlab="Cell rank", ylab="Genes per cell (log)",
       log='y', main= paste('Count detection rank', as.character(Project(seurat), sep = ' - ')))
  
  sorted[[i]] <- seurat
}
endo_sorted[["percent.mt"]] <- PercentageFeatureSet(endo_sorted, pattern = "^mt-", assay = 'RNA')

VlnPlot(object = unsorted[[1]], features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), group.by = 'orig.ident')
VlnPlot(object = unsorted[[2]], features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), group.by = 'orig.ident')
VlnPlot(object = unsorted[[3]], features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), group.by = 'orig.ident')
VlnPlot(object = unsorted[[4]], features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), group.by = 'orig.ident')
VlnPlot(object = unsorted[[5]], features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), group.by = 'orig.ident')
VlnPlot(object = unsorted[[6]], features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), group.by = 'orig.ident')
VlnPlot(object = unsorted[[7]], features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), group.by = 'orig.ident')
VlnPlot(object = unsorted[[8]], features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), group.by = 'orig.ident')


# Let's first try no filtering at all
dat_unsorted <- merge(unsorted[[1]], unsorted[-1])
dat_sorted <- merge(sorted[[1]], unsorted[-1])
dat_sorted <- merge(dat_sorted, endo_sorted)
dat <- merge(dat_unsorted, dat_sorted)

VlnPlot(object = dat, features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), group.by = 'orig.ident')

######################################################################
#Add some meta-data
######################################################################


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
  sets[[i]] <- FindVariableFeatures(sets[[i]], selection.method = 'disp', nfeatures = 1000, num.bin = 50)
  sets[[i]] <- ScaleData(sets[[i]], features = rownames(sets[[i]]))
}

feat <- SelectIntegrationFeatures(object.list = sets, nfeatures = 1000)
anchors <- FindIntegrationAnchors(object.list = sets, anchor.features = feat, k.anchor = 15, k.filter = 50, k.score = 20,
                                  normalization.method = 'LogNormalize', 
                                  nn.method = 'annoy', reduction = 'cca')
dat <- IntegrateData(anchorset = anchors, features.to.integrate = genes, k.weight = 150)

dat <- FindVariableFeatures(dat, nfeatures = 3000, selection.method = 'disp')
dat <- ScaleData(dat)
dat <- RunPCA(dat, npcs = 100)
ElbowPlot(dat, ndims = 100)

dat <- FindNeighbors(object = dat, features = VariableFeatures(dat))
dat <- FindClusters(object = dat, features = VariableFeatures(dat), algorithm = 2)
dat$gene_clusters <- dat$seurat_clusters

dat <- FindNeighbors(object = dat, dims = 1:50)
dat <- FindClusters(object = dat, algorithm = 2)
dat$pca_clusters <- dat$seurat_clusters

dat <- RunUMAP(object = dat, dims = 1:50)
UMAPPlot(dat, group.by = 'Sample') 

UMAPPlot(dat, group.by = 'Genotype') 
UMAPPlot(dat, group.by = 'Sorting') 

UMAPPlot(dat, group.by = 'pca_clusters', pt.size = 0.5) 
UMAPPlot(dat, group.by = 'gene_clusters', pt.size = 0.5) 

FeaturePlot(dat, features = 'Foxp3', order = T, min.cutoff = 0, max.cutoff = 1, pt.size = 1)
FeaturePlot(dat, features = 'Ikzf2', order = T)
FeaturePlot(dat, features = 'Ccr2', order = T, min.cutoff = 0, max.cutoff = 1)
FeaturePlot(dat, features = 'Cd209a', order = T, min.cutoff = 0, max.cutoff = 1)
FeaturePlot(dat, features = 'Fabp5', order = T, min.cutoff = 0, max.cutoff = 2)
FeaturePlot(dat, features = 'Myh11', order = T, min.cutoff = 0, max.cutoff = 2)
FeaturePlot(dat, features = 'Mki67', order = T, min.cutoff = 0, max.cutoff = 2)

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
as.h5Seurat(dat, 'Alencar.h5Seurat', overwrite = T)
Convert('Alencar.h5Seurat', dest = 'h5ad', overwrite = T, assay = 'RNA')


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

