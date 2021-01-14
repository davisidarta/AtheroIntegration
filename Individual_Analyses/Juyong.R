######################################################################
#Analysis of Kim, Juyong Brian et al data
######################################################################
#Load packages 
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(reticulate)
library(plotly)
library(cerebroApp)
setwd("~/Documents/Bioinfo/Vascular/Juyong/")

#Load data
counts <- as.matrix(read.delim('data/GSE150768_scRNAseq_rawmatrix.txt.gz', sep=' ', header = T, row.names = 1))
meta <- read.delim('data/GSE150768_scRNAseq_metadata.txt.gz', sep = ' ')
dat <- CreateSeuratObject(counts, meta.data = meta )

######################################################################
#QC and filtering 
######################################################################
dat$Sample <- dat$orig.ident
dats <- SplitObject(dat, split.by = 'Sample')

for(i in 1:length(dats)){
    seurat <- dats[[i]]
    
    seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-", assay = 'RNA')
    counts_per_cell <- Matrix::colSums(seurat)
    counts_per_gene <- Matrix::rowSums(seurat)
    genes_per_cell <- Matrix::colSums(seurat@assays$RNA@counts >0)
    
    hist(log10(counts_per_cell+1),main='Counts per cell',col='chartreuse4')
    hist(log10(genes_per_cell+1), main='Genes per cell', col='chartreuse4')
    
    plot(counts_per_cell, genes_per_cell, type = 'p', log='xy', pch=19, col="firebrick3", xlab="Counts per cell (log)", ylab="Genes per cell (log)", title(paste('Counts vs Genes per Cell', as.character((names(dats)[[i]]), sep = ' - '))))
    plot(sort(genes_per_cell), pch=19, col="firebrick3", xlab="Cell rank", ylab="Genes per cell (log)", 
         log='y', main= paste('Gene detection rank', as.character((names(dats)[[i]]), sep = ' - ')))
    plot(sort(counts_per_cell), pch=19, col="firebrick3", xlab="Cell rank", ylab="Genes per cell (log)",
         log='y', main= paste('Count detection rank', as.character((names(dats)[[i]]), sep = ' - ')))
    
    dats[[i]] <- seurat
  }  
}

# Let's try to preserve as much cells as possible during filtering first
# tdT-WT 1 #### Doubt in filtering
VlnPlot(object = dats[[1]], features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'))
dats[[1]] <- subset(dats[[1]], subset = 
                      nFeature_RNA > 1500 & nFeature_RNA < 3500 &
                      nCount_RNA > 6000 & nCount_RNA < 8500 &
                      percent.mt < 10)

# tdT-WT 2 #### Doubt in filtering
VlnPlot(object = dats[[2]], features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'))
dats[[2]] <- subset(dats[[2]], subset = 
                      nFeature_RNA > 1500 & nFeature_RNA < 3500 &
                      nCount_RNA > 6000 & nCount_RNA < 8250 &
                      percent.mt < 10)

# tdTAKO1 [TDT-A120RL_S4_L008] #### Doubt in filtering
VlnPlot(object = dats[[3]], features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'))
dats[[3]] <- subset(dats[[3]], subset = 
                      nFeature_RNA > 1500 & nFeature_RNA < 3500 &
                      nCount_RNA > 6250 & nCount_RNA < 8500 &
                      percent.mt < 10)

# tdTAKO2 [A305N-TDT_S6_L003] #### Doubt in filtering
VlnPlot(object = dats[[4]], features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'))
dats[[4]] <- subset(dats[[4]], subset = 
                      nFeature_RNA > 1500 & nFeature_RNA < 3500 &
                      nCount_RNA > 6250 & nCount_RNA < 8500 &
                      percent.mt < 10)

# tdTAKO3 [A305R-TDT_S5_L003] #### Doubt in filtering
VlnPlot(object = dats[[5]], features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'))
dats[[5]] <- subset(dats[[5]], subset = 
                      nFeature_RNA > 1500 & nFeature_RNA < 3500 &
                      nCount_RNA > 6250 & nCount_RNA < 8500 &
                      percent.mt < 10)

dat <- merge(dats[[1]], dats[-1])

######################################################################
# Let's compare some workflows
######################################################################
# Default workflow, no batch-correction
dat <- NormalizeData(dat)
dat <- FindVariableFeatures(dat, selection.method = 'disp', nfeatures = 5000)
dat <- ScaleData(dat, features = rownames(dat))
#
dat <- RunPCA(dat, npcs = 100)
ElbowPlot(dat, ndims = 100)
#
dat <- FindNeighbors(object = dat, features = VariableFeatures(dat))
dat <- FindClusters(object = dat, features = VariableFeatures(dat), algorithm = 2)
dat$gene_clusters <- dat$seurat_clusters
#
dat <- FindNeighbors(object = dat, dims = 1:30)
dat <- FindClusters(object = dat, algorithm = 2)
dat$pca_clusters <- dat$seurat_clusters
#
dat <- RunUMAP(object = dat, dims = 1:30)
UMAPPlot(dat, group.by = 'Sample') 
UMAPPlot(dat, group.by = 'cond') 
UMAPPlot(dat, group.by = 'celltype', label = T, repel = T, label.size = 7) 
FeaturePlot(dat, features = c('Cnn1', 'Lum', 'Cyp1b1', 'Ahr'), order = T, min.cutoff = 0)

# Default workflow, no batch-correction
###############################################################################
# Run dbMAP on Default workflow, no batch-correction
###############################################################################
library(reticulate)
np <- reticulate::import("numpy")
pd <- reticulate::import("pandas")
sp <- reticulate::import("scipy")
dbmap <- reticulate::import('dbmap')

data <- t(dat@assays$RNA@data[VariableFeatures(dat),])
a <- r_to_py(data)
b <- a$tocsr()
diff <- dbmap$diffusion$Diffusor(n_components = as.integer(80), n_neighbors = as.integer(30),
                                 transitions = as.logical(F),
                                 norm = as.logical(F), ann_dist = as.character('lp'), p = as.numeric(0.5),
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

dat <- FindNeighbors(dat, reduction = 'sc', dims = 1:ncol(dat@reductions$sc@cell.embeddings), annoy.metric = 'cosine', graph.name = 'dbgraph')
dat <- FindClusters(dat, resolution = 1, graph.name = 'dbgraph', algorithm = 2)

dat <- RunUMAP(dat, reduction = 'sc', dims = 1:ncol(dat@reductions$sc@cell.embeddings), n.neighbors = 30, init = 'spectral',
               min.dist = 0.8, spread = 1.2, learning.rate = 1.5, n.epochs = 300, reduction.key = 'dbMAP_', reduction.name = 'dbmap')

dat <- RunUMAP(dat, reduction = 'sc', n.components = 3, dims = 1:ncol(dat@reductions$sc@cell.embeddings),
               min.dist = 0.85, spread = 1.2, n.epochs = 1000, reduction.key = 'dbMAP_3D_', reduction.name = 'dbmap3d')

DimPlot(dat, reduction = 'dbmap', group.by = 'gene_clusters', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'pca_clusters', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'seurat_clusters', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'Sample', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'cond', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'celltype', pt.size = 1, label = T, repel = T, label.size = 7)
FeaturePlot(dat, reduction = 'dbmap', features = c('Cnn1', 'Lum', 'Cyp1b1', 'Ahr'), order = T, min.cutoff = 0)
###############################################################################
#
# CCA anchoring
genes <- rownames(dat)
sets <- SplitObject(dat, split.by = 'Sample')

for(i in 1:length(sets)){
  sets[[i]] <- NormalizeData(sets[[i]])
  sets[[i]] <- FindVariableFeatures(sets[[i]], selection.method = 'disp', nfeatures = 2000)
  sets[[i]] <- ScaleData(sets[[i]], features = rownames(sets[[i]]))
}

feat <- SelectIntegrationFeatures(object.list = sets, nfeatures = 2000)
anchors <- FindIntegrationAnchors(object.list = sets, anchor.features = feat, nn.method = 'annoy')
dat <- IntegrateData(anchorset = anchors, features.to.integrate = genes)
dat <- ScaleData(dat)

dat <- RunPCA(dat, npcs = 100)
ElbowPlot(dat, ndims = 100)

dat <- FindNeighbors(object = dat, features = VariableFeatures(dat))
dat <- FindClusters(object = dat, features = VariableFeatures(dat), algorithm = 2,
                    resolution = 1.5)
dat$gene_clusters <- dat$seurat_clusters

dat <- FindNeighbors(object = dat, dims = 1:30)
dat <- FindClusters(object = dat, algorithm = 2, resolution = 1.5)
dat$pca_clusters <- dat$seurat_clusters

dat <- RunUMAP(object = dat, dims = 1:30)
UMAPPlot(dat, group.by = 'Sample', pt.size = 0.5) 
UMAPPlot(dat, group.by = 'pca_clusters', pt.size = 0.5) 
UMAPPlot(dat, group.by = 'gene_clusters', pt.size = 0.5) 
UMAPPlot(dat, group.by = 'celltype', label = T, repel = T, label.size = 7) 
FeaturePlot(dat, features = c('Cnn1', 'Lum', 'Cyp1b1', 'Ahr'), order = T, min.cutoff = 0)

###############################################################################
# Run dbMAP on the integrated data
###############################################################################
library(reticulate)
np <- reticulate::import("numpy")
pd <- reticulate::import("pandas")
sp <- reticulate::import("scipy")
dbmap <- reticulate::import('dbmap')

data <- t(dat@assays$integrated@data[VariableFeatures(dat),])
a <- r_to_py(data)
b <- a$tocsr()
diff <- dbmap$diffusion$Diffusor(n_components = as.integer(60), n_neighbors = as.integer(20),
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
dat <- FindClusters(dat, resolution = 0.4, graph.name = 'dbgraph', algorithm = 2)

dat <- RunUMAP(dat, reduction = 'sc', dims = 1:ncol(dat@reductions$sc@cell.embeddings), n.neighbors = 20, init = 'spectral',
               min.dist = 0.8, spread = 1.5, learning.rate = 1.5, n.epochs = 2000, reduction.key = 'dbMAP_', reduction.name = 'dbmap')

dat <- RunUMAP(dat, reduction = 'sc', n.components = 3, dims = 1:ncol(dat@reductions$sc@cell.embeddings),
               min.dist = 0.8, spread = 1.2, n.epochs = 6000, reduction.key = 'dbMAP_3D_', reduction.name = 'dbmap3d')

DimPlot(dat, reduction = 'dbmap', group.by = 'seurat_clusters', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'gene_clusters', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'Group', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'celltype', label = T, repel = T, label.size = 7, pt.size = 0.5) 

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
dat <- CellCycleScoring(dat, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

FeaturePlot(dat, reduction = 'dbmap', features = 'Mki67', pt.size = 0.5, min.cutoff = 0, max.cutoff = 1, order = T) 


# Convert to h5Seurat, h5ad
library(SeuratDisk)
as.h5Seurat(dat, 'Juyong.h5Seurat', overwrite = T)
Convert('Juyong.h5Seurat', dest = 'h5ad', overwrite = T, assay = 'RNA')


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
                                              "pca_clusters",  'seurat_clusters', 'celltype', 'cond', 'Sample'
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
dat <- getMarkerGenes(dat, organism = 'mm', groups = c('seurat_clusters', 'celltype', 'cond', 'Sample'), min_pct = 0.3, assay = 'integrated', test = 'wilcox')
dat <- getEnrichedPathways(dat, 
                           databases = c('Transcription_Factor_PPIs', 'Chromossome_Location', 'Mouse_Gene_Atlas', 
                                         'MGI_Mammalian_Phenotype_Level_3', 'MGI_Mammalian_Phenotyoe_Level_4_2019', 'MGI_Mammalian_Phenotype_2017',
                                         'OMIM_Disease', 'OMIM_Expanded', 'MsigDB_Computational', 'UK_Biobank_GWAS_v1',
                                         'KEGG_2019_Mouse', 'PheWeb_2019', 'WikiPathways_2019_Mouse', 'GWAS_Catalog_2019',
                                         'GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018',
                                         'Panther_2016', 'BioCarta_2016', 'Reactome_2016', 'Kegg_2016'))
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('seurat_clusters', 'celltype', 'cond', 'Sample'),
                                        GMT_file =  'h.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('seurat_clusters', 'celltype', 'cond', 'Sample'),
                                        GMT_file =  'c2.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('seurat_clusters', 'celltype', 'cond', 'Sample'),
                                        GMT_file =  'c3.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('seurat_clusters', 'celltype', 'cond', 'Sample'),
                                        GMT_file =  'c5.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'integrated', 
                                        groups = c('seurat_clusters', 'celltype', 'cond', 'Sample'),
                                        GMT_file =  'c3.tft.v7.1.symbols.gmt')

cerebro <- exportFromSeurat(dat, file = 'Juyong.crb', experiment_name = 'Juyong et al', organism = 'mm',
                            cell_cycle = 'Phase', groups = c('seurat_clusters', 'celltype', 'cond', 'Sample'),
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

######################################################################
# Find markers
######################################################################
# Seurat cluster markers
library(dplyr)
Idents(dat) <- 'pca_clusters'
seurat.markers <- FindAllMarkers(dat, only.pos = F, min.pct = 0.25,
                                 assay = 'RNA', logfc.threshold = 0.25
                                  )
seurat.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top10 <- seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(dat, features = top10$gene) + NoLegend()

seurat.markers$ID <- rownames(seurat.markers) # For GAM and GSEA shiny applications
seurat.markers$log2FC <- seurat.markers$avg_logFC
write.csv(seurat.markers, file = 'KO.csv', row.names = F)
seurat.markers <- seurat.markers[-2]
write.table(seurat.markers, file = 'KO.tsv', sep = '\t', row.names = F)



# KO vs WT
Idents(dat) <- 'cond'
KO <- FindMarkers(dat, ident.1 = 'KO',
                  assay = 'integrated', test.use = 'LR')
KO$ID <- rownames(KO) # For GAM and GSEA shiny applications
KO$log2FC <- KO$avg_logFC
write.csv(KO, file = 'KO.csv', row.names = F)
KO <- KO[-2]
write.table(KO, file = 'KO.tsv', sep = '\t', row.names = F)
KO = KO[order(KO$log2FC),]
DoHeatmap(dat, features = rownames(KO), group.by = 'Sample')

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
