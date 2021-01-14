######################################################################
#Analysis of Pan et al human data (will use as human reference)
######################################################################

#Load packages 
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(reticulate)
library(cerebroApp)
setwd("~/Documents/Bioinfo/Ebru/Pan/")

#Load data
# SMC-derived (Myh11-CreERT2+ = ZsGreen+ = cells that came from SMCs)
human1 <- as.matrix(read.delim('data/Human/GSM4705589_RPE004_matrix.txt.gz', sep = '\t', header = T, row.names = 1))
human1 <- CreateSeuratObject(human1, project = 'Subject 1')

human2 <- as.matrix(read.delim('data/Human/GSM4705590_RPE005_matrix.txt.gz', sep = '\t', header = T, row.names = 1))
human2 <- CreateSeuratObject(human2, project = 'Subject 2')

human3 <- as.matrix(read.delim('data/Human/GSM4705591_RPE006_matrix.txt.gz', sep = '\t', header = T, row.names = 1))
human3 <- CreateSeuratObject(human3, project = 'Subject 3')


######################################################################
#QC and filtering 
######################################################################
dats <- list(human1, human2, human3)

seurat <- list()

# This is the filtered public-released data. Let's check how filtered it is.
for(i in 1:length(dats)){
  seurat <- dats[[i]]
      seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-", assay = 'RNA')
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

VlnPlot(dats[[1]], features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'))
VlnPlot(dats[[2]], features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'))
VlnPlot(dats[[3]], features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt'))

# Actually looks quite good! Great job, Dr. Pan and colleagues !!!
# We won't filter this any further
# Merge data sets

dat <- merge(dats[[1]], dats[-1])

VlnPlot(object = dat, features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), group.by = 'orig.ident')
dat$Subject <- dat$orig.ident


######################################################################
#Default workflow, integrate with CCA anchoring and cluster with graph Louvain
######################################################################
genes <- rownames(dat)
sets <- SplitObject(dat, split.by = 'Subject')

for(i in 1:length(sets)){
  sets[[i]] <- NormalizeData(sets[[i]])
  sets[[i]] <- FindVariableFeatures(sets[[i]], selection.method = 'disp', nfeatures = 3000, num.bin = 100)
  sets[[i]] <- ScaleData(sets[[i]], features = rownames(sets[[i]]))
}

feat <- SelectIntegrationFeatures(object.list = sets, nfeatures = 3000)
anchors <- FindIntegrationAnchors(object.list = sets, anchor.features = feat, reference = 1, 
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

dat <- FindNeighbors(object = dat, dims = 1:30)
dat <- FindClusters(object = dat, algorithm = 2, resolution = 1.5)
dat$pca_clusters <- dat$seurat_clusters


dat <- RunUMAP(object = dat, dims = 1:30)
UMAPPlot(dat, group.by = 'Subject') 

UMAPPlot(dat, group.by = 'gene_clusters') 
UMAPPlot(dat, group.by = 'pca_clusters', pt.size = 0.5) 

FeaturePlot(dat, features = 'MKI67', order = T, min.cutoff = 0, max.cutoff = 1, pt.size = 1)
FeaturePlot(dat, features = 'Ikzf2', order = T)
FeaturePlot(dat, features = 'Ccr2', order = T, min.cutoff = 2.5)
FeaturePlot(dat, features = 'Cd209a', order = T, min.cutoff = 0, max.cutoff = 1)
FeaturePlot(dat, features = 'FABP5', order = T, min.cutoff = 0, max.cutoff = 4)

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
diff <- dbmap$diffusion$Diffusor(n_components = as.integer(80), n_neighbors = as.integer(15),
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
# Post-diffusion processing (structure clustering, MAP)
##################################################################################

dat <- FindNeighbors(dat, reduction = 'sc', dims = 1:ncol(dat@reductions$sc@cell.embeddings), annoy.metric = 'cosine', graph.name = 'dbgraph')
dat <- FindClusters(dat, resolution = 1, graph.name = 'dbgraph', algorithm = 2)

dat <- RunUMAP(dat, reduction = 'sc', dims = 1:ncol(dat@reductions$sc@cell.embeddings), n.neighbors = 15, init = 'spectral',
               min.dist = 0.8, spread = 1.5, learning.rate = 1.5, n.epochs = 3000, reduction.key = 'dbMAP_', reduction.name = 'dbmap')

dat <- RunUMAP(dat, reduction = 'sc', n.components = 3, dims = 1:ncol(dat@reductions$sc@cell.embeddings),
               min.dist = 0.8, spread = 0.8, learning.rate = 1.5, n.epochs = 800, reduction.key = 'dbMAP_3D_', reduction.name = 'dbmap3d')

DimPlot(dat, reduction = 'dbmap', group.by = 'seurat_clusters', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'gene_clusters', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'Subject', pt.size = 1)

FeaturePlot(dat, reduction = 'dbmap', features = 'MKI67', pt.size = 0.5, min.cutoff = 0, max.cutoff = 1, order = T) 

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
dat <- CellCycleScoring(dat, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)


plot.data <- FetchData(object = dat, vars = c("dbMAP3D_1", "dbMAP3D_2", "dbMAP3D_3", 
                                              "pca_clusters",  'gene_clusters', 'seurat_clusters', 'Subject', 'pred.hcl'
))
plot.data$label <- dat$pred.hcl
plot_ly(data = plot.data, 
        x = ~dbMAP3D_1, y = ~dbMAP3D_2, z = ~dbMAP3D_3, 
        color = ~pred.hcl, 
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 3, width=0.5),
        text=~pred.hcl,
        hoverinfo="text", plot_bgcolor = 'black') 



# Convert to h5Seurat, h5ad
library(SeuratDisk)
as.h5Seurat(dat, 'Pan_Human.H5Seurat', overwrite = T)
Convert('Pan_Human.H5Seurat', dest = 'h5ad', overwrite = T, assay = 'RNA')


######################################################################
#Plot and explore for annotation
########################################################### ###########
####################################################################
# Annotation
####################################################################
hcl_aorta <- readRDS("~/Documents/Bioinfo/Ebru/HCL/HCL_Aorta.Rds")
transfer.anchors <- FindTransferAnchors(reference = hcl_aorta, query = dat, reduction = 'cca',
                                        dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = hcl_aorta$Celltype, weight.reduction = 'cca',
                            dims = 1:30)
dat <- AddMetaData(dat, metadata = predictions)
dat$pred.hcl <- dat$predicted.id
DimPlot(dat, reduction = 'dbmap', group.by = 'pred.hcl', label = T, repel = T, label.size = 5, pt.size = 0.5)

DimPlot(dat, reduction = 'dbmap', group.by = 'Group', pt.size = 0.5)
DimPlot(dat, reduction = 'dbmap', group.by = 'pca_clusters', pt.size = 0.5)
DimPlot(dat, reduction = 'dbmap', group.by = 'seurat_clusters', pt.size = 0.5)

# General markers
FeaturePlot(dat, reduction = 'dbmap', features = 'CD209', pt.size = 0.5, min.cutoff = 0, max.cutoff = 1, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'TREM2', pt.size = 0.5, min.cutoff = 0, max.cutoff = 1, order = T) 


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
dat <- addPercentMtRibo(dat, organism = 'hg', gene_nomenclature = 'name', assay = 'RNA')
dat <- getMarkerGenes(dat, organism = 'hg', groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'Subject', 'pred.hcl'), min_pct = 0.3, assay = 'integrated', test = 'wilcox')

dat <- getEnrichedPathways(dat, 
                           databases = c('Transcription_Factor_PPIs', 'Chromossome_Location', 'Mouse_Gene_Atlas', 
                                         'MGI_Mammalian_Phenotype_Level_3', 'MGI_Mammalian_Phenotyoe_Level_4_2019', 'MGI_Mammalian_Phenotype_2017',
                                         'OMIM_Disease', 'OMIM_Expanded', 'MsigDB_Computational', 'UK_Biobank_GWAS_v1',
                                         'KEGG_2019_Mouse', 'PheWeb_2019', 'WikiPathways_2019_Mouse', 'GWAS_Catalog_2019',
                                         'GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018',
                                         'Panther_2016', 'BioCarta_2016', 'Reactome_2016', 'Kegg_2016'))
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'RNA', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'Subject', 'pred.hcl'),
                                        GMT_file =  'h.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'RNA', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'Subject', 'pred.hcl'),
                                        GMT_file =  'c2.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'RNA', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'Subject', 'pred.hcl'),
                                        GMT_file =  'c3.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'RNA', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'Subject', 'pred.hcl'),
                                        GMT_file =  'c5.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'RNA', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'Subject', 'pred.hcl'),
                                        GMT_file =  'c3.tft.v7.1.symbols.gmt')

cerebro <- exportFromSeurat(dat, file = 'Pan_Human.crb', experiment_name = 'Pan et al Human Carothid atherosclerosis', organism = 'hg', assay = 'integrated', add_all_meta_data = T,
                            cell_cycle = 'Phase',  groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'Subject', 'pred.hcl'),
                            nUMI = 'nCount_RNA', nGene = 'nFeature_RNA')
launchCerebro()

