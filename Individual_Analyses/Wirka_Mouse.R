######################################################################
#Analysis of Wirka et al mouse data 
######################################################################

#Load packages 
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(reticulate)
library(cerebroApp)
setwd("~/Documents/Bioinfo/Vascular/Wirka/")

#Load data
data <- as.matrix(read.delim('data/GSE131776_mouse_scRNAseq_wirka_et_al_GEO.txt', sep = '\t', header = T, row.names = 1))
dat <- CreateSeuratObject(data, project = 'Wirka et al. - Mouse scRNAseq Data')

cit_rna <- as.matrix(read.delim('data/GSE131777_mouse_16wk_citeseq_RNA_wirka_et_al_GEO.txt', sep = '\t', header = T, row.names = 1))
cit <- CreateSeuratObject(cit_rna, project = 'Wirka et al. - Mouse CITEseq Data')

cit_adt <- as.matrix(read.delim('data/GSE131777_mouse_16wk_citeseq_RNA_wirka_et_al_GEO.txt', sep = '\t', header = T, row.names = 1))
cit[['ADT']] <- CreateAssayObject(cit_rna)

######################################################################
#QC and filtering 
######################################################################
# scRNAseq
dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^mt-", assay = 'RNA')

counts_per_cell <- Matrix::colSums(dat)
counts_per_gene <- Matrix::rowSums(dat)
genes_per_cell <- Matrix::colSums(dat@assays$RNA@counts >0)

hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')

plot(counts_per_cell, genes_per_cell, log='xy', col='wheat', title('counts vs genes per cell'))
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
plot(sort(counts_per_cell), xlab='cell', log='y', main='counts per cell (ordered)')

VlnPlot(object = dat, features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), pt.size = 0.1)

# Prefiltered data. Just remove the tips
dat <- subset(dat, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 &
                nCount_RNA > 500 & nCount_RNA < 15000)


# CITEseq
cit[["percent.mt"]] <- PercentageFeatureSet(cit, pattern = "^mt-", assay = 'RNA')

counts_per_cell <- Matrix::colSums(cit)
counts_per_gene <- Matrix::rowSums(cit)
genes_per_cell <- Matrix::colSums(cit@assays$RNA@counts >0)

hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')

plot(counts_per_cell, genes_per_cell, log='xy', col='wheat', title('counts vs genes per cell'))
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
plot(sort(counts_per_cell), xlab='cell', log='y', main='counts per cell (ordered)')

VlnPlot(object = cit, features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), pt.size = 0.1)
VlnPlot(object = cit, features = c("nCount_ADT", 'nFeature_ADT'), pt.size = 0.1)

# Prefiltered data. Looks quite good and won't need to filter further.

#Add some metadata
l_scrnaseq <- rep('10X scRNAseq', times = ncol(dat))
l_cite <- rep('CITE-seq', times = ncol(cit))
label <- c(l_scrnaseq, l_cite)

dat <- merge(dat, cit)
dat <- AddMetaData(dat, label, col.name = 'Sample')

VlnPlot(object = dat, features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), 
        pt.size = 0.01, group.by = 'Sample')

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
                                  normalization.method = 'LogNormalize', dims = 1:50,  reference =2, # Use CITEseq as ref
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
FeaturePlot(dat, features = 'Mki67', order = T, min.cutoff = 0, max.cutoff = 1, pt.size = 1)
FeaturePlot(dat, features = 'Trem2', order = T, min.cutoff = 0, max.cutoff = 3, pt.size = 1)
FeaturePlot(dat, features = 'Cd209a', order = T, min.cutoff = 0, max.cutoff = 3, pt.size = 1)
FeaturePlot(dat, features = 'Myh11', order = T, min.cutoff = 0, max.cutoff = 3, pt.size = 1)
FeaturePlot(dat, features = 'Cd68', order = T, min.cutoff = 0, max.cutoff = 2, pt.size = 1)

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
diff <- dbmap$diffusion$Diffusor(n_components = as.integer(150), n_neighbors = as.integer(15),
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
               min.dist = 0.8, spread = 1.5, learning.rate = 1.2, n.epochs = 5000, reduction.key = 'dbMAP_', reduction.name = 'dbmap')

dat <- RunUMAP(dat, reduction = 'sc', n.components = 3, dims = 1:ncol(dat@reductions$sc@cell.embeddings),
               min.dist = 0.8, spread = 0.8, learning.rate = 1.5, n.epochs = 800, reduction.key = 'dbMAP_3D_', reduction.name = 'dbmap3d')

DimPlot(dat, reduction = 'dbmap', group.by = 'seurat_clusters', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'gene_clusters', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'Sample', pt.size = 0.5) 

FeaturePlot(dat, reduction = 'dbmap',features = 'Mki67', order = T, min.cutoff = 0, max.cutoff = 1, pt.size = 1)
FeaturePlot(dat, reduction = 'dbmap',features = 'Trem2', order = T, min.cutoff = 0, max.cutoff = 3, pt.size = 1)
FeaturePlot(dat, reduction = 'dbmap',features = 'Cd209a', order = T, min.cutoff = 0, max.cutoff = 3, pt.size = 1)
FeaturePlot(dat, reduction = 'dbmap',features = 'Myh11', order = T, min.cutoff = 0, max.cutoff = 3, pt.size = 1)
FeaturePlot(dat, reduction = 'dbmap',features = 'Cd68', order = T, min.cutoff = 0, max.cutoff = 2, pt.size = 1)

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
as.h5Seurat(dat, 'Wirka_Mouse.h5Seurat', overwrite = T)
Convert('Wirka_Mouse.h5Seurat', dest = 'h5ad', overwrite = T, assay = 'integrated')


######################################################################
#Plot and explore for annotation
########################################################### ###########
####################################################################
# Annotation
####################################################################
kal <- readRDS("~/Documents/Bioinfo/Ebru/Kalluri/Kalluri.Rds")
transfer.anchors <- FindTransferAnchors(reference = kal, query = dat, reduction = 'cca', reference.assay = 'RNA',
                                        query.assay = 'integrated',
                                        dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = kal$Sub.Cluster, weight.reduction = 'cca',
                            dims = 1:30)
dat <- AddMetaData(dat, metadata = predictions)
dat$pred.kal <- dat$predicted.id
DimPlot(dat, reduction = 'dbmap', group.by = 'pred.kal', label = T, repel = T, label.size = 5, pt.size = 0.5)

DimPlot(dat, reduction = 'dbmap', group.by = 'Group', pt.size = 0.5)
DimPlot(dat, reduction = 'dbmap', group.by = 'pca_clusters', pt.size = 0.5)
DimPlot(dat, reduction = 'dbmap', group.by = 'seurat_clusters', pt.size = 0.5)

# General markers
FeaturePlot(dat, reduction = 'dbmap', features = 'CD209', pt.size = 0.5, min.cutoff = 0, max.cutoff = 1, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'TREM2', pt.size = 0.5, min.cutoff = 0, max.cutoff = 1, order = T) 
FeaturePlot(dat, reduction = 'dbmap', features = 'MYH11', pt.size = 0.5, min.cutoff = 0, max.cutoff = 1, order = T) 


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
dat <- getMarkerGenes(dat, organism = 'mm', groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'pred.kal', 'Sample'), min_pct = 0.3, assay = 'RNA', test = 'wilcox')

dat <- getEnrichedPathways(dat, 
                           databases = c('Transcription_Factor_PPIs', 'Chromossome_Location', 'Mouse_Gene_Atlas', 
                                         'MGI_Mammalian_Phenotype_Level_3', 'MGI_Mammalian_Phenotyoe_Level_4_2019', 'MGI_Mammalian_Phenotype_2017',
                                         'OMIM_Disease', 'OMIM_Expanded', 'MsigDB_Computational', 'UK_Biobank_GWAS_v1',
                                         'KEGG_2019_Mouse', 'PheWeb_2019', 'WikiPathways_2019_Mouse', 'GWAS_Catalog_2019',
                                         'GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018',
                                         'Panther_2016', 'BioCarta_2016', 'Reactome_2016', 'Kegg_2016'))
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'RNA', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'pred.kal', 'Sample'),
                                        GMT_file =  'h.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'RNA', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'pred.kal', 'Sample'),
                                        GMT_file =  'c2.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'RNA', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'pred.kal', 'Sample'),
                                        GMT_file =  'c3.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'RNA', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'pred.kal', 'Sample'),
                                        GMT_file =  'c5.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'RNA', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'pred.kal', 'Sample'),
                                        GMT_file =  'c3.tft.v7.1.symbols.gmt')

cerebro <- exportFromSeurat(dat, file = 'Wirka_Mouse.crb', experiment_name = 'Wirka et al. mouse atherosclerosis', organism = 'mm', assay = 'RNA', add_all_meta_data = T,
                            cell_cycle = 'Phase',  groups = c('Sample','pca_clusters', 'gene_clusters', 'seurat_clusters', 'Sample'),
                            nUMI = 'nCount_RNA', nGene = 'nFeature_RNA')
launchCerebro()

