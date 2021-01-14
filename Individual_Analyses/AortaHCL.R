######################################################################
#Analysis of aortic human tissue scRNASeq data - Human Cell Landscape
######################################################################

#Load packages 
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(reticulate)
library(cerebroApp)
setwd("~/Documents/Bioinfo/Vascular/HCL/")

#Load data
counts <- read.delim('data/AdultArtery1.rawdge.txt', sep = ' ', row.names = 1, header = T)
cells <- gsub(".*_1.","",colnames(counts))
colnames(counts) <- cells

meta <-read.csv('data/Adult-Artery1_rmbatchAnno.csv')
rownames(meta) <- meta$Cell_barcode

dat <- CreateSeuratObject(counts = counts, meta.data = meta, project = 'Human Artery HCL')

######################################################################
#QC and filtering 
######################################################################

dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^MT-", assay = 'RNA')

counts_per_cell <- Matrix::colSums(dat)
counts_per_gene <- Matrix::rowSums(dat)
genes_per_cell <- Matrix::colSums(dat@assays$RNA@counts >0)

hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')

plot(counts_per_cell, genes_per_cell, log='xy', col='wheat', title('counts vs genes per cell'))
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')
plot(sort(counts_per_cell), xlab='cell', log='y', main='counts per cell (ordered)')

VlnPlot(object = dat, features = c("nCount_RNA", 'nFeature_RNA', 'percent.mt'), pt.size = 0.01)

# Filter a little bit on top
dat <- subset(dat, subset = nFeature_RNA < 1500 &
              nCount_RNA < 2000 & percent.mt < 20)

######################################################################
# Default Seurat workflow 
######################################################################

dat <- NormalizeData(dat)
dat <- FindVariableFeatures(dat, selection.method = 'disp', nfeatures = 5000)
dat <- ScaleData(dat, features = rownames(dat))

dat <- RunPCA(dat, npcs = 100)
ElbowPlot(dat, ndims = 100)

dat <- FindNeighbors(object = dat, dims = 1:30)
dat <- FindClusters(object = dat, algorithm = 2, resolution = 1.5)
dat$pca_clusters <- dat$seurat_clusters

dat <- FindNeighbors(object = dat, features = VariableFeatures(dat))
dat <- FindClusters(object = dat, features = VariableFeatures(dat), algorithm = 2,
                    resolution = 1.5)
dat$gene_clusters <- dat$seurat_clusters

dat <- RunUMAP(object = dat, dims = 1:30)
UMAPPlot(dat, group.by = 'pca_clusters', pt.size = 0.5) 
UMAPPlot(dat, group.by = 'Celltype', pt.size = 0.5) 

FeaturePlot(dat, features = 'MKI67', order = T, min.cutoff = 0, max.cutoff = 0.1, pt.size = 1)
FeaturePlot(dat, features = 'CD209', order = T, min.cutoff = 0, max.cutoff = 0.1, pt.size = 1)
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

data <- t(dat@assays$RNA@data[VariableFeatures(dat),])
a <- r_to_py(data)
b <- a$tocsr()
diff <- dbmap$diffusion$Diffusor(n_components = as.integer(100), n_neighbors = as.integer(15),
                                 transitions = as.logical(F),
                                 norm = as.logical(F), ann_dist = as.character('cosine'),
                                 n_jobs = as.integer(8), kernel_use = as.character('decay_adaptive')
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

dat <- RunUMAP(dat, reduction = 'sc', dims = 1:ncol(dat@reductions$sc@cell.embeddings), n.neighbors = 15, init = 'spectral',
               min.dist = 0.65, spread = 1.5, learning.rate = 1.5, n.epochs = 800, reduction.key = 'dbMAP_', reduction.name = 'dbmap')

dat <- RunUMAP(dat, reduction = 'sc', n.components = 3, dims = 1:ncol(dat@reductions$sc@cell.embeddings),
               min.dist = 0.45, spread = 0.8, n.epochs = 800, reduction.key = 'dbMAP_3D', reduction.name = 'dbmap3d')

DimPlot(dat, reduction = 'dbmap', group.by = 'Celltype', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'gene_clusters', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'Group', pt.size = 1)
DimPlot(dat, reduction = 'dbmap', group.by = 'Phase', pt.size = 1)


plot.data <- FetchData(object = dat, vars = c("dbMAP3D_1", "dbMAP3D_2", "dbMAP3D_3", 
                                              "Celltype", "pca_clusters", 'gene_clusters', 'seurat_clusters', 
                                              'Mki67'))
plot.data$label <- dat$Celltype
plot_ly(data = plot.data, 
        x = ~dbMAP3D_1, y = ~dbMAP3D_2, z = ~dbMAP3D_3, 
        color = ~Celltype, 
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 3, width=0.5),
        text= ~Celltype, 
        hoverinfo="text") 
plot <- DimPlot(dat, reduction = "dbmap") + NoLegend()
cycling <- subset(dat, subset = Mki67 > 2)
LabelPoints(plot = plot, points = colnames(cycling), repel = TRUE)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
dat <- CellCycleScoring(dat, s.features = s.genes, g2m.features = g2m.genes, set.ident = F, nbin=10)

# Convert to h5Seurat, h5ad
saveRDS(dat, 'HCL_Aorta.Rds')
library(SeuratDisk)

dat@assays$RNA@data <- dat@assays$RNA@counts
as.h5Seurat(dat, 'HCL_Aorta.h5Seurat', overwrite = T)
Convert('HCL_Aorta.h5Seurat', dest = 'h5ad', overwrite = T, assay = 'RNA')

library(cerebroApp)
dat <- addPercentMtRibo(dat, organism = 'hg', gene_nomenclature = 'name', assay = 'RNA')
dat <- getMarkerGenes(dat, organism = 'hg', groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'Celltype', 'Phase'), min_pct = 0.3, assay = 'RNA', test = 'wilcox')

dat <- getEnrichedPathways(dat, 
                           databases = c('Transcription_Factor_PPIs', 'Chromossome_Location', 'Mouse_Gene_Atlas', 
                                         'MGI_Mammalian_Phenotype_Level_3', 'MGI_Mammalian_Phenotyoe_Level_4_2019', 'MGI_Mammalian_Phenotype_2017',
                                         'OMIM_Disease', 'OMIM_Expanded', 'MsigDB_Computational', 'UK_Biobank_GWAS_v1',
                                         'KEGG_2019_Mouse', 'PheWeb_2019', 'WikiPathways_2019_Mouse', 'GWAS_Catalog_2019',
                                         'GO_Biological_Process_2018', 'GO_Cellular_Component_2018', 'GO_Molecular_Function_2018',
                                         'Panther_2016', 'BioCarta_2016', 'Reactome_2016', 'Kegg_2016'))
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'RNA', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'Celltype', 'Phase'),
                                        GMT_file =  'h.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'RNA', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'Celltype', 'Phase'),
                                        GMT_file =  'c2.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'RNA', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'Celltype', 'Phase'),
                                        GMT_file =  'c3.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'RNA', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'Celltype', 'Phase'),
                                        GMT_file =  'c5.all.v7.1.symbols.gmt')
dat <- performGeneSetEnrichmentAnalysis(dat, assay = 'RNA', 
                                        groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'Celltype', 'Phase'),
                                        GMT_file =  'c3.tft.v7.1.symbols.gmt')

cerebro <- exportFromSeurat(dat, file = 'HCL_Aorta.crb', experiment_name = 'Human Cell Landscape Aorta', organism = 'hg', assay = 'RNA', add_all_meta_data = T,
                            cell_cycle = 'Phase',  groups = c('pca_clusters', 'gene_clusters', 'seurat_clusters', 'Celltype', 'Phase'),
                            nUMI = 'nCount_RNA', nGene = 'nFeature_RNA')
launchCerebro()
