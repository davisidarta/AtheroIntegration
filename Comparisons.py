# Converting individual datasets from Seurat to AnnData

# Import packages
import scanpy as sc
import pandas as pd
from scipy.sparse import csr_matrix
import dbmap as dm
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px

# Set WD and initial settings
wd = '/home/davi/Documents/Bioinfo/Vascular/'

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=300, frameon=False, figsize=(6, 6), facecolor='white', fontsize=6,
                              vector_friendly=True)

results_file = wd + 'MouseIntegration.h5ad'

###
# Load AnnDatas
###
win = sc.read_h5ad(wd + 'Winkels/Winkels.h5ad')
win = win.raw.to_adata()
#
coc = sc.read_h5ad(wd + 'Cochain/Cochain.h5ad')
coc = coc.raw.to_adata()
#
lin = sc.read_h5ad(wd + 'Lin/Lin.h5ad')
lin = lin.raw.to_adata()
#
sha = sc.read_h5ad(wd + 'Sharma/Sharma.h5ad')
sha = sha.raw.to_adata()
#
kim = sc.read_h5ad(wd + 'Kim/Kim.h5ad')
kim = kim.raw.to_adata()
#
al = sc.read_h5ad(wd + 'Alencar/Alencar.h5ad')
al = al.raw.to_adata()
#
wl = sc.read_h5ad(wd + 'Williams/Williams.h5ad')
wl = wl.raw.to_adata()
#
pan = sc.read_h5ad(wd + 'Pan/Pan_Mouse.h5ad')
pan = pan.raw.to_adata()
#
wirk = sc.read_h5ad(wd + 'Wirka/Wirka_Mouse.h5ad')
wirk = wirk.raw.to_adata()
#
den = sc.read_h5ad(wd + 'Deng/Deng.h5ad')
den = den.raw.to_adata()
#
juy = sc.read_h5ad(wd + 'Juyong/Juyong.h5ad')
juy = juy.raw.to_adata()
#
dob = sc.read_h5ad(wd + 'Dobnikar/Dobnikar.h5ad')
dob = dob.raw.to_adata()

adatas = [win, coc, lin, sha, kim, al, wl, pan, wirk, den, juy, dob]

adata = adatas[0].concatenate(adatas[1], adatas[2], adatas[3], adatas[4], adatas[5], adatas[6], adatas[7],
                              adatas[8], adatas[9], adatas[10], adatas[11],
                              batch_categories=['Winkels', 'Cochain', 'Lin', 'Sharma', 'Kim', 'Alencar', 'Williams',
                                                'Pan', 'Wirka', 'Deng', 'Juyong Kim', 'Dobnikar'],
                              batch_key='Study', join='outer')

sc.pp.normalize_per_cell(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.3, max_mean=80, min_disp=0.8, batch_key='Study', n_bins=100, flavor='cell_ranger')
sc.pl.highly_variable_genes(adata)

print("Highly variable genes intersection: %d" % sum(adata.var.highly_variable_intersection))

print("Number of batches where gene is variable:")
print(adata.var.highly_variable_nbatches.value_counts())

var_genes_batch = adata.var.highly_variable_nbatches > 0

var_select = adata.var.highly_variable_nbatches > 2
var_genes = var_select.index[var_select]
len(var_genes)
adata.var['var_genes'] = var_select

adata.raw = adata

adata = adata[:, adata.var.var_genes]

data = csr_matrix(adata.X)

diff = dm.diffusion.Diffusor(ann_dist='cosine', n_jobs=10, n_neighbors=50, n_components=500, norm=False,
                             transitions=False, kernel_use='simple').fit(data)
mms = diff.transform(data)
mms = np.array(mms)
res = diff.return_dict()
plt.plot(range(0, len(res['EigenValues'])), res['EigenValues'])
adata.obsm['X_adapmap'] = mms

# Run BBKNN on the diffusion space
sc.external.pp.bbknn(adata, batch_key='Study', n_pcs=500, metric='manhattan')

# UMAP of the shared diffusion space
sc.tl.umap(adata, min_dist=0.7, spread=1.2, alpha=1.5, maxiter=800, init_pos='spectral', n_components=2)
sc.pl.umap(adata, color=['Study', 'Group', 'pred.coc'], size=50)

adata = coc.concatenate(kim, lin, sha, win, batch_categories=['Cochain', 'Kim', 'Lin', 'Sharma', 'Winkels'],
                        batch_key='Study', join='outer')

# Run Harmony on the diffusion space
sc.external.pp.harmony_integrate(adata, key='Group', basis='X_adapmap', adjusted_basis='X_adapmap_harmony')

sc.pp.neighbors(adata, use_rep='X_adapmap_harmony', n_neighbors=15, metric='cosine')
sc.tl.umap(adata, min_dist=0.6, spread=1.5, alpha=1.5, maxiter=300)
sc.pl.umap(adata, color=['Study', 'Group', 'pred.coc'], size=50)

import scanorama

corrected = scanorama.correct_scanpy(adatas)

sc.tl.umap(adata_merge, min_dist=0.7, spread=1.2, alpha=1.5, maxiter=800, init_pos='spectral', n_components=2)

adata_merge.obsm['X_dbmap_dbbknn'] = adata_merge.obsm['X_umap']
sc.pl.embedding(adata_merge, basis='dbmap', color=['Study', 'pred.coc'], size=50)
sc.pl.embedding(adata_merge, basis='dbmap', color=['Cd209a', 'Mki67', 'Trem2'], size=50)

adata_merge.write(results_file)


























