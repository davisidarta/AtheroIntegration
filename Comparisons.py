# diffusion-based batch-balanced k-nearest-neighbors
# for integration of non-linear single-cell data

# Import some libraries
import scanpy as sc
import pandas as pd
from scipy.sparse import csr_matrix
import dbmap as dm
import numpy as np

# Set WD and initial settings
wd = '/home/davi/Documents/Bioinfo/Vascular/'  # Disha, this should be changed to your working directory (HPCC?)

sc.settings.verbosity = 1            
sc.logging.print_versions()

results_file = wd + 'AtheroIntegration.h5ad'

# Load AnnDatas
adata1 = sc.read_h5ad(wd + 'Data1/Data1.h5ad')
adat1 = adata1.raw.to_adata()
# .
# .
# .
# All datasets should be individually filtered, analyzed and saved prior to integration

# Create an array of adatas
adatas = [adata1, adata2, adata3] # and so forth

# Run anisotropic diffusion maps for each adata
for adata in adatas:
  adata = adata[:, adata.var.highly_variable]
  data = csr_matrix(adata.X)
  diff = dm.diffusion.Diffusor(ann_dist='cosine', n_jobs=10, n_neighbors=15, n_components=100, norm=False, transitions=False, kernel_use='simple').fit(data)
  mms = diff.transform(data)
  mms = np.array(mms)
  res = diff.return_dict()
  plt.plot(range(0, len(res['EigenValues'])), res['EigenValues'])
  adata.obsm['X_adapmap'] = mms
  adata.obsm['X_pca_true'] = adata.obsm['X_pca'] # Let's trick bbKNN, as the scanpy external API does not allow choosing another low-dimensional representation.
  adata.obsm['X_pca'] = adata.obsm['X_adapmap'] 

# Concatenate all adatas into a single adata
all_adatas = coc.concatenate(adatas[0], adatas[1], adatas[2], batch_categories=['Adata1', 'Adata2', 'Adata3'], batch_key = 'Study', join ='outer')

# Run BBKNN on the diffusion space
sc.external.pp.bbknn(all_adatas, batch_key='Study')  

# Run Harmony on the diffusion space
sc.external.pp.harmony_integrate(all_adatas, key='Group', basis='X_adapmap', adjusted_basis='X_adapmap_harmony')
































