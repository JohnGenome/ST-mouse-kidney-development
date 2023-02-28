## Backgrounds ##
# RNA velocity: the first time derivative of the spliced mRNA abundance
# RNA velocity is determined by the balance between production of spliced mRNA from unspliced mRNA, and the mRNA degradation

import scvelo as scv
import numpy as np
import pandas as pd
import os
from scipy.sparse import csr_matrix

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization
scv.settings.set_figure_params('scvelo', dpi_save=2400)

## 1 Read your data
adata = scv.read("/path/to/mapped_RCI4O.loom", cache=True)
KNN_mat = pd.read_csv('/path/to/1_KNN_matrix.Npns.filtered.csv',index_col=0)
KNN_df = pd.read_csv('/path/to/1_KNN_metadata.Npns.filtered.csv',index_col=0)
KNN_df = KNN_df[~(KNN_df.UMAP1.isna())] # remove NaN values

KNN_mat = KNN_mat.loc[KNN_df.barcodeBA]
cell_list1 = KNN_mat.columns.values
cell_list2 = KNN_mat.index.values
KNN_mat = np.mat(KNN_mat.values)
np.sum(cell_list2 == KNN_df.barcodeBA) # 1564

## 2 Rearrange data
adata = adata[cell_list1,]
bdata = adata.copy()
bdata = bdata[cell_list2,]
# Transfer cluster id --> Categories
bdata.obs["Ann_v5"] = KNN_df.Ann_v5.astype('str').to_list()
bdata.obs["Ann_v5"] = bdata.obs["Ann_v5"].astype('category')
bdata.obs["Ann_v5"].cat.set_categories(['8','7','24','6','25','12','10','2','27','20','28','14','29'],inplace=True)

bdata.obs["Ann_KNN"] = KNN_df.Top_Ann.astype('str').to_list()
bdata.obs["Ann_KNN"] = bdata.obs["Ann_KNN"].astype('category')
bdata.obs["Ann_KNN"].cat.set_categories(['8','7','24','6','25','12','10','2','27','20','28','14','29'],inplace=True)
np.sum(bdata.obs.Ann_KNN == bdata.obs.Ann_v5) # 1481

# np.sum(bdata.X) # 1133046 same as spliced
# np.sum(bdata.layers['matrix']) # 3227691 sum of ambiguous, spliced, unspliced
# np.sum(bdata.layers['ambiguous']) # 183289
# np.sum(bdata.layers['spliced']) # 1133046
# np.sum(bdata.layers['unspliced']) # 1911356
scv.pl.proportions(bdata,groupby = "Ann_v5",save = "0_spliced&unspliced_counts.pdf")

## 3 Merge into metacell
## np.dot works on int dtype is very slow, about 30 times than float
# https://github.com/numba/numba/issues/5391
bdata.X = np.dot(KNN_mat,adata.X.todense());bdata.X = csr_matrix(bdata.X, dtype=np.float32)
# bdata.layers["ambiguous"] = np.dot(KNN_mat,adata.layers["ambiguous"].todense());bdata.layers["ambiguous"] = csr_matrix(bdata.layers["ambiguous"], dtype=np.float32)
# bdata.layers["matrix"] = np.dot(KNN_mat,adata.layers["matrix"].todense());bdata.layers["matrix"] = csr_matrix(bdata.layers["matrix"], dtype=np.float32)
bdata.layers["spliced"] = np.dot(KNN_mat,adata.layers["spliced"].astype(np.float32).todense());bdata.layers["spliced"] = csr_matrix(bdata.layers["spliced"], dtype=np.float32)
bdata.layers["unspliced"] = np.dot(KNN_mat,adata.layers["unspliced"].astype(np.float32).todense());bdata.layers["unspliced"] = csr_matrix(bdata.layers["unspliced"], dtype=np.float32)
bdata.write('data/E16_Npns.SeuratUMAP.h5ad')
bdata = scv.read('data/E16_Npns.SeuratUMAP.h5ad')

## 4 Basic preprocessing
scv.pp.filter_and_normalize(bdata, min_counts=200, n_top_genes=2000) # (1566, 55401) -> (1566, 10013) -> (1566, 3000)
scv.pp.moments(bdata, n_pcs=30, n_neighbors=30)
# Velocity Tools
scv.tl.recover_dynamics(bdata)
scv.tl.velocity(bdata, mode='dynamical') # min_r2=1e-10,min_likelihood=1e-10
scv.tl.velocity_graph(bdata)

## 5 Visualization
# Ref:https://github.com/theislab/scvelo/issues/884
KNN_df = pd.read_csv('../Out/1_KNN_metadata.Npns.filtered.csv',index_col=0)
KNN_df = KNN_df[~(KNN_df.UMAP1.isna())] # remove NaN values
bdata.obsm['X_umap'] = KNN_df[["UMAP1","UMAP2"]].values
# switch LoH color
scv.pl.velocity_embedding(bdata, basis = 'umap',color = "Ann_KNN",save = "0_velocity_embedding.39genes.n_nei30.pdf", arrow_length=5, arrow_size=3, dpi=2400,figsize = (4,4.5),palette=("#000080","#A3A500","#4169E1","#00BAE0","#FF9999","#FF0000","#EA8331","#AD0000","#FF62BC","#39B600","#FFD032","#FFFF00","#A2FF00"))
scv.pl.velocity_embedding_stream(bdata, basis = 'umap',color = "Ann_KNN",save = "1_velocity_embedding_stream.39genes.n_nei30.png",size=40,dpi=2400,linewidth=.5,arrow_size=.4,density=3,palette=("#000080","#A3A500","#4169E1","#00BAE0","#FF9999","#FF0000","#EA8331","#39B600","#AD0000","#FF62BC","#FFD032","#FFFF00","#A2FF00"))

bdata.obs.to_csv('data/E16_Npns.metadata.csv')


