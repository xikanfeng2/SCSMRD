######## pesudotime analysis using scanpy 
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns 
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
adata=sc.read_csv('muscle_cell_scanpy.csv')
sc.pl.highest_expr_genes(adata, n_top=20, save=True)
sc.pp.recipe_zheng17(adata,5000)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=4, n_pcs=20)
sc.tl.leiden(adata)
sc.tl.paga(adata, groups='leiden')
# optional figures show the network graph
sc.pl.paga(adata, threshold=0.03, show=False,save='paga.pdf')
sc.tl.draw_graph(adata, init_pos='paga',layout= 'lgl')
sc.tl.draw_graph(adata, init_pos='paga',layout= 'kk')
sc.pl.draw_graph(adata, color='leiden', legend_loc='on data',save='graph.pdf')
sc.tl.draw_graph(adata, init_pos='paga')
###save cell types and select original accodring to the biological information
df=pd.DataFrame(adata.obs['leiden'])
df.to_csv('mouse-scanpy-clust.csv')
adata.uns['iroot'] = np.flatnonzero(adata.obs['leiden']  == '4')[0]
sc.tl.dpt(adata)
df2=pd.DataFrame(adata.obs['dpt_pseudotime'])
df2.to_csv('cell_pesudotime.csv')
celltype=pd.read_csv('muscle_cell_ann.csv',index_col=0,header=None)
adata.obs['leiden_anno']=celltype
sc.pl.draw_graph(adata, color=['leiden_anno', 'dpt_pseudotime'], legend_loc='on data',save=True,legend_fontsize=6,edges=True)
plt.figure(figsize=(30,30))
