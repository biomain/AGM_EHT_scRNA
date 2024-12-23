import sys
import decoupler as dc
import scanpy as sc
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import pypath
import tkinter as tk
#import scimap as sm
#import scanpy.io.anndata as anndata
sc.settings.set_figure_params(dpi=200, frameon=False)
sc.set_figure_params(dpi=200)
sc.set_figure_params(figsize=(4, 4))

#adata = home/cs/python/decopular/AGMCELLTYPE_anndata.h5ad
#adata = "//wsl.localhost/Ubuntu/home/cs/python/decopular/AGM45_anndata.h5ad"
#print (adata)


#adata = ad.read_h5ad("sys.argv[1]")


adata = ad.read_h5ad("AGMF_anndata.h5ad")
print(adata)


#bdata = adata[adata.obs.namedclusters == ("zero", "1", "2","3","4","5","6","7","8","9","10","11","12","13","14","HPC","16","17","18","HE","ERY","21","MK","23","24","25","26","27","28")]

#adata = sm.hl.dropFeatures(adata, drop_groups=['29'], groups_column='namedclusters')

net = dc.get_collectri(organism='human', split_complexes=False)
net


dc.run_ulm(
    mat=adata,
    net=net,
    source='source',
    target='target',
    weight='weight',
    verbose=True,
	use_raw=False
)



#print(bdata.obsm['ulm_estimate'])

#adata.obsm['collectri_ulm_estimate'] = adata.obsm['ulm_estimate'].copy()
#adata.obsm['collectri_ulm_pvals'] = adata.obsm['ulm_pvals'].copy()
#print(adata)

#ad.AnnData.write(adata, "ulmAGM45.h5ad")

#adata = ad.read_h5ad("ulmAGM45.h5ad")


acts = dc.get_acts(adata, obsm_key='ulm_estimate')
#print(acts)


df = dc.rank_sources_groups(acts, groupby='fc', reference='rest', method='t-test_overestim_var')
print(df)
df.to_csv('TF_fc_vs_all.csv')
n_markers = 10
source_markers = df.groupby('group').head(n_markers).groupby('group')['names'].apply(lambda x: list(x)).to_dict()

df = dc.rank_sources_groups(acts, groupby='fg', reference='rest', method='t-test_overestim_var')
print(df)
df.to_csv('TF_fg_vs_all.csv')

dc.plot_network(
    net=net,
    n_sources=['RUNX1', 'PRDM16', 'SPI1','KLF1','FOXC1'],
    n_targets=15,
    node_size=100,
    s_cmap='white',
    t_cmap='white',
    c_pos_w='darkgreen',
    c_neg_w='darkred',
    figsize=(10, 10)
)






#print(source_markers)

#CAN ALSO TURN ON 


plt.switch_backend('TkAgg')

sc.pl.matrixplot(acts, source_markers,'namedclusters', dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='RdBu_r')


sc.pl.matrixplot(acts, source_markers,'fg', dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='RdBu_r')

df = dc.rank_sources_groups(acts, groupby='cell.typep', reference='rest', method='t-test_overestim_var')
print(df)

n_markers = 10
source_markers = df.groupby('group').head(n_markers).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
print(source_markers)

plt.switch_backend('TkAgg')

sc.pl.matrixplot(acts, source_markers,'cell.typep', dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='RdBu_r')

#desired_order = ["1","2","3","7","10","HE","HPC","MK","ERY","6","8","13","14","18","zero","4","5","9","11","12","16","17","21","23","24","25","27","28"]
#sc.pl.matrixplot(adata, source_markers, 'namedclusters', groupby= desired_order, dendrogram=True, standard_scale='var', colorbar_title='Z-scaled scores', cmap='RdBu_r')
#sc.pl.matrixplot(adata, source_markers, 'namedclusters', groupby= 'cell.type', dendrogram=True, standard_scale='var', colorbar_title='Z-scaled scores', cmap='RdBu_r')

#CAN TURN ON IT DOES CELLTYPE AND CLUSTER
#sc.pl.matrixplot(acts, source_markers,'namedclusters', dendrogram=True, standard_scale='var',
             #    colorbar_title='Z-scaled scores', cmap='RdBu_r')

# groupby("1","2","3","7","10","HE","HPC","MK","ERY","6","8","13","14","18","zero","4","5","9","11","12","16","17","21","23","24",,"25","27","28"))



#df = dc.rank_sources_groups(acts, groupby='namedclusters', reference='re>#print(df)

#n_markers = 3
#source_markers = df.groupby('group').head(n_markers).groupby('group')['n>source_markers

#sc.pl.matrixplot(acts, source_markers,'namedclusters' ,dendrogram=True, colorbar_title='Z-scaled scores', cmap='RdBu_r')
