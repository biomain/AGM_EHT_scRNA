#Transcription factor identification from gene expression using decoupler. -CTS 2024
#Despite the amazing oppertunity to have an "R" pun this is a python package meant to look at anndata objects. 
#There are several different functions, but this script covers how to take an adata object and use the 
#gene expression to determine the active transcritpion factors for each group of cells.

#Despite my earlier claims of missed oppertunities for R puns, there is actually an R version of decoupler, but it 
#lacks some of the functionallity and also takes significantly longer and much more memory. 

#Aside from downloading the following packages, you do need to convert a seurat object into
#an anndata file. This can be done easily with the "Seurat Converter" script. 

#Also the code will stop running if you have a graph open. If you close it then it will resume. 
#This is a slight issue when running directly from command lines such as terminal or Ubuntu. 


#pip install decoupler

import sys
import decoupler as dc
import scanpy as sc
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import pypath
import tkinter as tk

sc.settings.set_figure_params(dpi=200, frameon=False)
sc.set_figure_params(dpi=200)
sc.set_figure_params(figsize=(4, 4))


# Load a seurat object in the form of anndata object. There are multiple options for this and conversion 
#can be done with an R script.

#adata = home/cs/python/decopular/AGMLIV_anndata.h5ad
#print (adata)

#adata = ad.read_h5ad("sys.argv[1]")
#print(adata)

adata = ad.read_h5ad("AGM45_anndata.h5ad")
print(adata)


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


#If you want to save the results of the ulu to the adata object then you can run this and it will add the ulu
# to a copy of the original adata object. 

#adata.obsm['collectri_ulm_estimate'] = adata.obsm['ulm_estimate'].copy()
#adata.obsm['collectri_ulm_pvals'] = adata.obsm['ulm_pvals'].copy()
#print(adata)
#ad.AnnData.write(adata, "ulmAGM45.h5ad")
#adata = ad.read_h5ad("ulmAGM45.h5ad")



#Acts is a portion of the adata object for faster use. 
acts = dc.get_acts(adata, obsm_key='ulm_estimate')
#print(acts)


#This creates a long form dataframe with pvalues, fold overall changes, and other stats for each TF in
#each group.
#The groupby argument may change depending on how your seurat object metadata is named. Edit this to change the level
# that you want to see compared. works for clusters, cell types, and anything else as long as its in the meta.data

df = dc.rank_sources_groups(acts, groupby='cell.type', reference='rest', method='t-test_overestim_var')
#print(df)
# sets the amount of TFs per group you want to get. Also shows key TFs per group

n_markers = 3
source_markers = df.groupby('group').head(n_markers).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
#print(source_markers)


#This creates a heatmap of the TFs and if there is overlap in the groups. If dendrogram is false, then 
# the heatmap alligns by the order of the groups in the adata file rather then based on phylo tree
plt.switch_backend('TkAgg')

sc.pl.matrixplot(acts, source_markers,'cell.type', dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='RdBu_r')


#Now that we have our canidates, we can look at them individually and graphically 


#This will create umap and violin plots of the activity of a specefic transcription factor. 
sc.pl.umap(acts, color=['RUNX1', 'cell.type'], cmap='RdBu_r', vcenter=0)
sc.pl.violin(acts, keys=['RUNX1'], groupby='cell.type', rotation=90)



#This creates a network of the transcription factors and their targets. It also shows interactions between the 
#selected TFs in terms of targets. For instance this one covers tgfb signaling. Also there seems to be a slight 
#issue with this graph actually appearing but if you just have another graph coded after it then they both appear
#in separate windows

dc.plot_network(
    net=net,
    n_sources=['ZEB1', 'SMAD2', 'SMAD3'],
    n_targets=15,
    node_size=100,
    s_cmap='white',
    t_cmap='white',
    c_pos_w='darkgreen',
    c_neg_w='darkred',
    figsize=(10, 10)
)
sc.pl.umap(acts, color=['RUNX1', 'cell.type'], cmap='RdBu_r', vcenter=0)



