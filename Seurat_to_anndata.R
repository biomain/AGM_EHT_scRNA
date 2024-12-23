######################## DECOUPLER ##############################
#This requires a bit of set up, a command line terminal, and python 
#installed on your computer, but this is the best option for determining 
#active or overactive transcription factors.

#First we create an anndata object and then we move on to the python 
#script "Tfact.py". 

library(scCustomize)
library(reticulate)
library(Seurat)

#You need anndata and a virtual environment. This should create a simple env
#and install anndata
virtualenv_create("r-reticulate")
py_install("anndata")

#Alternatively, you if you have an existing virtual environment, you can use it 
#here
reticulate::use_virtualenv("python3")
virtualenv_install("python3", "anndata")
adata <- import("anndata")

#next we load in the Seurat object if not already there, and we convert it 
#to an anndata object and save it. This can take a few minutes. 
# Your Seurat object is the first argument. 

#readRDS("path/to/your/seurat.RDS")

as.anndata(x = AGM, file_path = "~/", file_name = "AGMF_anndata.h5ad")
#From here you can run the anndata object through decoupler. 