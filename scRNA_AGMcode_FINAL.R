# WELCOME TO THE SCRIPT THAT CONTAINS THE MOST OF THE BIOINFORMATICS 
# CODE NEEDED TO REPLICATE THE RESULTS FOUND IN 
"Integrated local and systemic communication factors regulate nascent hematopoietic
progenitor escape during developmental hematopoiesis"
#Most of the code is in R, however some of it is in python is found in 
#an alternate script. 
#Most of this is just modifications from the seurat and cellchat 
# guides on their website or github. 

#LIBRARYIES 
#There are quite a few invovled with this work. Load all of these. 


library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
library(dplyr)
library(data.table)
library(RColorBrewer)
library(SingleR)



#OBJECTS. 
#The most important part of this script is the creation of the objects. 
#All of our data comes from the 10X scRNA data produced in the Calvanese 2022 
#Paper and is found at GSE162950. 

#AGM The first thing we need to do is to make the AGM dataset. This is the main 
#object we work with. It contains scRNA from embronic AGMs from week 4.5, 5, 5, and 6
#post conception. Yes there are duplicates of week 5. 

#The main plan is to download the datasets and merge them together. 
#There is a slight problem in the form of the batch effect, but we will reduce 
#That with harmony intergration. You will need to set your own paths 

# Load the AGMdataset
AGM_4.5_week_data <- Read10X("PATH/fastqs_SRA/fullQscores_fastqs/aorta-gonad-mesonephros4-5week/outs/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
AGM_4.5_week <- CreateSeuratObject(counts = AGM_4.5_week_data, project = "AGM_4_5week", min.cells = 3)


AGM_5_week_data <- Read10X("PATH/fastqs_SRA/fullQscores_fastqs/aorta-gonad-mesonephros-5-week/outs/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
AGM_5_week <- CreateSeuratObject(counts = AGM_5_week_data, project = "AGM_5week", min.cells = 3)

#THIS IS NOT ACTUALLY A DIFFERENT TIME POINT. But is is labeled as 5.5 and 5b at times. For our purposes it is not relevent though. 
AGM_5.5_week_data <- Read10X("PATH/fastqs_SRA/fullQscores_fastqs/aorta-gonad-mesonephros-5-5-week/outs/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
AGM_5.5_week <- CreateSeuratObject(counts = AGM_5.5_week_data, project = "AGM_5_5week", min.cells = 3)

AGM_6_week_data <- Read10X("PATH/fastqs_SRA/fullQscores_fastqs/aorta-gonad-mesonephros-6-week/outs/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
AGM_6_week <- CreateSeuratObject(counts = AGM_6_week_data, project = "AGM_6week", min.cells = 3)

#merge all
AGM_merged_all <- merge(AGM_4.5_week, y = c(AGM_5_week, AGM_5.5_week, AGM_6_week))

AGM_merged_all[["percent.mt"]] <- PercentageFeatureSet(AGM_merged_all, pattern = "^MT-")
AGM_merged_filt <- subset(AGM_merged_all, subset = nFeature_RNA > 100 & percent.mt < 10)
AGM_merged_filt$orig.ident <- factor(x = AGM_merged_filt$orig.ident, levels = c("AGM_4_5week", "AGM_5week", "AGM_5_5week", "AGM_6week"))

AGM_merged_filt <- subset(AGM_merged_filt, subset = nCount_RNA > 500)
saveRDS(AGM_merged_filt, file = "AGM_merged_filt.rds")
AGM_merged_filt <- NormalizeData(AGM_merged_filt) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
AGM_merged_filt <- RunHarmony(AGM_merged_filt, group.by.vars = "orig.ident")
AGM_merged_filt <- RunUMAP(AGM_merged_filt, reduction = "harmony", dims = 1:30)
AGM_merged_filt <- FindNeighbors(AGM_merged_filt, reduction = "harmony", dims = 1:30) %>% FindClusters()

#Alternativly, In the github, I have uploaded the AGM dataset if you would like to skip a few steps and 
#50 gigabite scRNA 10X files. 
#CLusters at "fc" and groups are "fg"

#What are we looking at SiingleR time 

#SingleR to identify groups of cells. You provide it a scRNA dataset with annotations
# and then it matches those to your seurat dataset. In this case we actually need to convert our seurat 
# object to a SingleCellExperiment object. 

# We will be using blueprint and HPCA as our scRNA datasets for annotation. 
# Now this is slightly depreciated, and while it works for the version we are using, 
# it may not work for yours so feel free to use the celldex package to download 
# the datasets. 

#Going forward I will refer to it as the AGM.
AGM <-AGM_merged_filt

#You can also start with the github seurat object here. 

bp.ref <- BlueprintEncodeData()
hp.ref <- HumanPrimaryCellAtlasData()
###bpalt.ref <- celldex::BlueprintEncodeData()### ALTERNATE IF DEPRECIATED 

options(ggrepel.max.overlaps = 100)

AGM.sce <- as.SingleCellExperiment(AGM)
#Now this does have a chance of not working in newer seurat objects especially if 
#they have "layers". Consult the SingleR/Seurat documentation if this is the case for you. 
# More of an issue for non harmony intergrated datasets. It is possible to revert seurat objects 
# to earlier versions and then run the as.SingleCellExperiment function. 

#For the most part the .sce object is only temproary for singleR to be run. 
#afterwards we can delete it and go back to the classsic seurat object. 

#These next 3 lines can and will take a bit. I am running this on a laptop with 32 gigs of RAM 
# and it takes about 10 minutes for the hpref.main code to run. 

bpref.main <- SingleR(test = AGM.sce,assay.type.test = 1,ref = bp.ref,labels = bp.ref$label.main)
bpref.fine <- SingleR(test = AGM.sce,assay.type.test = 1,ref = bp.ref,labels = bp.ref$label.fine)
hpref.main <- SingleR(test = AGM.sce,assay.type.test = 1,ref = hp.ref,labels = hp.ref$label.main)

AGM@meta.data$bp.main   <- bpref.main$labels
AGM@meta.data$bp.fine   <- bpref.fine$labels
AGM@meta.data$hp.main   <- hpref.main$labels

#We have now added the Single R annotations to the AGM seurat object. 
# you can treat them like any other metadata slot such as original ident 
# or seurat clusters.

DimPlot(AGM, label = T, repel= T, label.size = 3, group.by = "bp.fine") + NoLegend()
DimPlot(AGM, label = T, repel= T, label.size = 3, group.by = "hp.main") + NoLegend()


#Now as you have likely noticed, the dimplot with the annotations is a bit 
# "disorganized". There are quite a few different labels, so to simplify things 
# and actually determine what our cells are we can use Hoverlocator. 
#SingleR graphs can be quite granular so using hoverlocator can make it clear which 
#clusters are which cell types. Additionally, the gene labels will indicate original 
# identity within AGM or Liver. for the AGM cells 1_1 1_2 and so on correspond to 
# week 4.5,5, 5.5, and 6. You can switch between grouping by cluster or ident 
# and showing cluster or ident on the hoverlocator. You can also show multiple 
#variables by using the code for gg4. There are 2 main annotation databases used for this 
#being the Human primary cell atlas and the blueprint data. The former is better 
#for more general annotations while bp is superior for immune cells. We used both 
# to annotate the cells. As with any cell annotation there is some ambiguity 
# and should be taken with a slight grain of salt. Also the fact that this is 
# embryonic tissue complicates things. Cell annotation via scRNA expression 
# is not more art than science, but it certianly has some "art" to it. 




#You can change any of these around. The group.by will refer to the coloring 
# of the cells on the dimplot. You can pick any metadata varible. 
# The hover locator part will list the varibles you see when you mouse over 
# a cell on the dimplot. You can add multiple metadata values. 
# Also if you are running this on a High powered cluster or in command line
# I would recomend swithing over to Rstudio of some type as it tends to crash on 
#command line 


gg1 <- DimPlot(AGM, label = T, repel= T, label.size = 3, group.by = "seurat_clusters") + NoLegend()
HoverLocator(gg1, information = FetchData(AGM,vars = "seurat_clusters"))

gg2 <- DimPlot(AGM, label = T, repel= T, label.size = 3, group.by = "bp.main") + NoLegend()
HoverLocator(gg2, information = FetchData(AGM,vars = "bp.main"))

gg3 <- DimPlot(AGM, label = T, repel= T, label.size = 3, group.by = "seurat_clusters") + NoLegend()
HoverLocator(gg3, information = FetchData(AGM,vars = "hp.main"))

gg4 <- DimPlot(AGM, label = T, repel= T, label.size = 3, group.by = "seurat_clusters") + NoLegend()
HoverLocator(gg4, information = FetchData(AGM,vars = c("bp.main", "seurat_clusters")))


#A bit more modular of the code. 
groupings <- "bp.fine"
vars = c("fg", "bp.fine")
gg5 <- DimPlot(AGM, label = T, repel= T, label.size = 3, group.by = groupings) + NoLegend()
HoverLocator(gg5, information = FetchData(AGM,vars = vars))



# Okay we have intepreted the identties of our clusters. 
# Once identified we need to group the cells and prepare for cellchat. 
# while one could simply memorize all the clusters numbers by staring at dimplots 
# for hours on end, it is a bit more practical to label the clusters. 


#If you are using the github seruat object, then you can skip this part. 

#If you are attempting to use this code for your own project, I recomend 
# that you name your object prior to making it a cellchat. 
#You can alter the names after the fact but it is a hassle.

"ALSO CELLCHAT HATES THE NUMBER 0, if you have any cluster with the name
0 then it will not run properly."



AGM <- SetIdent(AGM, value = seurat_clusters)
#AGM <- RenameIdents(AGM, "19" = "HE")# and so on 
AGM <-RenameIdents(AGM, "0" = "Stromal1",
             "1" = "Endo1",
             "2" = "Endo2",
             "3" = "Endo3",
             "4" = "Fibro1",
             "5" = "Fibro2",
             "6" = "Macro1",
             "7" = "Endo4",
             "8" = "Lymph",
             "9" = "Fibro3",
             "10" = "Endo5",
             "11" = "Peri",
             "12" = "Fibro4",
             "13" = "Gran",
             "14" = "Macro2",
             "15" = "HSC",
             "16" = "Fibro5",
             "17" = "Stromal2",
             "18" = "Macro3",
             "19" = "HE",
             "20" = "ERY",
             "21" = "Kidney1",
             "22" = "MK",
             "23" = "Fibro6",
             "24" = "Liver",
             "25" = "Kidney2",
             "26" = "Kidney3",
             "27" = "Fibro7",
             "28" = "Stromal3"
)
AGM@meta.data$fc <- AGM@meta.data$seurat_clusters
AGM <- SetIdent(AGM, value = fc)

#For the group level just do repeate the process but without the numbers of the clusters 
# so it wouldnt be Endo 1, Endo 2 and Endo3 but rahter Endo, Endo, Endo 

#CELLCHAT TIME. 
#Time to chat about CellChat
# Cellchat is our tool for examining the ligand receptor interactions 
#between the liver and the AGM. It generally runs the same each time, but you 
#can make changes to a few parameters. Database and group.by do have some variability
# the entire cellchat database is quite large and if it needs to look at every type 
# of interaction it may find irrelevent interactions, waste computing power/time, 
# and obsure interactions. If you are looking for long distance signalling 
# then cell cell contact or ECM interactions arent going to be relevent or realisitic. 
# you can set the program to look at all interactions or only a subset.

# Additionally, you can group cells by any meta label in the object. 
# you can group them by cell type or cluster. We choose to do both and made 2 objects
AGM <- SetIdent(AGM, value = "fc")
cellchat <- createCellChat(object = AGM, group.by = "ident", assay = "RNA")
CellChatDB <- CellChatDB.human
###########################################################################
#This gives you the option of every possible dataset or, just 1 of the 4 types. 
# you can also do 2 or 3 if that is desired
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 
#If you want it all use the code below instead.
#CellChatDB.use <- CellChatDB
#########################################################################
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#############################################################
#IF YOU ARE ON CELLCHAT 2.1+ THEN USE SMOOTHDATA. 
#THE PROJECTDATAFUNCTION HAS BEEN REMOVED IN NEWER VERSIONS.
cellchat <- smoothData(cellchat, adj = PPI.human)
#cellchat <- projectData(cellchat, PPI.human)
##################################################
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
#df.net <- subsetCommunication(cellchat, signaling = c("CXCL"))
#pathways.show <- c("CXCL")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat)
saveRDS(cellchat, "AGM_FC_ALL.rds")

AGM <- SetIdent(AGM, value = "fg")
cellchat <- createCellChat(object = AGM, group.by = "ident", assay = "RNA")
CellChatDB <- CellChatDB.human
###########################################################################
#This gives you the option of every possible dataset or, just 1 of the 4 types. 
# you can also do 2 or 3 if that is desired
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 
#If you want it all use the code below instead.
#CellChatDB.use <- CellChatDB
#########################################################################
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#############################################################
#IF YOU ARE ON CELLCHAT 2.1+ THEN USE SMOOTHDATA. 
#THE PROJECTDATAFUNCTION HAS BEEN REMOVED IN NEWER VERSIONS.
cellchat <- smoothData(cellchat, adj = PPI.human)
#cellchat <- projectData(cellchat, PPI.human)
##################################################
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
#df.net <- subsetCommunication(cellchat, signaling = c("CXCL"))
#pathways.show <- c("CXCL")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat)
saveRDS(cellchat, "AGM_FG_ALL.rds")

#CELLCHAT ANALYSIS 
#So after about an hour to 5 depending on your computer, you will now have a 
#cellchat object. What can we do with it and what can we learn from it. 


#If graphing isnt your style then dataframes are the way to go with cellchat.
#You can specefy the targets and the recpients. If you are using a named cluster 
#such as the HE, then you need to set it up such that targets.use = c("HE", any other targets you want)
df.ccLRI<-subsetCommunication(cellchat, slot.name = "net")
df.ccPATHS<-subsetCommunication(cellchat, slot.name = "netP")

#GSEA and Marker Genes 

#We have our clusters and cell types, but how can we find out what they are doing 
#and what makes them unique. Well, there is Seurat Find Markers and FGSEA
#to determine those facts. FindMarkers finds genes of differential expression between 
#two metadata idents and fgsea determines if those markers match with a geneset 
#such as inflammation.

library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(presto)
library(Seurat)

#Download the GSEA librarys from msigdb 

h1_db <- msigdbr(species = "Homo sapiens", category = c("H"))
h1_db <- h1_db %>% split(x = .$gene_symbol, f = .$gs_name)
c2_db <- msigdbr(species = "Homo sapiens", category = c("C2"))
c2_db <- c2_db %>% split(x = .$gene_symbol, f = .$gs_name)
c3_db <- msigdbr(species = "Homo sapiens", category = c("C3"))
c3_db <- c3_db %>% split(x = .$gene_symbol, f = .$gs_name)
c5_db <- msigdbr(species = "Homo sapiens", category = c("C5"))
c5_db <- c5_db %>% split(x = .$gene_symbol, f = .$gs_name)

#Find your markers for your cluster of interest. If you set ident2 to be blank 
#then it will do ident 1 verses all other clusters. 

#AGMHEvsEndo <- FindMarkers(AGM, group.by = "fc",ident.1 = c("HE"), ident.2 = c("Endo1","Endo2","Endo3","Endo4","Endo5"))
#AGMHEvsHSC <- FindMarkers(AGM, group.by = "fc",ident.1 = c("HE"), ident.2 = c("HPSC"))
##AGMHEvsHSCMK <- FindMarkers(AGM, group.by = "fc",ident.1 = c("HE"), ident.2 = c("HPSC","MK"))
#AGMHSCvsEndo <- FindMarkers(AGM, group.by = "fc",ident.1 = c("HPSC"), ident.2 = c("Endo1","Endo2","Endo3","Endo4","Endo5"))
#AGMHSCvsEndoHE <- FindMarkers(AGM, group.by = "fc",ident.1 = c("HPSC"), ident.2 = c("Endo1","Endo2","Endo3","Endo4","Endo5", "HE"))


#This is where we actually do the GSEA. Select your marker/comparison and 
#adjust the names to fit the GSEA library you are using. 
#Currently this will show the Gene sets that are up-regulated in the HE compared
#to the endothelial cells. To change the database adjust the h1_db to be another database


Markers_AGM_HEvsENDO <- FindMarkers(AGM, c("HE"), c("Endo"))
exm <-rownames_to_column(Markers_AGM_HEvsENDO,"features")
respls <- exm %>%
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(features, avg_log2FC)
ranks <- deframe(respls)
AGMHEvsEndo_h1<- fgsea(pathways = h1_db,
                       stats    = ranks,
                       minSize  = 5,
                       maxSize  = 500, nproc=1)
Markers_AGM_HEvsENDO <- FindMarkers(AGM, c("HE"), c("Endo"))
exm <-rownames_to_column(Markers_AGM_HEvsENDO,"features")
respls <- exm %>%
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(features, avg_log2FC)
ranks <- deframe(respls)
AGMHEvsEndo_c2<- fgsea(pathways = c2_db,
                       stats    = ranks,
                       minSize  = 5,
                       maxSize  = 500, nproc=1)

Markers_AGM_HEvsHSC <- FindMarkers(AGM, ident.1 = c("HE"), ident.2 = c("HSC"))
exm <-rownames_to_column(Markers_AGM_HEvsHSC,"features")
respls <- exm %>%
  arrange(desc(avg_log2FC)) %>% 
  dplyr::select(features, avg_log2FC)
ranks <- deframe(respls)
AGMHEvsHSC_h1<- fgsea(pathways = h1_db,
                       stats    = ranks,
                       minSize  = 5,
                       maxSize  = 500, nproc=1)


#GSEA Figures. 
#What if you want to map multiple pathways but they are in different datasets 
#List the pathways you want. This for instance would be Hallmark and C2. 
#This was used for 2D

Hpaths = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_PI3K_AKT_MTOR_SIGNALING","HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_ANGIOGENESIS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
C2paths= c("AZARE_STAT3_TARGETS", "BENPORATH_PRC2_TARGETS", "KONDO_EZH2_TARGETS", "LU_TUMOR_ANGIOGENESIS_UP","CERVERA_SDHB_TARGETS_2")

#Place your pathway lists here and create the F1` pathway list. The f1 is the combined. 
H1_selected <- AGMHEvsEndo_h1 %>%
  filter(pathway %in% Hpaths)
C2_selected <- AGMHEvsEndo_c2 %>%
  filter(pathway %in% C2paths)
HEvEndoF1 <- rbind(H1_selected, C2_selected, fill = TRUE)


ggplot(HEvEndoF1,
       aes(x = NES, y = pathway)) + 
  geom_point(aes(color = padj , size = 8)) +
  theme_bw() +
  scale_colour_gradient(limits=c(0, .05), high = "blue", low = "red")  +
  scale_size_area(limits=c(0,10), breaks=c(1,5,10)) + #max_size = 0.41, 
  ylab(NULL) +
  xlim(-4, 4) + scale_y_discrete(limits = HEvEndoF1$pathway) + 
  theme(axis.text.x = element_text(size=10, color="black"),
        axis.text.y = element_text(size=6, color="black")) +
  ggtitle("HEvEndo") 


#At this point you may be asking, "Where is the liver" 
#Well we will be making a new Seurat object with both the Liver 
#cells and select AGM cells. Now you can make it so it isnt select, 
#but it gets rather computationally intensive and can make figures 
#far less readable. We will be taking the Endo, HE, HSC, and MK cells from the 
#AGM and adding them into a Seurat object containing week 4.5, week 5b, and week 6 
# Fetal(technically embryonic) Liver. There is no liver from 5a. Additionally, 
#5A and 5B are the same time just different replicates. 

#We will not be doing integration with this. Yes this does lead to the batch effect, 
#But at the same time, integration leads to the loss of distinction between critical cells
#like the HE. This does add a grain of salt, but non integration is the best option
#given the data and methods cellchat uses. It is a catch 22.

load("/Users/CS/Documents/RObjects/EDFig2_5/seurat_object.Rdata")
LIV45 = UpdateSeuratObject(object = sample)

load("/Users/CS/Documents/RObjects/EDFig3_12/seurat_object.Rdata")
LIV5 = UpdateSeuratObject(object = sample)

load("/Users/CS/Documents/RObjects/EDFig3_13/seurat_object.Rdata")
LIV6 = UpdateSeuratObject(object = sample)


AGMless <- subset(AGM, idents= c("Endo", "HE","Mk","HSC"))



metadata <- LIV45@meta.data
metadata$orig.ident <- "LIV45"
LIV45@meta.data <- metadata

metadata <- LIV5@meta.data
metadata$orig.ident <- "LIV5"
LIV5@meta.data <- metadata

metadata <- LIV6@meta.data
metadata$orig.ident <- "LIV6"
LIV6@meta.data <- metadata


metadata <- AGMless@meta.data
metadata$orig.ident <- "AGML"
AGMless@meta.data <- metadata

#This can crash weaker computers so save before hand. 
LIV45 <- subset(LIV45, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mito < .1)
LIV5 <- subset(LIV5, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mito < .1)
LIV6 <- subset(LIV6, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mito < .1)

CAL <- merge(AGMless, y = c(LIV5, LIV6, LIV45), add.cell.ids = c("AGM", "LIV5","LIV6", "LIV45"))
CAL <- NormalizeData(object = CAL)
CAL <- FindVariableFeatures(object = CAL)
CAL <- ScaleData(object = CAL)
CAL <- RunPCA(object = CAL)
CAL <- FindNeighbors(object = CAL, dims = 1:30)
CAL <- FindClusters(object = CAL)
CAL <- RunUMAP(object = CAL, dims = 1:30)
DimPlot(object = CAL, reduction = "umap")

CAL <- SetIdent(CAL, value = "fg")
cellchat <- createCellChat(object = CAL, group.by = "ident", assay = "RNA")
CellChatDB <- CellChatDB.human
###########################################################################
#This gives you the option of every possible dataset or, just 1 of the 4 types. 
# you can also do 2 or 3 if that is desired
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 
#If you want it all use the code below instead.
#CellChatDB.use <- CellChatDB
#########################################################################
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#############################################################
#IF YOU ARE ON CELLCHAT 2.1+ THEN USE SMOOTHDATA. 
#THE PROJECTDATAFUNCTION HAS BEEN REMOVED IN NEWER VERSIONS.
cellchat <- smoothData(cellchat, adj = PPI.human)
#cellchat <- projectData(cellchat, PPI.human)
##################################################
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
#df.net <- subsetCommunication(cellchat, signaling = c("CXCL"))
#pathways.show <- c("CXCL")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat)
saveRDS(cellchat, "CAL_FG_ALL.rds")







#FIGURES
#A picture is worth 1000 words. A figure is probably more  

#I usually set the pathways just for modular code. 
pathways.show <- "CXCL"

#ScatterPlots 
netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show, do.label = F)

#Heatmaps 
netVisual_heatmap(cellchat,signaling = pathways.show ,remove.isolate = T, targets.use = "HE")

#ViolinPlots # The additional lines are to control the label order
VlnPlot(AGM, "PROCR", idents = c("HE","HSC","Endo"), pt.size = 0) + scale_x_discrete(limits = c("Endo","HE","HSC")) + NoLegend() + theme(axis.text.y = element_text(size=20, color="black"))

#Chord Plots 
netVisual_chord_cell(cellchatAL, "TNF", targets.use = "HE")

#Now if you would like the exact code for the exact figures. They can be found below

#Figure1 
#1A
DimPlot(AGM, group.by = "fc")
#1B
redo=c("AFP","APOA1","ALB","NPHS2","NPHS1","GATA3","REN","LUM","CRABP1","COL1A1","POSTN","PTN","PDGFRA","C1QA","CD14","LYVE1","RNASE2","LYZ","CD69","IL7R","CD7","HLF","GFI1","MYB","MYCN","ITGA2B","SELP","GP9","HBE1","HBZ","GYPA","GJA5","CXCR4","CDH5","KDR","TIE1","APLNR","NRP2")
DotPlot(AGM, features=redo, group.by = "fgdot" ,cols=c("grey90","red3")) +
  coord_flip() + theme(axis.text.x=element_text(angle=45, hjust=1))

#1C, 1F, 1I. Arrows and recoloring was done in Biorender 
netAnalysis_signalingRole_scatter(cellchatAGMG, signaling = "TGFb", do.label = F)
netAnalysis_signalingRole_scatter(cellchatAGMG, signaling = "CXCL", do.label = F)
netAnalysis_signalingRole_scatter(cellchatAGMG, signaling = "ADGRG", do.label = F)

# 1D, 1G, 1J
VlnPlot(AGM,features = c("TGFB1","TGFBR1","TGFB2", "TGFBR2","ACVR1") ,group.by = "fg", pt.size = 0, stack = T, flip = T, idents = c("HE","HSC", "Endo","Mk","Peri", "Macro", "Gran")) + scale_x_discrete(limits = c("Gran","Peri","Mk","Endo","HE","HSC")) + NoLegend() + theme(axis.text.y = element_text(size=20, color="black"))
VlnPlot(AGM,features = c("CXCL12", "CXCR4") ,group.by = "fg", pt.size = 0, stack = T, flip = T, idents = c("HE","HSC", "Endo", "Stromal","Fibro","Peri")) + scale_x_discrete(limits = c("Stromal","Fibro","Peri","Endo","HE","HSC")) + NoLegend() + theme(axis.text.y = element_text(size=20, color="black"))
VlnPlot(AGM,features = c("KDR", "GATA2","F2R","ADGRG6", "COL4A1", "COL4A2", "COL4A5") ,group.by = "fg", pt.size = 0, stack = T, flip = T, idents = c("HE","HSC", "Endo", "Fibro","Peri")) + scale_x_discrete(limits = c("Fibro","Peri","Endo","HE","HSC")) + NoLegend() + theme(axis.text.y = element_text(size=20, color="black"))

#1E 1H 1K
netVisual_bubble(cellchatAGMG, signaling = "TGFb", targets.use = "HE", remove.isolate = T, sources.use = c("Endo","MK","Peri","HE","Gran")) + scale_size_continuous(range = c(1,12))
netVisual_bubble(cellchatAGMG, signaling = "CXCL", targets.use = "HE", remove.isolate = T) + scale_size_continuous(range = c(1,12))
netVisual_bubble(cellchatAGMG, signaling = "ADGRG", targets.use = "HE", remove.isolate = T, sources.use = c("Endo","Fibro","Peri","HE"))

#Figure 2 Some of this is in python unfortunately and found in that script. 

#2A
netVisual_chord_cell(cellchatAGMG, "TNF", targets.use = "HE")
#2B
VlnPlot(AGM,features = c("TNF", "TNFRSF1A","TNFRSF1B") ,group.by = "fg", pt.size = 0, stack = T, flip = T, idents = c("HE","HSC", "Endo", "Gran","Macro")) + scale_x_discrete(limits = c("Gran","Macro","Endo","HE","HSC")) + NoLegend() + theme(axis.text.y = element_text(size=20, color="black"))
#2C
netVisual_bubble(cellchatAL, signaling = "TNF", remove.isolate = T, targets.use = c("HE")) + scale_size_continuous(range = c(1,12))
#2D
H1_selected <- AGMHEvsEndo_h1 %>%
  filter(pathway %in% Hpaths)
C2_selected <- AGMHEvsEndo_c2 %>%
  filter(pathway %in% C2paths)
HEvEndoF1 <- rbind(H1_selected, C2_selected, fill = TRUE)


ggplot(HEvEndoF1,
       aes(x = NES, y = pathway)) + 
  geom_point(aes(color = padj , size = 8)) +
  theme_bw() +
  scale_colour_gradient(limits=c(0, .05), high = "blue", low = "red")  +
  scale_size_area(limits=c(0,10), breaks=c(1,5,10)) + #max_size = 0.41, 
  ylab(NULL) +
  xlim(-4, 4) + scale_y_discrete(limits = HEvEndoF1$pathway) + 
  theme(axis.text.x = element_text(size=10, color="black"),
        axis.text.y = element_text(size=6, color="black")) +
  ggtitle("HEvEndo") 

H1_selected <- AGMHEvsHSC_h1 %>%
  filter(pathway %in% Hpaths)
C2_selected <- AGMHEvsHSC_c2 %>%
  filter(pathway %in% C2paths)
HEvHSCF1 <- rbind(H1_selected, C2_selected, fill = TRUE)


ggplot(HEvHSCF1,
       aes(x = NES, y = pathway)) + 
  geom_point(aes(color = padj , size = 8)) +
  theme_bw() +
  scale_colour_gradient(limits=c(0, .05), high = "blue", low = "red")  +
  scale_size_area(limits=c(0,10), breaks=c(1,5,10)) + #max_size = 0.41, 
  ylab(NULL) +
  xlim(-4, 4) + scale_y_discrete(limits = HEvHSCF1$pathway) + 
  theme(axis.text.x = element_text(size=10, color="black"),
        axis.text.y = element_text(size=6, color="black")) +
  ggtitle("HEvHSC") 

#2E Python

#2F
netVisual_chord_cell(cellchatAGMG, "SPP1", targets.use = "HE")

#2G
VlnPlot(AGM,features = c("SPP1", "ITGAV","ITGA5","ITGB1") ,group.by = "fg", pt.size = 0, stack = T, flip = T, idents = c("HE","HSC", "Endo", "Stromal","Macro")) + scale_x_discrete(limits = c("Macro","Endo","HE","HSC")) + NoLegend() + theme(axis.text.y = element_text(size=20, color="black"))

#2H
VlnPlot(AGM,features = c("MMP2","MMP9","MMP14","MMP16","MMP28") ,group.by = "fg", pt.size = 0, stack = T, flip = T, idents = c("HE","HSC", "Endo")) + scale_x_discrete(limits = c("Endo","HE","HSC")) + NoLegend() + theme(axis.text.y = element_text(size=20, color="black"))


#Figure 3

#3A
DimPlot(AL)
#3B
redoLiver = c("SELP","ITGA2B","HBZ","CD38","CD14","C1QA","RNASE2","LYZ","MECOM","MYB","CXCR4","GJA5","CDH5","ESAM","KDR","DCN","BGN","COL1A1", "CCBE1", "COLEC10","ACTA2","TAGLN","MYL9","DES","ALB","TF","TTR",
              "HNF4A")
DotPlot(CAL, features = redoLiver, cluster.idents = T, group.by = "fg")+ RotatedAxis()+ coord_flip()  + NoLegend()+ scale_y_discrete(limits = c("Hept", "Afib", "Stel", "Endo","HE","HSC","Gran","Macro","Mono","GMP","Ery","Mk","Mk Liv"))

#3C, #3F, 3I 
netVisual_heatmap(cellchatCAL, signaling = "PROC")
netVisual_heatmap(cellchatCAL, signaling = "PARs")
netVisual_heatmap(cellchatCALC, signaling = "ANGPTL")

#3D,#3G,3J
netVisual_bubble(cellchatCAL, signaling = "PROC", remove.isolate = T)
netVisual_bubble(cellchatCAL, signaling = "PARs", remove.isolate = T)
netVisual_bubble(cellchatCAL, signaling = "ANGPTL", remove.isolate = T)

#3E, 3H, 3K 
VlnPlot(AGM, "PROCR", idents = c("HE","HSC","Endo"), pt.size = 0) + scale_x_discrete(limits = c("Endo","HE","HSC")) + NoLegend() + theme(axis.text.y = element_text(size=20, color="black"))
VlnPlot(AGM, c("F2R","PARD3"), idents = c("HE","HSC","Endo"), pt.size = 0) + scale_x_discrete(limits = c("Endo","HE","HSC")) + NoLegend() + theme(axis.text.y = element_text(size=20, color="black"))
VlnPlot(AGM, c("CDH5","ITGA5","CLDN5"), idents = c("HE","HSC","Endo"), pt.size = 0) + scale_x_discrete(limits = c("Endo","HE","HSC")) + NoLegend() + theme(axis.text.y = element_text(size=20, color="black"))





