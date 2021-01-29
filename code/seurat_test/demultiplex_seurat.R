#!/usr/bin/env Rscript

#USE: 
# Rscript --vanilla demultiplex_seurat.R \
# "path/to/filtered_feature_bc_matrix_RNA" \
# "path/to/filtered_feature_bc_matrix_HT" \
# "path/to/output_folder"

library(Seurat)

#--- args from comandline ----------

#args[1] = path/to/filtered_feature_bc_matrix_RNA
#args[2] = path/to/filtered_feature_bc_matrix_HT 
#args[3] = path/to/output_folder

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("provide paths to RNA and HTO feature matrix", call.=FALSE)
}else if(length(args)==2){
  args[3] = getwd()
}

#--- DATA ----------------------

#read umi matrix 
rna.data = Read10X(data.dir = args[1])

#read hashtag matrix 
hto.data = Read10X(data.dir = args[2])

#--- CELL BARCODES -------------

#select cell barcodes found in both gex and ht matrix  
joint = intersect(colnames(rna.data), colnames(hto.data))

#subset on joint cell barcodes
rna.data = rna.data[, joint]
hto.data = hto.data[, joint]

#--- SEURAT OBJECT -----------

#create seurat object with RNA assay 
data = CreateSeuratObject(counts = rna.data)

#add hashtag data as a new assay 
data[["HTO"]] = CreateAssayObject(counts = hto.data)

#normalize HTO data (centered log-ratio (CLR) transformation) 
data = NormalizeData(data, assay = "HTO", normalization.method = "CLR")

#--- DEMULTIPLEX --------------

#samples can be found in data@meta.data

data = HTODemux(data, assay = "HTO", positive.quantile = 0.99)


