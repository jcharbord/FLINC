## script extracting count matrix from seurat object of single-nucleus, and writing in .h5 file
setwd("/Users/jianwu/Documents/20250505_liver/FLINC/code/deconvolution")

library(Seurat)
library(SeuratObject)

## load metadata listing out sample names
metadata <- read.csv("../../data/st_analysis_pipeline/metadata_spatial_ffpe0.2.csv")
sample.list <- metadata$linkID_spatial_ffpe

## path to the directory where count matrices for each sample reside
data_dir <- "../../data/st_analysis_pipeline/"

## iterate through all folders in sample.list to load count matrices and images, and merge all into a single object
object.list <- lapply(sample.list, function(sample_id, data_dir){
  object <- Load10X_Spatial(data.dir = paste(data_dir, sample_id, sep = "/"))
},data_dir)

rm(data_dir, sample.list)

for (index in 1:length(object.list)) {
  object.list[[index]] <- RenameCells(object.list[[index]], add.cell.id = index)
}

object <- merge(object.list[[1]], object.list[[2]])
for(index in 3:length(object.list)){
  object <- merge(object, object.list[[index]])
}

## write .h5 file of ST data (ffpe)
library(rhdf5)

file.name <- "ffpe"

write_h5 <- function(object, file.name){
  
  h5createFile(paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"))
  
  h5createGroup(paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"),"matrix")
  
  barcodes <- colnames(object)
  
  h5write(barcodes, paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"),"matrix/barcodes")
  
  data <- object@assays[["Spatial"]]@data@x
  
  h5write(data, paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"),"matrix/data")
  
  indices <- object@assays[["Spatial"]]@data@i
  
  h5write(indices, paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"),"matrix/indices")
  
  shape <- object@assays[["Spatial"]]@data@Dim
  
  h5write(shape, paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"),"matrix/shape")
  
  indptr <- object@assays[["Spatial"]]@data@p
  
  h5write(indptr, paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"),"matrix/indptr")
  
  h5createGroup(paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"),"matrix/features")
  
  id <- object@assays[["Spatial"]]@data@Dimnames[[1]]
  
  h5write(id, paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"),"matrix/features/id")
  
  h5write(id, paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"),"matrix/features/name")
  
}

write_h5(object, file.name)



