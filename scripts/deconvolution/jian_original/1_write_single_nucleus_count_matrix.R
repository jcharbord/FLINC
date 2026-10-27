## script extracting count matrix from seurat object of single-nucleus, and writing in .h5 file
setwd("/Users/jianwu/Documents/20250505_liver/FLINC/code/deconvolution")

library(Seurat)
library(SeuratObject)

## load single-nucleus seurat object and extract only single-nucleus "cells" and features intersecting with ST seurat objects
object <- readRDS(file = "../../data/deconvolution/Integrated_with_Subannotations.rds")
object@images <- list()

common.feature <-  read.table(file = '../../data/deconvolution/feature.common.spatial.single.nucleus.txt', header = FALSE)
index.common.feature <- which((rownames(object) %in% common.feature$V1))

index.single_nuc <- which(object@meta.data[["assay_type"]]%in%"single_nuc")
file.name <- "single_nucleus"

object <- object[index.common.feature,index.single_nuc]

## write out cell-type info file, supervising downstream deconvolution
barcodes <- colnames(object)
label <- as.vector(object@meta.data[["Manual_Annotation"]])
object.label <- data.frame(label, row.names = barcodes, fix.empty.names = TRUE)
write.table(label, file = paste(".../../data/deconvolution/Human_liver_", file.name, "_celltype.txt", sep = ""), append = FALSE, sep = "\t", col.names = FALSE, row.names = barcodes, quote = FALSE)

## write .h5 file of single-nucleus
library(rhdf5)

write_h5 <- function(object, file.name){
  
  h5createFile(paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"))
  
  h5createGroup(paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"),"matrix")
  
  barcodes <- colnames(object)
  
  h5write(barcodes, paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"),"matrix/barcodes")
  
  data <- object@assays[["RNA"]]@data@x
  
  h5write(data, paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"),"matrix/data")
  
  indices <- object@assays[["RNA"]]@data@i
  
  h5write(indices, paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"),"matrix/indices")
  
  shape <- object@assays[["RNA"]]@data@Dim
  
  h5write(shape, paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"),"matrix/shape")
  
  indptr <- object@assays[["RNA"]]@data@p
  
  h5write(indptr, paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"),"matrix/indptr")
  
  h5createGroup(paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"),"matrix/features")
  
  id <- object@assays[["RNA"]]@data@Dimnames[[1]]
  
  h5write(id, paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"),"matrix/features/id")
  
  h5write(id, paste(".../../data/deconvolution/Human_liver", file.name ,"gene_count.h5", sep = "_"),"matrix/features/name")
  
}

write_h5(object, file.name)

