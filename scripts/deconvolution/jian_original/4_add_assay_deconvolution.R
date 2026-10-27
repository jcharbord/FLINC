## script adding deconvolution result of celltype fraction as an assay to ST seurat object
setwd("/Users/jianwu/Documents/20250505_liver/FLINC/code/deconvolution")

library(Seurat)

## load seurat object
object <- readRDS("../../data/st_analysis_pipeline/object.ffpe.rds")

## load deconvolution results
data <- read.table("../../data/deconvolution/single_nucleus_spot_celltype_frac.txt")

## add as assay
object@assays[["singlecell"]] <- CreateAssayObject(counts = as.sparse(t(data)))
object@assays[["singlecell"]]@key <- "singlecell_"

## save updated seurat object
saveRDS(object, file = "../../data/st_analysis_pipeline/object.ffpe.rds")
