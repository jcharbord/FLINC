## shell script inplementing deconvolution via STRIDE


conda activate stride


./bin/STRIDE deconvolve --sc-count ../../data/deconvolution/Human_liver_single_nucleus_gene_count.h5 \
--sc-celltype ../../data/deconvolution/Human_liver_single_nucleus_celltype.txt \
--st-count ../../data/deconvolution/Human_liver_ffpe_gene_count.h5 \
--outdir ../../data/deconvolution --outprefix single_nucleus --normalize

