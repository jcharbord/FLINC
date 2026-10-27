# Source history

This folder contains historical source files retained for traceability.

They document how the clustering, annotation and downstream characterization
analysis was prepared from the earlier FLINC workflow.

## Files

### rebuild_clustering_source_from_historical_analysis.R

Technical preparation script used while adapting the earlier analytical source.

### clustering_annotation_characterization_source.Rmd

Human-readable R Markdown source containing clustering, subclustering,
annotation and downstream biological characterization.

### clustering_annotation_characterization_extracted.R

Plain R code extracted from the R Markdown source using `knitr::purl()`.

The reader-facing analysis script:

`../analysis/02_cluster_annotate_and_characterize_spatial_domains.R`

uses the R Markdown source in this folder internally, extracts the R code and
executes it directly.

Therefore, users normally run the script in `scripts/analysis/`; they do not
need to manually run the files in this folder.

These files are preserved to maintain technical traceability between the
earlier analysis workflow and the validated reference analysis.
