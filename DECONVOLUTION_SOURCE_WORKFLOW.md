# Deconvolution source workflow

## Scope

This document records the provenance and validation of the spot-level cell-type
composition matrix used by the FLINC spatial-transcriptomics workflow. It is a
reproducibility input to the principal analysis: the principal workflow imports
the validated composition matrix; it does not silently rerun deconvolution.

## Canonical composition matrix

`/path/to/single_nucleus_spot_celltype_frac.txt` (external local input)

The matrix contains 26,056 Visium spots and 12 cell-type fraction columns. Each
row is one spatial spot; its fractions are relative cell-type compositions.
The canonical SHA-256 is
`aa2e9aef96fb51afeb5488a61ff374a03eb6ea9a96b724e96f04ce3cdd24815e`.

It is represented in the analysis Seurat object as the `singlecell` assay (key
`singlecell_`). These values are compositions/fractions, not sequencing counts,
even though `CreateAssayObject()` is used as the storage container.

## Attribution

The original STRIDE deconvolution workflow and its source scripts were provided
by Jian Wu. This public release documents the source workflow and adds
input-integrity validation for reproducible downstream use.

## Source-code provenance

The upstream source workflow by Jian Wu is available within this package at
`scripts/deconvolution/jian_original/`. The packaged copy was retrieved from
the original repository's `main` branch, `Deconvolution/`, without changing its
contents. The old `code/deconvolution/` folder is not required in this checkout. See
[`scripts/deconvolution/README.md`](scripts/deconvolution/README.md) for the
original-path and expression-slot checks needed before a fresh rerun:

| Step | Source file | Operation |
| --- | --- | --- |
| 1 | `1_write_single_nucleus_count_matrix.R` | Creates the single-nucleus reference count matrix and cell-type labels. |
| 2 | `2_write_ST_count_matrix.R` | Creates the Visium spatial count matrix. |
| 3 | `3_deconvolution_stride.sh` | Runs `STRIDE deconvolve` using those reference and spatial matrices. |
| 4 | `4_add_assay_deconvolution.R` | Imports the resulting fractions into the Seurat `singlecell` assay. |

The single-nucleus reference derives from the public human-liver single-nucleus
dataset reported by Alonso-Peña *et al.* (GEO: GSE185477), as recorded by the
source scripts.

## Reproducibility boundary

The released workflow validates and imports the recorded STRIDE output. A fresh
STRIDE calculation is not claimed by this release because the original STRIDE
runtime and reference-object environment were not frozen as a portable software
environment. The exact upstream export, STRIDE-command, and import scripts are
nevertheless included so that a full rerun can be performed when those inputs
and environment are available.

If spatial QC changes the retained spot set, STRIDE must be rerun using the
matching spatial count matrix; reusing this composition matrix after a changed
spot set would not be valid.

## Validation requirements

Before import or downstream analysis, confirm:

1. SHA-256 matches the manifest;
2. the table has 26,056 spot rows and 12 numeric fraction columns;
3. spot identifiers are unique and exactly align with the Seurat object;
4. all values are finite and lie in [0, 1] within numerical tolerance; and
5. each spot's fractions sum to one within tolerance (`1e-6`).

Run the validator from the repository root, supplying the actual local matrix path:

```bash
Rscript scripts/deconvolution/01_validate_deconvolution_input.R \
  --fractions /path/to/single_nucleus_spot_celltype_frac.txt \
  --output-dir reproducibility_records/Deconvolution
shasum -a 256 -c DECONVOLUTION_SHA256.txt
```

The validator is intentionally independent of Seurat so it can check the
published input before an R/Seurat environment is installed.

The source manifest checks files in this repository. The external fraction
matrix must be checked separately against the canonical hash above; it is not
included in the repository.
