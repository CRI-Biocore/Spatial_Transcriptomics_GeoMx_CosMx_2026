## GeoMx R package installer
## Run on login node. (Within conda_envs folder)
##   conda activate ./geomx_env
##   Rscript geomx_Rpackages_install.R

options(repos = c(CRAN = "https://cloud.r-project.org"))

## ── Step 1: GeomxTools and its Bioconductor dependencies ─────────────────────
## Note: XVector, GenomeInfoDb, Biostrings installed via conda (geomx_env.yaml)
message("Installing GeomxTools chain...")
## Note: ComplexHeatmap, SpatialExperiment, magick, png installed via conda (geomx_env.yaml)
BiocManager::install(c(
  "NanoStringNCTools",
  "GeomxTools",
  "BiocGenerics",
  "SummarizedExperiment",
  "SingleCellExperiment",
  "SpatialDecon",
  "limma",
  "edgeR",
  "PCAtools"
), update = FALSE, ask = FALSE)

## ── Step 2: CRAN packages ─────────────────────────────────────────────────────
message("Installing CRAN packages...")
## Note: png, igraph, data.table, reticulate, networkD3, umap installed via conda (geomx_env.yaml)
install.packages(c(
  "tidyverse",
  "readxl",
  "ggrepel",
  "ggalluvial",
  "patchwork",
  "cowplot",
  "circlize",
  "future",
  "future.apply"
))

message("All packages installed.")