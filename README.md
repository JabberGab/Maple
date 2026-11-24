
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Maple

<!-- badges: start -->

[![R-CMD-check](https://github.com/JabberGab/Maple/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JabberGab/Maple/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Installation

``` r
remotes::install_github("JabberGab/Maple")
```

## Quick start

Maple is a framework for stability calculation and multi-resolution
clustering for scRNA-seq data. The algorithm executes multiple
partitioning through Leiden clustering at different resolution values,
builds a stability tree and finds the optimal frontier for every
lineage. The workflow is divided in three main steps:

First Maple() function is called to obtain the observation-level
stability scores for the chosen resolution grid.

``` r
library(Maple)
sce <- Maple(sce, dimred = "PCA", dims = 30, grid = seq(0.1, 2, by = 0.1))
plotDiagnostics(sce)
```

The stability tree is then built using grow_Maple().

``` r
sce <- grow_Maple(sce)
plotTree(sce)
```

Finally, the tree is pruned to obtain a stable frontier of clusters,
maximizing the resolution for every lineage.

``` r
sce <- cut_Maple(sce)
plotPruning(sce)
```
