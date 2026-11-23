
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Maple

<!-- badges: start -->

[![R-CMD-check](https://github.com/JabberGab/Maple/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/JabberGab/Maple/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Installation

``` r
remotes::install_github("JabberGab/Maple")
#> Using GitHub PAT from the git credential store.
#> Downloading GitHub repo JabberGab/Maple@HEAD
#> cpp11        (0.5.2    -> 0.5.5   ) [CRAN]
#> RcppArmad... (15.2.2-1 -> 15.8.0-1) [CRAN]
#> Installing 2 packages: cpp11, RcppArmadillo
#> Installing packages into '/private/var/folders/q3/mt9rsrrj07j9f8myj1gb1hz00000gn/T/RtmpjNcQCA/temp_libpath505478f66dc6'
#> (as 'lib' is unspecified)
#> 
#> The downloaded binary packages are in
#>  /var/folders/q3/mt9rsrrj07j9f8myj1gb1hz00000gn/T//RtmpV5n2wY/downloaded_packages
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/private/var/folders/q3/mt9rsrrj07j9f8myj1gb1hz00000gn/T/RtmpV5n2wY/remotes532f2ffb9ea7/JabberGab-Maple-1eaa039/DESCRIPTION’ ... OK
#> * preparing ‘Maple’:
#> * checking DESCRIPTION meta-information ... OK
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * building ‘Maple_0.1.0.tar.gz’
#> Installing package into '/private/var/folders/q3/mt9rsrrj07j9f8myj1gb1hz00000gn/T/RtmpjNcQCA/temp_libpath505478f66dc6'
#> (as 'lib' is unspecified)
```

## Quick start

``` r
library(Maple)
sce <- Maple(sce, dimred = "PCA", dims = 30, grid = seq(0.1, 2, by = 0.1))
sce <- grow_Maple(sce)
sce <- cut_Maple(sce)
```
