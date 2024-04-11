
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cyDefine

`cyDefine` robustly identifies canonical cell types as well as
unassigned (and potentially novel) cells in various single-cell
cytometry datasets.

## Article

A preprint of `cyDefine` is available
[here](https://www.youtube.com/watch?v=dQw4w9WgXcQ&t=43s).

Please cite with `citation("cyDefine")`

## Install from GitHub

``` r
# To ensure Rstudio looks up BioConductor packages run:
setRepositories(ind = c(1:6, 8))
# Install cyCombine dependency
remotes::install_github("biosurf/cyCombine")

# Install cyDefine
remotes::install_github("claradrachmann/cyDefine", build_vignettes = TRUE)
```

## Documentation

To view documentation for the version of this package installed in your
system, start R and enter:

``` r
browseVignettes("cyDefine")
```

## Usage

Please see the vignettes for a detailed description of usage.

Here is a quick run-through of the functions:
