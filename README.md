
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cyDefine

`cyDefine` robustly identifies canonical cell types as well as
unassigned (and potentially novel) cells in various single-cell
cytometry datasets.

## Article

The manuscript for `cyDefine` is in preparation.

Please cite with `citation("cyDefine")`.

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

To view the documentation for using `cyDefine`, start R and enter:

``` r
browseVignettes("cyDefine")
```

<!-- ## Usage -->
<!-- Please see the vignettes for a detailed description of usage. -->
<!-- Here is a quick run-through of the main functionalities. -->

## Report issues

If you have any issues or questions regarding the use of `cyDefine`,
please do not hesitate to raise an issue on GitHub. In this way, others
may also benefit from the answers and discussions.
