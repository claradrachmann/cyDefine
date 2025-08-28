
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cyDefine Benchmark

This repository contains the code to rerun the cyDefine benchmark. The
benchmark compares the phenotype transfer tools `cyDefine`, `CyAnno`,
`Spectre`, and `CyTOF Linear Classifier` on the datasets `Levine13`,
`Levine32`, `Samusik`, and `POISED`.

| Dataset  | Source      |
|----------|-------------|
| Levine13 | HDCytoData  |
| Levine32 | HDCytoData  |
| Samusik  | HDCytoData  |
| POISED   | FR-FCM-Z2V9 |

Download the POISED dataset manually, unzip `ProcessedData.zip`, and put
it in `data/POISED`. The remaining datasets will be downloaded
automatically.

## Recreate figures

All results have been compiled in this repository. Figures can therefore
be recreated with:

``` sh
Rscript R/Figures.R
```

To rerun the entire benchmark, read on.

## Dependencies

The `CyTOF Linear Classifier` functions have been extracted and added in
`R/method_CLC.R` for a simpler implementation.

### Benchmark dependecies

Besides the R packages below, CyAnno runs in a custom-built mamba/conda
environment, so mamba or conda is required as well.

``` r
install.packages(c("dplyr", "ggplot2", "readr", "caret", "argparse", "MASS", "patchwork", "ggpubr"))
BiocManager::install(c("HDCytoData"))
```

``` sh
git clone https://github.com/claradrachmann/cyDefine.git
cd cyDefine
git checkout benchmark
```

### cyDefine dependecies

Install cyDefine and its dependencies

``` r
# To ensure Rstudio looks up BioConductor packages run:
setRepositories(ind = 1:2)
# If you are correcting cytometry data, install the following Bioconductor packages:
BiocManager::install(c("flowCore", "Biobase", "sva"))
# Install cyCombine dependency
remotes::install_github("biosurf/cyCombine", ref = "dev")

# Install cyDefine
remotes::install_github("claradrachmann/cyDefine")
```

### Spectre dependecies

Install Spectre

``` r
remotes::install_github(repo = "immunedynamics/spectre")
```

### CyAnno dependecies

Clone CyAnno

``` sh
git clone -o ../CyAnno https://github.com/abbioinfo/CyAnno.git
# Add benchmark-specific updates to the CyAnno script
cp CyAnno_changes/* ../CyAnno/

# Build environment
mamba env create -f CyAnno_changes/environment.yml
mamba activate cyanno
```

## Run the benchmark

The benchmark is orchestrated by a simple shell script. See its content
to run each module individually.

``` sh
bash src/00_run_all.sh
```
