#!/bin/bash


# Prepare datasets
bash src/01_prepare_data.sh

# Run cyDefine
bash src/02_run_cyDefine.sh
# Run CLC
bash src/02_run_CLC.sh
# Run Spectre
bash src/02_run_Spectre.sh
# Run CyAnno
mamba activate cyanno
cp CyAnno_changes/* ../CyAnno/
bash src/02_run_CyAnno.sh

# Evaluate
bash src/03_run_F1.sh

# Compile results
Rscript R/compile_results.R

# Create Figures
Rscript R/Figures.R

