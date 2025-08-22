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

exit

# Activate CyAnno environment
mamba activate cyanno

# POISED
## Data
Rscript R/data_POISED.R

## cyDefine
Rscript R/method_cyDefine.R \
  --name "POISED" \
  --train_unassigned true \
  --mtry 10 \
  --seed 145 \
  --threads 39
Rscript R/method_cyDefine.R \
  --name "POISED" \
  --train_unassigned false \
  --mtry 10 \
  --seed 145 \
  --threads 38

## F1
Rscript R/metric_f1.R \
  --dir "data/POISED" \
  --name "POISED" \
  --method "cyDefine" \
  --suffix "_wo_unassigned" \
  --true_col "celltype" \
  --predicted_col "model_prediction"
  # --model_prediction "predicted_celltype"

Rscript R/metric_f1.R \
  --dir "data/POISED" \
  --name "POISED" \
  --method "cyDefine" \
  --suffix "_wo_unassigned" \
  --true_col "celltype" \
  --predicted_col "predicted_celltype"
  # --model_prediction "predicted_celltype"

Rscript R/metric_f1.R \
  --dir "data/POISED" \
  --name "POISED" \
  --method "cyDefine" \
  --suffix "_w_unassigned" \
  --true_col "celltype" \
  --predicted_col "predicted_celltype"

## CLC
Rscript R/method_CLC.R \
  --name "POISED" \
  --unassigned true

## F1
Rscript R/metric_f1.R \
  --dir "data/POISED_Z2V9" \
  --name "POISED" \
  --method "CLC" \
  --true_col "celltype" \
  --predicted_col "predicted_celltype"

## CyAnno
Rscript R/method_CyAnno.R \
  --dir "data/POISED_Z2V9" \
  --name "POISED"

python CyAnno/CyAnno_general.py \
  --dir "data/POISED_Z2V9/CyAnno_data"

# Levine13
## Data
Rscript R/data_Levine13.R

## cyDefine
Rscript R/method_cyDefine.R \
  --dir "data/Levine13" \
  --name "Levine13" \
  --train_unassigned true \
  --mtry 4 \
  --seed 145
Rscript R/method_cyDefine.R \
  --dir "data/Levine13" \
  --name "Levine13" \
  --train_unassigned false \
  --mtry 4 \
  --seed 145

## F1
Rscript R/metric_f1.R \
  --dir "data/Levine13" \
  --name "Levine13" \
  --method "cyDefine" \
  --suffix "_wo_unassigned" \
  --true_col "celltype" \
  --predicted_col "predicted_celltype"
  # --predicted_celltype "model_prediction"

Rscript R/metric_f1.R \
  --dir "data/Levine13" \
  --name "Levine13" \
  --method "cyDefine" \
  --suffix "_w_unassigned" \
  --true_col "celltype" \
  --predicted_col "predicted_celltype"


## CLC
Rscript R/method_CLC.R \
  --dir "data/Levine13" \
  --name "Levine13"

## F1
Rscript R/metric_f1.R \
  --dir "data/Levine13" \
  --name "Levine13" \
  --method "CLC" \
  --true_col "celltype" \
  --predicted_col "predicted_celltype"

## CyAnno
Rscript R/method_CyAnno.R \
  --dir "data/Levine13" \
  --name "Levine13"
mamba activate cyanno
python CyAnno/CyAnno_general.py \
  --dir "data/Levine13/CyAnno_data"

# Levine32
## Data
Rscript R/data_Levine32.R

## cyDefine
Rscript R/method_cyDefine.R \
  --dir "data/Levine32" \
  --name "Levine32" \
  --train_unassigned true \
  --mtry 10 \
  --seed 145
Rscript R/method_cyDefine.R \
  --dir "data/Levine32" \
  --name "Levine32" \
  --train_unassigned false \
  --mtry 10 \
  --seed 145

## F1
Rscript R/metric_f1.R \
  --name "Levine32" \
  --method "cyDefine" \
  --suffix "_wo_unassigned" \
  --true_col "celltype" \
  --predicted_col "predicted_celltype"
  # --predicted_celltype "model_prediction"

Rscript R/metric_f1.R \
  --dir "data/Levine32" \
  --name "Levine32" \
  --method "cyDefine" \
  --suffix "_w_unassigned" \
  --true_col "celltype" \
  --predicted_col "predicted_celltype"

## CLC
Rscript R/method_CLC.R \
  --dir "data/Levine32" \
  --name "Levine32"

## F1
Rscript R/metric_f1.R \
  --dir "data/Levine32" \
  --name "Levine32" \
  --method "CLC" \
  --true_col "celltype" \
  --predicted_col "predicted_celltype"

## CyAnno
Rscript R/method_CyAnno.R \
  --dir "data/Levine32" \
  --name "Levine32"

python CyAnno/CyAnno_general.py \
  --dir "data/Levine32/CyAnno_data"

# Samusik
## Data
Rscript R/data_Samusik.R

## cyDefine
Rscript R/method_cyDefine.R \
  --dir "data/Samusik" \
  --name "Samusik" \
  --train_unassigned true \
  --mtry 10 \
  --seed 145
Rscript R/method_cyDefine.R \
  --dir "data/Samusik" \
  --name "Samusik" \
  --train_unassigned false \
  --mtry 10 \
  --seed 145

## F1
Rscript R/metric_f1.R \
  --dir "data/Samusik" \
  --name "Samusik" \
  --method "cyDefine" \
  --suffix "_wo_unassigned" \
  --true_col "celltype" \
  --predicted_col "predicted_celltype"
  # --predicted_celltype "model_prediction"

Rscript R/metric_f1.R \
  --dir "data/Samusik" \
  --name "Samusik" \
  --method "cyDefine" \
  --suffix "_w_unassigned" \
  --true_col "celltype" \
  --predicted_col "predicted_celltype"

## CLC
Rscript R/method_CLC.R \
  --dir "data/Samusik" \
  --name "Samusik"

## F1
Rscript R/metric_f1.R \
  --dir "data/Samusik" \
  --name "Samusik" \
  --method "CLC" \
  --true_col "celltype" \
  --predicted_col "predicted_celltype"

# CyAnno
## Activate CyAnno environment
mamba activate cyanno

data="POISED"
Rscript R/method_CyAnno.R \
  --name ${data}
SECONDS=0
python CyAnno/CyAnno_general.py \
  --dir "data/${data}/CyAnno_data"
echo "Runtime (seconds): $SECONDS" > data/${data}/${data}_CyAnno_w_unassigned_runtime.txt

Rscript R/method_CyAnno.R \
  --name "Levine13"
python CyAnno/CyAnno_general.py \
  --dir "data/Levine13/CyAnno_data"

## F1
Rscript R/metric_f1.R \
  --dir "data/Samusik" \
  --name "Samusik" \
  --method "CyAnno" \
  --true_col "celltype" \
  --predicted_col "predicted_celltype"

## Spectre
Rscript R/method_Spectre.R \
  --name "Samusik"
## Spectre
Rscript R/method_Spectre.R \
  --name "Levine13"
## Spectre
Rscript R/method_Spectre.R \
  --name "Levine32"
## Spectre
Rscript R/method_Spectre.R \
  --name "POISED"

## F1
Rscript R/metric_f1.R \
  --name "Samusik" \
  --method "Spectre"
Rscript R/metric_f1.R \
  --name "Levine13" \
  --method "Spectre"
