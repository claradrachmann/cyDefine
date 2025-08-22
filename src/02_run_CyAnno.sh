#!/bin/bash

## Activate CyAnno environment
mamba activate cyanno

datasets=("Levine13" "Levine32" "Samusik" "POISED")
echo "Running CyAnno"
for dataset in "${datasets[@]}"; do
    echo "Processing dataset: $dataset"
    echo "with unassigned"
    Rscript R/method_CyAnno.R \
      --name ${dataset} \
      --unassigned true \
      --force true
    SECONDS=0
    python ../CyAnno/CyAnno_general.py \
      --dir "data/${dataset}/CyAnno_data" \
      --unassigned true
    echo "Runtime (seconds): $SECONDS" > data/${dataset}/${dataset}_CyAnno_w_unassigned_runtime.txt
    echo "without unassigned"
    Rscript R/method_CyAnno.R \
      --name ${dataset} \
      --unassigned false \
      --force true
    SECONDS=0
    python ../CyAnno/CyAnno_general.py \
      --dir "data/${dataset}/CyAnno_data" \
      --unassigned false
    echo "Runtime (seconds): $SECONDS" > data/${dataset}/${dataset}_CyAnno_w_unassigned_runtime.txt
done
