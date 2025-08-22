#!/bin/bash

datasets=("Levine13" "Levine32" "Samusik" "POISED")
echo "Running cyDefine"
for dataset in "${datasets[@]}"; do
    echo "Processing dataset: $dataset"
    echo "With unassigned"
    Rscript R/method_cyDefine.R \
        --name "$dataset" \
        --unassigned true \
        --seed 145 \
        --threads 39
    echo "Without unassigned"
    Rscript R/method_cyDefine.R \
        --name "$dataset" \
        --unassigned false \
        --seed 145 \
        --threads 39
done
