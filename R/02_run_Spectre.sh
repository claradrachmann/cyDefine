#!/bin/bash

datasets=("Levine13" "Levine32" "Samusik" "POISED")
echo "Running Spectre"
for dataset in "${datasets[@]}"; do
    echo "Processing dataset: $dataset"
    echo "With unassigned"
    Rscript R/method_Spectre.R \
        --name "$dataset" \
        --unassigned true
    echo "Without unassigned"
    Rscript R/method_Spectre.R \
        --name "$dataset" \
        --unassigned false
done
