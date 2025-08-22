#!/bin/bash

datasets=("Levine13" "Levine32" "Samusik" "POISED")
echo "Running CyTOF Linear Classifier"
for dataset in "${datasets[@]}"; do
    echo "Processing dataset: $dataset"
    echo "With unassigned"
    Rscript R/method_CLC.R \
        --name "$dataset" \
        --unassigned true
    echo "Without unassigned"
    Rscript R/method_CLC.R \
        --name "$dataset" \
        --unassigned false
done
