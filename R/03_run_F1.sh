#!/bin/bash

methods=("cyDefine" "CLC" "Spectre" "CyAnno")
datasets=("Levine13" "Levine32" "Samusik" "POISED")

for method in "${methods[@]}"; do
    echo "Computing F1 $method"
    for dataset in "${datasets[@]}"; do
        echo "Processing dataset: $dataset"
        echo "With unassigned"
        Rscript R/metric_f1.R \
            --name "$dataset" \
            --method "$method" \
            --suffix "_w_unassigned"
        echo "Without unassigned"
        Rscript R/metric_f1.R \
            --name "$dataset" \
            --method "$method" \
            --suffix "_wo_unassigned"
    done
done
