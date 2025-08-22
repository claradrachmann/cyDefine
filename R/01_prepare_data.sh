#!/bin/bash

datasets=("Levine13" "Levine32" "Samusik" "POISED")
echo "Running cyDefine"
for dataset in "${datasets[@]}"; do
    echo "Preprocessing dataset: $dataset"
    ## Data
    Rscript "R/data_${dataset}.R"
done
