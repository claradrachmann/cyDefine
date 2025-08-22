# Dependencies
library(readr)
library(dplyr)
library(HDCytoData)
library(caret)
library(cyCombine)


# Get data from HDCytoData
system("mkdir -p data/Levine13")
levine13 <- HDCytoData::Levine_13dim_SE()
markers <- colnames(levine13)

# Object to data.frame
df <- SummarizedExperiment::assay(levine13, "exprs") |>
  as.data.frame() |>
  dplyr::mutate(celltype = as.character(levine13@elementMetadata$population_id)) |>
  transform_asinh(markers = markers)

# Set seed for reproducibility
set.seed(123)

# Create an index for stratified splitting
train_index <- createDataPartition(df$celltype, p = 0.5, list = FALSE)

# Split the data into reference and query
reference <- df[train_index, ]
query  <- df[-train_index, ]


# Store output
saveRDS(list(
  "reference" = reference,
  "query" = query,
  "markers" = markers,
  unassigned_name = "unassigned"
), file = "data/Levine13/preprocessed.rds")
