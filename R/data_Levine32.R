# Dependencies
library(readr)
library(dplyr)
library(HDCytoData)
library(cyCombine)


# Get data from HDCytoData
system("mkdir -p data/Levine32")
levine32 <- HDCytoData::Levine_32dim_SE()

markers <- c("CD45RA", "CD133", "CD19", "CD22", "CD11b", "CD4", "CD8", "CD34",
             "Flt3", "CD20", "CXCR4", "CD235ab", "CD45", "CD123", "CD321",
             "CD14", "CD33", "CD47", "CD11c", "CD7", "CD15", "CD16", "CD44",
             "CD38", "CD13", "CD3", "CD61", "CD117", "CD49d", "HLA-DR",
             "CD64", "CD41")

# Object to data.frame
df <- SummarizedExperiment::assay(levine32, "exprs") |>
  as.data.frame() |>
  dplyr::mutate(celltype = as.character(levine32@elementMetadata$population_id),
                sample = levine32@elementMetadata$patient_id) |>
  transform_asinh(markers = markers)

# Split the data into reference and query
reference <- dplyr::filter(df, sample == "H1")
query  <- dplyr::filter(df, sample == "H2")


# Store output
saveRDS(list(
  "reference" = reference,
  "query" = query,
  "markers" = markers,
  unassigned_name = "unassigned"
), file = "data/Levine32/preprocessed.rds")
