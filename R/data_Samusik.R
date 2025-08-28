# Dependencies
library(readr)
library(dplyr)
library(HDCytoData)
library(cyCombine)


# Get data from HDCytoData
system("mkdir -p data/Samusik")
Samusik <- HDCytoData::Samusik_all_SE()

markers <- c("CD45.2", "IgD", "CD11c", "CD3",
             "CD23", "CD34", "CD115", "CD19","CD8",
             "CD4", "CD11b", "CD27", "CD16_32", "SiglecF",
             "Foxp3", "CD5", "FceR1a", "TCRgd", "CCR7", "Sca1",
             "CD49b", "cKit", "CD150", "CD25", "TCRb", "CD43", "CD64",
             "CD138", "CD103", "IgM", "CD44", "MHCII")

# Object to data.frame
df <- SummarizedExperiment::assay(Samusik, "exprs") |>
  as.data.frame() |>
  dplyr::mutate(celltype = as.character(Samusik@elementMetadata$population_id),
                sample = Samusik@elementMetadata$sample_id) |>
  transform_asinh(markers = markers)

# Select reference sample
sample_counts <- table(df$sample)
sample_max <- names(sample_counts[sample_counts == max(sample_counts)])

# Split data into reference and query
reference <- dplyr::filter(df, sample == sample_max)
query <- dplyr::filter(df, sample != sample_max)


# Filter populations <10 cells
small_pop <- table(reference$celltype)
small_pop <- small_pop[small_pop < 10]

if (length(small_pop) > 0) {
  message("Removing celltypes: ", paste(names(small_pop)), " from the reference")
  reference <- reference |>
    dplyr::filter(!celltype %in% names(small_pop))
  query <- query |>
    mutate(celltype = case_when(
      celltype %in% names(small_pop) ~ "unassigned",
      TRUE ~ celltype
    ))
}

# Store output
saveRDS(list(
  "reference" = reference,
  "query" = query,
  "markers" = markers,
  unassigned_name = "unassigned"
), file = "data/Samusik/preprocessed.rds")

