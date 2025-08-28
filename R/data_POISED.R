# Dependencies
library(readr)
library(dplyr)
library(cyCombine)

dir <- "data/POISED"

if (!dir.exists(dir)) {
  # Get data from FlowRepository
  library(FlowRepositoryR)
  x <- flowRep.get("FR-FCM-Z2V9")
  download(x, dirpath = dir)
}
if (!dir.exists(file.path(dir, "LiveLabelledCells"))) {
  # Unzip processed files
  system(paste("unzip", file.path(dir, "ProcessedData.zip"), "-d", dir))
}

# Load data
data_dir <- file.path(dir, "LiveLabelledCells")

files <- list.files(data_dir, pattern = "csv", full.names = TRUE)

df <- readr::read_csv(files, id = "sample", col_select = -1, show_col_types = FALSE)

markers <- colnames(df)[!colnames(df) %in% c("labels", "condition", "sample")]

df_mod <- df |>
  dplyr::mutate(condition = sample |> stringr::str_extract("(Un|Pea)Stim"),
                sample = sample |> stringr::str_extract("P\\d+")) |>
  dplyr::rename(celltype = labels) |>
  cyCombine::transform_asinh(markers = markers)


# Select reference sample
sample_counts <- table(df_mod$sample)
sample_max <- names(sample_counts[sample_counts == max(sample_counts)])

# Split data into reference and query
reference <- dplyr::filter(df_mod, sample == sample_max) |>
  tidyr::unite(col = "sample", sample, condition)
query <- dplyr::filter(df_mod, sample != sample_max) |>
  tidyr::unite(col = "sample", sample, condition)

# Filter populations <10 cells
small_pop <- table(reference$celltype)
small_pop <- small_pop[small_pop < 10]

if (length(small_pop) > 0) {
  message("Removing celltypes: ", paste(names(small_pop)), " from the reference")
  reference <- reference |>
    dplyr::filter(!celltype %in% names(small_pop))
  query <- query |>
    mutate(celltype = case_when(
      celltype %in% names(small_pop) ~ "Unknown",
      TRUE ~ celltype
    ))
}

# Store output
saveRDS(list(
  "reference" = reference,
  "query" = query,
  "markers" = markers,
  unassigned_name = "Unknown"
), file = file.path(dir, "preprocessed.rds"))

