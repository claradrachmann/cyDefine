# Dependencies
library(readr)
library(dplyr)
library(cyCombine)

if (!dir.exists("data/MultiCenter")) {
  # Get data from FlowRepository
  library(FlowRepositoryR)
  x <- flowRep.get("FR-FCM-ZYTT")
  download(x, dirpath = "data/MultiCenter")
}

# Load data
data_dir <- "data/MultiCenter"
celltype_col <- "labels"

files <- list.files(data_dir, pattern = "^Multi.*csv$", full.names = TRUE)

df <- readr::read_csv(files, id = "sample", col_select = -1)

markers <- colnames(df)[!colnames(df) %in% c("labels", "sample")]

df <- df |>
  dplyr::mutate(condition = sample |> stringr::str_extract("(Un|Pea)Stim"),
                sample = sample |> stringr::str_extract("P\\d*_(Un|Pea)Stim")) |>
  dplyr::rename(celltype = labels) |>
  cyCombine::transform_asinh(markers = markers)

# Select reference sample
sample_counts <- table(df$sample)
sample_max <- names(sample_counts[sample_counts == max(sample_counts)])

# Split data into reference and query
reference <- dplyr::filter(df, sample == sample_max)
query <- dplyr::filter(df, sample != sample_max)


# Store output
saveRDS(list(
  "reference" = reference,
  "query" = query,
  "markers" = markers,
  unassigned_name = "Unknown"
), file = "data/MultiCenter_ZYTT/preprocessed.rds")
