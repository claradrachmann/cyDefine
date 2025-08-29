library(argparse)
library(dplyr)
library(cyCombine)
library(cyDefine)
library(ggplot2)
library(patchwork)

# Define argument parser
parser <- ArgumentParser(description="Run cyDefine with universal reference")

# Add arguments
# parser$add_argument("--output_dir", "-o", dest="output_dir", type="character", help="output directory where files will be saved")
# parser$add_argument("--dir", "-d", dest="dir", type="character", help="data directory")
parser$add_argument("--name", "-n", dest="name", type="character", help="name of the dataset")
# parser$add_argument("--data", dest="data", type="character", help="input file #1")
parser$add_argument("--mtry", dest="mtry", type="numeric", help="mtry", default = NULL)
parser$add_argument("--seed", dest="seed", type="numeric", help="seed", default = 145)
parser$add_argument("--threads", dest="threads", type="numeric", help="threads", default = 4)
# parser$add_argument("--metadata", dest="metadata", type="character", help="input file #2")

# Parse command-line arguments
opt <- parser$parse_args()

# output_dir <- opt$output_dir
name <- opt$name
dir <- file.path("data", name)
input <- readRDS(file.path(dir, "preprocessed.rds"))
mtry <- opt$mtry
seed <- opt$seed
threads <- opt$threads
# unassigned <- opt$unassigned
# if (!unassigned) input$reference <- input$reference[input$reference$celltype != input$unassigned_name, ]

query <- bind_rows(input$reference, input$query)

if (file.exists("data/pbmc_reference.rds")) {
  reference <- readRDS("data/pbmc_reference.rds")
} else {
  reference <- cyDefine::get_reference("pbmc", store = TRUE, path = "data")
}

if (name == "Levine13") {
  map_specific_from <- NULL
  map_specific_to = NULL
} else if (name == "Levine32") {
  map_specific_from <- c("Flt3")
  map_specific_to = c("CD135")
} else if (name == "Samusik") {
  map_specific_from <- c("CD16_32", "cKit", "TCRb", "MHCII", "CD45.2")
  map_specific_to = c("CD16", "CD117", "abTCR", "HLA-DR", "CD45")
} else if (name == "POISED") {
  map_specific_from <- c("integrin")
  map_specific_to = c("Integrin-7")
} else {
  map_specific_from <- NULL
  map_specific_to = NULL
}


query_mapped <- map_marker_names(
  query,
  query_markers = input$markers,
  ref_markers = cyDefine::pbmc_markers,
  using_pbmc = TRUE,
  map_specific_from = map_specific_from,
  map_specific_to = map_specific_to)


if (is.null(mtry)) mtry <- floor(length(input$markers)/3)
# metadata <- opt$metadata

t <- system.time({
message("Classifying ", name, " using cyDefine")

## Run cyDefine
classified <- cyDefine(
  reference = reference,
  query = query_mapped,
  exclude_redundant = TRUE,
  markers = cyCombine::get_markers(query_mapped),
  adapt_reference = TRUE,
  using_pbmc = TRUE,
  batch_correct = TRUE,
  norm_method = "scale",
  xdim = c(1,5),
  ydim = 5,
  identify_unassigned = TRUE,
  train_on_unassigned = FALSE,
  unassigned_name = input$unassigned_name,
  seed = seed,
  num.threads = threads,
  num.trees = 600,
  mtry = mtry,
  pb = TRUE
)
})

classified$query <- dplyr::filter(classified$query, celltype != "unassigned")

p1 <- plot_umap(
  classified$reference,
  classified$query,
  markers = cyCombine::get_markers(query_mapped))
# p1
p2 <- plot_umap(
  classified$query,
  col = "celltype",
  title = paste(name, "-", "Celltype"),
  markers = cyCombine::get_markers(query_mapped)) +
  plot_umap(
    classified$query,
    col = "predicted_celltype",
    title = paste(name, "-", "Predicted celltype"),
    markers = cyCombine::get_markers(query_mapped))
# p2

# Store results
ggsave(paste0("figs/", name, "_universal.png"), p1)
ggsave(paste0("figs/", name, "_true-v-predicted.png"), p2, width = 15)
# saveRDS(output, paste0(dir, "/", name, "_cyDefine", ".rds"))

# cat("Runtime (elapsed): ", t["elapsed"], file = paste0(dir, "/", name, "_cyDefine", suffix, "_runtime.txt"))


# message("Runtime: ", t["elapsed"])
