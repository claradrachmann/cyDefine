library(argparse)
library(dplyr)
library(lubridate)

# Define argument parser
parser <- ArgumentParser(description="Run CyAnno")

# Add arguments
# parser$add_argument("--output_dir", "-o", dest="output_dir", type="character", help="output directory where files will be saved")
# parser$add_argument("--dir", "-d", dest="dir", type="character", help="data directory")
parser$add_argument("--name", "-n", dest="name", type="character", help="name of the dataset")
parser$add_argument("--unassigned", dest="unassigned", type="logical", help="include unassigned", default=FALSE)
# Parse command-line arguments
opt <- parser$parse_args()

# output_dir <- opt$output_dir
name <- opt$name
dir <- file.path("data", name)
input <- readRDS(file.path(dir, "preprocessed.rds"))
unassigned <- opt$unassigned
out <- ifelse(unassigned, "_w_unassigned", "_wo_unassigned")

# if (unassigned) input$reference$celltype[input$reference$celltype == input$unassigned_name] <- "Unknown"
# if (!unassigned) input$reference <- input$reference[input$reference$celltype != input$unassigned_name, ]
input$query$celltype[input$query$celltype == input$unassigned_name] <- "Unknown"

input$query$id <- seq_len(nrow(input$query))

# Select newest results
output_dir <- list.dirs(dir, recursive = F)
output_dir <- output_dir[stringr::str_detect(output_dir, paste0("CyAnno_results", out))]
dates <- stringr::str_remove_all(output_dir, "^.*_unassigned_") |>
  dmy_hm()
output_dir <- output_dir[which(max(dates) == dates)]

result_files <- list.files(output_dir, pattern = "csv$", full.names = TRUE)
multiple_files <- length(result_files) > 1
for (file in result_files) {
  result <- read.csv(file)
  result_sample <- basename(file) |>
    stringr::str_remove_all("^sample_") |>
    stringr::str_remove_all("_labelled_expr.csv")


  if (multiple_files) {
    input$query[input$query$sample == result_sample, "predicted_celltype"] <- result$labels
  } else {
    input$query$predicted_celltype <- result$labels
  }
}
input$query <- input$query |>
  dplyr::arrange(id)

output <- data.frame(
  "celltype" = input$query$celltype,
  "predicted_celltype" = input$query$predicted_celltype
)

# Store results
suffix <- ifelse(unassigned, "_w_unassigned", "_wo_unassigned")
saveRDS(output, paste0(dir, "/", name, "_CyAnno", suffix, ".rds"))
