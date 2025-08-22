library(argparse)
library(dplyr)
source("../CyAnno/helper.R")

# Define argument parser
parser <- ArgumentParser(description="Run CyAnno")

# Add arguments
# parser$add_argument("--output_dir", "-o", dest="output_dir", type="character", help="output directory where files will be saved")
# parser$add_argument("--dir", "-d", dest="dir", type="character", help="data directory")
parser$add_argument("--name", "-n", dest="name", type="character", help="name of the dataset")
parser$add_argument("--force", "-f", dest="force", type="logical", help="rm old folder", default = FALSE)
parser$add_argument("--unassigned", dest="unassigned", type="logical", help="include unassigned", default=FALSE)
# Parse command-line arguments
opt <- parser$parse_args()

# output_dir <- opt$output_dir
name <- opt$name
dir <- file.path("data", name)
input <- readRDS(file.path(dir, "preprocessed.rds"))
unassigned <- opt$unassigned
out <- ifelse(unassigned, "_w_unassigned", "_wo_unassigned")
output_dir = file.path(dir, paste0("CyAnno_data", out))

if (dir.exists(output_dir) & opt$force) system(paste("rm -rf", output_dir))

if (unassigned) input$reference$celltype[input$reference$celltype == input$unassigned_name] <- "Unknown"
if (!unassigned) input$reference <- input$reference[input$reference$celltype != input$unassigned_name, ]
input$query$celltype[input$query$celltype == input$unassigned_name] <- "Unknown"

r <- format_data_for_cyanno(input, output_dir)
