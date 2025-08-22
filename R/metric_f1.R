library(argparse)
library(cyDefine)

# Define argument parser
parser <- ArgumentParser(description="Compute F1 score")

# Add arguments
# parser$add_argument("--output_dir", "-o", dest="output_dir", type="character", help="output directory where files will be saved"
# parser$add_argument("--dir", "-d", dest="dir", type="character", help="data directory")
parser$add_argument("--name", "-n", dest="name", type="character", help="name of the dataset")
parser$add_argument("--method", dest="method", type="character", help="name of the method")
parser$add_argument("--data", dest="data", type="character", help="input file")
parser$add_argument("--suffix", dest="suffix", type="character", help="analysis suffix", default = "_wo_unassigned")
parser$add_argument("--true_col", dest="true_col", type="character", help="True celltype column in Reference", default="celltype")
parser$add_argument("--predicted_col", dest="predicted_col", type="character", help="Predicted celltype column", default="predicted_celltype")

# Parse command-line arguments
opt <- parser$parse_args()

name <- opt$name
dir <- file.path("data", name, name)
method <- opt$method
suffix <- opt$suffix
input <- readRDS(paste0(dir, "_", method, suffix, ".rds"))
true_col <- opt$true_col
predicted_col <- opt$predicted_col

message("Computing F1 for dataset:", name, " for method:", method )

f1 <- cyDefine:::compute_f1(input[[predicted_col]], input[[true_col]])

message("Mean F1: ", mean(f1))
message("Median F1: ", median(f1))
saveRDS(f1, paste0(dir, "_", method, suffix, "_F1.rds"))
