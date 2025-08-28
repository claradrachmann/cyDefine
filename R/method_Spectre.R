library(argparse)
library(Spectre)
library(data.table)

# Define argument parser
parser <- ArgumentParser(description="Run Spectre")

# Add arguments
parser$add_argument("--name", "-n", dest="name", type="character", help="name of the dataset")
parser$add_argument("--seed", dest="seed", type="numeric", help="seed", default = 145)
parser$add_argument("--unassigned", dest="unassigned", type="logical", help="include unassigned", default=FALSE)

# Parse command-line arguments
opt <- parser$parse_args()

# output_dir <- opt$output_dir
name <- opt$name
dir <- file.path("data", name)
input <- readRDS(file.path(dir, "preprocessed.rds"))
unassigned <- opt$unassigned
if (!unassigned) input$reference <- input$reference[input$reference$celltype != input$unassigned_name, ]
reference <- as.data.table(input$reference)
query <- as.data.table(input$query)
seed <- opt$seed

t <- system.time({

message("Classifying ", name, " using Spectre")

knn.stats <- train.knn.classifier(
  dat = as.data.table(reference),
  use.cols = input$markers,
  label.col = "celltype",
  # method = "CV", num.folds = 10,
  seed = seed,
  min.num.neighbours = 3,
  max.num.neighbours = 15)

print(knn.stats)
k.max <- knn.stats$k[which(knn.stats$accuracy == max(knn.stats$accuracy, na.rm = TRUE))]
message("Best K: ", k.max)

classified <- Spectre::run.knn.classifier(
  train.dat = reference,
  unlabelled.dat = query,
  use.cols = input$markers,
  label.col = "celltype",
  seed = seed,
  num.neighbours = k.max)

})

output <- data.frame(
  "celltype" = classified$celltype,
  "predicted_celltype" = classified$Prediction
)

# Store results
suffix <- ifelse(unassigned, "_w_unassigned", "_wo_unassigned")
saveRDS(output, paste0(dir, "/", name, "_Spectre", suffix, ".rds"))

cat("Runtime (elapsed): ", t["elapsed"], file = paste0(dir, "/", name, "_Spectre", suffix, "_runtime.txt"))
message("Runtime: ", t["elapsed"])
