library(argparse)
library(MASS)

# Define argument parser
parser <- ArgumentParser(description="Run CyTOF Linear Classifier")

# Add arguments
# parser$add_argument("--output_dir", "-o", dest="output_dir", type="character", help="output directory where files will be saved")
parser$add_argument("--name", "-n", dest="name", type="character", help="name of the dataset")
parser$add_argument("--unassigned", dest="unassigned", type="logical", help="include unassigned", default=FALSE)

# Parse command-line arguments
opt <- parser$parse_args()

# output_dir <- opt$output_dir
name <- opt$name
dir <- file.path("data", name)
input <- readRDS(file.path(dir, "preprocessed.rds"))
unassigned <- opt$unassigned
if (!unassigned) input$reference <- input$reference[input$reference$celltype != input$unassigned_name, ]

RejectionThreshold <- 0.4

t <- system.time({

# Train model
LDAclassifier <- MASS::lda(input$reference[, input$markers], as.factor(input$reference[, "celltype"]))

Model <- list(LDAclassifier = LDAclassifier,Transformation = FALSE, markers = input$markers)


# Predict celltypes
Testing.data <- input$query[,Model$markers]
Cell.types = list()

Predictions <- predict(Model$LDAclassifier, Testing.data)
Post.max <- apply(Predictions$posterior, 1, max)
if (!unassigned) Predictions$class <- factor(Predictions$class,levels = c(levels(Predictions$class),input$unassigned_name))
Predictions$class[Post.max < RejectionThreshold] <- input$unassigned_name
Predictions <- list(as.character(Predictions$class))

})

output <- data.frame(
  "celltype" = input$query$celltype,
  "predicted_celltype" = Predictions[[1]]
)

suffix <- ifelse(unassigned, "_w_unassigned", "_wo_unassigned")
saveRDS(output, paste0(dir, "/", name, "_CLC", suffix, ".rds"))

cat("Runtime (elapsed): ", t["elapsed"], file = paste0(dir, "/", name, "_CLC", suffix, "_runtime.txt"))
message("Runtime: ", t["elapsed"])
