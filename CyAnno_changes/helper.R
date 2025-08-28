format_data_for_cyanno <- function(input, output_dir = file.path(dir, "CyAnno_data")) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Extract components from input
  reference <- input$reference
  query <- input$query
  markers <- input$markers

  # Create subdirectories
  handgated_dir <- file.path(output_dir, "handgated_files")
  livecells_dir <- file.path(output_dir, "livecells_files")
  query_dir <- file.path(output_dir, "query_files")

  dir.create(handgated_dir, showWarnings = FALSE)
  dir.create(livecells_dir, showWarnings = FALSE)
  dir.create(query_dir, showWarnings = FALSE)

  # Assume we have sample information - if not present, create artificial samples
  if (!"sample" %in% colnames(reference)) {
    # Create artificial sample IDs - you may need to adjust this logic
    reference$sample <- "sample_1"
    message("No 'sample' column found in reference data. Using 'sample_1' for all cells.")
  } else {
    reference$sample <- paste0("sample_", reference$sample)
  }

  if (!"sample" %in% colnames(query)) {
    query$sample <- "sample_2"
    message("No 'sample' column found in query data. Using 'sample_2' for all cells.")
  } else {
    query$sample <- paste0("sample_", query$sample)
  }

  # 1. Create hand-gated files (one file per cell type per sample)
  handgated_info <- data.frame()

  # Get unique combinations of sample and cell type
  sample_celltype_combinations <- reference %>%
    select(sample, celltype) %>%
    distinct()

  for (i in seq_len(nrow(sample_celltype_combinations))) {
    sample_id <- sample_celltype_combinations$sample[i]
    celltype_id <- sample_celltype_combinations$celltype[i]

    # Filter cells for this sample and cell type
    cells_subset <- reference %>%
      dplyr::filter(sample == sample_id, celltype == celltype_id) %>%
      select(all_of(markers))

    # Create filename
    filename <- paste0("handgated_", sample_id, "_", gsub("[^A-Za-z0-9]", "_", celltype_id), ".csv")
    filepath <- file.path(handgated_dir, filename)

    # Write CSV file
    write.csv(cells_subset, filepath, row.names = FALSE)

    # Add to handgated info
    handgated_info <- rbind(handgated_info, data.frame(
      file_path = filepath,
      celltype = celltype_id,
      sample_id = sample_id,
      stringsAsFactors = FALSE
    ))
  }

  # Write handgated.csv
  write.csv(handgated_info, file.path(output_dir, "handgated.csv"), row.names = FALSE)

  # 2. Create live cells training files (all cells per sample)
  livecells_training_info <- data.frame()

  # Get unique training samples
  training_samples <- unique(reference$sample)

  for (sample_id in training_samples) {
    # Get all cells from this sample (with only lineage markers)
    sample_cells <- reference %>%
      dplyr::filter(sample == sample_id) %>%
      select(all_of(markers))

    # Create filename
    filename <- paste0("livecells_training_", sample_id, ".csv")
    filepath <- file.path(livecells_dir, filename)

    # Write CSV file
    write.csv(sample_cells, filepath, row.names = FALSE)

    # Add to live cells info
    livecells_training_info <- rbind(livecells_training_info, data.frame(
      file_path = filepath,
      sample_id = sample_id,
      stringsAsFactors = FALSE
    ))
  }

  # Write LivecellsTraining.csv
  write.csv(livecells_training_info, file.path(output_dir, "LivecellsTraining.csv"), row.names = FALSE)

  # 3. Create query/unlabeled files
  livecells_query_info <- data.frame()

  # Get unique query samples
  query_samples <- unique(query$sample)

  for (sample_id in query_samples) {
    # Get all cells from this query sample (with only lineage markers)
    sample_cells <- query %>%
      dplyr::filter(sample == sample_id)  |>
      select(all_of(markers))

    # Create filename
    filename <- paste0("query_", sample_id, ".csv")
    filepath <- file.path(query_dir, filename)

    # Write CSV file
    write.csv(sample_cells, filepath, row.names = FALSE)

    # Add to query info
    livecells_query_info <- rbind(livecells_query_info, data.frame(
      file_path = filepath,
      sample_id = sample_id,
      stringsAsFactors = FALSE
    ))
  }

  # Write Livecells.csv (query)
  write.csv(livecells_query_info, file.path(output_dir, "Livecells.csv"), row.names = FALSE)

  # 4. Write lineage markers file
  writeLines(markers, file.path(output_dir, "lineage_markers.txt"))

  # Create summary
  summary_info <- list(
    n_training_samples = length(training_samples),
    n_query_samples = length(query_samples),
    n_celltypes = length(unique(reference$celltype)),
    n_markers = length(markers),
    celltypes = unique(reference$celltype),
    markers = markers
  )

  # Write summary as JSON for reference
  # jsonlite::write_json(summary_info, file.path(output_dir, "data_summary.json"), pretty = TRUE)

  # Print summary
  cat("CyAnno data formatting completed!\n")
  cat("Output directory:", output_dir, "\n")
  cat("Training samples:", length(training_samples), "\n")
  cat("Query samples:", length(query_samples), "\n")
  cat("Cell types:", paste(unique(reference$celltype), collapse = ", "), "\n")
  cat("Lineage markers:", paste(markers, collapse = ", "), "\n")
  cat("\nFiles created:\n")
  cat("- handgated.csv (", nrow(handgated_info), " entries)\n")
  cat("- LivecellsTraining.csv (", nrow(livecells_training_info), " entries)\n")
  cat("- Livecells.csv (", nrow(livecells_query_info), " entries)\n")
  cat("- lineage_markers.txt\n")
  cat("- Individual CSV files in subdirectories\n")

  return(list(
    handgated_info = handgated_info,
    livecells_training_info = livecells_training_info,
    livecells_query_info = livecells_query_info,
    summary = summary_info
  ))
}

