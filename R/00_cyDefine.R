
#' Run the full cyDefine pipeline
#' @inheritParams adapt_reference
#' @inheritParams batch_correct
#' @inheritParams classify_cells
#' @inheritParams identify_unassigned
#' @param reference Tibble of reference data (cells in rows, markers in
#' columns) if you want to apply a custom reference. If NULL, the Seurat
#' PBMC atlas will be adapted and used as reference.
#' @param batch_correct Boolean indicating whether or not you want to perform
#' batch correction via cyCombine
#' @param norm_method Normalization method for cyCombine batch correction.
#' Should be either 'rank', 'scale', or 'qnorm'. Default is 'scale'.
#' @param identify_unassigned Boolean indicating whether or not you want to
#' identify unassigned cells after classifying the cells.
#' @param using_pbmc Boolean indicating whether the Seurat PBMC atlas is
#' being supplied as reference. This info is needed to enable automated
#' screening of common marker names in `map_marker_names()` and tree-based
#' merging of cell type labels.
#'
#' @param adapt_reference Boolean indicating whether the provided reference
#' should be adapted to match the query marker panel. Default: FALSE if
#' `using_pbmc` is FALSE, otherwise TRUE.
#'
#' @param save_adapted_reference Default: NULL. Save the adapted reference in a
#'  specified object name or as "adapted_reference.rds" in a specified path.
#'  Usage: FAlSE/NULL - not stored.
#'  TRUE/path - stored as "adapted_reference.rds" at "./" or path.
#'  path/to/filename.rds - stored as "filename.rds" at path/to/.
#' @return Tibble of query data with added columns: "model_prediction"
#' indicating the canonical cell type predicted by the model, and
#' "predicted_celltype" indicating the final predicted cell type after
#' identifying unassigned cells, i.e. will either be "unassigned" or the same
#' as "model_prediction".
#' @export
#'
cyDefine <- function(
    reference,
    query,
    markers,
    using_pbmc = FALSE,
    adapt_reference = ifelse(
      using_pbmc,
      TRUE,
      FALSE),
    min_f1 = 0.7,
    batch_correct = TRUE,
    ref.batch = NULL,
    xdim = 6, ydim = 6,
    identify_unassigned = TRUE,
    norm_method = ifelse(
      using_pbmc,
      "rank",
      "scale"),
    covar = NULL,
    load_model = NULL,
    save_model = NULL,
    mtry = 22,
    splitrule = "gini",
    min.node.size = 1,
    num.trees = 300,
    num.threads = 4,
    save_adapted_reference = NULL,
    unassigned_name = "unassigned",
    train_on_unassigned = ifelse(
      using_pbmc,
      FALSE,
      TRUE),
    seed = 332,
    verbose = TRUE) {

  if (inherits(reference, "character")) {
    if (reference == "pbmc" | reference == "pbmc_reference") {
      reference <- getReference("pbmc", verbose = verbose)
    }
  }

  if (adapt_reference) {

    if (verbose) message("Adapting reference to query marker panel")

    reference <- adapt_reference(
      reference = reference,
      markers = markers,
      num.threads = num.threads,
      mtry = mtry,
      min_f1 = min_f1,
      using_pbmc = using_pbmc,
      verbose = verbose
      )
    # Store adapted reference
    if (!is(save_adapted_reference, "NULL")) {
      if (is(save_adapted_reference, "logical")){
        if (save_adapted_reference) saveRDS(reference, "adapted_reference.rds")
      } else {
        if (grepl("rds", tolower(save_adapted_reference))) {
          saveRDS(reference, save_adapted_reference)
        } else {
          saveRDS(reference, file.path(save_adapted_reference, "adapted_reference.rds"))
        }
      }
    }
  }

  if (batch_correct) {

    # Batch correction via cyCombine
    corrected <- batch_correct(
      reference = reference,
      query = query,
      markers = markers,
      xdim = xdim,
      ydim = ydim,
      norm_method = norm_method,
      ref.batch = ref.batch,
      seed = seed,
      verbose = verbose
      )
    reference <- corrected$reference
    query <- corrected$query

  }

  # Canonical cell type assignment
  query <- classify_cells(
    reference = reference,
    query = query,
    markers = markers,
    unassigned_name = unassigned_name,
    load_model = load_model,
    save_model = save_model,
    mtry = mtry,
    splitrule = splitrule,
    min.node.size = min.node.size,
    num.trees = num.trees,
    num.threads = num.threads,
    seed = seed,
    verbose = verbose
    )


  # Identification of unassigned cells
  if (identify_unassigned) {

    query <- identify_unassigned(
      query = query,
      reference = reference,
      markers = markers,
      mtry = mtry,
      num.threads = num.threads,
      unassigned_name = unassigned_name,
      train_on_unassigned = train_on_unassigned,
      seed = seed,
      verbose = verbose
      )
  }

  return (list("query" = query, "reference" = reference))
}


