
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
#' Should be either 'rank', 'scale' or 'qnorm'. Default is 'scale' unless
#' reference is NULL, then 'rank', assuming heavy batch effects between Seurat
#' atlas and query.
#' @param identify_unassigned Boolean indicating whether or not you want to
#' identify unassigned cells after classifying the cells.
#' @param using_seurat Boolean indicating whether the Seurat PBMC atlas is
#' being supplied as reference. This info is needed to enable automated
#' screening of common marker names in `map_marker_names()` and tree-based
#' merging of cell type labels.
#'
#' @param adapt_reference Boolean indicating whether the provided reference
#' should be adapted to match the query marker panel. Default: FALSE if
#' `using_seurat` is FALSE, otherwise TRUE.
#'
#' @return Tibble of query data with added columns: "model_prediction"
#' indicating the canonical cell type predicted by the model, and
#' "predicted_celltype" indicating the final predicted cell type after
#' identifying unassigned cells, i.e. will either be "unassigned" or the same
#' as "model_prediction".
#' @export
#'
cyDefine <- function(reference,
                     query,
                     markers,
                     using_seurat = FALSE,
                     adapt_reference = ifelse(using_seurat,
                                              TRUE,
                                              FALSE),
                     batch_correct = TRUE,
                     identify_unassigned = TRUE,
                     norm_method = ifelse(using_seurat,
                                          "rank",
                                          "scale"),
                     covar = NULL,
                     load_model = NULL,
                     save_model = NULL,
                     param_grid = expand.grid(
                       mtry = as.integer(seq(
                         from = round(0.25*length(markers)),
                         to = round(0.9*length(markers)),
                         length.out = 8)),
                       splitrule = "gini",
                       min.node.size = seq(1, 6, 2),
                       num.trees = seq(200, 500, 100)),
                     n_cv_folds = 5,
                     # save_adapted_reference    OBS: ADD THIS ARGUMENT!
                     unassigned_name = "unassigned",
                     train_on_unassigned = ifelse(using_seurat,
                                                  FALSE,
                                                  TRUE),
                     n_cores = 2,
                     seed = 332,
                     verbose = TRUE) {
# TODO: Remove 0 columns
  if (adapt_reference) {

    if (verbose) message("Adapting reference to query marker panel")

    reference <- adapt_reference(reference = reference,
                                 query = query,
                                 markers = markers,
                                 using_seurat = using_seurat,
                                 verbose = verbose)
  }

  if (batch_correct) {

    # Batch correction via cyCombine
    corrected <- batch_correct(reference = reference,
                               query = query,
                               markers = markers,
                               norm_method = norm_method,
                               seed = seed,
                               verbose = verbose)
    reference <- corrected$reference
    query <- corrected$query

  }

  # Canonical cell type assignment
  query <- classify_cells(reference = reference,
                          query = query,
                          markers = markers,
                          unassigned_name = unassigned_name,
                          load_model = load_model,
                          save_model = save_model,
                          param_grid = param_grid,
                          n_cv_folds = n_cv_folds,
                          n_cores = n_cores,
                          seed = seed,
                          verbose = verbose)


  # Identification of unassigned cells
  if (identify_unassigned) {

    query <- identify_unassigned(query = query,
                                 reference = reference,
                                 markers = markers,
                                 unassigned_name = unassigned_name,
                                 train_on_unassigned = train_on_unassigned,
                                 seed = seed,
                                 verbose = verbose)
  }

  return (query)
}


