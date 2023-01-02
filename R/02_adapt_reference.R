
#' Change query marker names to match reference marker names
#'
#' @param using_seurat Boolean indicating whether or not the Seurat PBMC atlas is
#' used as reference
#' @param ref_markers Character vector of reference markers
#' @param map_specific_from Optional: used together with `map_specific_to`.
#' Character vector of specific marker(s) of the query to manually map to
#' specific marker(s) of the reference (given by `map_specific_to`).
#' @param map_specific_to Optional: used together with `map_specific_from`.
#' Markers of the reference that the markers given by `map_specific_from` should
#' be mapped to.
#' @param verbose Verbosity
#' @param query Tibble of query data (cells in rows, markers in columns)
#' @param query_markers Character vector of query markers
#'
#' @return Tibble of query data modified to make marker names match reference
#' markers
#' @family adapt
#'
#' @examples
#'
#' map_marker_names(using_seurat = TRUE, query = example_query,
#' query_markers = example_markers, map_specific_from = c("CXCR5", "TCRgd"),
#' map_specific_to = c("CD185","gdTCR"))
#'
#' @export
map_marker_names <- function(query,
                             ref_markers,
                             query_markers,
                             using_seurat = TRUE,
                             map_specific_from = NULL,
                             map_specific_to = NULL,
                             verbose = TRUE) {

  # check required columns
  check_colnames(avail_colnames = colnames(query),
                 req_colnames = c(query_markers))

  if (!is.null(map_specific_from)) {

    if (verbose) {message("Mapping specified markers to reference markers:\n",
                          paste(map_specific_from,
                                "to",
                                map_specific_to,
                                collapse = '\n'),
                          sep='\n')}

    query <- query %>%
      dplyr::rename_at(dplyr::vars(dplyr::all_of(map_specific_from)),
                ~ map_specific_to)

    # update data markers
    query_markers <- c(query_markers[query_markers %!in% map_specific_from],
                     map_specific_to)
  }

  avail_markers_in_ref <- query_markers[query_markers %in% ref_markers]

  # further investigation if not all markers are found in reference
  if (length(avail_markers_in_ref) != length(query_markers)) {

    not_in_ref <- query_markers[query_markers %!in% ref_markers]

    if (using_seurat) {

      # common marker names
      common_names <- list("CXCR5"="CD185", "CCR4"="CD194", "CCR6"="CD196",
                           "TIM3"="CD366", "TIM-3"="CD366", "PDL1"="CD274",
                           "PD-L1"="CD274", "PD1"="CD279", "PD-1"="CD279",
                           "CTLA4"="CD152", "TCRgd"="gdTCR", "TCRab"="abTCR",
                           "CCR9"="CD199","CDw199"="CD199", "ICOS"="CD278",
                           "LAG3"="CD223", "CD40L"="CD154", "CD40LG"="CD154")

      old_names <- intersect(names(common_names), not_in_ref)
      new_names <- as.character(common_names[old_names])

      # update not_in_ref
      not_in_ref <- not_in_ref[not_in_ref %!in% old_names]
    }
    else (old_names <- c())

    # if still un-mapped markers, check minor differences in names
    for (name in not_in_ref) {

      # get index of reference marker matching name
      # ignore non-alphabetic and non-numeric signs + case
      ref_marker_idx <- which(toupper(stringr::str_replace_all(string = ref_markers,
                                                               pattern = "[^[:alnum:]]",
                                                               replacement = "")) ==
                                toupper(stringr::str_replace_all(string = name,
                                                                 pattern = "[^[:alnum:]]",
                                                                 replacement = "")))
      if (length(ref_marker_idx) == 1) {
        old_names <- c(old_names, name)
        new_names <- c(new_names, ref_markers[ref_marker_idx])
      }
    }

    if (!is.null(old_names)) {
      # rename markers
      if (verbose) {message(paste("Renaming",
                                  old_names,
                                  "to",
                                  new_names,
                                  sep = " ",
                                  collapse = '\n'))}

      query <- query %>%
        dplyr::rename_at(dplyr::all_of(old_names), ~ new_names)
    }

    # remove markers not found in reference
    not_in_ref <- not_in_ref[not_in_ref %!in% old_names]
    if (length(not_in_ref) > 0) {
      warning("The following markers were not detected to be in the reference and were excluded from the data:\n  ",
              paste(not_in_ref, collapse = ', '), "\n  ",
              "Markers of the reference can be found in 'seurat_markers'. ",
              "If needed, markers can be manually mapped to reference markers using the 'map_specfic_from' and 'map_specfic_to' arguments")

      query <- query %>%
        dplyr::select(-dplyr::all_of(not_in_ref))
    }
  }
  return (query)
}




#' Merge groups of similar cell populations
#'
#' @param populations_to_merge Tibble ???
#' @param reference Tibble of reference data (cells in rows, markers in columns)
#' @param name_unassigned_if_merged Cell types of the labeling tree, where if merged, the merged population should rather be named 'unassigned'
#' @family adapt
merge_populations <- function(populations_to_merge,
                              reference,
                              using_seurat = TRUE) {

  # get labels of clusters of merged populations
  populations_to_merge <- populations_to_merge %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(popu = popu,
                     merged_label = get_merged_label(popu,
                                                     using_seurat = using_seurat),
                     .groups = "drop")

  # modify cell type labels in reference to correspond to clusters
  reference <- dplyr::left_join(reference,
                                populations_to_merge %>%
                                  dplyr::select(popu,
                                                merged_label),
                                by = c("celltype" = "popu")) %>%
    dplyr::mutate(celltype = dplyr::coalesce(merged_label,
                                             celltype)) %>%
    dplyr::select(-merged_label)

  return (reference)
}



#' Get label of a merged population
#'
#' @param populations Character vector of individual cell types in merged population
#' @inheritParams merge_populations
#'
#' @return Label of merged population
#' @family adapt
get_merged_label <- function(populations, using_seurat = TRUE) {

  if (using_seurat) {

    check_package("reshape2")

    levels <- c(stringr::str_c("L",
                               1:(ncol(labeling_tree)-1)), "leaf")
    i <- 1

    while (i < length(levels)-1) {
      parent_popus <- unique(na.omit(labeling_tree[[levels[i]]]))

      for (parent in parent_popus) {

        leaf_popus <- dplyr::filter(labeling_tree,
                                    get(levels[i]) == parent) %>%
          dplyr::pull(leaf)

        # if all popus of parent node are to be merged
        if (length(setdiff(leaf_popus, populations)) == 0) {
          populations <- c(populations[(populations %!in% leaf_popus)], parent)
        }
      }
      i <- i+1
    }

  }

  # concatenate remaining populations
  label <- paste(populations,
                 collapse = ' / ')

  return (label)
}




#' Identify indistinguishable cell type populations based on query markers
#'
#' @param population_pairs Data frame of two character columns with one cell
#' type in each, and rows representing all pairs of cell types to be tested.
#' @param reference
#' @param min_f1
#' @param seed
#' @param verbose
#'
#' @return
#'
#' @examples
identify_similar_populations <- function(population_pairs,
                                         reference,
                                         markers,
                                         min_f1 = 0.6,
                                         seed = 332,
                                         verbose = TRUE) {

  check_colnames(colnames(reference), c(markers, "celltype"))

  # check number of cells of each cell type
  celltype_counts <- summary(factor(reference$celltype))

  if (length(celltype_counts[celltype_counts < 5]) > 0) {

    if (verbose) {message("OBS: Less than 5 cells of type(s):\n",
                          paste(names(celltype_counts[celltype_counts < 5]),
                                collapse = "\n"),
                          "\nare present in the reference. These will not be ",
                          "considered for population merging.")}

    reference <- reference %>%
      dplyr::group_by(celltype) %>%
      dplyr::filter(dplyr::n() >= 5) %>%
      dplyr::ungroup()
  }


  if (verbose) {message("Running GBM for pairs of subpopulations")}

  set.seed(seed)
  op <- pbapply::pboptions(type = "timer", char = "=")
  gbm_perf <- pbapply::pbapply(population_pairs,
                    1,
                    function(popu_pair) {
                      tmp_ref <- reference %>%
                        dplyr::filter(celltype %in% popu_pair)

                      X <- tmp_ref %>%
                        dplyr::select(dplyr::all_of(markers))
                      y <- tmp_ref %>%
                        dplyr::pull(celltype)

                      # stratified data partition
                      train_idx <- caret::createDataPartition(y,
                                                              list = FALSE,
                                                              p = 0.5)
                      X_train <- as.matrix(X[train_idx,])
                      X_test <- as.matrix(X[-train_idx,])
                      y_train <- factor(y[train_idx])
                      y_test <- factor(y[-train_idx])

                      model_weights <- tmp_ref[train_idx,] %>%
                        dplyr::group_by(celltype) %>%
                        dplyr::mutate(weight = nrow(tmp_ref[train_idx,])/dplyr::n()) %>%
                        dplyr::ungroup() %>%
                        dplyr::mutate(weight = weight/sum(weight)) %>%
                        dplyr::pull(weight)

                      model_gbm <- caret::train(X_train,
                                                y_train,
                                                method = 'gbm',
                                                weights = model_weights,
                                                preProcess = c('center','scale'),
                                                trControl = caret::trainControl(method = "none"),
                                                tuneGrid = expand.grid(n.trees = 100,
                                                                       interaction.depth = 3,
                                                                       n.minobsinnode = 10,
                                                                       shrinkage = 0.05),
                                                verbose = FALSE)

                      y_pred <- predict(object = model_gbm,
                                        newdata = X_test)

                      # ensure same levels for pred and obs
                      y_pred <- factor(y_pred,
                                       levels = levels(y_test))

                      # compute performance
                      perf <- crfsuite::crf_evaluation(y_pred,
                                                       y_test)

                      macro_f1 <- perf$overall["f1_mean"]

                      return (macro_f1)
                    })

  colnames(population_pairs) <- c("popu1", "popu2")

  # get pairs of similar populations
  similar_population_pairs <- dplyr::tibble(population_pairs,
                                            "F1" = gbm_perf) %>%
    dplyr::filter(F1 < min_f1) %>%
    dplyr::select(popu1, popu2)

  population_clusters <- igraph::components(
    igraph::graph_from_data_frame(similar_population_pairs))$membership

  similar_populations <- dplyr::tibble(popu = names(population_clusters),
                                       cluster = population_clusters)

  return (similar_populations)
}



#' Adapt Seurat PBMC reference to query marker panel
#'
#' @param reference Tibble of reference data (cells in rows, markers in columns)
#' @param markers Character vector of available markers present in both query
#' and reference.
#' @param min_f1 Minimum F1 score required for two cell types to be separated
#' @param exclude_celltypes Only relevant if `using_seurat = TRUE`. Character
#' vector of cell types of the Seurat reference to NOT consider during cell
#' type classification. Defaults to non-PBMCs.
#' @inheritParams map_marker_names
#' @inheritParams merge_populations
#'
#' @return Tibble of adapted reference
#' @export
#'
#' @examples
adapt_reference <- function(reference,
                            query,
                            markers,
                            using_seurat = TRUE,
                            exclude_celltypes = c("Doublet",
                                                  "Platelet",
                                                  "Eryth",
                                                  "HSPC"),
                            celltype_col = ifelse(using_seurat,
                                                  "celltype.l2",
                                                  "celltype"),
                            min_f1 = 0.7,
                            seed = 332,
                            verbose = TRUE) {

  # remove non-available markers and excluded cell types
  reference <- reference %>%
    dplyr::rename("celltype" = !!celltype_col) %>%
    dplyr::filter(celltype %!in% exclude_celltypes)

  # check merged against remaining for as long as new merges are performed
  # start by having all cell types as queries
  query_popus <- unique(reference$celltype)
  remaining_popus <- query_popus

  set.seed(seed)

  i <- 1
  while (length(query_popus) > 0) {

    if (verbose) {message("# ------- Population merging - Round ", i, " ------- #")}

    # find similar populations for queries against all remaining populations
    population_pairs <- expand.grid(query_popus,
                                    remaining_popus,
                                    stringsAsFactors = FALSE) %>%
      dplyr::filter(Var1 > Var2)


    similar_populations <- identify_similar_populations(
      population_pairs = population_pairs,
      reference = reference,
      markers = markers,
      min_f1 = min_f1)

    # merge groups of similar populations
    if (nrow(similar_populations) > 0) {

      if (verbose) {message("Merging groups of similar populations")}

      reference <- merge_populations(populations_to_merge = similar_populations,
                                     reference = reference,
                                     using_seurat = using_seurat)

      if (verbose) {message("Merging ", similar_populations)}
    }

    # update queries and remaining populations
    query_popus <- setdiff(unique(reference$celltype),
                           remaining_popus)

    remaining_popus <- unique(reference$celltype)

    i <- i + 1
  }

  if (verbose) {message("\nReference computed!")}

  return (reference)
}






