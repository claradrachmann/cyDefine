
#' Compute Mahalanobis distance from each cell to a given population
#'
#' @param cells Tibble of cells (cells in rows, markers in columns) from each of which the Mahalanobis distance to the population should be computed. ???
#' @param population Tibble of a population of cells (cells in rows, markers in columns) to which the Mahalanobis distance should be computed. ???
#' @param weights Numeric vector of marker weights for computing the weighted Mahalanobis distance. ???
#'
#' @return Numeric vector of Mahalanobis distance to population for each cell.
#'
mahalanobis_distances <- function(cells,
                                  population,
                                  weights = NULL) {

  check_package("matrixcalc")

  if (nrow(population) < 3) {
    warning("Less than 3 observations are available for population ",
            population)
    return (rep(NA, nrow(cells)))
  }

  # covariance and mean of population
  covar <- stats::cov(population)
  mu <- apply(population, 2, mean)

  if (!(matrixcalc::is.positive.definite(covar))) {
    warning("Decomposition impossible for population ", population,
            " as the covariance matrix is not positive definite")
    return (rep(NA, nrow(cells)))
  }

  # Mahalanobis distance from each cell to population
  if (is.null(weights)) {
    check_package("Rfast")
    distances <- apply(cells,
                       1,
                       Rfast::mahala,
                       mu = mu,
                       sigma = covar)
  }

  else {
    distances <- apply(as.matrix(cells),
                       1,
                       wmahalanobis,
                       center = mu,
                       cov = covar,
                       weight = weights)
  }
  return(distances)
}



#' Use median absolute deviation to find the upper boundary for outlier detection
#'
#' @param distances Numeric vector of Mahalanobis distances
#' @param MAD_factor Coefficient used to define the boundary. Generally, 3 is
#' considered very conservative, 2.5 is moderately conservative, and 2 is poorly
#' conservative.
#'
#' @return Upper limit
#'
MAD_max_distance <- function(distances, MAD_factor = 3) {

  # remove outliers and find maximum distance
  MAD_score <- stats::mad(distances, na.rm = TRUE)
  upper_thres <- stats::median(distances, na.rm = TRUE) + MAD_factor*MAD_score

  return (upper_thres)
}





#' Identify unassigned cells
#' @param query Already classified... ???
#' @inheritParams classify_cells
#' @param tune_params_on_unassigned ???
#' @param pct_expl_var ???
#' @param MAD_factor ???
#'
#' @return A tibble of unassigned cells
#' @export
#'
identify_unassigned <- function(reference,
                                query,
                                markers,
                                train_on_unassigned = TRUE,
                                unassigned_name = "unassigned",
                                pct_expl_var = 0.95,
                                MAD_factor = 2.5,
                                seed = 332,
                                verbose = TRUE) {

  check_colnames(colnames(reference), c("celltype", markers))
  check_colnames(colnames(query), c("model_prediction", markers))

  # keep track of ids
  if ("id" %!in% colnames(query)) {query$id <- 1:nrow(query)}

  # if (verbose) {message("Identifying unassigned cells using ", n_components,
  #                       " principal components and MAD-factor ", MAD_factor)}

  if (train_on_unassigned) {

    # check dimensions
    n_unassigned <- reference %>%
      dplyr::filter(celltype == !!unassigned_name) %>%
      nrow()

    if (n_unassigned < 20) {
      stop("Too few cells labeled '", unassigned_name, "' are present in the reference (",
           n_unassigned, ") to train on these cells. ",
           "Please use modify the 'unassigned_name' or use the unsupervised analysis type.")
    }

    op <- pbapply::pboptions(type = "timer", char = "=")
    preds <- pbapply::pbsapply(unique(query$model_prediction),

                    function(popu) {
                      if (verbose) {message("Filtering unassigned cells from ", popu)}

                      # cell type specific query and reference
                      celltype_query <- dplyr::filter(query,
                                                      model_prediction == popu)

                      celltype_ref <- dplyr::filter(reference,
                                                    celltype == popu | celltype == !!unassigned_name)

                      # check dimension
                      n_popu <- celltype_ref %>%
                        dplyr::filter(celltype == !!popu) %>%
                        nrow()

                      if (n_popu < 30) {
                        warning("Too few cells are available for celltype '", popu,
                                "' for modelling, thus no cells predicted to belong to ", popu,
                                " will be identified as 'unassigned'")

                        return (dplyr::tibble("id" = celltype_query$id,
                                              "predicted_celltype" = !!unassigned_name))
                      }

                      y_pred <- classify_cells(reference = celltype_ref,
                                               query = celltype_query,
                                               markers = markers,
                                               unassigned_name = FALSE,
                                               return_pred = TRUE,
                                               n_cv_folds = 2,
                                               n_trees = 50,
                                               verbose = FALSE)


                      return (dplyr::tibble("id" = celltype_query$id,
                                            "predicted_celltype" = y_pred))
                    })

    preds <- dplyr::as_tibble(t(preds)) %>%
      tidyr::unnest(c(id,
                      predicted_celltype))

    query <- dplyr::left_join(query,
                              preds,
                              by = "id")
  }

  else {
    op <- pbapply::pboptions(type = "timer", char = "=")
    distances <- pbapply::pbsapply(unique(query$model_prediction),

                        function(popu) {

                          # cell type specific query and reference
                          celltype_query <- dplyr::filter(query,
                                                          model_prediction == popu)

                          celltype_ref <- dplyr::filter(reference,
                                                        celltype == popu)

                          # check dimensions
                          if (nrow(celltype_ref) <= ncol(celltype_ref)) {
                            warning("Fewer observations than number of markers are available for celltype '", popu, "'\n",
                                    "Mahalanobis distance cannot be computed, thus no cells predicted to belong to ", popu, " will be identified as 'unassigned'")

                            return (dplyr::tibble("id" = celltype_query$id,
                                                  "max_distance" = rep(Inf, nrow(celltype_query)),
                                                  "distance" = rep(-Inf, nrow(celltype_query))))
                          }

                          celltype_all <- dplyr::bind_rows(celltype_query,
                                                           celltype_ref)

                          # compute cell type specific PCA
                          pca_embed <- stats::prcomp(x = celltype_all %>%
                                                       dplyr::select(dplyr::all_of(markers)),
                                                     retx = TRUE,
                                                     center = TRUE,
                                                     scale. = TRUE)

                          # find no of PCs
                          cum_expl_var <- summary(pca_embed)$importance["Cumulative Proportion",]
                          n_components <- which(cum_expl_var >= pct_expl_var)[1]


                          pca_all <- dplyr::as_tibble(pca_embed$x[,1:n_components])

                          # split back into data and reference
                          query_pcs <- pca_all[1:nrow(celltype_query),]
                          ref_pcs <- pca_all[(nrow(celltype_query)+1):nrow(pca_all),]


                          # compute weighted Mahalanobis distance from each reference cell to its own population
                          weights <- summary(pca_embed)$importance["Proportion of Variance",][1:n_components]
                          # weights <- rep(1, n_components)

                          ref_distances <- mahalanobis_distances(cells = ref_pcs,
                                                                 population = ref_pcs,
                                                                 weights = diag(weights, length(weights)))

                          # get max distance, using MAD to sort out outliers
                          max_distance <- MAD_max_distance(ref_distances,
                                                           MAD_factor = MAD_factor)

                          # get Mahalanobis distances to predicted populations for all cells
                          query_distances <- mahalanobis_distances(cells = query_pcs,
                                                                   population = ref_pcs,
                                                                   weights = diag(weights, length(weights)))

                          return (dplyr::tibble("id" = celltype_query$id,
                                                "max_distance" = max_distance,
                                                "distance" = query_distances))
                        })

    distances <- dplyr::as_tibble(t(distances)) %>%
      tidyr::unnest(c(id,
                      max_distance,
                      distance))

    query <- dplyr::left_join(query,
                              distances,
                              by = "id")

    # change predicted population to "unassigned" for cells with distance > max_distance to their predicted celltype
    query <- query %>%
      dplyr::mutate(predicted_celltype = ifelse(distance > max_distance,
                                                "unassigned",
                                                as.character(model_prediction)))
  }

  return (query %>% dplyr::arrange(id))
}


