
#' Summary function for using macro-average F1-score as metric in caret::trainControl
#'
#' @inheritParams caret::defaultSummary
#'
#' @return Macro-average F1 score
#'
macroF1_summary <- function(data, lev, model = NULL) {

  # multiclass classification performance
  perf <- crfsuite::crf_evaluation(factor(data$pred, levels=lev),
                                   factor(data$obs, levels=lev))

  # make F1s of NA count as 0
  f1 <- perf$bylabel$f1 %>%
    tidyr::replace_na(0) %>%
    mean()

  names(f1) <- "macroF1"

  return (f1)
}



#' Classify canonical cell types
#'
#' @param reference Tibble of reference data (cells in rows, markers in columns)
#' @param query Tibble of query data (cells in rows, markers in columns)
#' @param markers Character vector of available markers
#' @param unassigned_name Name used for unassigned cells
#' @param train_on_unassigned Boolean indicating whether unassigned cells should be included in model training.
#' Recommended when reference and query are samples stemming from the same experiment and unassigned
#' cells of the query are assumed to be representative of those found in the reference.
#' @param load_model Optional: Path to an rda file of a previously trained model
#' @param save_model Optional: Path to save an rda file of the trained model
#' @param param_grid Optional: Modify parameter grid to train on
#' @param return_pred Boolean indicating if only predictions should be returned as a character vector
#' @param n_cv_folds Number of cross-validation folds used for hyperparameter tuning
#' @param n_trees Number of trees in the random forest
#' @param seed Random seed
#' @param n_cores The number of cores to use for parallel execution. Obs:
#' parallel execution only performed if the 'doParallel' package is installed!
#' @param verbose Verbosity
#'
#' @return Tibble of query data with an added column of the predicted cell type, 'model_prediction'
#' @export
#'
classify_cells <- function(reference,
                           query,
                           markers,
                           unassigned_name = "unassigned",
                           load_model = NULL,
                           save_model = NULL,
                           param_grid = expand.grid(
                             mtry = as.integer(seq(
                               from = round(0.25*length(markers)),
                               to = round(0.9*length(markers)),
                               length.out = 8)),
                             min.node.size = seq(1, 6, 2),
                             splitrule = "gini"),
                           return_pred = FALSE,
                           n_cv_folds = 5,
                           n_trees = 100,
                           seed = 332,
                           n_cores = 2,
                           verbose = TRUE) {

  check_colnames(colnames(reference), c("celltype", markers))
  check_colnames(colnames(query), c(markers))

  # remove unassigned cells from reference prior to classification
  reference <- reference %>%
    dplyr::filter(celltype != !!unassigned_name)

  # keep track of ids
  if ("id" %!in% colnames(query)) {query$id <- 1:nrow(query)}

  if (!is.null(load_model)) {
    if (verbose) {message("Loading saved model: ", load_model)}
    load(load_model)
  }

  else {

    # stratified cross-validation
    set.seed(seed)
    train_control <- caret::trainControl(method = "cv",
                                         number = n_cv_folds,
                                         summaryFunction = macroF1_summary,
                                         preProcOptions = c("center", "scale"),
                                         allowParallel = TRUE,
                                         verboseIter = verbose)
    if(is.null(param_grid)){
      param_grid <- expand.grid(
        mtry = as.integer(seq(
          from = round(0.25*length(markers)),
          to = round(0.9*length(markers)),
          length.out = 8)),
        min.node.size = seq(1, 6, 2),
        splitrule = "gini"
        )
    }


    if (verbose) {message("Training random forest using ", n_cv_folds,
                          "-fold cross-validation for tuning hyperparameters")}

    model_weights <- reference %>%
      dplyr::group_by(celltype) %>%
      dplyr::mutate(weight = nrow(reference)/dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(weight = weight/sum(weight)) %>%
      dplyr::pull(weight)


    # enable parallel computation
    if (check_package("doParallel", required = FALSE)) {

      doParallel::registerDoParallel(cores = n_cores)
      if (verbose) {message("Number of parallel workers: ",
                            foreach::getDoParWorkers())}
    }

    else {
      if (verbose) {message("Package doParallel is not installed. ",
                            "Continuing without parallellization")}
    }

    rf_model <- caret::train(x = as.data.frame(reference[, markers]),
                             y = factor(reference$celltype),
                             weights = model_weights,
                             method = 'ranger',
                             trControl = train_control,
                             tuneGrid = param_grid,
                             metric = "macroF1",
                             num.trees = n_trees,
                             importance = "impurity",
                             oob.error = FALSE)

    # stop parallel computing
    if (check_package("doParallel", required = FALSE)) {
      doParallel::stopImplicitCluster()
    }

    if (verbose) {message("Optimal model parameters selected to be:\n  ",
                          paste(names(rf_model$bestTune),
                                rf_model$bestTune,
                                sep = ": ",
                                collapse = "\n  "))}

    # save model
    if (!is.null(save_model)) {save(rf_model, file = save_model)}
  }

  pred <- predict(object = rf_model,
                  newdata = query[, markers])

  if (return_pred) {return (pred)}

  # add model predictions to query
  query <- dplyr::bind_cols(query, "model_prediction" = pred)

  return (query %>% dplyr::arrange(id))
}

