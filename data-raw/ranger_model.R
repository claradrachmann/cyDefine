# Modify the ranger function so we can include num.trees as an argument in parameter tuning

# Code adapted from:
# https://rstudio-pubs-static.s3.amazonaws.com/519785_adeb1a4592584939986b1b23103d213e.html

# Get the existing method for the ranger() function
ranger_model = caret::getModelInfo(model = "ranger")
# Change the parameters that may be tuned to include num.trees
ranger_model$ranger$parameters = data.frame("parameter" = c("mtry", "splitrule", "min.node.size", "num.trees"),
                                           "class" = c("numeric", "character", "numeric", "numeric"),
                                           "label" = c("Selected Predictors",
                                                       "Number of Trees",
                                                       "Minimal node size",
                                                       "Splitting rule"))

# Edit the model fit function to include a num.trees argument in the call to the ranger function()
ranger_model$ranger$fit = function (x, y, wts, param, lev, last, classProbs,
                                   ...) {
  if ((!is.data.frame(x)) || dplyr::is.tbl(x))
    x <- as.data.frame(x)
  x$.outcome <- y
  if (!is.null(wts)) {
    out <- ranger::ranger(dependent.variable.name = ".outcome",
                          data = x, mtry = param$mtry, num.trees = param$num.trees,
                          splitrule = param$splitrule, min.node.size = param$min.node.size,
                          probability = classProbs, case.weights = wts, ...)
  }
  else {
    out <- ranger::ranger(dependent.variable.name = ".outcome",
                          data = x, mtry = param$mtry, num.trees = param$num.trees,
                          splitrule = param$splitrule, min.node.size = param$min.node.size,
                          write.forest = TRUE, probability = classProbs, ...)
  }
  if (!last)
    out$y <- y
  out
}

usethis::use_data(ranger_model, overwrite = TRUE)



