
#' Helper function to get spaceRAT scaffolds
#'
#' This function acquires a given prebuilt scaffold from the scaffolds
#' available on Zenodo. To prevent repeated downloads, `store` can be set to
#' `TRUE`. This will create a "scaffolds" folder, where the scaffolds are
#' stored.
#' @importFrom utils data
#' @param name Name of scaffold to get
#' @param store Whether or not to store the scaffold for future use
#' @param path Path to folder for storing scaffolds
#' @param timeout Setting download timeout
#' @param verbose Verbosity
#' @return A spaceRAT scaffold
#' @export
get_reference <- function(
    name = "pbmc",
    store = FALSE,
    path = "data",
    timeout = 120,
    verbose = TRUE){

  check_package("zen4R")
  # Extract zenodo location
  zenodo <- get_zenodo(tolower(name))

  if (verbose) message("Retrieving reference..")

  # If scaffold is already downloaded
  if (file.exists(file.path(path, zenodo$filename))){
    reference <- readRDS(file.path(path, zenodo$filename))
    return(reference)
  }
  if(store){
    # Create dir
    message("Putting '", name,"_reference' in folder '", path, "'")
    system2(c("mkdir", "-p", path))
  } else{
    path <- "./"
    if (verbose) message("Reference will be removed after download.")
  }

  # Download scaffold
  zen4R::download_zenodo(
    doi = zenodo$doi,
    files = zenodo$filename,
    path = path,
    quiet = store && verbose,
    timeout = timeout)

  reference <- readRDS(file.path(path, zenodo$filename))
  if(!store) {
    system2(c("rm", file.path(path, zenodo$filename)))
    message("Reference deleted.")
    }
  return(reference)
}




#' Extract Zenodo information about a reference
#'
#' @param name Name of reference
#'
#' @return A data.frame with Zenodo location of scaffold
#' @noRd
get_zenodo <- function(name = "pbmc"){
  # Get info of scaffolds
  allReferences <- get_all_references()

  if(length(name) > 1) stop(
    "Only one reference can be requested at a time. ",
    "The available are:\n",
    paste(allReferences$fullName, collapse=", "))

  # Remove "_scaffold" from name
  referenceName <- gsub("_reference", "", name)

  # Check name matches
  in_name <- any(referenceName %in% allReferences$name)
  if(!in_name) stop(
    name, " is not an available reference. The available are:\n",
    paste(allReferences$fullName, collapse=", "))

  # Subset references
  zenodo <- allReferences[allReferences$name == referenceName,]
  return(zenodo)
}




#' Get data frame of all references uploaded to Zenodo
#'
#' @param doi The DOI of the Zenodo entry
#'
#' @return A data.frame of scaffolds
#' @noRd
#'
#' @examples
#' all_references <- get_all_references()
get_all_references <- function(doi = "10.5281/zenodo.15394344"){
  zen <- suppressMessages(zen4R::get_zenodo(doi))
  allReferences <- do.call(rbind, lapply(zen$files, function(x) {
    data.frame(
      name = gsub("\\_reference\\..*", "", x$filename),
      # version = gsub(".*\\.v(.*)\\_.*", "\\1", x$filename),
      doi = doi,
      fullName = gsub("\\.rds", "", x$filename),
      filename = x$filename,
      stringsAsFactors = FALSE
    )
  }))
  allReferences <- allReferences[zen$versions$is_latest,]
  return(allReferences)
}
