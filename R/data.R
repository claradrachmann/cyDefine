#' Markers of the Seurat PBMC atlas
#'
#' @format A character vector of length 241
#'
"pbmc_markers"


#' The modified Seurat PBMC atlas
#'
#' Load the modified PBMC reference with
#' `pbmc_reference <- get_reference("pbmc")`
#' See `?get_reference` for options to store the reference locally for subsequent faster load times.
#'
#' @format A tibble of 161,764 rows and 244 variables
#' @source \url{https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat}
#'
"pbmc_reference"


#' An example of a CyTOF reference, already prepared by cyDefine as seen in
#' ´data-raw/example_reference.R´
#'
#' @format A tibble of 80,405 rows and 46 columns
#' @source \url{http://flowrepository.org/id/FR-FCM-Z2YR}
#'
"example_reference"


#' An example of a CyTOF query, already prepared by cyDefine as seen in
#' ´data-raw/example_query.R´
#'
#' @format A tibble of 79,179 rows and 45 columns
#' @source \url{http://flowrepository.org/id/FR-FCM-Z2YR}
#'
"example_query"


#' The marker panel of `example_reference` and `example_query`, i.e., the panel
#' of the Z2YR CyTOF dataset.
#'
#' @format A character vector of length 42
#'
"example_markers"
