
#' Prepare a directory of .fcs files
#'
#' This is a wrapper function of cyCombine's data preparation functions, with
#' some extra added functionality, taking you from a directory of FCS files or
#' a flowset to a tidy tibble.
#'
#' @param flowset Optional: Prepare a flowset instead of a directory of fcs files
#' @inheritParams cyCombine::compile_fcs
#' @inheritParams cyCombine::convert_flowset
#' @inheritParams cyCombine::transform_asinh
#' @param transform If TRUE, the data will be transformed; if FALSE, it will not.
#' @param extract_filename_regex Optional: Use, if there are details that you
#' want to keep (e.g. sample, batch or cell type information) saved in the
#' filenames. Can be used together with or in place of metadata. Should be a
#' string with a regex with groups capturing the information of interest.
#' Example: "Helios2_(Plate\\d+)_(Sample\\d+)_" could extract plate and well
#' from files named something like "Helios2_Plate21_Sample1_ctrl.fcs".
#' @param extract_filename_into Optional: Only if extract_filename_regex is
#' given. A character vector of names corresponding to the capture groups in
#' extract_filename_regex. These names will represent the column names of the
#' resulting data frame. Example: (matching the example above)
#' extract_filename_into = c("batch", "sample").
#'
#' @return Tibble of data (cells in rows, markers in columns)
#'
#' @examples
#' data_dir <- system.file("extdata", package = "cyDefine", mustWork = TRUE)
#'
#' uncorrected <- data_dir %>%
#'   prepare_data(extract_filename_regex = "Helios2_(Plate\\d+)_(Sample\\d+)_",
#'   extract_filename_into = c("batch", "sample"),
#'   markers = markers)
#'
#' \dontrun{
#' uncorrected <- data_dir %>%
#'   prepare_data(metadata = "metadata.csv",
#'   markers = markers,
#'   filename_col = "FCS_name",
#'   batch_ids = "Batch",
#'   condition = "condition",
#'   down_sample = TRUE,
#'   sample_size = 100000)
#' }
#'
#' @export
prepare_data <- function(data_dir = NULL,
                         flowset = NULL,
                         markers = NULL,
                         pattern = "\\.fcs",
                         extract_filename_regex = NULL,
                         extract_filename_into = NULL,
                         metadata = NULL,
                         filename_col = "filename",
                         sample_ids = NULL,
                         batch_ids = NULL,
                         condition = NULL,
                         anchor = NULL,
                         down_sample = FALSE,
                         sample_size = 500000,
                         sampling_type = "random",
                         seed = 473,
                         panel = NULL,
                         panel_channel = "fcs_colname",
                         panel_antigen = "antigen",
                         transform = TRUE,
                         cofactor = 5,
                         derand = TRUE,
                         compensate = FALSE,
                         .keep = FALSE,
                         clean_colnames = TRUE,
                         verbose = TRUE){

  # Stop if no data is given
  if(is.null(data_dir) & is.null(flowset)) stop("No data given.")

  # directory of FCS files
  if(!is.null(data_dir)){
    # Remove slash at end of data_dir
    if(data_dir %>% endsWith("/")) data_dir <- data_dir %>% stringr::str_sub(end = -2)

    if (verbose) {message("Preparing FCS files in directory ", data_dir)}

    # Compile directory to flowset
    if(is.null(flowset)){
      flowset <- data_dir %>%
        cyCombine::compile_fcs(pattern = pattern)
    }

    # Compensate for spectral overlap
    if (compensate) {
      if (verbose) message("Compensating for spectral overlap between fluorescence channels")
      comp <- flowCore::fsApply(flowset,
                                function(x) flowCore::spillover(x)$SPILL,
                                simplify = FALSE)
      flowset <- flowCore::compensate(flowset, comp)
    }

    # Look for metadata in data_dir
    if(!is.null(metadata)){
      if("data.frame" %!in% class(metadata)){
        if(!file.exists(file.path(metadata)) & file.exists(file.path(data_dir, metadata))) metadata <- file.path(data_dir, metadata)
      }
    }
  }

  # Convert flowset to dataframe
  if (verbose) message("Converting flowset to data frame")
  fcs_data <- flowset %>%
    cyCombine::convert_flowset(metadata = metadata,
                               filename_col = filename_col,
                               sample_ids = sample_ids,
                               batch_ids = batch_ids,
                               condition = condition,
                               anchor = anchor,
                               down_sample = down_sample,
                               sample_size = sample_size,
                               sampling_type = sampling_type,
                               seed = seed,
                               panel = panel,
                               panel_channel = panel_channel,
                               panel_antigen = panel_antigen,
                               clean_colnames = clean_colnames) %>%
    # Transform dataset with asinh
    purrr::when(transform ~ cyCombine::transform_asinh(., markers = markers,
                                                       cofactor = cofactor,
                                                       derand = derand,
                                                       .keep = .keep),
                ~ .)

  # Extract relevant information from file names (saved in 'sample' column)
  if (!is.null(extract_filename_regex)) {
    fcs_data <- fcs_data %>%
      tidyr::extract(col = sample,
                     into = extract_filename_into,
                     regex = extract_filename_regex)
  }

  if (verbose) message("Done!")
  return(fcs_data)
}

