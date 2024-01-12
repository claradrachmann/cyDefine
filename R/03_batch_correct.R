#' Perform batch correction via cyCombine
#'
#' @param reference Tibble of reference data (cells in rows, markers in columns)
#' @param query Tibble of query data (cells in rows, markers in columns)
#' @param markers Character vector of available markers
#' @inheritParams cyCombine::batch_correct
#' @param verbose Verbosity
#'
#' @return Named list of tibbles of batch corrected reference and query data
#' @export
#'
batch_correct <- function(reference,
                          query,
                          markers,
                          covar = NULL,
                          norm_method = "scale",
                          seed = 332,
                          verbose = TRUE) {
  check_colnames(colnames(reference), c("celltype", markers))
  check_colnames(colnames(query), c(markers))


  # add batch and sample columns if not provided
  if ("sample" %!in% colnames(reference)) {
    if (verbose) {
      message("Sample information is not provided for reference - assuming one sample.")
    }
    reference <- reference %>%
      dplyr::mutate(sample = "ref_sample")
  }

  if ("sample" %!in% colnames(query)) {
    if (verbose) {
      message("Sample information is not provided for query - assuming one sample.")
    }
    query <- query %>%
      dplyr::mutate(sample = "dat_sample")
  }

  if ("batch" %!in% colnames(reference)) {
    if (verbose) {
      message("Batch information is not provided for reference - assuming one batch.")
    }
    reference <- reference %>%
      dplyr::mutate(batch = "ref_batch")
  }

  if ("batch" %!in% colnames(query)) {
    if (verbose) {
      message("Batch information is not provided for query - assuming one batch.")
    }
    query <- query %>%
      dplyr::mutate(batch = "dat_batch")
  }

  # check uniqueness of samples
  if (length(intersect(
    unique(reference$sample),
    unique(query$sample)
  )) > 0) {
    warning(
      "Overlapping sample ID(s) found between reference and query. ",
      "Assuming that these represent different samples. Adding '_ref' ",
      "and '_query', respectively, to the end of overlapping sample ID(s)"
    )

    for (s in intersect(
      unique(reference$sample),
      unique(query$sample)
    )) {
      reference <- reference %>%
        dplyr::mutate(sample = replace(
          sample,
          sample == s,
          paste0(s, "_ref")
        ))
      query <- query %>%
        dplyr::mutate(sample = replace(
          sample,
          sample == s,
          paste0(s, "_query")
        ))
    }
  }

  if (length(intersect(
    unique(reference$batch),
    unique(query$batch)
  )) > 0 &
    verbose) {
    message(
      "Overlapping batch IDs found between reference and query. ",
      "If originating from different batches, make the IDs unique and ",
      "rerun batch correction."
    )
  }


  # combine query and reference
  uncorrected <- dplyr::full_join(reference,
    query,
    by = intersect(
      colnames(reference),
      colnames(query)
    )
  )
  uncorrected$id <- 1:nrow(uncorrected)

  # run batch correction
  corrected <- cyCombine::batch_correct(
    df = uncorrected,
    markers = markers,
    norm_method = norm_method, # "rank" is recommended when heavy batch effects
    covar = covar,
    rlen = 10, # higher values recommended if 10 does not appear to perform well
    seed = seed
  )

  corrected_ref <- corrected %>%
    dplyr::filter(sample %in% unique(reference$sample)) %>%
    dplyr::select(colnames(reference))
  corrected_query <- corrected %>%
    dplyr::filter(sample %in% unique(query$sample)) %>%
    dplyr::select(colnames(query))

  return(list(
    "reference" = corrected_ref,
    "query" = corrected_query
  ))
}
