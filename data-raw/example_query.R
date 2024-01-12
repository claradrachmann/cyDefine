example_query <- prepare_data(
  data_dir = system.file("extdata", package = "cyDefine", mustWork = TRUE),
  pattern = "Plate4_Sample1",
  markers = markers,
  transform = TRUE,
  cofactor = 5,
  derand = TRUE,
  compensate = FALSE,
  extract_filename_regex = "^\\d{2}Feb18_Helios2_(Plate\\d+)_(Sample\\d+)",
  extract_filename_into = c("batch", "sample")
)

usethis::use_data(example_query, overwrite = TRUE)
