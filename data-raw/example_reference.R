
example_reference <- prepare_data(
  data_dir = system.file("extdata", package = "cyDefine", mustWork = TRUE),
  pattern = "Plate3_Sample1",
  markers = markers,
  transform = TRUE,
  cofactor = 5,
  derand = TRUE,
  compensate = FALSE,
  extract_filename_regex = "\\d{2}Feb18_Helios2_(Plate\\d+)_(Sample\\d+)_HIMCctrl_(.+)$",
  extract_filename_into = c("batch", "sample", "celltype"))

usethis::use_data(example_reference, overwrite = TRUE)
