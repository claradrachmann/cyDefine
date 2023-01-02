seurat_markers <- seurat_reference %>%
  dplyr::select(-dplyr::starts_with("celltype")) %>%
  colnames()

usethis::use_data(seurat_markers, overwrite = TRUE)
