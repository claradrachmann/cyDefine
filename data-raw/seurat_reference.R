library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(ggplot2)

# The Seurat PBMC reference is obtained from https://atlas.fredhutch.org/nygc/multimodal-pbmc/
system("wget https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat")

seu_obj <- LoadH5Seurat("/home/projects/cytograted/data/aml/scrna-seq/pbmc_multimodal.h5seurat",
                        assays = "ADT",
                        reductions = FALSE,
                        graphs = FALSE,
                        images = FALSE
)

# Remove raw data
system("rm pbmc_multimodal.h5seurat")

# --------------- cyCombine ---------------

seu_uncorrected <- seu_obj@assays$ADT$

  # --------------- Integration ---------------

# split by orig.ident which is sample of origin (patient + time point)
seu_obj_list <- SplitObject(seu_obj, split.by = "orig.ident")

seu_obj_list <- lapply(
  X = seu_obj_list,
  FUN = function(x) {
    # CLR transform per cell
    x <- NormalizeData(x, normalization.method = "CLR", margin = 2)
    x <- FindVariableFeatures(x, selection.method = "vst")
  }
)

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seu_obj_list)
anchors <- FindIntegrationAnchors(
  object.list = seu_obj_list,
  anchor.features = features
)

# integrate data
seurat_integrated <- IntegrateData(anchorset = anchors)

seurat_integrated <- readRDS("/home/projects/cytograted/data/seurat_integrated.rds")


# --------------- Make tibble of normalized data and metadata --------------- #

# combine marker expression and cell population label
seurat_reference <- seurat_integrated@assays$integrated@data %>%
  as.matrix() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  as_tibble()

seurat_metadat <- seurat_integrated@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  as_tibble() %>%
  select(cell_id, starts_with("celltype"))

seurat_reference <- left_join(seurat_reference, seurat_metadat, by = "cell_id")

## Fix duplicated markers ----

if (FALSE) {
  simple_markers <- colnames(seurat_reference) |>
    stringr::str_remove("-[12]$")
  duplicated_markers <- simple_markers[duplicated(simple_markers) & !(simple_markers %in% c("Notch", "TCR", "Rat-IgG1", "CD3"))]
  duplicated_location <- which(simple_markers %in% duplicated_markers)
  # duplicated_markers <- c(duplicated_markers, "CD3", "TCR")
  seurat_reference[duplicated_markers] <- NA
  seurat_reference <- seurat_reference |>
    # rowwise() |>
    mutate(across(all_of(duplicated_markers),
                  ~ pmax(
                    get(paste0(cur_column(), "-1")),
                    get(paste0(cur_column(), "-2"))
                  ),
                  .names = "{col}"))
  # mutate(across(all_of(duplicated_markers),
  #               ~ rowMeans(select(seurat_reference,
  #                 paste0(cur_column(), "-1"),
  #                 paste0(cur_column(), "-2"))
  #               ),
  #               .names = "{col}"))
  seurat_reference[, duplicated_location] <- NULL

  seurat_reference <- seurat_reference |>
    mutate(`CD3-1`= NULL) |>
    rename(
      "CD3" = "CD3-2",
      "gdTCR" = "TCR-1",
      "abTCR" = "TCR-2"
    )

  # dup_ref <- lapply(setNames(duplicated_markers, duplicated_markers), function(marker) {
  #   dup_max <- apply(seurat_reference[, paste0(marker, c("-1", "-2"))], 1, max)
  #   return(dup_max)
  # })
  # seurat_reference[, duplicated_location] <- NULL
  # seurat_reference <- cbind(seurat_reference, dup_ref)

  usethis::use_data(seurat_reference, overwrite = TRUE)


  # ------------ Modify naming convention of markers present in two versions ------------ #
} else {
  # check markers with several versions
  multi_version_markers <- c("CD3", "CD4", "CD56", "CD11b", "CD38", "CD44", "CD26", "CD275", "CD45", "CD133", "CD138", "TCR")
  pccs <- as_tibble(matrix(ncol = length(multi_version_markers), nrow = 1))
  colnames(pccs) <- multi_version_markers

  # compute correlation between marker versions
  for (mar in multi_version_markers) {
    pccs[[mar]] <- cor(
      seurat_reference[[str_c(mar, "-1")]],
      seurat_reference[[str_c(mar, "-2")]]
    )
  }

  # add columns with mean for multi-version markers with high PCC
  seurat_reference[multi_version_markers[pccs > 0.5]] <- NA
  seurat_reference <- seurat_reference |>
    mutate(across(all_of(multi_version_markers[pccs > 0.5]),
                  ~ pmax(
                    get(paste0(cur_column(), "-1")),
                    get(paste0(cur_column(), "-2"))
                  ),
                  .names = "{col}"))

  for (mar in multi_version_markers[pccs > 0.5]) {
    seurat_reference <- seurat_reference %>%
      mutate("{mar}" := rowMeans(select(
        seurat_reference,
        str_c(mar, "-1"),
        str_c(mar, "-2")
      )))
  }

  # for rest of multi-version markers, use the marker seeming most informative
  # (evaluated manually)
  seurat_reference$CD26 <- seurat_reference$`CD26-2`
  seurat_reference$CD275 <- seurat_reference$`CD275-2`
  seurat_reference$CD45 <- seurat_reference$`CD45-2`
  seurat_reference$CD133 <- seurat_reference$`CD133-1`
  seurat_reference$CD138 <- seurat_reference$`CD138-1`
  seurat_reference$CD3 <- seurat_reference$`CD3-2`


  multi_version_markers_suffix <- unlist(lapply(multi_version_markers, function(marker) {
    c(paste0(marker, "-1"), paste0(marker, "-2"))
  }))



  # also rename CD3-1 and CD3-2 + TCR-1 and TCR-2 to more informative names
  seurat_reference <- seurat_reference %>%
    rename(
      "CD3D" = "CD3-1",
      "CD3E" = "CD3-2",
      "gdTCR" = "TCR-1",
      "abTCR" = "TCR-2"
    )

  seurat_reference <- select(seurat_reference, -any_of(multi_version_markers_suffix))

  usethis::use_data(seurat_reference, overwrite = TRUE)
}
