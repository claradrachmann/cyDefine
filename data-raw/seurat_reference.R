library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(ggplot2)

# The Seurat PBMC reference is obtained from https://atlas.fredhutch.org/nygc/multimodal-pbmc/
system("wget -o ../data/ https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat")

seu_reference <- LoadH5Seurat("../data/pbmc_multimodal.h5seurat", #"/home/projects/cytograted/data/aml/scrna-seq/pbmc_multimodal.h5seurat"
                        assays = "ADT",
                        reductions = FALSE,
                        graphs = FALSE,
                        images = FALSE
)

# Remove raw data
system("rm pbmc_multimodal.h5seurat")

# --------------- cyCombine ---------------

devtools::load_all("../cyCombine/")

SeuratObject::DefaultAssay(seu_reference) <- "ADT"
# seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst")
# markers <- VariableFeatures(seu_obj)
seu_reference$batch <- seu_reference$orig.ident
seu_reference <- normalize_seurat(seu_reference, norm_method = "CLR")
# seu_reference <- normalize_seurat(seu_reference, norm_method = "scale")
seu_reference <- create_som_seurat(seu_reference, markers = VariableFeatures(seu_reference), xdim = 4, ydim = 3, cluster_method = "kohonen", layer = "scale.data")
seu_reference <- correct_data_seurat(seu_reference, mc.cores = 4, pb = T, markers = VariableFeatures(seu_reference))

seu_reference <- Seurat::ScaleData(seu_reference)

seu_reference <- Seurat::RunPCA(seu_reference, features = rownames(seu_reference))
seu_reference <- Seurat::RunUMAP(seu_reference, dims = 1:10) #, features = markers)
# seu_reference <- Seurat::RunUMAP(seu_reference, features = markers)
#
devtools::load_all("../cyDefine/")
colors <- cyDefine:::get_distinct_colors(unique(seu_reference$celltype.l2), FALSE)

Seurat::DimPlot(seu_reference, reduction = "umap", group.by = "batch")
Seurat::DimPlot(seu_reference, reduction = "umap", group.by = "celltype.l2", cols = colors)
Seurat::DimPlot(seu_reference, reduction = "umap", group.by = "Labels")






pbmc_reference <- SeuratObject::LayerData(seu_reference, "data") %>%
  as.matrix() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  as_tibble()

pbmc_metadata <- seu_reference[[]] %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  as_tibble() %>%
  select(cell_id, starts_with("celltype"))

pbmc_reference <- left_join(pbmc_reference, pbmc_metadata, by = "cell_id")



simple_markers <- colnames(pbmc_reference) |>
  stringr::str_remove("-[12]$")
duplicated_markers <- simple_markers[duplicated(simple_markers) & !(simple_markers %in% c("Notch", "TCR", "Rat-IgG1"))]
duplicated_location <- which(simple_markers %in% duplicated_markers)
# duplicated_markers <- c(duplicated_markers, "CD3", "TCR")
pbmc_reference[duplicated_markers] <- NA
pbmc_reference <- pbmc_reference |>
  # rowwise() |>
  mutate(across(all_of(duplicated_markers),
                ~ pmax(
                  get(paste0(cur_column(), "-1")),
                  get(paste0(cur_column(), "-2"))
                ),
                .names = "{col}"))
# mutate(across(all_of(duplicated_markers),
#               ~ rowMeans(select(pbmc_reference,
#                 paste0(cur_column(), "-1"),
#                 paste0(cur_column(), "-2"))
#               ),
#               .names = "{col}"))
pbmc_reference[, duplicated_location] <- NULL

pbmc_reference <- pbmc_reference |>
  # mutate(`CD3-1`= NULL) |>
  rename(
    # "CD3" = "CD3-2",
    "gdTCR" = "TCR-1",
    "abTCR" = "TCR-2"
  )
# pbmc_reference$CD3E <- pbmc_reference$CD3
# dup_ref <- lapply(setNames(duplicated_markers, duplicated_markers), function(marker) {
#   dup_max <- apply(pbmc_reference[, paste0(marker, c("-1", "-2"))], 1, max)
#   return(dup_max)
# })
# pbmc_reference[, duplicated_location] <- NULL
# pbmc_reference <- cbind(pbmc_reference, dup_ref)

usethis::use_data(pbmc_reference, overwrite = TRUE)
# saveRDS(pbmc_reference, file = "~/mnt/cr2/people/s153398/cytograted/data/seurat_reference_250422.RDS")


if (FALSE) { # Legacy


seurat_integrated |>
  Seurat::ScaleData() |>
  Seurat::RunPCA() |>
  Seurat::RunUMAP(dims = 1:10) |>
  Seurat::DimPlot(reduction = "umap", group.by = "orig.ident")


SeuratObject::LayerData(seu_reference, "data")["CD8", ] |>
  hist()
uncorrected$CD8 |> hist()


data <- SeuratObject::LayerData(seu_reference, "counts")[markers, ]
md <- seu_reference[[]] |>
  as_tibble(rownames = "cell_id")
non_markers <- union(cyCombine::non_markers, colnames(md))
uncorrected <- data |>
  t() |>
  as_tibble(rownames = "cell_id") |>
  left_join(md, by = "cell_id") |>
  normalize(norm_method = "CLR")

labels <- uncorrected |>
  normalize(norm_method = "scale") |>
  create_som(ydim = 4, xdim = 4)

corrected <- batch_correct(uncorrected, label = labels)

seu_reference[["cyCombine2"]] <- SeuratObject::CreateAssayObject(data = corrected |> dplyr::select(!any_of(non_markers), cell_id) |> tibble::column_to_rownames("cell_id") |> t() |> as.matrix(), key = "cycombine2_")
SeuratObject::DefaultAssay(seu_reference) <- "cyCombine2"
seu_reference <- Seurat::ScaleData(seu_reference)
seu_reference <- Seurat::RunPCA(seu_reference, features = markers)
seu_reference <- Seurat::RunUMAP(seu_reference, dims = 1:10) #, features = markers)
# seu_reference <- Seurat::RunUMAP(seu_reference, features = markers)
Seurat::DimPlot(seu_reference, reduction = "umap", group.by = "batch")
Seurat::DimPlot(seu_reference, reduction = "umap", group.by = "celltype.l2")
Seurat::DimPlot(seu_reference, reduction = "umap", group.by = "Labels")

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

seurat_integrated <- readRDS("~/mnt/cr2/people/s153398/cytograted/data/seurat_integrated.rds")


# --------------- Make tibble of normalized data and metadata --------------- #

# combine marker expression and cell population label
seurat_reference <- seurat_integrated@assays$integrated@data %>%
  as.matrix() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  as_tibble()

seurat_metadata <- seurat_integrated@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  as_tibble() %>%
  select(cell_id, starts_with("celltype"))

seurat_reference <- left_join(seurat_reference, seurat_metadata, by = "cell_id")
}
## Fix duplicated markers ----

if (FALSE) {
  simple_markers <- colnames(seurat_reference) |>
    stringr::str_remove("-[12]$")
  duplicated_markers <- simple_markers[duplicated(simple_markers) & !(simple_markers %in% c("Notch", "TCR", "Rat-IgG1"))]
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
    # mutate(`CD3-1`= NULL) |>
    rename(
      # "CD3" = "CD3-2",
      "gdTCR" = "TCR-1",
      "abTCR" = "TCR-2"
    )
  seurat_reference$CD3E <- seurat_reference$CD3
  # dup_ref <- lapply(setNames(duplicated_markers, duplicated_markers), function(marker) {
  #   dup_max <- apply(seurat_reference[, paste0(marker, c("-1", "-2"))], 1, max)
  #   return(dup_max)
  # })
  # seurat_reference[, duplicated_location] <- NULL
  # seurat_reference <- cbind(seurat_reference, dup_ref)

  usethis::use_data(seurat_reference, overwrite = TRUE)


  # ------------ Modify naming convention of markers present in two versions ------------ #
} else if (FALSE) {
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
