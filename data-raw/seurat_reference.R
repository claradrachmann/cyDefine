library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(ggplot2)

# multi.h5seurat is obtained from https://atlas.fredhutch.org/nygc/multimodal-pbmc/
seu_obj <- LoadH5Seurat("multi.h5seurat",
                        assays = "ADT",
                        reductions = FALSE,
                        graphs = FALSE,
                        images = FALSE)



# --------------- integration --------------- #

# split by orig.ident which is sample of origin (patient + time point)
seu_obj_list <- SplitObject(seu_obj, split.by = "orig.ident")

seu_obj_list <- lapply(X = seu_obj_list,
                       FUN = function(x) {
                         # CLR transform per cell
                         x <- NormalizeData(x, normalization.method = "CLR", margin = 2)
                         x <- FindVariableFeatures(x,  selection.method = "vst")
                       })

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seu_obj_list)
anchors <- FindIntegrationAnchors(object.list = seu_obj_list,
                                  anchor.features = features)

# integrate data
seurat_integrated <- IntegrateData(anchorset = anchors)


saveRDS(seurat_integrated, file = "/home/projects/cytograted/data/seurat_integrated.rds")

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




# ------------ Modify naming convention of markers present in two versions ------------ #

# check markers with several versions
multi_version_markers <- c("CD3", "CD4", "CD56", "CD11b", "CD38", "CD44", "CD26", "CD275", "CD45", "CD133", "CD138", "TCR")
pccs <- as_tibble(matrix(ncol = length(multi_version_markers), nrow = 1))
colnames(pccs) <- multi_version_markers

# compute correlation between marker versions
for (mar in multi_version_markers) {
  pccs[[mar]] <- cor(seurat_reference[[str_c(mar, "-1")]],
                     seurat_reference[[str_c(mar, "-2")]])
}

# add columns with mean for multi-version markers with high PCC
for (mar in multi_version_markers[pccs > 0.5]) {
  seurat_reference <- seurat_reference %>%
    mutate("{mar}" := rowMeans(select(seurat_reference,
                                      str_c(mar, "-1"),
                                      str_c(mar, "-2"))))
}

# for rest of multi-version markers, use the marker seeming most informative
# (evaluated manually)
seurat_reference$CD26 <- seurat_reference$`CD26-2`
seurat_reference$CD275 <- seurat_reference$`CD275-2`
seurat_reference$CD45 <- seurat_reference$`CD45-2`
seurat_reference$CD133 <- seurat_reference$`CD133-1`
seurat_reference$CD138 <- seurat_reference$`CD138-1`

# also rename CD3-1 and CD3-2 + TCR-1 and TCR-2 to more informative names
seurat_reference <- seurat_reference %>%
  rename("CD3E" = "CD3-1",
         "CD3D" = "CD3-2",
         "gdTCR" = "TCR-1",
         "abTCR" = "TCR-2")
# %>%
  # select(-cell_id)



usethis::use_data(seurat_reference, overwrite = TRUE)
