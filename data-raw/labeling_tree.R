
# nested list of hierarchical cell type levels
labeling_tree <- list(
  "Myeloid cell" = list("Monocyte" = c("CD14 Mono",
                                       "CD16 Mono"),
                        "DC" = list("pDC" = c("ASDC",
                                              "pDC"),
                                    "cDC" = c("cDC1",
                                              "cDC2"))),
  "Lymphoid cell" = list("B cell" = c("B naive",
                                      "B intermediate",
                                      "B memory",
                                      "Plasmablast"),
                         "T cell" = list("abT" = list("CD4 T cell" = list("CD4 Naive",
                                                                          "CD4 CTL",
                                                                          "CD4 Proliferating",
                                                                          "Treg",
                                                                          "CD4 Memory" = c("CD4 TCM",
                                                                                           "CD4 TEM")),
                                                      "CD8 T cell" = list("CD8 Naive",
                                                                          "CD8 Proliferating",
                                                                          "CD8 Memory" = c("CD8 TCM",
                                                                                           "CD8 TEM")),
                                                      "dnT",
                                                      "MAIT"),
                                         "gdT"),
                         "ILC" = list("NK" = c("NK",
                                               "NK_CD56bright",
                                               "NK Proliferating"),
                                      "ILC")))


labeling_tree <- reshape2::melt(labeling_tree,
                                value.name = "leaf") %>%
  dplyr::as_tibble()  %>%
  dplyr::na_if("")


usethis::use_data(labeling_tree, overwrite = TRUE)
