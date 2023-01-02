
#' Get distinct colors for a range of cell types
#'
#' @param populations Character vector of populations to get colors for
#' @param add_unassigned Boolean indicating whether the color 'black' should be
#' added for unassigned cells
#' @return Named list with values HEX colors and populations as names
#' @export
#'
#' @examples
get_distinct_colors <- function(populations, add_unassigned = TRUE) {

  # unassigned cells will be assigned to black
  populations <- populations[populations != "unassigned"]

  colors <- c("#DC050C", "#7BAFDE", "#FDB462","#33A02C", "#FB8072",
              "#1965B0", "#882E72", "#B2DF8A", "#B17BA6", "#A6761D",
              "#E7298A", "#55A1B1", "#E6AB02", "#7570B3", "#ba5ce3",
              "#FAE174", "#B5651D", "#E78AC3", "#aeae5c", "#FF7F00",
              "#56ff0d", "#0C0C91", "#BEAED4", "#1e90ff", "#aa8282",
              "#0AC8D4", "#808000", "#7800FA", "#00FAFA", "#641400",
              "#8DD3C7", "#666666", "#999999", "#d4b7b7", "#8600bf",
              "#00bfff", "#ffff00", "#D4E1C8", "#D470C8", "#64C870",
              "#64C80C", "#0C00FA", "#FA00FA", "#707A00")

  colors <- colors[1:length(populations)]
  names(colors) <- populations

  # add unassigned cells
  if (add_unassigned) {colors <- c(colors, "unassigned" = "#000000")}

  return(colors)
}




#' Visualize reference and query by UMAPs
#'
#' @param reference Tibble of reference data (cells in rows, markers in columns)
#' @param query Tibble of query data (cells in rows, markers in columns)
#' @param markers Character vector of available markers
#' @param colors Optional: Named list of colors per cell type. Values should be
#' HEX colors and names should be the unique cell types
#' @param build_umap_on Character describing what data should be applied for
#' performing the UMAP embedding. Options are "both" (default), "reference",
#' and "query". "Both" implies that the UMAP embedding is done on both the
#' reference and query data.
#' @param shuffle Boolean indicating whether rows should be shuffled before plotting
#' @param down_sample Boolean indicating whether cells should be down-sampled
#' @param sample_n Number of cells to sample from each of reference and query,
#' if down_sample = TRUE
#' @param return_data Boolean indicating if the data used for plotting should be
#' returned rather than a ggplot2 object
#' @param seed Random seed
#' @param verbose Verbosity
#'
#' @return ggplot2 plot
#' @export
#'
#' @examples
plot_umap <- function(reference,
                      query,
                      markers,
                      colors = NULL,
                      build_umap_on = "both",
                      shuffle = TRUE,
                      down_sample = TRUE,
                      sample_n = 100000,
                      return_data = FALSE,
                      seed = 332,
                      verbose = TRUE) {

  # check required packages
  check_package("ggplot2")
  check_package("uwot")
  check_package("patchwork")

  check_colnames(colnames(reference), c("celltype", markers))
  check_colnames(colnames(query), c("predicted_celltype", markers))

  if (verbose) {message("Visualizing reference and query by UMAP")}
  set.seed(seed)

  # shuffle data
  if (shuffle) {
    reference <- reference[sample(nrow(reference)),]
    query <- query[sample(nrow(query)),]
  }

  if (is.null(colors)) {colors <- get_distinct_colors(unique(reference$celltype))}

  if (down_sample & (sample_n < nrow(reference) | sample_n < nrow(query))) {
    if (verbose) {
      message("Down-sampling reference and query to ",
              sample_n,
              " cells each")
    }
    reference <- reference %>%
      dplyr::slice_sample(n = min(sample_n, nrow(reference)))
    query <- query %>%
      dplyr::slice_sample(n = min(sample_n, nrow(reference)))
  }

  if (build_umap_on == "both") {

    if (verbose) {"Computing UMAP embedding of all cells of reference and query"}

    # UMAP embedding of both reference and query
    full_umap <- dplyr::bind_rows(reference %>%
                                    dplyr::select(dplyr::all_of(markers)),
                                  query %>%
                                    dplyr::select(dplyr::all_of(markers))) %>%
      uwot::umap(n_neighbors = 15,
                 min_dist = 0.2,
                 metric = "euclidean",
                 ret_model = FALSE)

    ref_umap <- full_umap[1:nrow(reference),]
    query_umap <- full_umap[(nrow(reference)+1):nrow(full_umap),]

  }

  else if (build_umap_on == "reference") {

    if (verbose) {"Computing UMAP embedding only of reference cells. Be aware that this can hide potential novel populations in query!"}

    # UMAP embedding of reference
    ref_umap_embed <- reference %>%
      dplyr::select(dplyr::all_of(markers)) %>%
      uwot::umap(n_neighbors = 15,
                 min_dist = 0.2,
                 metric = "euclidean",
                 ret_model = TRUE)

    ref_umap <- ref_umap_embed$embedding

    if (verbose) {"Projecting query cells onto reference UMAP embedding"}

    # projection of new data
    query_umap <- query %>%
      dplyr::select(dplyr::all_of(markers)) %>%
      uwot::umap_transform(model = ref_umap_embed)
  }


  if (return_data) {return (list("reference" = ref_umap,
                                 "query" = query_umap))}

  # avoid too long cell type labels
  names(colors) <- stringr::str_wrap(names(colors), 20)

  # get number of columns in legend
  n_legend_cols <- dplyr::case_when(length(colors) > 45 ~ 4,
                                    length(colors) > 30 ~ 3,
                                    length(colors) > 15 ~ 2,
                                    TRUE ~ 1)

  # plot reference colored by cell type
  ref_embedding <- ref_umap %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(celltype = stringr::str_wrap(reference$celltype, 20)) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$V1,
                                 y = .data$V2)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$celltype),
                        size = 0.5,
                        alpha = 0.5) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2,
                                                                      alpha = 1),
                                                  ncol = n_legend_cols)) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Reference") +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2") +
    ggplot2::scale_colour_manual("Cell type",
                                 values = colors)

  # plot query colored by predicted cell type on top of reference in grey
  query_embedding <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = ref_umap[,1],
                                     y = ref_umap[,2]),
                        size = 0.5,
                        alpha = 1,
                        color = "grey87") +
    ggplot2::geom_point(ggplot2::aes(x = query_umap[,1],
                                     y = query_umap[,2],
                                     color = stringr::str_wrap(query$predicted_celltype, 20)),
                        size = 0.5,
                        alpha = 0.5,
                        show.legend = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Query - predicted cell types") +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2") +
    ggplot2::scale_colour_manual("Cell type",
                                 values = colors)

  p <- ref_embedding +
    query_embedding + ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                                     axis.text.y = ggplot2::element_blank(),
                                     axis.ticks.y = ggplot2::element_blank()) +
    patchwork::plot_layout(guides = 'collect')

  return (p)
}



#' Make bar plot of predicted cell type abundances
#'
#' @inheritParams plot_umap
#' @param predicted_populations Character vector of (predicted) populations to plot abundance for
#'
#' @return
#' @export
#'
#' @examples
plot_abundance <- function(predicted_populations,
                           colors = NULL,
                           return_data = FALSE,
                           verbose = TRUE) {

  check_package("ggplot2")

  if (verbose) {message("Visualizing abundance of predicted cell types")}

  if (is.null(colors)) {colors <- get_distinct_colors(unique(predicted_populations))}

  freqs <- table(predicted_populations) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(prop = n/sum(n),
                  label = paste0(round(100*prop,1), "%"),
                  col = colors[predicted_populations])

  if (return_data) {return(freqs)}

  p <- freqs %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$predicted_populations,
                                 y = .data$prop,
                                 fill = .data$predicted_populations)) +
    ggplot2::geom_bar(stat = "identity",
                      width = 0.9) +
    ggplot2::geom_text(data = freqs,
                       ggplot2::aes(
                         # x = .data$predicted_populations,
                         # y = .data$prop,
                         label = paste0(round(100*.data$prop, 1), "%")),
                       position = ggplot2::position_dodge(width = 1),
                       vjust = -0.5,
                       hjust = 0.5,
                       size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                       vjust = 0.5,
                                                       hjust = 1,
                                                       size = 10),
                   legend.position = "none",
                   panel.grid.major.x = ggplot2::element_blank()) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::ggtitle("Cell type abundance") +
    ggplot2::xlab("Predicted cell types") +
    ggplot2::ylab("Proportion of cells")

  return (p)
}



#' Plot heatmap of marker expressions per population
#'
#' @inheritParams plot_umap
#' @param data Tibble with data (cells in rows, markers in columns)
#' @param population_col Name of the cell type column in 'data' that you want to
#' plot expressions for
#' @param markers_to_plot Character vector of the markers that you want to
#' include in the heatmap
#' @param title Title of the plot
#'
#' @return
#' @export
#'
#' @examples

plot_heatmap <- function(data,
                         population_col = "predicted_celltype",
                         markers_to_plot,
                         title = "Average marker expression of predicted cell types\n",
                         return_data = FALSE,
                         verbose = TRUE) {

  check_package("pheatmap")
  check_colnames(colnames(data), c(population_col, markers_to_plot))

  if (verbose) {message("Visualizing expression of predicted cell types")}

  data <- data %>%
    dplyr::group_by_at(population_col) %>%
    dplyr::summarise_at(dplyr::all_of(markers_to_plot), mean) %>%
    dplyr::mutate_at(dplyr::all_of(markers_to_plot), scale)

  plot_dat <- data %>%
    dplyr::select(dplyr::all_of(markers_to_plot)) %>%
    as.matrix() %>%
    t()

  colnames(plot_dat) <- data[[population_col]]

  if (return_data) {return(plot_dat)}


  p <- pheatmap::pheatmap(mat = plot_dat,
                          scale = "none",
                          main = title,
                          cluster_rows = TRUE,
                          cluster_cols = TRUE,

                          #breaks = my.breaks,
                          color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name =
                                                                    "RdBu")))(100),
                          legend = TRUE,
                          border_color = "grey",
                          treeheight_row = 30,
                          treeheight_col = 30,
                          fontsize = 10,
                          fontsize_row = 8,
                          fontsize_col = 8,
                          angle_col = 90)

  return (p)
}





