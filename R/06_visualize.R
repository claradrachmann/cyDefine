#' Get distinct colors for a range of cell types
#'
#' @param populations Character vector of populations to get colors for
#' @param add_unassigned Boolean indicating whether the color 'black' should be
#' added for unassigned cells
#' @return Named list with values HEX colors and populations as names
#' @export
#'
get_distinct_colors <- function(populations, add_unassigned = TRUE) {
  # unassigned cells will be assigned to black
  populations <- populations[populations != "unassigned"]

  colors <- c(
    "#DC050C", "#7BAFDE", "#FDB462", "#33A02C", "#FB8072",
    "#1965B0", "#882E72", "#B2DF8A", "#B17BA6", "#A6761D",
    "#E7298A", "#55A1B1", "#E6AB02", "#7570B3", "#ba5ce3",
    "#FAE174", "#B5651D", "#E78AC3", "#aeae5c", "#FF7F00",
    "#56ff0d", "#AA0A0A", "#BEAED4", "#1e90ff", "#aa8282",
    "#0AC8D4", "#808000", "#7800FA", "#00FAFA", "#641A00",
    "#8DD3C7", "#666666", "#999999", "#d4b7b7", "#8600bf",
    "#00bfff", "#ffff00", "#D4E1C8", "#D470C8", "#64C870",
    "#64C80C", "#0C00FA", "#FA00FA", "#707A00"
  )

  colors <- colors[1:length(populations)]
  names(colors) <- populations

  # add unassigned cells
  if (add_unassigned) {
    colors <- c(colors, "unassigned" = "#000000")
  }

  return(colors)
}

expand_colors <- function(adapted_celltypes, original_celltypes = NULL, colors = NULL) {

  unique_celltypes <- sort(unique(adapted_celltypes))
  if (is.null(colors)) {
    return(get_distinct_colors(unique_celltypes))
  }

  if (!all(unique_celltypes %in% names(colors))) {
    new_celltypes <- unique_celltypes[!unique_celltypes %in% names(colors)]
    cell_counts <- table(original_celltypes)
    colors[new_celltypes] <- sapply(new_celltypes, function(celltype) {
      sub_cells <- unique(original_celltypes[adapted_celltypes == celltype])
      sub_cell_counts <- cell_counts[sub_cells]
      major_celltype <- names(sub_cell_counts[sub_cell_counts == max(sub_cell_counts)])
      colors[major_celltype]
    })
  }
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
plot_umap <- function(reference,
                      query,
                      markers,
                      colors = NULL,
                      query_color_col = "predicted_celltype",
                      build_umap_on = "both",
                      shuffle = TRUE,
                      down_sample = TRUE,
                      sample_n = 10000,
                      return_data = FALSE,
                      seed = 332,
                      verbose = TRUE) {
  # check required packages
  check_package("ggplot2")
  check_package("uwot")
  check_package("patchwork")

  check_colnames(colnames(reference), c("celltype", markers))
  check_colnames(colnames(query), c(query_color_col, markers))

  if (verbose) {
    message("Visualizing reference and query by UMAP")
  }
  set.seed(seed)

  # shuffle data
  if (shuffle) {
    reference <- reference[sample(nrow(reference)), ]
    query <- query[sample(nrow(query)), ]
  }

  if (is.null(colors)) {
    colors <- get_distinct_colors(unique(reference$celltype))
  }

  if (down_sample & (sample_n < nrow(reference) | sample_n < nrow(query))) {
    if (verbose) {
      message(
        "Down-sampling reference and query to ",
        sample_n,
        " cells each"
      )
    }
    reference <- reference %>%
      dplyr::slice_sample(n = min(sample_n, nrow(reference)))
    query <- query %>%
      dplyr::slice_sample(n = min(sample_n, nrow(query)))
  }

  if (build_umap_on == "both") {
    if (verbose) {
      "Computing UMAP embedding of all cells of reference and query"
    }

    # UMAP embedding of both reference and query
    full_umap <- dplyr::bind_rows(
      reference %>%
        dplyr::select(dplyr::all_of(markers)),
      query %>%
        dplyr::select(dplyr::all_of(markers))
    ) %>%
      uwot::umap(
        n_neighbors = 15,
        min_dist = 0.2,
        metric = "euclidean",
        ret_model = FALSE
      )

    ref_umap <- full_umap[1:nrow(reference), ]
    query_umap <- full_umap[(nrow(reference) + 1):nrow(full_umap), ]
  } else if (build_umap_on == "reference") {
    if (verbose) {
      "Computing UMAP embedding only of reference cells. Be aware that this can hide potential novel populations in query!"
    }

    # UMAP embedding of reference
    ref_umap_embed <- reference %>%
      dplyr::select(dplyr::all_of(markers)) %>%
      uwot::umap(
        n_neighbors = 15,
        min_dist = 0.2,
        metric = "euclidean",
        ret_model = TRUE
      )

    ref_umap <- ref_umap_embed$embedding

    if (verbose) {
      "Projecting query cells onto reference UMAP embedding"
    }

    # projection of new data
    query_umap <- query %>%
      dplyr::select(dplyr::all_of(markers)) %>%
      uwot::umap_transform(model = ref_umap_embed)
  }


  if (return_data) {
    colnames(ref_umap) <- c("UMAP1", "UMAP2")
    colnames(query_umap) <- c("UMAP1", "UMAP2")

    return(list(
      "reference" = dplyr::bind_cols(reference, ref_umap),
      "query" = dplyr::bind_cols(query, query_umap)
    ))
  }

  # avoid too long cell type labels
  names(colors) <- stringr::str_wrap(names(colors), 20)

  # get number of columns in legend
  n_legend_cols <- dplyr::case_when(
    length(colors) > 45 ~ 4,
    length(colors) > 30 ~ 3,
    length(colors) > 15 ~ 2,
    TRUE ~ 1
  )

  # plot reference colored by cell type
  ref_embedding <- ref_umap %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(celltype = stringr::str_wrap(reference$celltype, 20)) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = .data$V1,
      y = .data$V2
    )) +
    ggplot2::geom_point(ggplot2::aes(color = .data$celltype),
      size = 0.5,
      alpha = 0.5
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(
        size = 2,
        alpha = 1
      ),
      ncol = n_legend_cols
    )) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Reference") +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2") +
    ggplot2::scale_colour_manual("Cell type",
      values = colors
    )

  # plot query colored by predicted cell type on top of reference in grey
  query_embedding <- ggplot2::ggplot() +
    ggplot2::geom_point(
      ggplot2::aes(
        x = ref_umap[, 1],
        y = ref_umap[, 2]
      ),
      size = 0.5,
      alpha = 1,
      color = "grey87"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        x = query_umap[, 1],
        y = query_umap[, 2],
        color = stringr::str_wrap(query[[query_color_col]], 20)
      ),
      size = 0.5,
      alpha = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Query - predicted cell types") +
    ggplot2::xlab("UMAP1") +
    ggplot2::ylab("UMAP2") +
    ggplot2::scale_colour_manual("Cell type",
      values = colors
    )

  p <- ref_embedding +
    query_embedding + ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    ) +
    patchwork::plot_layout(guides = "collect")

  return(p)
}



#' Make bar plot of predicted cell type abundances
#'
#' @inheritParams plot_umap
#' @param predicted_populations Character vector of (predicted) populations to plot abundance for
#'
#' @return A ggplot2 bar plot of abundances
#' @export
#'
plot_abundance <- function(predicted_populations,
                           colors = NULL,
                           return_data = FALSE,
                           verbose = TRUE) {
  check_package("ggplot2")

  if (verbose) {
    message("Visualizing abundance of predicted cell types")
  }

  if (is.null(colors)) {
    colors <- get_distinct_colors(unique(predicted_populations))
  }

  freqs <- table(predicted_populations) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(
      prop = n / sum(n),
      label = paste0(round(100 * prop, 1), "%"),
      col = colors[predicted_populations]
    )

  if (return_data) {
    return(freqs)
  }

  p <- freqs %>%
    ggplot2::ggplot(ggplot2::aes(
      x = .data$predicted_populations,
      y = .data$prop,
      fill = .data$predicted_populations
    )) +
    ggplot2::geom_bar(
      stat = "identity",
      width = 0.9
    ) +
    ggplot2::geom_text(
      data = freqs,
      ggplot2::aes(
        # x = .data$predicted_populations,
        # y = .data$prop,
        label = paste0(round(100 * .data$prop, 1), "%")
      ),
      position = ggplot2::position_dodge(width = 1),
      vjust = -0.5,
      hjust = 0.5,
      size = 3
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 10
      ),
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank()
    ) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::ggtitle("Cell type abundance") +
    ggplot2::xlab("Predicted cell types") +
    ggplot2::ylab("Proportion of cells")

  return(p)
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
#' @return a pheatmap heatmap
#' @export
#'
plot_heatmap <- function(data,
                         population_col = "predicted_celltype",
                         markers_to_plot,
                         title = "Average marker expression of predicted cell types\n",
                         return_data = FALSE,
                         verbose = TRUE) {
  check_package("pheatmap")
  check_colnames(colnames(data), c(population_col, markers_to_plot))

  if (verbose) {
    message("Visualizing expression of predicted cell types")
  }

  data <- data %>%
    dplyr::group_by_at(population_col) %>%
    dplyr::summarise_at(dplyr::all_of(markers_to_plot), mean) %>%
    dplyr::mutate_at(dplyr::all_of(markers_to_plot), scale)

  plot_dat <- data %>%
    dplyr::select(dplyr::all_of(markers_to_plot)) %>%
    as.matrix() %>%
    t()

  colnames(plot_dat) <- data[[population_col]]

  if (return_data) {
    return(plot_dat)
  }


  p <- pheatmap::pheatmap(
    mat = plot_dat,
    scale = "none",
    main = title,
    cluster_rows = TRUE,
    cluster_cols = TRUE,

    # breaks = my.breaks,
    color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(
      n = 7, name =
        "RdBu"
    )))(100),
    legend = TRUE,
    border_color = "grey",
    treeheight_row = 30,
    treeheight_col = 30,
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 8,
    angle_col = 90
  )

  return(p)
}

#' Plot a Diagram from Merged Cell Types
#'
#' This function generates a diagram that visualizes the merging of cell types based on
#' a provided list or data frame. The diagram illustrates how original cell types (nodes) merge into
#' final cell types.
#'
#' @param input A `list` or `data.frame`. If a `data.frame`, it should have been processed using
#'   `adapt_reference()` and must contain columns `cluster`, `celltype`, and `celltype_original`.
#'   If a `list`, it should be structured such that each key represents a merged cell type and
#'   each value is a vector of original cell types.
#' @param colors A named vector where the names correspond to cell types (both original and merged),
#'   and the values are color codes to be used for filling the nodes in the diagram.
#' @param fontcolor_nodes A named vector (optional) that specifies font colors for specific nodes.
#'   If a node's font color is not provided in this vector, `default_fontcolor` will be used.
#' @param default_fontcolor A `character` string specifying the default font color to be used
#'   for nodes that are not included in `fontcolor_nodes`. Default is `"black"`.
#'
#' @return A rendered diagram produced by `DiagrammeR::grViz()` that visually represents the merging
#'   of cell types. The diagram includes nodes for both original and merged cell types, with edges
#'   showing the relationship between them.
#'
#' @details The function first checks if the input is a data frame or a list. If it is a data frame,
#'   it calls `get_merge_list()` to convert it into the required list structure. The function then
#'   generates a DOT language script to define the diagram, including node styling, edge creation,
#'   and graph layout. Finally, it uses `DiagrammeR::grViz()` to render the flowchart.
#'
#' @examples
#' \dontrun{
#' # Assuming `adapted_reference` is a preprocessed tibble from `adapt_reference()`
#' fontcolor_nodes <- c("NK" = "white", "NK Proliferating / NK" = "white")
#'
#' plot_diagram(adapted_reference, fontcolor_nodes)
#' }
#'
#' @seealso `adapt_reference()`, `get_merge_list()`
#'
#' @export
plot_diagram <- function(input, colors = NULL, fontcolor_nodes = NULL, default_fontcolor = "black") {

  cyDefine:::check_package("DiagrammeR")

  if ("data.frame" %in% class(input)) {
    merge_list <- get_merge_list(input)
  } else if (class(input) == "list") {
    merge_list <- input
  }


  if (all(c("celltype", "celltype_original") %in% names(input))) {
    colors <- expand_colors(input$celltype, input$celltype_original, colors)
  } else if (is.null(colors)) {
    colors <- cyDefine::get_distinct_colors(unique(unlist(merge_list)))
  } else if (!all(names(merge_list) %in% names(colors))) {
    colors[names(merge_list)] <- sapply(names(merge_list), function(merged) {
      colors[stringr::str_split(merged, " /")[[1]][1]]
    })
  }

  rm(input)
  # Initialize the diagram
  diagram_code <- "
    digraph flowchart {
      # Global node styles
      node [fontname = Helvetica, style = filled, color = black, shape = box];
  "

  # To track the incrementing value for node names
  node_counter <- 0

  # Function to get a unique node name by appending an incrementing value
  get_unique_node_name <- function(node_name) {
    node_counter <<- node_counter + 1
    return(paste0(node_name, "_", node_counter))
  }

  # Loop through the merge list to create nodes and edges
  for (merge_into in names(merge_list)) {
    merged_nodes <- merge_list[[merge_into]]

    # Set the font color based on the input vector or default
    fontcolor <- ifelse(merge_into %in% names(fontcolor_nodes),
                        fontcolor_nodes[merge_into], default_fontcolor)

    # Create a unique identifier for the merge_into node
    merge_into_unique <- get_unique_node_name(merge_into)

    # Define the merged node
    diagram_code <- paste0(diagram_code, "
      \"", merge_into_unique, "\" [shape=ellipse, fillcolor='", colors[merge_into], "', fontcolor='", fontcolor, "', label = \"", merge_into, "\"];
    ")

    # Define the individual nodes and edges
    for (node in merged_nodes) {
      fontcolor <- ifelse(node %in% names(fontcolor_nodes),
                          fontcolor_nodes[node], default_fontcolor)

      # Create a unique identifier for the node
      node_unique <- get_unique_node_name(node)

      diagram_code <- paste0(diagram_code, "
        \"", node_unique, "\" [shape=box, fillcolor='", colors[node], "', fontcolor='", fontcolor, "', label = \"", node, "\"];
        \"", node_unique, "\" -> \"", merge_into_unique, "\";
      ")
    }
  }

  # Close the diagram definition
  diagram_code <- paste0(diagram_code, "
      # Set global graph attributes
      graph [layout = dot, rankdir = TB];
    }
  ")

  # Render the diagram
  DiagrammeR::grViz(diagram_code)
}

#' Get Merge List for create_diagram
#'
#' This function generates a list mapping merged cell types to their original cell types
#' based on the clustering results from an adapted reference dataset. The function assumes
#' that the dataset has been preprocessed using `adapt_reference()`.
#'
#' @param adapted_reference A `data.frame` or `tibble` that contains the adapted reference data.
#'   The data should include the columns `cluster`, `celltype`, and `celltype_original`.
#'   These columns are expected to be present in the dataset after running `adapt_reference()`.
#'
#' @return A named list where each element represents a merged cell type. The names of the list
#'   correspond to the unique cell types in the `celltype` column, and the values are lists of
#'   the original cell types (`celltype_original`) that were merged into the final cell type.
#'
#' @details The function filters the provided reference data to include only unique combinations
#'   of `cluster`, `celltype`, and `celltype_original`. It then groups the data by `cluster` and
#'   creates a list where each unique `celltype` is associated with its corresponding original
#'   cell types. The list is returned in a format suitable for use in generating flowcharts or
#'   diagrams that visualize the merging of cell types.
#'
#' @examples
#' \dontrun{
#' # Assuming `adapted_reference` is a preprocessed tibble from `adapt_reference()`
#' merge_list <- get_merge_list(adapted_reference)
#' print(merge_list)
#' }
#'
get_merge_list <- function(adapted_reference) {
  stopifnot(
    "Please run adapt_reference() before creating a diagram." =
      all(c("cluster", "celltype", "celltype_original") %in%
            colnames(adapted_reference)))

  merge_list <- adapted_reference |>
    dplyr::distinct(cluster, celltype, celltype_original) |>
    dplyr::filter(!is.na(cluster)) |>
    dplyr::group_by(cluster) |>
    dplyr::reframe(
      celltype = unique(celltype),
      celltype_original = list(celltype_original)) %>%
    dplyr::select(-cluster) |>
    tibble::deframe()
  return(merge_list)
}

#' Create an Alluvial Plot
#'
#' This function creates an alluvial plot to compare true and predicted clusterings.
#'
#' @param df A data frame containing the classified data.
#' @param true A character string indicating the column name for true labels.
#' @param predicted A character string indicating the column name for predicted labels.
#' @param color A named vector of colors for the true labels.
#' @param color_by A character string indicating which variable to use for filling the alluvia ("true" or "predicted").
#' @param n An integer indicating the sample size.
#' @param seed An integer indicating the seed value for random sampling.
#'
#' @return A ggplot object representing the alluvial plot.
#'
#' @examples
#' # Example usage
#' plot_alluvium(df = classified, true = "celltype", predicted = "model_prediction",
#'               color = celltype_colors, color_by = "true", n = 500, seed = 1286)
#'
#' @importFrom dplyr sample_n
#' @import ggplot2
#' @import ggalluvial
#'
#' @export
plot_alluvium <- function(df, true = "celltype", predicted = "model_prediction", color = NULL, color_by = c("true", "predicted"), n = 500, seed = 1286) {
  color_by <- match.arg(color_by)
  color_by <- if (color_by == "true") true else predicted
  stopifnot("Please set 'true' as a column in the input data." = true %in% colnames(df))
  stopifnot("Please set 'predicted' as a column in the input data." = predicted %in% colnames(df))
  if (is(color, "NULL")) {
    color <- cyDefine::get_distinct_colors(unique(df[[color_by]]), add_unassigned = F)
  }

  # Set the seed for reproducibility
  set.seed(seed)

  if (n < nrow(df)) {
    df <- dplyr::sample_n(df, size = n)
  }

  # Create the alluvial plot
  plot <- ggplot(df, aes_string(axis1 = true, axis2 = predicted)) +
    geom_alluvium(aes_string(fill = color_by), alpha = .7, reverse = TRUE) +
    geom_stratum(width = 1/4) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("True", "Predicted"), expand = c(0.15, 0.05)) +
    scale_fill_manual(values = sort(color)) +
    theme_void() +
    theme(
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14, face = "bold"),
      legend.position = "none"
    ) +
    ggtitle("True vs Predicted")

  return(plot)
}
