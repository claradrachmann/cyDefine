library(tidyverse)
devtools::load_all()

tof_downsample_density <-
  function(
    tof_tibble,
    group_cols = NULL,
    density_cols = where(tof_is_numeric),
    target_num_cells,
    target_prop_cells,
    target_percentile = 0.03,
    outlier_percentile = 0.01,
    distance_function = c("euclidean", "cosine"),
    density_estimation_method = c("mean_distance", "sum_distance", "spade"),
    ...
  ) {

    # check distance_function
    distance_function <- rlang::arg_match(distance_function)

    #initialize result
    result <-
      tof_tibble %>%
      dplyr::mutate(..cell_id = 1:nrow(tof_tibble))

    # extract group and density column names
    group_names <-
      tidyselect::eval_select(
        expr = rlang::enquo(group_cols),
        data = tof_tibble
      ) %>%
      names()
    print(group_names)

    density_names <-
      tidyselect::eval_select(
        expr = rlang::enquo(density_cols),
        data = tof_tibble
      ) %>%
      names()
    print(density_names)


    # nest data needed to compute densities for each group
    nested_data <-
      result %>%
      dplyr::select(
        .data$..cell_id,
        dplyr::any_of(group_names),
        dplyr::any_of(density_names)
      ) %>%
      tidyr::nest(cell_ids = .data$..cell_id, data = {{density_cols}})

    # find local density estimates for each cell in all groups
    nested_data <-
      nested_data %>%
      dplyr::mutate(
        densities =
          purrr::map(
            .x = .data$data,
            .f = tidytof::tof_estimate_density,
            method = density_estimation_method,
            augment = FALSE,
            ...
          )
      )

    # save the local density estimates as one vector per group
    densities <- purrr::map(.x = nested_data$densities, .f = ~.x[[1]])
    nested_data <-
      nested_data %>%
      dplyr::select(.data$cell_ids, dplyr::any_of(group_names)) %>%
      dplyr::mutate(
        densities = densities
      )

    # if target_num_cells and target_prop_cells are not specified, directly use
    # target_percentile to do the downsampling
    if (missing(target_prop_cells) & missing(target_num_cells)) {

      # store a vector of ..cell_ids for which cells to keep from tof_tibble
      chosen_cells <-
        nested_data %>%
        dplyr::mutate(
          target_density =
            purrr::map_dbl(
              .x = densities,
              # protect against returning 0 cells when there are many degenerate
              # densities by computing the target percentile after removing
              # the cells with 0 density (which will be removed during the
              # final filtering step anyway)
              .f = ~ quantile(subset(.x, .x > 0), probs = target_percentile),
            )
        ) %>%
        tidyr::unnest(cols = c(.data$cell_ids, .data$densities)) %>%
        dplyr::group_by(dplyr::across({{group_cols}})) %>%
        dplyr::arrange(.data$densities) %>%
        dplyr::mutate(
          rank = 1:dplyr::n(),
          percentile = rank / dplyr::n()
        ) %>%
        dplyr::mutate(
          sample_prob = (.data$target_density / .data$densities)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(sample_value = stats::runif(n = dplyr::n())) %>%
        dplyr::filter(
          .data$percentile > outlier_percentile,
          .data$sample_prob > .data$sample_value
        ) %>%
        dplyr::pull(.data$..cell_id)

      # if either target_num_cells or target_prop_cells are specified, find
      # a threshold that approximates the requested number or proportion of
      # cells
    } else {
      # if target_num_cells is not constant, compute it for each group
      # using target_prop_cells
      if (missing(target_num_cells)) {
        target_num_cells_vector <-
          purrr::map_int(
            .x = nested_data$cell_ids,
            .f =
              ~ as.integer(ceiling(nrow(.x) * target_prop_cells))
          )
        # if target_num_cells is specified, shoot for that constant number
        # in each group
      } else {
        target_num_cells_vector <-
          rep(x = target_num_cells, times = nrow(nested_data))
      }

      chosen_cells <-
        nested_data %>%
        dplyr::mutate(
          num_cells = target_num_cells_vector,
          # remove cells at a lower percentile than the outlier percentile
          # and also any cells whose local densities are exactly 0
          cells_to_remove =
            purrr::map(
              .x = .data$densities,
              .f = ~
                dplyr::tibble(
                  densities = .x,
                  ..cell_ids = 1:length(.x)
                ) %>%
                dplyr::arrange(.data$densities) %>%
                dplyr::mutate(
                  rank = 1:dplyr::n(),
                  percentile = rank / dplyr::n(),
                ) %>%
                dplyr::filter(
                  .data$percentile <= outlier_percentile | .data$densities == 0
                ) %>%
                dplyr::pull(.data$..cell_ids)
            ),
          densities =
            purrr::map2(
              .x = .data$densities,
              .y = .data$cells_to_remove,
              .f = ~ .x[-.y]
            ),
          cell_ids =
            purrr::map2(
              .x = .data$cell_ids,
              .y = .data$cells_to_remove,
              .f = ~ .x[-.y, ]
            )
        ) %>%
        # use tof_spade_downsampling to find which cells should be retained
        # after downsampling
        dplyr::transmute(
          sampled_cells =
            purrr::pmap(
              .l = list(.data$densities, .data$cell_ids, .data$num_cells),
              .f = tof_spade_downsampling
            )
        ) %>%
        dplyr::pull(.data$sampled_cells) %>%
        c(recursive = TRUE)

    }

    # filter only selected cells out of the original tof_tibble
    result <-
      result %>%
      dplyr::filter(.data$..cell_id %in% chosen_cells) %>%
      dplyr::select(-.data$..cell_id)

    return(result)
  }


# Density dependent down-sampling via tidytof
# performed as proposed in https://pubmed.ncbi.nlm.nih.gov/21964415/

seurat_down <- tidytof::tof_downsample_density(
  tof_tibble = seurat_reference,
  group_cols = celltype.l2,
  density_cols = all_of(seurat_markers),
  target_num_cells = 1000,
  # The local density percentile (i.e. a value between 0 and 1)
  # below which cells should be considered outliers (and discarded)
  outlier_percentile = 0.01,
  distance_function = "euclidean",
  density_estimation_method = "spade"
)













-KnnDensity <- function(k, min, max, n, nn.ids.df, nn.dists.df,
                       numcluster = numcluster, #table.lengths
                       table.breaks, offset) {
  all.densities.list <- list()
  for(i in 1:nrow(nn.ids.df)) {
    ## Alternatives for calculating density
    ## Dist to kth nearest neighbor:
    #density.by.knn <- abs(nn.dists.df[i,k])
    ## Transformed metric from X.shift paper:
    #density.by.Xshift <- (1/(n*(longest.dist^d))) * (sum(c(1:k)^d)/sum(compile.dists))^d
    ## Sum of distances from 1:kth nearest neighbor:
    #density.by.knn <- sum(abs(nn.dists.df[i,1:k]))
    ## Mean of distances from 1:kth nearest neighbor:
    density.by.knn <- sum(abs(nn.dists.df[i,1:k]))/k
    all.densities.list[[paste(i,".nn", sep='')]] <- density.by.knn
  }
  densities.df <- rbind(matrix(unlist(all.densities.list), byrow=T))
  normalized.densities <- BBmisc::normalize(densities.df, method = "range", range = c(min, max), margin = 2, on.constant = "quiet")
  if (offset == 0) {
    names(normalized.densities) <- (1:numcluster)
  } else {
    names(normalized.densities) <- offset:table.breaks[n + 2]
  }
  return(normalized.densities)
}




# density <- seurat_reference %>%
#   spade::SPADE.density()


nns <- RANN::nn2(data = seurat_reference[1:1000, seurat_markers], searchtype="priority", eps=0.1)
temp_nnids.df <- as.data.frame(nns$nn.idx)
temp_nndists.df <- as.data.frame(nns$nn.dists)
nn.ids.df <- temp_nnids.df[,2:length(temp_nnids.df)]
nn.dists.df <- temp_nndists.df[,2:length(temp_nndists.df)]
# numcluster <- nrow(clusters)
#TODO What to set for k?
density <- KnnDensity(k=15, min = 1, max = 100, n=0, nn.ids.df = nn.ids.df,
                      nn.dists.df = nn.dists.df)



if (max(density) == 0.0)
  warning(paste(infilename,"has degenerate densities, possibly due to many identical observations",sep=" "))

boundary <- quantile(density, 0.01, names = FALSE)


