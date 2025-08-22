library(dplyr)
library(ggplot2)

# Load all F1 score files recursively from data folder
f1_files <- list.files("data", pattern = "_F1\\.rds$", recursive = TRUE, full.names = TRUE)
f1_files <- f1_files[stringr::str_detect(f1_files, "unassigned")]
# Create a list to store results
f1_results <- list()

# Load each file and extract information from filename
for (file in f1_files) {
  # Read the RDS file
  f1_data <- readRDS(file)
  if (length(f1_data) == 0) {
    message(file)
    next}

  # Extract basename without extension
  basename <- tools::file_path_sans_ext(basename(file))

  # Parse filename to extract dataset, method, and other info
  # Assuming format: dataset_method_suffix_F1.rds
  parts <- strsplit(basename, "_")[[1]]

  # Remove "F1" from the end
  parts <- parts[-length(parts)]
  print(parts)
  # Extract dataset (first part) and method (combine remaining parts)
  dataset <- parts[1]
  method <- parts[2]
  unassigned <- paste(parts[-1:-2], collapse = "_")

  # Store results with metadata
  f1_results[[basename]] <- data.frame(
    celltype = names(f1_data),
    f1_score = f1_data,
    dataset = dataset,
    method = method,
    unassigned = unassigned,
    file = basename, row.names = NULL
  )
}

# Combine all results into one dataframe
df_f1 <- dplyr::bind_rows(f1_results)

# Clean up method names to match the plot
df_f1 <- df_f1 |>
  mutate(
    method = case_when(
      str_detect(method, "CLC") ~ "CyTOF Linear\nClassifier",
      TRUE ~ method
    ),
    method = factor(method, levels = c("CyTOF Linear\nClassifier", "Spectre", "CyAnno", "cyDefine"))
  )

# Figure 3c Boxplot ----
fig3c <- df_f1 |>
  # dplyr::filter(stringr::str_detect(file, "wo_unassigned", negate = TRUE)) |>
  ggplot(aes(x = dataset, y = f1_score, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               outlier.size = 0,
               outlier.alpha = 0.7) +
  facet_wrap(~unassigned) +
  geom_point(shape = 21, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c(
    "cyDefine" = "#A8D8A8",
    "CyTOF Linear\nClassifier" = "#F4A460",
    "CyAnno" = "#FA8072",
    "Spectre" = "#87CEEB"
  )) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.25)) +
  labs(
    title = "F1 score per cell type\n(including unassigned cells)",
    x = "Dataset",
    y = "F1 score per cell type",
    fill = ""
  ) +
  # theme_classic() +
  ggpubr::theme_pubr() +
  theme(
    plot.title = element_text(hjust = 0, size = 11),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right",
    panel.grid.major.y = element_line(color = "grey90", size = 0.3),
    panel.grid.minor.y = element_line(color = "grey95", size = 0.2)
  )

fig3c


# Heatmap ----
# Calculate average F1 scores by dataset and method
df_heatmap <- df_f1 |>
  group_by(dataset, method, unassigned) |>
  summarise(
    avg_f1 = mean(f1_score, na.rm = TRUE),
    .groups = "drop"
  ) |>
  group_by(dataset, unassigned) |>
  mutate(
    is_max = avg_f1 == max(avg_f1, na.rm = TRUE),
    # Round to 3 decimal places for display
    avg_f1 = round(avg_f1, 3),
    # Create text labels
    label = ifelse(avg_f1 >= 0.9,
                   paste0(avg_f1, "†"),
                   as.character(avg_f1)),
    # Add asterisk for missing data (if needed)
    label = ifelse(is.na(avg_f1), "*", label),
    is_max = ifelse(is.na(is_max), FALSE, is_max),
  )

# Create heatmaps for both conditions
create_heatmap <- function(data, condition_filter, title) {
  data |>
    filter(unassigned == condition_filter) |>
    ggplot(aes(x = method, y = dataset, fill = avg_f1)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = ifelse(.data$is_max, .data$label, ""), fontface = "bold"), size = 3.5) +
    geom_text(aes(label = ifelse(!.data$is_max, .data$label, "")), size = 3.5) +
    # geom_text(aes(label = label),
    #           color = "black",
    #           size = 3.5,
    #           fontface = "bold") +
    # scale_fill_gradient(low = "#e3f5e1", high = "#30a719", na.value = "white",
    scale_fill_gradient2(
      low = "#F5F5F5",
      mid = "lightgreen",
      high = "darkgreen",
      na.value = "white",
    midpoint = 0.4,
      limits = c(0.2, 1.0),
      name = "Average\nF1 score",
      breaks = seq(0.4, 1.0, 0.1)
    ) +
    labs(
      title = title,
      x = NULL,
      y = "Dataset"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11, hjust = 0),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(hjust = 1),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.4, "cm")
    ) +
    coord_fixed()
}

# Create heatmap for including unassigned cells
fig_heatmap_a <- create_heatmap(
  df_heatmap,
  "w_unassigned",
  "Average F1 score\n(including unassigned cells)"
)

# Create heatmap for excluding unassigned cells
fig_heatmap_b <- create_heatmap(
  df_heatmap,
  "wo_unassigned",
  "Average F1 score\n(excluding unassigned cells)"
)

print(fig_heatmap_a)
print(fig_heatmap_b)

# Combine heatmaps side by side
library(patchwork)
combined_heatmap <- fig_heatmap_a + fig_heatmap_b +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a")

print(combined_heatmap)

# Save plots
ggsave("f1_boxplot.png", fig3c, width = 12, height = 6, dpi = 300)
ggsave("f1_heatmap_combined.png", combined_heatmap, width = 14, height = 6, dpi = 300)
ggsave("f1_heatmap_with_unassigned.png", fig_heatmap_a, width = 7, height = 5, dpi = 300)
ggsave("f1_heatmap_without_unassigned.png", fig_heatmap_b, width = 7, height = 5, dpi = 300)
