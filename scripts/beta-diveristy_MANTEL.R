##### Community Composition: NMDS and Mantel Tests ##########################

library(tidyverse)
library(vegan)
library(patchwork)

set.seed(12345)

data_dir <- "data"
output_dir <- file.path("outputs", "community_composition_ranch_tests")
figure_dir <- file.path("figures", "community_composition_ranch_tests")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- read_csv(file.path(data_dir, "1.1_Metadata.csv"), show_col_types = FALSE)

required_metadata <- c(
  "Sample", "Core_ID", "Location", "Practice", "depth", "depth_num",
  "pH", "total_c", "maom_c", "pom_c", "total_n", "maom_n", "pom_n",
  "percent_clay", "percent_sand", "percent_silt", "Adaptive_Score",
  "mean_annual_temp", "mean_annual_precip", "PET_mm_year", "radiation"
)
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata) > 0) {
  stop("Missing metadata columns: ", paste(missing_metadata, collapse = ", "))
}

##### 1. Read and Align ASV Tables
read_asv_table <- function(path, marker, minimum_reads = 1) {
  table <- read_csv(path, show_col_types = FALSE)
  id_column <- intersect(c("ASV", "OTU"), names(table))[1]
  if (is.na(id_column)) stop("No ASV or OTU column found in ", path)
  
  asv_ids <- table[[id_column]]
  counts <- table %>%
    dplyr::select(-all_of(c(id_column, intersect("Taxon", names(table))))) %>%
    as.matrix() %>%
    t()
  colnames(counts) <- asv_ids
  
  storage.mode(counts) <- "numeric"
  common <- intersect(metadata$Sample, rownames(counts))
  meta <- metadata %>%
    filter(Sample %in% common) %>%
    arrange(match(Sample, common))
  counts <- counts[meta$Sample, , drop = FALSE]
  
  read_depth <- rowSums(counts, na.rm = TRUE)
  retained <- read_depth >= minimum_reads
  excluded <- tibble(
    Marker = marker,
    Sample = meta$Sample[!retained],
    Read_Depth = read_depth[!retained],
    Reason = paste0("Fewer than ", minimum_reads, " reads")
  )
  
  list(
    marker = marker,
    counts = counts[retained, , drop = FALSE],
    metadata = meta[retained, , drop = FALSE],
    excluded = excluded
  )
}

community_data <- list(
  `16S` = read_asv_table(
    file.path(data_dir, "1.2_16S_ASV_Table.csv"), "16S", minimum_reads = 1
  ),
  ITS = read_asv_table(
    file.path(data_dir, "1.3_ITS_ASV_Table.csv"), "ITS", minimum_reads = 1000
  )
)

##### 2. Transform Communities and Run NMDS
run_nmds <- function(data) {
  community <- if (data$marker == "16S") {
    decostand(data$counts, method = "total")
  } else {
    decostand(data$counts, method = "hellinger")
  }
  
  bray <- vegdist(community, method = "bray")
  ordination <- metaMDS(
    bray, k = 2, trymax = 200, autotransform = FALSE, trace = FALSE
  )
  
  coordinates <- scores(ordination, display = "sites") %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    left_join(data$metadata, by = "Sample")
  
  list(
    marker = data$marker,
    metadata = data$metadata,
    community = community,
    bray = bray,
    ordination = ordination,
    coordinates = coordinates,
    stress = ordination$stress
  )
}

nmds_results <- map(community_data, run_nmds)

##### 3. Core-Aware Mantel Tests
# Soil variables are tested with permutations restricted within Core_ID.
# Climate and management variables are tested across ranch-level mean
# communities using all 719 non-identity permutations of the six ranches.
# The observed ordering is included by mantel when calculating the p-value.
environment_variables <- tribble(
  ~Variable,          ~Column,                ~Scale,
  "MAT",              "mean_annual_temp",     "Ranch",
  "PET",              "PET_mm_year",          "Ranch",
  "pH",               "pH",                   "Sample",
  "MAP",              "mean_annual_precip",   "Ranch",
  "MAOM-N",           "maom_n",               "Sample",
  "TN",               "total_n",              "Sample",
  "MAOC",             "maom_c",               "Sample",
  "SOC",              "total_c",              "Sample",
  "Solar radiation",  "radiation",            "Ranch",
  "Sand",             "percent_sand",         "Sample",
  "Adaptive score",   "Adaptive_Score",       "Ranch",
  "Depth",            "depth_num",            "Sample",
  "Clay",             "percent_clay",         "Sample",
  "POM-N",            "pom_n",                "Sample",
  "POC",              "pom_c",                "Sample",
  "Silt",             "percent_silt",         "Sample"
)

run_mantel <- function(result, variable, column, scale, permutations = 9999) {
  meta <- result$metadata
  community <- result$community
  
  if (scale == "Ranch") {
    ranch_rows <- split(seq_len(nrow(meta)), meta$Location)
    ranch_community <- map_dfr(ranch_rows, ~as_tibble_row(
      colMeans(community[.x, , drop = FALSE])
    )) %>% as.matrix()
    ranch_values <- map_dbl(ranch_rows, function(rows) {
      values <- unique(na.omit(meta[[column]][rows]))
      if (length(values) != 1) stop(column, " is not constant within ranch")
      values
    })
    
    community_distance <- vegdist(ranch_community, method = "bray")
    environmental_distance <- dist(ranch_values)
    n <- length(ranch_rows)
    ranch_permutations <- permute::allPerms(
      n, control = permute::how(
        within = permute::Within(type = "free"), observed = FALSE
      )
    )
    test <- mantel(
      community_distance, environmental_distance,
      method = "spearman", permutations = ranch_permutations
    )
    statistic <- unname(test$statistic)
    p_value <- test$signif
  } else {
    keep <- is.finite(meta[[column]])
    meta <- meta[keep, , drop = FALSE]
    community <- community[keep, , drop = FALSE]
    community_distance <- vegdist(community, method = "bray")
    environmental_distance <- dist(meta[[column]])
    test <- mantel(
      community_distance, environmental_distance,
      method = "spearman", permutations = permutations,
      strata = meta$Core_ID
    )
    statistic <- unname(test$statistic)
    p_value <- test$signif
    n <- nrow(meta)
  }
  
  tibble(
    Marker = result$marker, Variable = variable, Column = column, Scale = scale,
    n = n, Mantel_Rho = statistic, p_value = p_value
  )
}

mantel_results <- imap_dfr(nmds_results, function(result, marker) {
  pmap_dfr(environment_variables, ~run_mantel(result, ..1, ..2, ..3))
}) %>%
  group_by(Marker) %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    Significance = case_when(
      p_adjusted <= 0.001 ~ "***",
      p_adjusted <= 0.01 ~ "**",
      p_adjusted <= 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ungroup()

##### 4. NMDS Figure
practice_depth_colors <- c(
  "Adaptive: 0-15 cm" = "#B85C00", "Adaptive: 15-30 cm" = "#F39C34",
  "Adaptive: 30-50 cm" = "#FDD9A6", "Conventional: 0-15 cm" = "#5B2A9D",
  "Conventional: 15-30 cm" = "#9687CB", "Conventional: 30-50 cm" = "#D7DBED"
)
location_levels <- sort(unique(metadata$Location))
location_shapes <- setNames(c(17, 15, 16, 18, 4, 7), location_levels)

# Omitted from the figure only; retained in all analyses.
display_exclusions_16S <- c(
  "WYCO_0923_05_PR_D3_B_C03_XT",
  "WYCO_0923_05_PR_D2_B_C03_XT",
  "WYCO_0923_015_AD_D3_B_C05_XT"
)

make_nmds_plot <- function(result, title) {
  plot_data <- result$coordinates %>%
    mutate(
      Depth_Label = recode(depth, D1 = "0-15 cm", D2 = "15-30 cm", D3 = "30-50 cm"),
      Practice_Depth = paste(Practice, Depth_Label, sep = ": ")
    )
  
  if (result$marker == "16S") {
    plot_data <- plot_data %>%
      filter(!Sample %in% display_exclusions_16S)
  }
  
  axis_span <- max(
    diff(range(plot_data$NMDS1)),
    diff(range(plot_data$NMDS2))
  ) * 1.08
  x_midpoint <- mean(range(plot_data$NMDS1))
  y_midpoint <- mean(range(plot_data$NMDS2))
  x_limits <- x_midpoint + c(-0.5, 0.5) * axis_span
  y_limits <- y_midpoint + c(-0.5, 0.5) * axis_span
  
  ggplot(plot_data, aes(NMDS1, NMDS2, color = Practice_Depth, shape = Location)) +
    geom_point(size = 3.2, alpha = 0.90) +
    annotate(
      "text", x = x_limits[2], y = y_limits[2],
      label = sprintf("Stress = %.3f", result$stress),
      hjust = 1.1, vjust = 1.3, size = 4
    ) +
    scale_color_manual(values = practice_depth_colors, name = NULL) +
    scale_shape_manual(values = location_shapes, name = "Ranch") +
    labs(title = title, x = "NMDS1", y = "NMDS2") +
    coord_fixed(xlim = x_limits, ylim = y_limits, ratio = 1) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(), plot.title = element_text(face = "bold", hjust = 0.5),
      legend.key.height = grid::unit(0.45, "cm")
    )
}

nmds_16S_plot <- make_nmds_plot(nmds_results[["16S"]], "Bacteria and Archaea")
nmds_ITS_plot <- make_nmds_plot(nmds_results[["ITS"]], "Fungi")
nmds_panel <- ((nmds_16S_plot / nmds_ITS_plot) + plot_layout(guides = "collect")) &
  theme(legend.position = "left")

##### 5. Mantel Butterfly Figure
variable_order <- environment_variables$Variable
axis_limit <- max(abs(mantel_results$Mantel_Rho), na.rm = TRUE) * 1.25

butterfly_data <- mantel_results %>%
  mutate(
    Variable = factor(Variable, levels = rev(variable_order)),
    Plot_Rho = if_else(Marker == "16S", -abs(Mantel_Rho), abs(Mantel_Rho)),
    Star_X = Plot_Rho + if_else(Marker == "16S", -0.018, 0.018),
    Star_Hjust = if_else(Marker == "16S", 1, 0)
  )

butterfly_plot <- ggplot(
  butterfly_data, aes(Plot_Rho, Variable, fill = Marker)
) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linewidth = 0.7) +
  geom_text(
    aes(x = Star_X, label = Significance, hjust = Star_Hjust),
    size = 5, fontface = "bold", show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("16S" = "#286765", "ITS" = "#8D9C82"),
    labels = c("16S" = "Bacteria and Archaea", "ITS" = "Fungi"),
    name = NULL
  ) +
  scale_x_continuous(
    limits = c(-axis_limit, axis_limit),
    labels = function(x) abs(x), expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(title = "Mantel tests", x = "Absolute Mantel rho", y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top", axis.line.y = element_blank(),
    axis.ticks.y = element_blank(), axis.text.y = element_text(color = "black")
  )

combined_figure <- nmds_panel | butterfly_plot
combined_figure <- combined_figure + plot_layout(widths = c(1.20, 1))

##### 6. Save Results and Figures
nmds_coordinates <- imap_dfr(nmds_results, ~mutate(.x$coordinates, Marker = .y))
nmds_stress <- imap_dfr(nmds_results, ~tibble(Marker = .y, Stress = .x$stress))
excluded_samples <- map_dfr(community_data, "excluded")

write_csv(nmds_coordinates, file.path(output_dir, "NMDS_coordinates.csv"))
write_csv(nmds_stress, file.path(output_dir, "NMDS_stress.csv"))
write_csv(mantel_results, file.path(output_dir, "Mantel_results.csv"))
write_csv(excluded_samples, file.path(output_dir, "excluded_low_read_samples.csv"))

save_figure <- function(plot, name, width, height) {
  ggsave(file.path(figure_dir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 300, bg = "white")
  ggsave(file.path(figure_dir, paste0(name, ".pdf")), plot,
         width = width, height = height)
}

save_figure(nmds_panel, "Figure_4A_NMDS", 8.5, 10)
save_figure(butterfly_plot, "Figure_4B_Mantel_butterfly", 7.5, 8)
save_figure(combined_figure, "Figure_4_NMDS_and_Mantel", 12, 10)

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
print(nmds_stress)
print(mantel_results, n = Inf)
message("Community-composition analysis complete. Results saved to: ", output_dir)





