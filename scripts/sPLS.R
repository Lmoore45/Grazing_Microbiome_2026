################################################################################
### Sparse partial least squares (sPLS) regression analysis
###
### Identifies depth-adjusted ASVs associated with MAOC using nested
### leave-one-ranch-out validation. Filtering, tuning, and ASV selection occur
### within training data only. Separate VIP bar charts and fold-selection
### heatmaps are produced for 16S and ITS.
################################################################################


##### Step 1: Load Libraries and Set Analysis Options
library(tidyverse)
library(mixOmics)
library(paletteer)
library(patchwork)

set.seed(42)

metadata_path <- file.path("data", "1.1_Metadata.csv")
asv_paths <- c(
  "16S" = file.path("data", "1.2_16S_ASV_Table.csv"),
  "ITS" = file.path("data", "1.3_ITS_ASV_Table.csv")
)

output_dir <- file.path("outputs", "sPLS")
figure_dir <- file.path("figures", "sPLS")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

minimum_prevalence <- 0.10
minimum_total_count <- 10
pseudocount <- 0.5
ncomp_grid <- 1:3
keepX_grid <- c(5, 10, 15, 25)
minimum_selection_frequency <- 0.50
number_of_taxa <- 10
depth_levels <- c("D1", "D2", "D3")


##### Step 2: Load and Check Metadata
metadata <- read_csv(metadata_path, show_col_types = FALSE) %>%
  mutate(
    Location = factor(Location),
    depth = factor(depth, levels = depth_levels)
  )

required_columns <- c(
  "Sample", "Location", "Triangle", "Core_ID", "depth", "maom_c"
)

stopifnot(
  all(required_columns %in% names(metadata)),
  !anyDuplicated(metadata$Sample),
  !anyDuplicated(metadata[c("Core_ID", "depth")])
)


##### Step 3: Define Data-Preparation Functions
# Read an ASV table and return samples as rows and ASVs as columns.
read_asv_data <- function(path) {
  asv_table <- read_csv(path, show_col_types = FALSE)
  id_column <- intersect(c("ASV", "OTU"), names(asv_table))
  stopifnot(length(id_column) == 1, "Taxon" %in% names(asv_table))
  
  asv_table <- asv_table %>%
    dplyr::rename(ASV = all_of(id_column))
  
  sample_columns <- setdiff(names(asv_table), c("ASV", "Taxon"))
  counts <- t(as.matrix(asv_table[, sample_columns]))
  storage.mode(counts) <- "numeric"
  rownames(counts) <- sample_columns
  colnames(counts) <- asv_table$ASV
  
  taxonomy <- asv_table %>%
    dplyr::select(ASV, Taxon) %>%
    separate(
      Taxon,
      into = c(
        "Domain", "Phylum", "Class", "Order",
        "Family", "Genus", "Species"
      ),
      sep = ";", fill = "right", extra = "drop", remove = FALSE
    ) %>%
    mutate(
      across(Domain:Species, ~ str_remove(.x, "^[^_]+__")),
      Phylum = replace_na(Phylum, "Unknown")
    )
  
  list(counts = counts, taxonomy = taxonomy)
}

# Filter ASVs using training samples, then apply CLR to the same retained ASVs
# in the training and held-out data.
prepare_clr <- function(train_counts, test_counts) {
  keep <- colMeans(train_counts > 0) >= minimum_prevalence &
    colSums(train_counts) >= minimum_total_count
  
  if (!any(keep)) stop("No ASVs passed the training-data filters.")
  
  clr <- function(x) {
    logged <- log(x[, keep, drop = FALSE] + pseudocount)
    sweep(logged, 1, rowMeans(logged), FUN = "-")
  }
  
  list(train = clr(train_counts), test = clr(test_counts))
}

# Remove the average categorical depth effect using training coefficients only.
# The fitted depth effect is then added back to held-out MAOC predictions.
adjust_for_depth <- function(X_train, X_test, y_train, depth_train, depth_test) {
  D_train <- model.matrix(~ factor(depth_train, levels = depth_levels))
  D_test <- model.matrix(~ factor(depth_test, levels = depth_levels))
  coefficients <- qr.solve(D_train, cbind(y_train, X_train))
  
  list(
    X_train = X_train - D_train %*% coefficients[, -1, drop = FALSE],
    X_test = X_test - D_test %*% coefficients[, -1, drop = FALSE],
    y_train = drop(y_train - D_train %*% coefficients[, 1]),
    y_depth_test = drop(D_test %*% coefficients[, 1])
  )
}

prepare_split <- function(counts, y, depth, train, test) {
  clr <- prepare_clr(
    counts[train, , drop = FALSE],
    counts[test, , drop = FALSE]
  )
  
  adjust_for_depth(
    clr$train, clr$test, y[train], depth[train], depth[test]
  )
}

predict_spls <- function(model, newdata, ncomp) {
  drop(predict(model, newdata)$predict[, 1, ncomp])
}

safe_cor <- function(x, y) {
  complete <- is.finite(x) & is.finite(y)
  if (sum(complete) < 3 || sd(x[complete]) == 0 || sd(y[complete]) == 0) {
    return(NA_real_)
  }
  cor(x[complete], y[complete])
}


##### Step 4: Define Grouped Tuning and Feature-Selection Functions
# Tune within each outer training set by leaving out one additional ranch.
tune_spls <- function(counts, y, depth, ranch) {
  expand_grid(ncomp = ncomp_grid, keepX = keepX_grid) %>%
    pmap_dfr(function(ncomp, keepX) {
      predictions <- map_dfr(sort(unique(ranch)), function(held_out) {
        test <- ranch == held_out
        train <- !test
        split <- prepare_split(counts, y, depth, train, test)
        
        ncomp_fit <- min(
          ncomp, nrow(split$X_train) - 1, ncol(split$X_train)
        )
        keepX_fit <- min(keepX, ncol(split$X_train))
        
        model <- spls(
          split$X_train,
          split$y_train,
          ncomp = ncomp_fit,
          keepX = rep(keepX_fit, ncomp_fit),
          mode = "regression",
          scale = TRUE
        )
        
        tibble(
          observed = y[test],
          predicted = split$y_depth_test +
            predict_spls(model, split$X_test, ncomp_fit)
        )
      })
      
      tibble(
        ncomp = ncomp,
        keepX = keepX,
        RMSE = sqrt(mean((predictions$observed - predictions$predicted)^2))
      )
    }) %>%
    arrange(RMSE, ncomp, keepX) %>%
    slice(1)
}

# Extract ASVs with nonzero loadings and record importance and direction in the
# current outer training set.
extract_selected_asvs <- function(model, X, y, ncomp) {
  loadings <- model$loadings$X[, seq_len(ncomp), drop = FALSE]
  selected <- rownames(loadings)[rowSums(abs(loadings)) > 0]
  
  vip <- mixOmics::vip(model)
  if (is.null(dim(vip))) {
    vip <- matrix(vip, ncol = 1, dimnames = list(colnames(X), "comp1"))
  }
  
  association <- map_dbl(selected, ~ safe_cor(X[, .x], y))
  
  tibble(
    ASV = selected,
    VIP = apply(vip[selected, seq_len(ncomp), drop = FALSE], 1, max),
    Training_Association = association,
    Direction = case_when(
      association > 0 ~ "Positive",
      association < 0 ~ "Negative",
      TRUE ~ NA_character_
    )
  )
}

# Sum relative abundance for ASVs selected within a training fold.
sum_relative_abundance <- function(counts, asvs) {
  if (length(asvs) == 0) return(rep(0, nrow(counts)))
  rowSums(counts[, asvs, drop = FALSE]) / rowSums(counts)
}


##### Step 5: Run Nested Leave-One-Ranch-Out sPLS for One Marker
run_spls <- function(marker, path) {
  asv_data <- read_asv_data(path)
  
  sample_data <- metadata %>%
    filter(
      Sample %in% rownames(asv_data$counts),
      !is.na(Location), !is.na(Core_ID), !is.na(depth),
      !is.na(maom_c), maom_c > 0
    )
  
  counts <- asv_data$counts[sample_data$Sample, , drop = FALSE]
  y <- log(sample_data$maom_c)
  ranch <- as.character(sample_data$Location)
  depth <- sample_data$depth
  n_outer_folds <- n_distinct(ranch)
  
  predictions <- list()
  selections <- list()
  fold_sums <- list()
  tuning_choices <- list()
  
  for (held_out in sort(unique(ranch))) {
    message(marker, ": holding out ", held_out)
    test <- ranch == held_out
    train <- !test
    
    tuning <- tune_spls(
      counts[train, , drop = FALSE], y[train], depth[train], ranch[train]
    )
    split <- prepare_split(counts, y, depth, train, test)
    
    ncomp <- min(
      tuning$ncomp, nrow(split$X_train) - 1, ncol(split$X_train)
    )
    keepX <- min(tuning$keepX, ncol(split$X_train))
    
    model <- spls(
      split$X_train,
      split$y_train,
      ncomp = ncomp,
      keepX = rep(keepX, ncomp),
      mode = "regression",
      scale = TRUE
    )
    
    selected <- extract_selected_asvs(
      model, split$X_train, split$y_train, ncomp
    ) %>%
      mutate(
        Held_Out_Location = held_out,
        Held_Out_Association = map_dbl(
          ASV,
          ~ safe_cor(split$X_test[, .x], y[test] - split$y_depth_test)
        ),
        Held_Out_R2 = Held_Out_Association^2
      )
    
    positive <- selected %>%
      filter(Direction == "Positive") %>%
      pull(ASV)
    negative <- selected %>%
      filter(Direction == "Negative") %>%
      pull(ASV)
    
    predictions[[held_out]] <- tibble(
      Marker = marker,
      Sample = sample_data$Sample[test],
      Core_ID = sample_data$Core_ID[test],
      depth = depth[test],
      Held_Out_Location = held_out,
      Observed_Log_MAOC = y[test],
      Predicted_Log_MAOC = split$y_depth_test +
        predict_spls(model, split$X_test, ncomp)
    )
    
    selections[[held_out]] <- selected
    
    # These fold-specific sums are safe inputs for the later PLS because the
    # held-out ranch did not contribute to ASV selection.
    fold_sums[[held_out]] <- tibble(
      Marker = marker,
      Held_Out_Location = held_out,
      Sample = sample_data$Sample,
      Data_Role = if_else(test, "Test", "Train"),
      Positive_ASV_Sum = sum_relative_abundance(counts, positive),
      Negative_ASV_Sum = sum_relative_abundance(counts, negative)
    )
    
    tuning_choices[[held_out]] <- tuning %>%
      mutate(Marker = marker, Held_Out_Location = held_out, .before = 1)
  }
  
  predictions <- bind_rows(predictions)
  selections <- bind_rows(selections)
  fold_sums <- bind_rows(fold_sums)
  tuning_choices <- bind_rows(tuning_choices)
  
  fold_selections <- selections %>%
    mutate(Marker = marker, .before = 1) %>%
    left_join(asv_data$taxonomy, by = "ASV")
  
  performance <- predictions %>%
    summarize(
      Marker = marker,
      n = n(),
      RMSE = sqrt(mean((Observed_Log_MAOC - Predicted_Log_MAOC)^2)),
      MAE = mean(abs(Observed_Log_MAOC - Predicted_Log_MAOC)),
      R2_CV = 1 - sum((Observed_Log_MAOC - Predicted_Log_MAOC)^2) /
        sum((Observed_Log_MAOC - mean(Observed_Log_MAOC))^2)
    )
  
  performance_by_location <- predictions %>%
    group_by(Marker, Held_Out_Location) %>%
    summarize(
      n = n(),
      RMSE = sqrt(mean((Observed_Log_MAOC - Predicted_Log_MAOC)^2)),
      MAE = mean(abs(Observed_Log_MAOC - Predicted_Log_MAOC)),
      .groups = "drop"
    )
  
  selection_stability <- selections %>%
    group_by(ASV) %>%
    summarize(
      Selection_Frequency = n_distinct(Held_Out_Location) / n_outer_folds,
      Mean_VIP = mean(VIP, na.rm = TRUE),
      Mean_Association = mean(Held_Out_Association, na.rm = TRUE),
      Mean_R2 = mean(Held_Out_R2, na.rm = TRUE),
      Positive_Frequency = mean(Direction == "Positive", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(asv_data$taxonomy, by = "ASV") %>%
    mutate(Marker = marker, .before = 1) %>%
    arrange(desc(Selection_Frequency), desc(Mean_R2), desc(Mean_VIP))
  
  list(
    performance = performance,
    performance_by_location = performance_by_location,
    predictions = predictions,
    selection_stability = selection_stability,
    fold_selections = fold_selections,
    fold_sums = fold_sums,
    tuning_choices = tuning_choices
  )
}


##### Step 6: Run Both Markers and Save Analysis Tables
spls_results <- imap(asv_paths, ~ run_spls(.y, .x))

combined_performance <- map_dfr(spls_results, "performance")
performance_by_location <- map_dfr(spls_results, "performance_by_location")
outer_predictions <- map_dfr(spls_results, "predictions")
selection_stability <- map_dfr(spls_results, "selection_stability")
fold_selections <- map_dfr(spls_results, "fold_selections")
fold_specific_sums <- map_dfr(spls_results, "fold_sums")
tuning_choices <- map_dfr(spls_results, "tuning_choices")

write_csv(combined_performance, file.path(output_dir, "sPLS_performance.csv"))
write_csv(
  performance_by_location,
  file.path(output_dir, "sPLS_performance_by_location.csv")
)
write_csv(
  outer_predictions,
  file.path(output_dir, "sPLS_outer_predictions.csv")
)
write_csv(
  selection_stability,
  file.path(output_dir, "sPLS_ASV_selection_stability.csv")
)
write_csv(
  fold_selections,
  file.path(output_dir, "sPLS_fold_level_ASV_selections.csv")
)
write_csv(
  fold_specific_sums,
  file.path(output_dir, "sPLS_fold_specific_ASV_sums_for_PLS.csv")
)
write_csv(
  tuning_choices,
  file.path(output_dir, "sPLS_tuning_choices.csv")
)


##### Step 7: Prepare Taxon Labels and Plot Separate VIP Bar Charts
# Each chart shows the ten highest-VIP ASVs selected in at least half of the
# outer folds. Genus is used when available; otherwise the finest available
# rank is shown. Repeated labels are numbered within each marker.
clean_taxon_rank <- function(x, prefix) {
  x <- str_remove(as.character(x), paste0("^", prefix, "__"))
  x[x %in% c("", "NA", "Unknown")] <- NA_character_
  x
}

prepare_vip_bar_data <- function(data) {
  data %>%
    filter(
      Selection_Frequency >= minimum_selection_frequency,
      is.finite(Mean_VIP)
    ) %>%
    arrange(desc(Mean_VIP), desc(Selection_Frequency)) %>%
    slice_head(n = number_of_taxa) %>%
    mutate(
      Phylum = clean_taxon_rank(Phylum, "p"),
      Family = clean_taxon_rank(Family, "f"),
      Genus = clean_taxon_rank(Genus, "g"),
      Phylum = replace_na(Phylum, "Unknown"),
      Taxon_base = case_when(
        !is.na(Genus) ~ paste0("g_", Genus),
        !is.na(Family) ~ paste0("f_", Family),
        Phylum != "Unknown" ~ paste0("p_", Phylum),
        TRUE ~ paste0("ASV_", str_sub(ASV, 1, 6))
      )
    ) %>%
    arrange(Taxon_base, ASV) %>%
    group_by(Taxon_base) %>%
    mutate(
      Taxon_label = if (n() > 1) {
        paste0(Taxon_base, " (ASV ", row_number(), ")")
      } else {
        Taxon_base
      }
    ) %>%
    ungroup() %>%
    arrange(Mean_VIP) %>%
    mutate(Taxon_label = factor(Taxon_label, levels = Taxon_label))
}

make_phylum_colors <- function(data) {
  phyla <- sort(unique(data$Phylum))
  colors <- colorRampPalette(
    as.character(paletteer_d("ltc::crbhits"))
  )(length(phyla))
  setNames(colors, phyla)
}

make_vip_bar_plot <- function(data, title) {
  phylum_colors <- make_phylum_colors(data)
  
  ggplot(data, aes(x = Mean_VIP, y = Taxon_label, fill = Phylum)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.4) +
    geom_vline(
      xintercept = 1, linetype = "dotted", linewidth = 0.7, color = "black"
    ) +
    scale_fill_manual(
      name = "Phylum", values = phylum_colors, breaks = names(phylum_colors)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = title,
      subtitle = "ASVs selected in at least 50% of leave-one-ranch-out folds",
      x = "Mean sPLS VIP score",
      y = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
      axis.text.y = element_text(color = "black"),
      legend.position = "right"
    )
}

vip_bar_data_16S <- prepare_vip_bar_data(
  spls_results[["16S"]]$selection_stability
)
vip_bar_data_ITS <- prepare_vip_bar_data(
  spls_results[["ITS"]]$selection_stability
)

vip_bar_plot_16S <- make_vip_bar_plot(
  vip_bar_data_16S, "Bacteria and Archaea (16S)"
)
vip_bar_plot_ITS <- make_vip_bar_plot(vip_bar_data_ITS, "Fungi (ITS)")

ggsave(
  file.path(figure_dir, "Figure_5A_sPLS_16S_VIP_bar.png"),
  vip_bar_plot_16S, width = 8, height = 6, dpi = 300, bg = "white"
)
ggsave(
  file.path(figure_dir, "Figure_5A_sPLS_16S_VIP_bar.pdf"),
  vip_bar_plot_16S, width = 8, height = 6
)
ggsave(
  file.path(figure_dir, "Figure_5B_sPLS_ITS_VIP_bar.png"),
  vip_bar_plot_ITS, width = 8, height = 6, dpi = 300, bg = "white"
)
ggsave(
  file.path(figure_dir, "Figure_5B_sPLS_ITS_VIP_bar.pdf"),
  vip_bar_plot_ITS, width = 8, height = 6
)

##### Step 8: Plot Selection Direction and Frequency Across Ranch Folds
# Rows match the VIP charts, and columns represent leave-one-ranch-out folds.
# Tile color shows whether each ASV had a positive or negative association with
# MAOC in the training data. The adjacent bar shows overall selection frequency.

prepare_selection_heatmap <- function(fold_data, bar_data) {
  location_levels <- unique(fold_data$Held_Out_Location)
  location_levels <- location_levels[
    order(readr::parse_number(location_levels))
  ]
  
  taxon_information <- bar_data %>%
    dplyr::mutate(
      Taxon_label = as.character(Taxon_label)
    ) %>%
    dplyr::select(
      ASV,
      Taxon_label,
      Selection_Frequency,
      Mean_VIP
    )
  
  taxon_levels <- taxon_information %>%
    dplyr::arrange(Mean_VIP) %>%
    dplyr::pull(Taxon_label)
  
  tidyr::expand_grid(
    ASV = taxon_information$ASV,
    Held_Out_Location = location_levels
  ) %>%
    dplyr::left_join(
      fold_data %>%
        dplyr::filter(
          ASV %in% taxon_information$ASV
        ) %>%
        dplyr::select(
          ASV,
          Held_Out_Location,
          Direction
        ),
      by = c("ASV", "Held_Out_Location")
    ) %>%
    dplyr::left_join(
      taxon_information,
      by = "ASV"
    ) %>%
    dplyr::mutate(
      Selection = dplyr::case_when(
        Direction == "Positive" ~ "Positive",
        Direction == "Negative" ~ "Negative",
        TRUE ~ "Not selected"
      ),
      Selection = factor(
        Selection,
        levels = c(
          "Positive",
          "Negative",
          "Not selected"
        )
      ),
      Held_Out_Location = factor(
        Held_Out_Location,
        levels = location_levels
      ),
      Taxon_label = factor(
        Taxon_label,
        levels = taxon_levels
      )
    )
}


make_selection_figure <- function(data, title) {
  
  ##### Selection Heatmap
  heatmap <- ggplot(
    data,
    aes(
      x = Held_Out_Location,
      y = Taxon_label,
      fill = Selection
    )
  ) +
    geom_tile(
      color = "white",
      linewidth = 0.8
    ) +
    scale_fill_manual(
      name = "Training association",
      values = c(
        "Positive" = "#2A9D8F",
        "Negative" = "#E76F51",
        "Not selected" = "grey92"
      ),
      breaks = c(
        "Positive",
        "Negative",
        "Not selected"
      ),
      drop = FALSE
    ) +
    scale_x_discrete(
      labels = function(x) {
        stringr::str_remove(x, "^Site ")
      }
    ) +
    scale_y_discrete(drop = FALSE) +
    labs(
      title = title,
      x = "Held-out ranch",
      y = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      axis.text = element_text(
        color = "black"
      ),
      legend.position = "bottom"
    )
  
  ##### Selection-Frequency Panel
  frequency_data <- data %>%
    dplyr::distinct(
      Taxon_label,
      Selection_Frequency
    )
  
  frequency <- ggplot(
    frequency_data,
    aes(
      x = Selection_Frequency,
      y = Taxon_label
    )
  ) +
    geom_col(
      width = 0.68,
      fill = "grey35"
    ) +
    geom_text(
      aes(
        label = scales::percent(
          Selection_Frequency,
          accuracy = 1
        )
      ),
      hjust = -0.15,
      size = 3.4
    ) +
    scale_x_continuous(
      limits = c(0, 1.18),
      breaks = c(0, 0.5, 1),
      labels = scales::percent_format(
        accuracy = 1
      ),
      expand = expansion(
        mult = c(0, 0)
      )
    ) +
    scale_y_discrete(drop = FALSE) +
    labs(
      x = "Selection frequency",
      y = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  ##### Combine Heatmap and Frequency Panel
  combined_figure <- heatmap + frequency +
    patchwork::plot_layout(
      widths = c(3.8, 1.4),
      guides = "collect"
    )
  
  combined_figure &
    theme(
      legend.position = "bottom"
    )
}


##### Prepare Separate 16S and ITS Data
selection_data_16S <- prepare_selection_heatmap(
  spls_results[["16S"]]$fold_selections,
  vip_bar_data_16S
)

selection_data_ITS <- prepare_selection_heatmap(
  spls_results[["ITS"]]$fold_selections,
  vip_bar_data_ITS
)


##### Create Separate 16S and ITS Figures
selection_figure_16S <- make_selection_figure(
  selection_data_16S,
  "Bacteria and Archaea (16S)"
)

selection_figure_ITS <- make_selection_figure(
  selection_data_ITS,
  "Fungi (ITS)"
)

print(selection_figure_16S)
print(selection_figure_ITS)


##### Save Figures
ggsave(
  file.path(
    figure_dir,
    "Figure_5C_sPLS_16S_selection_heatmap.png"
  ),
  selection_figure_16S,
  width = 10,
  height = 6.2,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(
    figure_dir,
    "Figure_5C_sPLS_16S_selection_heatmap.pdf"
  ),
  selection_figure_16S,
  width = 10,
  height = 6.2
)

ggsave(
  file.path(
    figure_dir,
    "Figure_5D_sPLS_ITS_selection_heatmap.png"
  ),
  selection_figure_ITS,
  width = 10,
  height = 6.2,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(
    figure_dir,
    "Figure_5D_sPLS_ITS_selection_heatmap.pdf"
  ),
  selection_figure_ITS,
  width = 10,
  height = 6.2
)

message(
  "Separate 16S and ITS selection heatmaps saved to: ",
  figure_dir
)

##### Step 9: Print Results and Record the R Session
print(vip_bar_plot_16S)
print(vip_bar_plot_ITS)
print(selection_figure_16S)
print(selection_figure_ITS)
print(combined_performance)
print(performance_by_location)

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

message("Analysis complete. Results saved to: ", output_dir)
message("Figures saved to: ", figure_dir)


vip_bar_data_ITS %>%
  dplyr::filter(
    stringr::str_detect(
      stringr::str_to_lower(as.character(Taxon_label)),
      "unidentified"
    )
  ) %>%
  dplyr::select(
    ASV, Taxon, Phylum, Class, Order,
    Family, Genus, Species, Mean_VIP,
    Selection_Frequency
  ) %>%
  print(width = Inf)


dplyr::bind_rows(
  vip_bar_data_ITS %>%
    dplyr::filter(
      stringr::str_detect(
        as.character(Taxon_label),
        "518303"
      )
    ) %>%
    dplyr::mutate(Marker = "ITS"),
  
  vip_bar_data_16S %>%
    dplyr::filter(
      stringr::str_detect(
        as.character(Taxon_label),
        "Actinomycetota"
      )
    ) %>%
    dplyr::mutate(Marker = "16S")
) %>%
  dplyr::select(
    Marker, ASV, Taxon, Phylum, Class,
    Order, Family, Genus, Species,
    Mean_VIP, Selection_Frequency
  ) %>%
  print(width = Inf)




vip_bar_data_ITS %>%
  dplyr::mutate(
    Taxon_label = as.character(Taxon_label)
  ) %>%
  dplyr::arrange(
    dplyr::desc(Mean_VIP)
  ) %>%
  dplyr::select(
    ASV,
    Taxon_label,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    Selection_Frequency,
    Mean_VIP,
    Mean_Association
  ) %>%
  print(n = Inf, width = Inf)






##### Summarize Final sPLS Results for Manuscript
minimum_selection_frequency <- 0.50

##### 1. Out-of-ranch performance
model_performance <- dplyr::bind_rows(
  spls_results[["16S"]]$performance,
  spls_results[["ITS"]]$performance
) %>%
  dplyr::select(dplyr::any_of(
    c("Marker", "n", "RMSE", "MAE", "R2_CV")
  ))

cat("\nOUT-OF-RANCH MODEL PERFORMANCE\n")
print(model_performance, n = Inf, width = Inf)

##### 2. Stable ASVs
stable_taxa <- dplyr::bind_rows(
  spls_results[["16S"]]$selection_stability,
  spls_results[["ITS"]]$selection_stability
) %>%
  dplyr::filter(Selection_Frequency >= minimum_selection_frequency) %>%
  dplyr::mutate(
    Taxon_Label = dplyr::case_when(
      !is.na(Genus) & !is.na(Species) ~
        paste0("g_", Genus, "; s_", Species),
      !is.na(Genus) ~ paste0("g_", Genus),
      !is.na(Family) ~ paste0("f_", Family),
      !is.na(Order) ~ paste0("o_", Order),
      !is.na(Class) ~ paste0("c_", Class),
      !is.na(Phylum) ~ paste0("p_", Phylum),
      TRUE ~ paste0("ASV_", stringr::str_sub(ASV, 1, 6))
    ),
    Direction = dplyr::if_else(
      Mean_Association >= 0, "Positive", "Negative"
    )
  ) %>%
  dplyr::arrange(
    Marker,
    dplyr::desc(Selection_Frequency),
    dplyr::desc(Mean_VIP)
  ) %>%
  dplyr::select(
    Marker, ASV, Taxon_Label,
    Phylum, Class, Order, Family, Genus, Species,
    Selection_Frequency, Mean_VIP, Mean_Association,
    Positive_Frequency, Direction, Mean_R2
  )

cat("\nSTABLE ASVs\n")
print(stable_taxa, n = Inf, width = Inf)

##### 3. Number of stable ASVs
stable_taxa_counts <- stable_taxa %>%
  dplyr::count(Marker, Direction, name = "Number_of_ASVs") %>%
  tidyr::complete(
    Marker,
    Direction = c("Positive", "Negative"),
    fill = list(Number_of_ASVs = 0)
  )

cat("\nNUMBER OF STABLE ASVs\n")
print(stable_taxa_counts, n = Inf, width = Inf)

##### 4. Phylum representation
stable_phyla <- stable_taxa %>%
  dplyr::mutate(Phylum = tidyr::replace_na(Phylum, "Unknown")) %>%
  dplyr::count(Marker, Phylum, name = "Number_of_ASVs") %>%
  dplyr::group_by(Marker) %>%
  dplyr::mutate(
    Percent = 100 * Number_of_ASVs / sum(Number_of_ASVs)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Marker, dplyr::desc(Number_of_ASVs))

cat("\nPHYLA REPRESENTED AMONG STABLE ASVs\n")
print(stable_phyla, n = Inf, width = Inf)

##### 5. Selection in individual held-out-ranch folds
fold_selections <- dplyr::bind_rows(
  spls_results[["16S"]]$fold_selections,
  spls_results[["ITS"]]$fold_selections
) %>%
  dplyr::filter(ASV %in% stable_taxa$ASV) %>%
  dplyr::arrange(Marker, ASV, Held_Out_Location)

cat("\nSELECTION BY HELD-OUT RANCH\n")
print(fold_selections, n = Inf, width = Inf)
