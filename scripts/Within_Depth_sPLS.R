################################################################################
### Within-depth sPLS sensitivity analysis
###
### Runs separate 16S and ITS sPLS models for D1, D2, and D3. Because each
### core contributes only one sample within a depth, this removes repeated
### depths within cores and prevents associations from being driven only by
### average differences among depths.
################################################################################


##### Step 1: Load Libraries and Set Analysis Options
library(tidyverse)
library(mixOmics)

set.seed(42)

metadata_path <- file.path("data", "1.1_Metadata.csv")
asv_paths <- c(
  "16S" = file.path("data", "1.2_16S_ASV_Table.csv"),
  "ITS" = file.path("data", "1.3_ITS_ASV_Table.csv")
)

output_dir <- file.path("outputs", "sPLS_within_depth")
figure_dir <- file.path("figures", "sPLS_within_depth")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

minimum_prevalence <- 0.10
minimum_total_count <- 10
pseudocount <- 0.5
ncomp_grid <- 1:3
keepX_grid <- c(5, 10, 15, 25)
minimum_selection_frequency <- 0.50
depth_levels <- c("D1", "D2", "D3")


##### Step 2: Load and Check Metadata
metadata <- read_csv(metadata_path, show_col_types = FALSE) %>%
  mutate(
    Location = factor(Location),
    depth = factor(depth, levels = depth_levels)
  )

stopifnot(
  all(c("Sample", "Location", "Core_ID", "depth", "maom_c") %in%
        names(metadata)),
  !anyDuplicated(metadata$Sample),
  !anyDuplicated(metadata[c("Core_ID", "depth")])
)


##### Step 3: Define Data-Preparation Functions
read_asv_data <- function(path) {
  table <- read_csv(path, show_col_types = FALSE)
  id <- intersect(c("ASV", "OTU"), names(table))
  stopifnot(length(id) == 1)
  
  table <- table %>% rename(ASV = all_of(id))
  sample_columns <- setdiff(names(table), c("ASV", "Taxon"))
  counts <- t(as.matrix(table[, sample_columns]))
  storage.mode(counts) <- "numeric"
  colnames(counts) <- table$ASV
  
  taxonomy <- table %>%
    dplyr::select(ASV, Taxon) %>%
    separate(
      Taxon,
      into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
      sep = ";", fill = "right", extra = "drop", remove = FALSE
    ) %>%
    mutate(
      across(Domain:Species, ~ str_remove(.x, "^[^_]+__")),
      Phylum = replace_na(Phylum, "Unknown")
    )
  
  list(counts = counts, taxonomy = taxonomy)
}

# ASV filtering is estimated from training samples only.
prepare_clr <- function(train_counts, test_counts) {
  keep <- colMeans(train_counts > 0) >= minimum_prevalence &
    colSums(train_counts) >= minimum_total_count
  
  transform_clr <- function(x) {
    logged <- log(x[, keep, drop = FALSE] + pseudocount)
    sweep(logged, 1, rowMeans(logged), "-")
  }
  
  train_clr <- transform_clr(train_counts)
  test_clr <- transform_clr(test_counts)
  
  variable_asvs <- apply(train_clr, 2, sd) > 0
  list(
    train = train_clr[, variable_asvs, drop = FALSE],
    test = test_clr[, variable_asvs, drop = FALSE]
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

mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}


##### Step 4: Define Grouped Tuning and Feature Selection
# Hyperparameters are tuned within the outer training data by leaving out one
# additional ranch at a time.
tune_spls <- function(counts, y, ranch) {
  expand_grid(ncomp = ncomp_grid, keepX = keepX_grid) %>%
    pmap_dfr(function(ncomp, keepX) {
      predictions <- map_dfr(unique(ranch), function(held_out) {
        test <- ranch == held_out
        train <- !test
        clr <- prepare_clr(
          counts[train, , drop = FALSE],
          counts[test, , drop = FALSE]
        )
        
        ncomp_fit <- min(ncomp, nrow(clr$train) - 1, ncol(clr$train))
        keepX_fit <- min(keepX, ncol(clr$train))
        
        model <- spls(
          clr$train, y[train],
          ncomp = ncomp_fit,
          keepX = rep(keepX_fit, ncomp_fit),
          mode = "regression",
          scale = TRUE
        )
        
        tibble(
          observed = y[test],
          predicted = predict_spls(model, clr$test, ncomp_fit)
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
    Training_Direction = if_else(association >= 0, "Positive", "Negative")
  )
}


##### Step 5: Run One Marker at One Depth
run_within_depth_spls <- function(marker, path, current_depth) {
  asv_data <- read_asv_data(path)
  
  sample_data <- metadata %>%
    filter(
      depth == current_depth,
      Sample %in% rownames(asv_data$counts),
      !is.na(Location), !is.na(Core_ID),
      !is.na(maom_c), maom_c > 0
    )
  
  # This is the key repeated-measures check: one observation per core.
  stopifnot(!anyDuplicated(sample_data$Core_ID))
  
  counts <- asv_data$counts[sample_data$Sample, , drop = FALSE]
  y <- log(sample_data$maom_c)
  ranch <- as.character(sample_data$Location)
  
  predictions <- list()
  selections <- list()
  tuning_choices <- list()
  
  for (held_out in sort(unique(ranch))) {
    message(marker, " ", current_depth, ": holding out ", held_out)
    test <- ranch == held_out
    train <- !test
    
    tuning <- tune_spls(
      counts[train, , drop = FALSE], y[train], ranch[train]
    )
    clr <- prepare_clr(
      counts[train, , drop = FALSE],
      counts[test, , drop = FALSE]
    )
    
    ncomp <- min(tuning$ncomp, nrow(clr$train) - 1, ncol(clr$train))
    keepX <- min(tuning$keepX, ncol(clr$train))
    
    model <- spls(
      clr$train, y[train],
      ncomp = ncomp,
      keepX = rep(keepX, ncomp),
      mode = "regression",
      scale = TRUE
    )
    
    selected <- extract_selected_asvs(model, clr$train, y[train], ncomp) %>%
      mutate(
        Marker = marker,
        Depth = current_depth,
        Held_Out_Location = held_out,
        Held_Out_Association = map_dbl(
          ASV, ~ safe_cor(clr$test[, .x], y[test])
        ),
        Held_Out_Direction = case_when(
          Held_Out_Association > 0 ~ "Positive",
          Held_Out_Association < 0 ~ "Negative",
          TRUE ~ NA_character_
        ),
        Direction_Agreement = Training_Direction == Held_Out_Direction,
        .before = 1
      )
    
    predictions[[held_out]] <- tibble(
      Marker = marker,
      Depth = current_depth,
      Sample = sample_data$Sample[test],
      Core_ID = sample_data$Core_ID[test],
      Held_Out_Location = held_out,
      Observed_Log_MAOC = y[test],
      Predicted_Log_MAOC = predict_spls(model, clr$test, ncomp)
    )
    
    selections[[held_out]] <- selected
    tuning_choices[[held_out]] <- tuning %>%
      mutate(
        Marker = marker,
        Depth = current_depth,
        Held_Out_Location = held_out,
        .before = 1
      )
  }
  
  predictions <- bind_rows(predictions)
  selections <- bind_rows(selections)
  
  performance <- predictions %>%
    summarize(
      Marker = marker,
      Depth = current_depth,
      n = n(),
      RMSE = sqrt(mean((Observed_Log_MAOC - Predicted_Log_MAOC)^2)),
      MAE = mean(abs(Observed_Log_MAOC - Predicted_Log_MAOC)),
      R2_CV = 1 - sum((Observed_Log_MAOC - Predicted_Log_MAOC)^2) /
        sum((Observed_Log_MAOC - mean(Observed_Log_MAOC))^2)
    )
  
  stability <- selections %>%
    group_by(Marker, Depth, ASV) %>%
    summarize(
      Selection_Frequency = n_distinct(Held_Out_Location) /
        n_distinct(ranch),
      Mean_VIP = mean_or_na(VIP),
      Mean_Held_Out_Association = mean_or_na(Held_Out_Association),
      Mean_Held_Out_R2 = mean_or_na(Held_Out_Association^2),
      Training_Positive_Frequency = mean_or_na(
        as.numeric(Training_Direction == "Positive")
      ),
      Held_Out_Positive_Frequency = mean_or_na(
        as.numeric(Held_Out_Direction == "Positive")
      ),
      Direction_Agreement = mean_or_na(as.numeric(Direction_Agreement)),
      .groups = "drop"
    ) %>%
    left_join(asv_data$taxonomy, by = "ASV") %>%
    arrange(desc(Selection_Frequency), desc(Direction_Agreement), desc(Mean_VIP))
  
  list(
    performance = performance,
    predictions = predictions,
    selections = selections,
    stability = stability,
    tuning = bind_rows(tuning_choices)
  )
}


##### Step 6: Run 16S and ITS Separately for D1, D2, and D3
analysis_grid <- crossing(
  Marker = names(asv_paths),
  Depth = depth_levels
)

results <- pmap(
  analysis_grid,
  ~ run_within_depth_spls(..1, asv_paths[[..1]], ..2)
)

performance <- map_dfr(results, "performance")
outer_predictions <- map_dfr(results, "predictions")
selected_by_fold <- map_dfr(results, "selections")
selection_stability <- map_dfr(results, "stability")
tuning_choices <- map_dfr(results, "tuning")


##### Step 7: Compare Stable ASVs Across Depths and With the Pooled Analysis
stable_asvs <- selection_stability %>%
  filter(Selection_Frequency >= minimum_selection_frequency)

cross_depth_consistency <- stable_asvs %>%
  group_by(Marker, ASV) %>%
  summarize(
    Number_of_Depths = n_distinct(Depth),
    Depths = paste(sort(unique(Depth)), collapse = ", "),
    Mean_Selection_Frequency = mean(Selection_Frequency),
    Mean_Direction_Agreement = mean(Direction_Agreement, na.rm = TRUE),
    Taxon = first(Taxon),
    Phylum = first(Phylum),
    Genus = first(Genus),
    .groups = "drop"
  ) %>%
  arrange(Marker, desc(Number_of_Depths), desc(Mean_Selection_Frequency))

pooled_file <- file.path("outputs", "sPLS", "sPLS_ASV_selection_stability.csv")

if (file.exists(pooled_file)) {
  pooled_stable <- read_csv(pooled_file, show_col_types = FALSE) %>%
    filter(Selection_Frequency >= minimum_selection_frequency) %>%
    dplyr::select(Marker, ASV, Pooled_Selection_Frequency = Selection_Frequency)
  
  pooled_overlap <- stable_asvs %>%
    inner_join(pooled_stable, by = c("Marker", "ASV")) %>%
    arrange(Marker, Depth, desc(Selection_Frequency))
} else {
  pooled_overlap <- tibble()
  warning("Pooled sPLS stability file was not found; overlap was not calculated.")
}


##### Step 8: Save Results and a Performance Summary Figure
write_csv(performance, file.path(output_dir, "within_depth_performance.csv"))
write_csv(
  outer_predictions,
  file.path(output_dir, "within_depth_outer_predictions.csv")
)
write_csv(
  selected_by_fold,
  file.path(output_dir, "within_depth_selected_ASVs_by_fold.csv")
)
write_csv(
  selection_stability,
  file.path(output_dir, "within_depth_ASV_selection_stability.csv")
)
write_csv(
  cross_depth_consistency,
  file.path(output_dir, "within_depth_cross_depth_consistency.csv")
)
write_csv(
  pooled_overlap,
  file.path(output_dir, "within_depth_pooled_ASV_overlap.csv")
)
write_csv(
  tuning_choices,
  file.path(output_dir, "within_depth_tuning_choices.csv")
)

performance_plot <- ggplot(
  performance,
  aes(x = Depth, y = R2_CV, color = Marker, group = Marker)
) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  labs(
    x = "Soil depth",
    y = expression("Leave-one-ranch-out " * R^2),
    color = "Marker"
  ) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

print(performance_plot)

ggsave(
  file.path(figure_dir, "sPLS_within_depth_performance.png"),
  performance_plot, width = 7, height = 5, dpi = 300
)
ggsave(
  file.path(figure_dir, "sPLS_within_depth_performance.pdf"),
  performance_plot, width = 7, height = 5
)

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

print(performance)
print(cross_depth_consistency, n = 30)
message("Analysis complete. Results saved to: ", output_dir)





cross_depth_consistency %>%
  count(Marker, Number_of_Depths)

stable_asvs %>%
  group_by(Marker, Depth) %>%
  summarize(
    Stable_ASVs = n(),
    Directionally_Consistent = sum(Direction_Agreement >= 0.67, na.rm = TRUE),
    .groups = "drop"
  )

pooled_overlap %>%
  count(Marker, Depth, name = "Pooled_ASV_Overlap")

pooled_overlap %>%
  filter(Depth == "D1") %>%
  dplyr::select(
    Marker, ASV, Selection_Frequency,
    Pooled_Selection_Frequency,
    Held_Out_Positive_Frequency,
    Direction_Agreement,
    Phylum, Genus
  ) %>%
  arrange(Marker, desc(Selection_Frequency)) %>%
  print(n = 50)


pooled_overlap %>%
  filter(Depth == "D1") %>%
  dplyr::select(
    Marker,
    ASV,
    Held_Out_Positive_Frequency,
    Direction_Agreement,
    Phylum,
    Genus
  ) %>%
  print(n = Inf, width = Inf)
