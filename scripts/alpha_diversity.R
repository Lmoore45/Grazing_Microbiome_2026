################################################################################
### Alpha diversity: 16S and ITS
### Accounts for soil depth, ranch, and repeated depths from the same core.
################################################################################

##### Step 1: Libraries and Paths
library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
library(ggh4x)

metadata_path <- file.path("data", "1.1_Metadata.csv")
asv_paths <- c(
  "16S" = file.path("data", "1.2_16S_ASV_Table.csv"),
  "ITS" = file.path("data", "1.3_ITS_ASV_Table.csv")
)
output_dir <- file.path("outputs", "alpha_diversity")
figure_dir <- file.path("figures", "alpha_diversity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


##### Step 2: Calculate Alpha Diversity
# Evenness is Pielou's J. Shannon is calculated from untransformed counts.
calculate_alpha <- function(path, marker) {
  asv <- read_csv(path, show_col_types = FALSE)
  stopifnot(all(c("ASV", "Taxon") %in% names(asv)), !anyDuplicated(asv$ASV))
  
  samples <- setdiff(names(asv), c("ASV", "Taxon"))
  counts <- t(as.matrix(asv[, samples, drop = FALSE]))
  storage.mode(counts) <- "numeric"
  rownames(counts) <- samples
  if (anyNA(counts) || any(counts < 0)) stop(marker, " contains invalid counts.")
  
  richness <- rowSums(counts > 0)
  shannon <- vegan::diversity(counts, index = "shannon")
  tibble(
    Sample = samples, Marker = marker, Read_Depth = rowSums(counts),
    Richness = richness,
    Evenness = if_else(richness > 1, shannon / log(richness), NA_real_),
    Shannon = shannon
  )
}

metadata <- read_csv(metadata_path, show_col_types = FALSE)
required <- c("Sample", "maom_c", "Practice", "depth", "Core_ID", "Location")
stopifnot(all(required %in% names(metadata)), !anyDuplicated(metadata$Sample))

alpha_all <- imap_dfr(asv_paths, calculate_alpha) %>%
  inner_join(metadata %>% dplyr::select(all_of(required)), by = "Sample")

# Empty sequencing libraries are technical failures, not zero-diversity samples.
excluded_libraries <- alpha_all %>% filter(Read_Depth == 0)
alpha_data <- alpha_all %>%
  filter(Read_Depth > 0) %>%
  mutate(
    Marker = factor(Marker, c("16S", "ITS")),
    depth = factor(depth, c("D1", "D2", "D3")),
    across(c(Practice, Location, Core_ID), factor)
  )

alpha_long <- alpha_data %>%
  pivot_longer(c(Richness, Evenness, Shannon),
               names_to = "Metric", values_to = "Diversity") %>%
  mutate(Metric = factor(Metric, c("Richness", "Evenness", "Shannon")))


##### Step 3: Model Helpers
prepare_data <- function(data, require_maoc = TRUE) {
  data <- data %>%
    filter(is.finite(Diversity)) %>%
    drop_na(depth, Location, Core_ID, Read_Depth)
  if (require_maoc) data <- data %>% filter(is.finite(maom_c), maom_c > 0)
  
  data %>% mutate(
    Diversity_z = as.numeric(scale(Diversity)),
    log_MAOC_z = if (require_maoc) as.numeric(scale(log(maom_c))) else NA_real_
  )
}

add_read_depth <- function(data) {
  adjusted <- sd(log1p(data$Read_Depth)) > 1e-8
  data$log_Reads_z <- if (adjusted) as.numeric(scale(log1p(data$Read_Depth))) else 0
  list(data = data, adjusted = adjusted)
}

extract_effect <- function(model, term = "Diversity_z") {
  coefficients <- coef(summary(model))
  estimate <- coefficients[term, "Estimate"]
  standard_error <- coefficients[term, "Std. Error"]
  p_column <- grep("^Pr\\(", colnames(coefficients), value = TRUE)
  interval <- if (inherits(model, "merMod")) {
    suppressMessages(confint(model, parm = term, method = "Wald"))
  } else {
    confint(model, parm = term)
  }
  tibble(
    Estimate = estimate, Standard_Error = standard_error,
    CI_Lower = interval[1], CI_Upper = interval[2],
    p_value = coefficients[term, p_column]
  )
}


##### Step 4: Alpha Diversity and MAOC
# Depth and ranch are fixed effects; Core_ID is a random intercept.
fit_maoc <- function(data) {
  prepared <- add_read_depth(prepare_data(data))
  data <- prepared$data
  formula <- if (prepared$adjusted) {
    log_MAOC_z ~ Diversity_z + depth + Location + log_Reads_z + (1 | Core_ID)
  } else {
    log_MAOC_z ~ Diversity_z + depth + Location + (1 | Core_ID)
  }
  model <- lmer(formula, data = data, REML = FALSE)
  extract_effect(model) %>% mutate(
    n = nrow(data), Cores = n_distinct(data$Core_ID),
    Read_Depth_Adjusted = prepared$adjusted,
    Singular_Model = isSingular(model), .before = 1
  )
}

maoc_results <- alpha_long %>%
  group_by(Marker, Metric) %>% group_modify(~ fit_maoc(.x)) %>%
  ungroup() %>% mutate(p_adjusted = p.adjust(p_value, "BH"))


##### Step 5: Grazing Practice and Alpha Diversity
fit_practice <- function(data) {
  prepared <- add_read_depth(prepare_data(data, require_maoc = FALSE) %>%
                               drop_na(Practice))
  data <- prepared$data
  formula <- if (prepared$adjusted) {
    Diversity_z ~ Practice + depth + Location + log_Reads_z + (1 | Core_ID)
  } else {
    Diversity_z ~ Practice + depth + Location + (1 | Core_ID)
  }
  model <- lmer(formula, data = data, REML = FALSE)
  test <- as.data.frame(anova(model, type = 2, ddf = "Satterthwaite"))
  tibble(
    n = nrow(data), Cores = n_distinct(data$Core_ID),
    F_value = test["Practice", "F value"],
    p_value = test["Practice", "Pr(>F)"],
    Read_Depth_Adjusted = prepared$adjusted,
    Singular_Model = isSingular(model)
  )
}

practice_results <- alpha_long %>%
  group_by(Marker, Metric) %>% group_modify(~ fit_practice(.x)) %>%
  ungroup() %>% mutate(p_adjusted = p.adjust(p_value, "BH"))


##### Step 6: Within-Depth Sensitivity Analysis
# Each depth-specific model contains one observation per core.
fit_within_depth <- function(data) {
  prepared <- add_read_depth(prepare_data(data))
  data <- prepared$data
  formula <- if (prepared$adjusted) {
    log_MAOC_z ~ Diversity_z + Location + log_Reads_z
  } else {
    log_MAOC_z ~ Diversity_z + Location
  }
  model <- lm(formula, data = data)
  extract_effect(model) %>% mutate(
    n = nrow(data), Cores = n_distinct(data$Core_ID),
    Read_Depth_Adjusted = prepared$adjusted, .before = 1
  )
}

within_depth_results <- alpha_long %>%
  group_by(Marker, Metric, depth) %>%
  group_modify(~ fit_within_depth(.x), .keep = TRUE) %>%
  ungroup() %>% mutate(p_adjusted = p.adjust(p_value, "BH"))

overall_summary <- alpha_data %>%
  group_by(Marker) %>%
  summarize(
    n = n(),
    Richness_Mean = mean(Richness), Richness_SD = sd(Richness),
    Richness_Min = min(Richness), Richness_Max = max(Richness),
    Shannon_Mean = mean(Shannon), Shannon_SD = sd(Shannon),
    Evenness_Mean = mean(Evenness, na.rm = TRUE),
    Evenness_SD = sd(Evenness, na.rm = TRUE), .groups = "drop"
  )

depth_summary <- alpha_data %>%
  group_by(Marker, depth) %>%
  summarize(
    n = n(), Richness_Mean = mean(Richness),
    Richness_SD = sd(Richness), .groups = "drop"
  )


##### Step 7: Figures
practice_colors <- c(
  "Adaptive" = "#be7245", "Prescriptive" = "#998ec3",
  "Conventional" = "#998ec3"
)
metric_labels <- c(
  "Richness" = "Richness", "Evenness" = "Pielou's evenness",
  "Shannon" = "Shannon diversity"
)
marker_labels <- c(
  "16S" = "Bacteria and Archaea (16S)", "ITS" = "Fungi (ITS)"
)
format_q <- function(x) if_else(
  x < 0.001, "FDR q < 0.001",
  paste0("FDR q = ", formatC(x, format = "f", digits = 3))
)
practice_labels <- practice_results %>%
  transmute(Marker, Metric, Label = format_q(p_adjusted))
maoc_labels <- maoc_results %>%
  transmute(Marker, Metric, Label = format_q(p_adjusted))

practice_plot <- ggplot(alpha_long, aes(depth, Diversity, fill = Practice)) +
  geom_boxplot(outlier.shape = NA) +
  geom_text(
    data = practice_labels,
    aes(x = "D3", y = Inf, label = Label), inherit.aes = FALSE,
    hjust = 1.05, vjust = 1.3, size = 3.2
  ) +
  facet_grid2(
    rows = vars(Marker), cols = vars(Metric),
    scales = "free_y", independent = "y",
    labeller = labeller(Marker = marker_labels, Metric = metric_labels)
  ) +
  scale_fill_manual(values = practice_colors) +
  labs(x = "Soil depth", y = "Alpha diversity", fill = "Practice") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.text.y.right = element_text(angle = 270)
  )

maoc_plot <- ggplot(alpha_long, aes(maom_c, Diversity)) +
  geom_point(aes(color = Practice, shape = Practice), size = 2.8, alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "gray30") +
  geom_text(
    data = maoc_labels,
    aes(x = Inf, y = Inf, label = Label), inherit.aes = FALSE,
    hjust = 1.1, vjust = 1.3, size = 3.2
  ) +
  facet_grid2(
    rows = vars(Marker), cols = vars(Metric),
    scales = "free_y", independent = "y",
    labeller = labeller(Marker = marker_labels, Metric = metric_labels)
  ) +
  scale_color_manual(values = practice_colors) +
  labs(x = expression("MAOC (mg C g"^{-1} * " soil)"), y = "Alpha diversity",
       color = "Practice", shape = "Practice") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.text.y.right = element_text(angle = 270)
  )

forest_plot <- maoc_results %>%
  mutate(
    Marker = recode(Marker, "16S" = "Bacteria and Archaea (16S)",
                    "ITS" = "Fungi (ITS)"),
    Metric = factor(Metric, rev(levels(alpha_long$Metric))),
    Label = format_q(p_adjusted)
  ) %>%
  ggplot(aes(Estimate, Metric)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray35") +
  geom_segment(aes(x = CI_Lower, xend = CI_Upper, yend = Metric), linewidth = 0.8) +
  geom_point(shape = 21, size = 3.2, fill = "#74a9b0") +
  geom_text(aes(x = CI_Upper, label = Label), hjust = -0.15, size = 3.4) +
  facet_wrap(~Marker) +
  scale_y_discrete(labels = metric_labels) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.38))) +
  labs(x = "Standardized association with MAOC", y = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())


##### Step 8: Save Results
write_csv(alpha_data, file.path(output_dir, "alpha_diversity_values.csv"))
write_csv(excluded_libraries, file.path(output_dir, "excluded_empty_libraries.csv"))
write_csv(overall_summary, file.path(output_dir, "alpha_diversity_summary.csv"))
write_csv(depth_summary, file.path(output_dir, "alpha_diversity_by_depth_summary.csv"))
write_csv(maoc_results, file.path(output_dir, "alpha_diversity_MAOC_mixed_models.csv"))
write_csv(practice_results, file.path(output_dir, "alpha_diversity_practice_mixed_models.csv"))
write_csv(within_depth_results, file.path(output_dir, "alpha_diversity_MAOC_within_depth.csv"))

save_figure <- function(plot, name, width, height) {
  ggsave(file.path(figure_dir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 300, bg = "white")
  ggsave(file.path(figure_dir, paste0(name, ".pdf")), plot,
         width = width, height = height)
}
save_figure(practice_plot, "Alpha_diversity_by_depth_and_practice", 12, 7)
save_figure(maoc_plot, "Alpha_diversity_vs_MAOC", 12, 7)
save_figure(forest_plot, "Alpha_diversity_MAOC_forest_plot", 9, 4.5)

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
print(overall_summary)
print(depth_summary)
print(maoc_results)
print(practice_results)
print(within_depth_results)
message("Alpha-diversity analysis complete. Results saved to: ", output_dir)

