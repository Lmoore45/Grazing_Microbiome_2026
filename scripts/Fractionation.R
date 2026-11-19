library(dplyr)
library(ggplot2)
library(tidyr)
library(ggthemes)
library(ggpubr)

# Load your metadata
metadata <- read.csv("data/1.1_Metadata.csv")

### ==== BAR CHART SECTION ====

# Calculate total SOC
metadata <- metadata %>%
  select(Sample, total_c, maom_c, pom_c, Location, depth, Practice) %>% 
  drop_na()

# Prepare long format for bar chart
summary_df <- metadata %>%
  select(Location, pom_c, maom_c, total_c) %>%
  drop_na() %>%
  pivot_longer(cols = c(pom_c, maom_c, total_c),
               names_to = "Fraction",
               values_to = "Value") %>%
  mutate(
    Fraction = recode(Fraction,
                      "pom_c" = "POC",
                      "maom_c" = "MAOC",
                      "total_c" = "SOC"),
    Location = factor(Location,
                      levels = c("Site 19", "Site 5", "Site 6", "Site 16", "Site 13", "Site 15")),
    Fraction = factor(Fraction, levels = c("SOC", "MAOC", "POC"))
  )

# Summary stats for bar chart
bar_summary <- summary_df %>%
  group_by(Location, Fraction) %>%
  summarise(
    mean_value = mean(Value, na.rm = TRUE),
    se = sd(Value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# Bar chart
bar_plot <- ggplot(bar_summary, aes(x = Fraction, y = mean_value, fill = Fraction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
  geom_errorbar(aes(ymin = mean_value - se, ymax = mean_value + se),
                width = 0.2, position = position_dodge(width = 0.9)) +
  facet_wrap(~ Location, nrow = 1) +
  labs(x = "Carbon Fraction", y = "Carbon Concentration (mg C/g soil)",
       title = "Distribution of SOC Fractions Within Each Site") +
  scale_fill_manual(values = c("POC" = "gray60", "MAOC" = "#A3A990", "SOC" = "#BB946C")) +
  theme_tufte(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(color = "black", fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )

print(bar_plot)


### ==== Summaries ====

fractions_long <- metadata %>%
  select(Sample, Location, depth, Practice, total_c, maom_c, pom_c) %>%
  drop_na(Location, depth, total_c, maom_c, pom_c) %>%
  pivot_longer(
    cols = c(total_c, maom_c, pom_c),
    names_to = "Fraction",
    values_to = "Value"
  ) %>%
  mutate(
    Fraction = recode(Fraction,
                      "total_c" = "SOC",
                      "maom_c"  = "MAOC",
                      "pom_c"   = "POC"),
    # set these if you want consistent ordering:
    Location = factor(Location, levels = c("Site 19","Site 5","Site 6","Site 16","Site 13","Site 15")),
    Fraction = factor(Fraction, levels = c("SOC","MAOC","POC")),
    # if depth is like D1/D2/D3, set order explicitly:
    depth = factor(depth, levels = sort(unique(depth)))
  )


site_depth_summary <- fractions_long %>%
  group_by(Location, depth, Fraction) %>%
  summarise(
    n = n(),
    mean = mean(Value),
    sd = sd(Value),
    se = sd(Value) / sqrt(n),
    min = min(Value),
    max = max(Value),
    range = max - min,
    .groups = "drop"
  )

site_depth_summary


site_avg_range_across_depths <- site_depth_summary %>%
  group_by(Location, Fraction) %>%
  summarise(
    mean_of_depth_means = mean(mean),
    mean_range_across_depths = mean(range),
    .groups = "drop"
  )

site_avg_range_across_depths

overall_summary <- fractions_long %>%
  group_by(Fraction) %>%
  summarise(
    n = n(),
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    se = sd(Value, na.rm = TRUE) / sqrt(n),
    min = min(Value, na.rm = TRUE),
    max = max(Value, na.rm = TRUE),
    range = max - min,
    cv_percent = 100 * sd / mean,   # coefficient of variation
    .groups = "drop"
  )

overall_summary

site_depth_summary <- fractions_long %>%
  group_by(Location, depth, Fraction) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    .groups = "drop"
  )

depth_range_within_site <- site_depth_summary %>%
  group_by(Location, Fraction) %>%
  summarise(
    max_depth_mean = max(mean),
    min_depth_mean = min(mean),
    depth_range = max_depth_mean - min_depth_mean,
    .groups = "drop"
  )

depth_range_within_site

site_depth_means <- fractions_long %>%
  group_by(Location, depth, Fraction) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    .groups = "drop"
  )

site_depth_means <- site_depth_means %>%
  mutate(depth = factor(depth, levels = c("D1","D2","D3")))

depth_wide <- site_depth_means %>%
  pivot_wider(
    names_from = depth,
    values_from = mean
  )

depth_differences <- depth_wide %>%
  mutate(
    D1_minus_D2 = D1 - D2,
    D1_minus_D3 = D1 - D3,
    D2_minus_D3 = D2 - D3
  )

depth_differences

depth_wide <- site_depth_means %>%
  pivot_wider(
    names_from = depth,
    values_from = mean
  )

depth_fold_changes <- depth_wide %>%
  mutate(
    fold_D1_D2 = D1 / D2,
    fold_D1_D3 = D1 / D3,
    fold_D2_D3 = D2 / D3
  )

depth_fold_changes

depth_fold_changes <- depth_wide %>%
  mutate(
    # ratios (add epsilon in case of zeros)
    fold_D1_D2 = (D1 + 1e-6) / (D2 + 1e-6),
    fold_D1_D3 = (D1 + 1e-6) / (D3 + 1e-6),
    fold_D2_D3 = (D2 + 1e-6) / (D3 + 1e-6),
    
    # log2 fold-changes (use same epsilon)
    log2FC_D1_D2 = log2((D1 + 1e-6) / (D2 + 1e-6)),
    log2FC_D1_D3 = log2((D1 + 1e-6) / (D3 + 1e-6)),
    log2FC_D2_D3 = log2((D2 + 1e-6) / (D3 + 1e-6))
  )

fold_summary <- depth_fold_changes %>%
  group_by(Fraction) %>%
  summarise(
    mean_fold_D1_D3 = mean(fold_D1_D3, na.rm = TRUE),
    sd_fold_D1_D3 = sd(fold_D1_D3, na.rm = TRUE),
    .groups = "drop"
  )

poc_percent_overall <- metadata %>%
  drop_na(total_c, pom_c) %>%
  summarise(
    total_SOC = sum(total_c),
    total_POC = sum(pom_c),
    percent_POC_of_SOC = 100 * total_POC / total_SOC
  )

poc_percent_overall

maoc_fold_difference <- metadata %>%
  drop_na(maom_c) %>%
  summarise(
    max_MAOC = max(maom_c),
    min_MAOC = min(maom_c),
    fold_high_to_low = max_MAOC / min_MAOC
  )

maoc_fold_difference

### ==== BOX PLOT SECTION (Supplemental) ====

# Subset and pivot longer for box plot
data_long <- metadata %>%
  select(depth, Practice, pom_c, maom_c, total_c) %>%
  pivot_longer(cols = c(pom_c, maom_c, total_c),
               names_to = "Fraction", values_to = "Value") %>%
  mutate(
    Fraction = recode(Fraction,
                      "pom_c" = "POC",
                      "maom_c" = "MAOC",
                      "total_c" = "SOC"),
    Fraction = factor(Fraction, levels = c("SOC", "MAOC", "POC")),
    Practice = factor(Practice, levels = c("Conventional", "Adaptive"))
  )

# Calculate p-values (but don't plot them)
p_values <- data_long %>%
  group_by(depth, Fraction) %>%
  summarise(p_value = t.test(Value ~ Practice)$p.value,
            .groups = "drop")

# Reorder the fraction levels
data_long$Fraction <- factor(data_long$Fraction, levels = c("MAOC", "POC", "SOC"))

# Updated box plot (no p-values shown)
box_plot <- ggplot(data_long, aes(x = Practice, y = Value, fill = Fraction)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), alpha = 0.6) +
  # Optional jitter line is left commented out
  # geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
  #             size = 1.5, alpha = 0.6, color = "black") +
  facet_grid(Fraction ~ depth, scales = "free_y", switch = "y") +
  labs(x = "Practice", y = "SOC Concentration (mg C/g soil)", fill = "Carbon Fraction") +
  theme_tufte(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(color = "black", size = 1),
    strip.placement = "outside"
  ) +
  scale_fill_manual(values = c("POC" = "gray60", "MAOC" = "#A3B18A", "SOC" = "#B35806"))

# Print the plot
print(box_plot)

