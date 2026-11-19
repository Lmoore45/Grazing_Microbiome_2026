# Build Figure S2 from one climate record per ranch.
# Run from the project root: Rscript climate.R
# Optionally pass another metadata CSV as the first argument.
# Results are saved relative to the working directory.

library(httr)
library(jsonlite)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
metadata_file <- if (length(args) > 0) args[1] else "data/1.1_Metadata.csv"
metadata <- read_csv(metadata_file, show_col_types = FALSE) %>%
  mutate(Practice = recode(Practice, Prescriptive = "Conventional"))

climate_labels <- c(
  mean_annual_precip = "Annual Precipitation (mm)",
  mean_annual_temp = "Mean Annual Temperature (°C)",
  PET_mm_day = "Potential Evapotranspiration (mm/day)",
  radiation = "Solar Radiation (MJ/m²/day)"
)
climate_columns <- names(climate_labels)

# Each ranch should have a single management category.
ranch_metadata <- metadata %>%
  distinct(Location, Practice)

if (anyNA(ranch_metadata$Location) || anyDuplicated(ranch_metadata$Location)) {
  stop("Each ranch must have a Location and a single management category.")
}
if (!setequal(ranch_metadata$Practice, c("Adaptive", "Conventional"))) {
  stop("Practice must identify adaptive and conventional ranches.")
}

fetch_power_data <- function(lat, lon, start_date, end_date) {
  base_url <- "https://power.larc.nasa.gov/api/temporal/daily/point"
  response <- GET(
    url = base_url,
    query = list(
      parameters = "T2M,T2M_MAX,T2M_MIN,PRECTOTCORR,ALLSKY_SFC_SW_DWN",
      community = "AG",
      longitude = lon,
      latitude = lat,
      start = gsub("-", "", start_date),
      end = gsub("-", "", end_date),
      format = "JSON"
    )
  )
  stop_for_status(response)
  data <- fromJSON(content(response, as = "text", encoding = "UTF-8"))$properties$parameter
  
  tibble(
    Date = as.Date(names(data$T2M), format = "%Y%m%d"),
    Temperature = unlist(data$T2M),
    Max_Temperature = unlist(data$T2M_MAX),
    Min_Temperature = unlist(data$T2M_MIN),
    Precipitation = unlist(data$PRECTOTCORR),
    Radiation = unlist(data$ALLSKY_SFC_SW_DWN)
  )
}

calculate_pet <- function(t_mean, t_max, t_min, solar_radiation) {
  0.0023 * solar_radiation * (t_mean + 17.8) * sqrt(t_max - t_min)
}

# Reuse the published climate values in metadata when available.
# Otherwise fetch the original 2023 daily climate data once per ranch.
if (all(climate_columns %in% names(metadata))) {
  climate_data <- metadata %>%
    select(Location, Practice, all_of(climate_columns)) %>%
    distinct()
} else {
  unique_locations <- metadata %>%
    mutate(across(c(Latitude, Longitude), as.numeric)) %>%
    group_by(Location) %>%
    summarize(
      Latitude = mean(Latitude, na.rm = TRUE),
      Longitude = mean(Longitude, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (any(!is.finite(unique_locations$Latitude)) ||
      any(!is.finite(unique_locations$Longitude))) {
    stop("Each ranch must have latitude and longitude to fetch climate data.")
  }
  
  daily_climate_data <- unique_locations %>%
    rowwise() %>%
    mutate(Climate_Data = list(fetch_power_data(
      Latitude, Longitude, "2023-01-01", "2023-12-31"
    ))) %>%
    ungroup() %>%
    unnest(Climate_Data)
  
  daily_variables <- c("Temperature", "Max_Temperature", "Min_Temperature",
                       "Precipitation", "Radiation")
  if (anyNA(daily_climate_data$Date) ||
      any(!is.finite(as.matrix(daily_climate_data[daily_variables]))) ||
      any(as.matrix(daily_climate_data[daily_variables]) == -999)) {
    stop("NASA POWER returned missing daily climate values.")
  }
  days_by_ranch <- daily_climate_data %>%
    group_by(Location) %>%
    summarize(days = n_distinct(Date), .groups = "drop")
  if (nrow(days_by_ranch) != nrow(ranch_metadata) || any(days_by_ranch$days != 365L)) {
    stop("Expected 365 days of climate data for each ranch in 2023.")
  }
  
  mean_annuals <- daily_climate_data %>%
    group_by(Location) %>%
    summarize(
      mean_annual_temp = mean(Temperature),
      annual_high_temp = mean(Max_Temperature),
      annual_low_temp = mean(Min_Temperature),
      mean_annual_precip = sum(Precipitation),
      radiation = mean(Radiation),
      .groups = "drop"
    ) %>%
    mutate(PET_mm_day = calculate_pet(
      mean_annual_temp, annual_high_temp, annual_low_temp, radiation
    ))
  
  climate_data <- ranch_metadata %>%
    left_join(mean_annuals, by = "Location") %>%
    select(Location, Practice, all_of(climate_columns))
}

# Check that collapsing the soil-sample rows yielded one complete record per ranch.
if (anyDuplicated(climate_data$Location)) {
  stop("Climate values differ among soil samples from the same ranch.")
}
if (any(!is.finite(as.matrix(climate_data[climate_columns])))) {
  stop("The ranch-level climate table contains missing or non-finite values.")
}

climate_long <- climate_data %>%
  mutate(Practice = factor(Practice, levels = c("Adaptive", "Conventional"))) %>%
  pivot_longer(all_of(climate_columns), names_to = "Climate_Variable", values_to = "Value") %>%
  mutate(Climate_Variable = factor(
    Climate_Variable, levels = climate_columns, labels = unname(climate_labels)
  ))

# Compare management categories using ranches as the independent units.
# The plotted comparisons are descriptive.
climate_tests <- climate_long %>%
  group_by(Climate_Variable) %>%
  summarize(
    p = wilcox.test(Value ~ Practice, exact = anyDuplicated(Value) == 0L)$p.value,
    .groups = "drop"
  )
print(climate_tests)

figure_s2 <- ggplot(climate_long, aes(x = Practice, y = Value, color = Practice)) +
  geom_boxplot(aes(fill = Practice), width = 0.55, alpha = 0.30, outlier.shape = NA) +
  geom_point(size = 2.8) +
  facet_wrap(~Climate_Variable, scales = "free_y", ncol = 2) +
  scale_color_manual(
    values = c("Adaptive" = "#BE7245", "Conventional" = "#542788"),
    aesthetics = c("colour", "fill")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(x = NULL, y = NULL)

if (interactive()) print(figure_s2)

output_dir <- file.path("outputs", "climate")
figure_dir <- file.path("figures", "climate")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(climate_data, file.path(output_dir, "ranch_climate.csv"))
write_csv(climate_tests, file.path(output_dir, "climate_tests.csv"))
ggsave(file.path(figure_dir, "Figure_S2.pdf"), plot = figure_s2,
       width = 8, height = 6, units = "in")
