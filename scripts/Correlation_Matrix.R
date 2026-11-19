# Load necessary libraries
library(dplyr)
library(tidyr)
library(corrplot)
library(ggplot2)
library(RColorBrewer)

# Read in the data
metadata <- read.csv("data/1.1_Metadata.csv")


# Ensure all columns are numeric and remove rows with NA values
metadata[metadata == "#N/A"] <- NA
metadata <- na.omit(metadata)
metadata <- metadata %>%
  mutate(
    Richness_16s = as.numeric(Richness_16s),
    Evenness_16s = as.numeric(Evenness_16s),
    Shannons_16s = as.numeric(Shannons_16s),
    Richness_ITS = as.numeric(Richness_ITS),
    Evenness_ITS = as.numeric(Evenness_ITS),
    Shannons_ITS = as.numeric(Shannons_ITS),
    Biomass = as.numeric(Biomass),
  )



# Select the columns you want to include in the correlation matrix
metadata_selected <- metadata %>%
  select(Adaptive_Score, depth_num, Biomass, `percent_clay`, `percent_sand`, `percent_silt`, pH, total_c, maom_c, pom_c, 
         total_n, maom_n, pom_n, c_n_ratio, Richness_16s, Evenness_16s, Shannons_16s, Richness_ITS, Evenness_ITS, Shannons_ITS, mean_annual_temp, mean_annual_precip,radiation,PET_mm_year, Latitude, Longitude)

# Create the correlation matrix
cor = cor(metadata_selected, use='pairwise.complete.obs')

# Basic correlation plot (upper triangle)
corrplot(cor, type = "upper")


#combine with a significance test
p.mat <- cor.mtest(metadata_selected, conf.level = .95)$p

# Define the color palette you want to use 
#col <- brewer.pal(n = 11, name = "RdBu") 

#rerun correlation plot including only significant results
corrplot(cor, type="upper", p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         method='square', tl.col= "black")
