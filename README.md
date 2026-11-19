# Identifying Candidate Microbial Indicators of Mineral-Associated Organic Carbon in Semi-Arid Rangelands

This repository contains the R scripts and input data for analyses of microbial communities, soil properties, climate, and grazing management across six ranches in Colorado and Wyoming. The analyses examine how bacterial, archaeal, and fungal communities relate to mineral-associated organic carbon (MAOC).

## Data

| File | Contents |
| --- | --- |
| [`1.1_Metadata.csv`](data/1.1_Metadata.csv) | Sample IDs, soil depths, grazing categories, soil and climate measurements, and microbial community summaries |
| [`1.2_16S_ASV_Table.csv`](data/1.2_16S_ASV_Table.csv) | Bacterial and archaeal 16S amplicon sequence variant (ASV) counts and taxonomy |
| [`1.3_ITS_ASV_Table.csv`](data/1.3_ITS_ASV_Table.csv) | Fungal ITS ASV counts and taxonomy |
| [`PCA_Figure3_samples.csv`](data/PCA_Figure3_samples.csv) | The 75 samples used in the Figure 3 principal component analysis |

Raw sequencing reads are available under [NCBI BioProject PRJNA1271444](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1271444). Processed ASV tables are archived on [Zenodo](https://doi.org/10.5281/zenodo.19700145).

Public coordinates are county centroids to protect the locations of working ranches. Use the climate measurements supplied in the metadata for the analyses. [`Climate_API_Pull.R`](scripts/Climate_API_Pull.R) demonstrates the NASA POWER data retrieval using the public coordinates.

## Setup

Open [`Manuscript_Code.Rproj`](Manuscript_Code.Rproj) in RStudio or run the scripts with `Rscript` from the repository root. Install the required packages in R:

```r
install.packages(c(
  "corrplot", "factoextra", "FactoMineR", "ggh4x", "ggpubr", "ggthemes",
  "httr", "jsonlite", "lme4", "lmerTest", "paletteer", "patchwork",
  "permute", "pls", "RColorBrewer", "scales", "tidyverse", "vegan"
), repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install("mixOmics")
```

## Analysis scripts

| Script | Analysis and outputs |
| --- | --- |
| [`Fractionation.R`](scripts/Fractionation.R) | Soil organic carbon, MAOC, and particulate organic carbon summaries and plots by ranch and soil depth |
| [`PCA.R`](scripts/PCA.R) | Principal component analysis and biplots for Figure 3 |
| [`beta-diveristy_MANTEL.R`](scripts/beta-diveristy_MANTEL.R) | 16S and ITS community composition, NMDS ordinations, and Mantel tests for Figure 4 |
| [`sPLS.R`](scripts/sPLS.R) | Depth-adjusted sparse partial least squares (sPLS) models for 16S and ITS ASVs, with leave-one-ranch-out validation and Figure 5 panels |
| [`Within_Depth_sPLS.R`](scripts/Within_Depth_sPLS.R) | sPLS sensitivity analysis within each soil depth |
| [`PLS.R`](scripts/PLS.R) | Partial least squares regression comparing microbial and environmental predictors of surface-soil MAOC. Produces the Figure 6 plot and `PLS_performance.csv` |
| [`Climate.R`](scripts/Climate.R) | Ranch-level climate comparisons and Figure S2 |
| [`Correlation_Matrix.R`](scripts/Correlation_Matrix.R) | Descriptive pairwise Pearson correlation plots for Figure S3 |
| [`alpha_diversity.R`](scripts/alpha_diversity.R) | 16S and ITS alpha diversity summaries, models, and Figure S5 plots |

For example, run the pooled sPLS before the within-depth analysis to include the comparison of selected ASVs:

```sh
Rscript scripts/sPLS.R
Rscript scripts/Within_Depth_sPLS.R
Rscript scripts/PLS.R
```

Scripts that save results write to `outputs/` and `figures/`. View plots from the fractionation, PCA, and correlation scripts in the RStudio Plots pane or in `Rplots.pdf` after a command-line run.

## License and contact

Code is available under the [MIT License](LICENSE). For questions about the data or analysis, contact Laura.Moore@colostate.edu or open an issue in this repository.
