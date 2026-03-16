# Microbial Communities Function as Ecological Indicators of Mineral-Associated Organic Carbon in Semi-Arid Rangelands
*Companion code repository for manuscript submission*

This repository supports the analyses and figures described in the manuscript, providing reproducible workflows for data processing, sPLS/PLS modeling, and supplementary figure generation using microbial and environmental data.

## Manuscript Abstract  
Mineral-associated organic carbon (MAOC) is the largest and typically most stable fraction of soil organic carbon (SOC), underpinning long-term carbon storage. Adaptive grazing management has been proposed as a strategy to enhance MAOC, yet its effectiveness remains inconsistently evaluated across rangeland ecosystems and climatic contexts. This uncertainty is particularly pronounced in semi-arid rangelands, where low precipitation and poorly constrained controls on SOC stabilization complicate interpretation of MAOC dynamics. Because microbial processing of organic inputs can contribute to MAOC formation, we tested whether insights from bacterial and fungal taxonomic information correspond with differences in MAOC across a grazing management and climatic gradient. To evaluate this possibility, we combined 16S rRNA and ITS amplicon sequencing with soil physicochemical measurements and climate data from six ranches in Colorado and Wyoming spanning an adaptive–conventional grazing gradient. We assessed microbial, environmental, management, and climatic variables associated with MAOC concentration and found that MAOC concentrations were higher at sites with adaptive grazing management. Supervised machine learning identified 23 bacterial/archaeal and 22 fungal ASVs associated with MAOC, including lineages within Actinomycetota, Verrucomicrobiota, and Ascomycota. Additionally, incorporating microbiome variables significantly improved model performance beyond soil physicochemical variables alone, increasing explained variance in MAOC from 71.1% to 84.4%. Notably, the relative abundances of MAOC-associated taxa ranked among the strongest contributors to MAOC variation and exceeded the explanatory value of broad diversity metrics or management. Together, these results indicate that microbial community composition captures structured ecological variation associated with MAOC not fully resolved by conventional soil and climatic measurements, highlighting the potential value of microbiome-informed approaches for understanding SOC stabilization in semi-arid rangelands.<img width="468" height="345" alt="image" src="https://github.com/user-attachments/assets/0f3306d1-16df-473c-bfae-9e232d293ef8" />




---
## Disclosure

Due to rancher privacy agreements, precise geographic coordinates used in the manuscript analyses are not included in this repository. The version of the code shared here uses county-level centroids as placeholders.
As a result, certain figures (e.g., climate figures, PLS, RDA) may differ slightly from those in the published manuscript. However, all primary findings and interpreations remain consistent. 

---

## Repository Structure and How to Use 

The code is organized by numbered sections. To fully reproduce the analyses, follow the instructions below for obtaining and using the appropriate input files and scripts.

### 1.0 – Data Preparation  
This section loads metadata and OTU tables.

You will need to download the following input files from the repository or dataset archive:
- `1.1_Metadata.csv` – Core metadata table  
- `1.2_16S_OTU_table.csv.zip` – Rarefied 16S OTU table (zipped)
- `1.3_ITS_OTU_table.txt` – Raw ITS OTU table  

---

### 2.0 – Climate Data  
This section uses an API script to pull climate variables for each site.

Due to rancher privacy agreements, precise site coordinates are not publicly provided. For reproducibility and transparency, the metadata included contains county centroid coordinates for each site, which differ from the actual coordinates used in the original analyses.
The script `2.0_climate_API_Pull.R` demonstrates the API call.

---

### 6.1, 6.2 – Sparse Partial Least Squares (sPLS)  
These scripts perform sPLS regression on both 16S and ITS data to identify microbial taxa most predictive of MAOC.

To run these scripts, be sure to also download and include:
- `6.0_VIP.R` – Custom function script used to calculate VIP scores from sPLS models.

---

### 7.0 – Partial Least Squares (PLS)  
Similar to sPLS but using all features (no sparsity constraint). Also requires:
- `6.0_VIP.R` – Shared VIP function used in both sPLS and PLS scripts.
- Note: This analysis originally used precise goegraphic coordinates, but the version shared here substitues county-level centroids to protect rancher privacy. As a result, minor differences may appear in spatially influenced model outputs, though overall patterns remain consisten. 

---

### S1.0, S2.1, S2.2 – Supplementary Figures  
These scripts generate additional exploratory figures and supporting plots used in the supplemental materials of the manuscript.

---


## License  
This repository is released under the [MIT License](LICENSE). You are free to use, modify, and share this code with attribution.

---

## Questions or Contributions?  
For inquiries related to the manuscript or code, please reach out to Laura.Moore@colostate.edu, or open an issue here in the repository.

