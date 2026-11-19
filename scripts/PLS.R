################################################################################
### Surface-soil MAOC: microbial and environmental PLS regression
###
### Compares alpha diversity, training-selected ASVs, and environmental variables
### using a common complete-case cohort at 0-15 cm. Validation evaluates:
###   1. Additional triangles within sampled ranches: 10 repeats of stratified CV.
###   2. An excluded ranch: six leave-one-ranch-out fits.
### ASV filtering, CLR transformation, selection and tuning are training-only,
### including repeated ASV selection inside PLS inner validation. No preselected
### ASV files are imported. Selection uses eligible marker-specific training
### samples; PLS fitting and evaluation use the common complete-case cohort.
###
### Run from the manuscript project root: source("PLS.R")
### Outputs (CSV): outputs/PLS/{PLS_performance, PLS_predictions,
###                           PLS_ASV_selections, PLS_VIP}.csv
### Figure 6 (PDF/PNG): figures/PLS/Figure_6.{pdf,png}
### Rerunning replaces these generated files; source data are unchanged.
###
### MAOC is natural-log transformed. Repeated predictions are not independent
### replicates; variation across repeats is descriptive, not a confidence interval.
### Figure 6 shows conditional full-model VIPs from ranch-held-out fits, retaining
### ASVs included in at least three folds. VIPs are not unique variance fractions.
### Nested validation: Varma & Simon (2006), doi:10.1186/1471-2105-7-91.
################################################################################

##### Step 1: Load Libraries and Define Settings
library(tidyverse)
library(mixOmics)
library(pls)

seed <- 42L
repeats <- 10L
set.seed(seed)
prevalence <- 0.10
total_reads <- 10
pseudocount <- 0.5
spls_components <- 1:3
keepX_grid <- c(5, 10, 15, 25)
taxa_per_marker <- 5L
max_pls_components <- 5L
desired_outer_folds <- 5L
desired_inner_folds <- 3L
markers <- c("16S", "ITS")
schemes <- c("Within ranch", "Ranch held out")

alpha_vars <- c("Richness_16s", "Evenness_16s", "Shannons_16s",
                "Richness_ITS", "Evenness_ITS", "Shannons_ITS")
baseline_vars <- c("pom_c", "c_n_ratio", "pH", "percent_clay", "Adaptive_Score",
                   "mean_annual_temp", "mean_annual_precip")
summary_vars <- c(alpha_vars, "Biomass")
models <- c("Alpha diversity only", "ASVs only", "Alpha diversity + ASVs",
            "Environmental baseline", "Environmental + ASVs",
            "Environmental + summaries", "Full integrated model")
environmental_models <- models[4:7]

asv_paths <- c("16S" = "data/1.2_16S_ASV_Table.csv",
               "ITS" = "data/1.3_ITS_ASV_Table.csv")
output_dir <- file.path("outputs", "PLS")
figure_dir <- file.path("figures", "PLS")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)


##### Step 2: Load and Align Surface-Soil Samples
read_asv_data <- function(path) {
  tab <- read_csv(path, show_col_types = FALSE, name_repair = "check_unique")
  id <- intersect(c("ASV", "OTU"), names(tab))
  stopifnot(length(id) == 1, !anyNA(tab[[id]]), !anyDuplicated(tab[[id]]))
  samples <- setdiff(names(tab), c(id, "Taxon"))
  counts <- t(as.matrix(tab[, samples]))
  storage.mode(counts) <- "numeric"
  colnames(counts) <- as.character(tab[[id]])
  stopifnot(all(is.finite(counts)), all(counts >= 0))
  taxonomy <- if ("Taxon" %in% names(tab)) tab$Taxon else rep(NA_character_, nrow(tab))
  list(counts = counts, taxonomy = setNames(taxonomy, tab[[id]]))
}

raw <- lapply(asv_paths, read_asv_data)
meta <- read_csv("data/1.1_Metadata.csv", show_col_types = FALSE)
required <- unique(c("Sample", "Location", "Core_ID", "depth", "maom_c",
                     baseline_vars, summary_vars))
stopifnot(all(c(required, "Triangle") %in% names(meta)), !anyDuplicated(meta$Sample))
numeric_vars <- unique(c("maom_c", baseline_vars, summary_vars))
stopifnot(all(vapply(meta[numeric_vars], is.numeric, logical(1))))
common_samples <- Reduce(intersect, lapply(raw, function(x) rownames(x$counts)))
dat <- meta %>%
  filter(depth == "D1", Sample %in% common_samples, maom_c > 0) %>%
  drop_na(all_of(required)) %>%
  mutate(across(c(Sample, Location, Core_ID, Triangle), as.character)) %>%
  as.data.frame()
stopifnot(!anyNA(dat$Triangle), all(nzchar(trimws(dat$Triangle))),
          all(is.finite(as.matrix(dat[numeric_vars]))))
dat$Triangle_Group <- as.character(interaction(dat$Location, dat$Triangle,
                                               drop = TRUE, lex.order = TRUE))
stopifnot(!anyDuplicated(dat$Core_ID), !anyDuplicated(dat$Triangle_Group))
X <- lapply(raw, function(x) x$counts[dat$Sample, , drop = FALSE])
stopifnot(all(vapply(X, function(x) all(rowSums(x) > 0), logical(1))))
y <- log(dat$maom_c)
ranch <- dat$Location
triangle <- dat$Triangle_Group
n <- nrow(dat)
counts_by_ranch <- table(ranch)
stopifnot(length(counts_by_ranch) == 6, min(counts_by_ranch) >= 2, sd(y) > 0)
outer_k <- as.integer(min(desired_outer_folds, min(counts_by_ranch)))
message("PLS cohort: ", n, " samples, six ranches; ", outer_k, " within-ranch folds.")
if (n != 25) warning("The cohort differs from the manuscript's 25 samples; check eligibility.")

# Selection eligibility requires only the marker data, MAOC and sampling IDs.
# Missing alpha diversity, DNA yield or environmental covariates does not remove
# a sample from this upstream selection pool. PLS eligibility above is unchanged.
selection_pools <- setNames(lapply(markers, function(m) {
  pm <- meta %>%
    filter(depth == "D1", Sample %in% rownames(raw[[m]]$counts),
           Location %in% unique(ranch), maom_c > 0) %>%
    drop_na(Sample, Location, Core_ID, maom_c) %>%
    mutate(across(c(Sample, Location, Core_ID, Triangle), as.character)) %>%
    as.data.frame()
  stopifnot(!anyNA(pm$Triangle), all(nzchar(trimws(pm$Triangle))),
            all(is.finite(pm$maom_c)))
  counts <- raw[[m]]$counts[pm$Sample, , drop = FALSE]
  positive <- rowSums(counts) > 0
  pm <- pm[positive, , drop = FALSE]
  counts <- counts[positive, , drop = FALSE]
  pm$Triangle_Group <- as.character(interaction(pm$Location, pm$Triangle,
                                                drop = TRUE, lex.order = TRUE))
  stopifnot(!anyDuplicated(pm$Core_ID), !anyDuplicated(pm$Triangle_Group),
            all(dat$Sample %in% pm$Sample))
  message(m, " selection pool: ", nrow(pm), " eligible surface samples; ",
          "PLS fitting/testing cohort: ", n, ".")
  list(sample = pm$Sample, ranch = pm$Location, triangle = pm$Triangle_Group,
       counts = counts, y = log(pm$maom_c))
}), markers)


##### Step 3: Define Training-Only Preprocessing and Validation Functions
subset_seed <- function(ids, salt = "", sample_ids = dat$Sample) {
  chars <- utf8ToInt(paste(c(as.character(seed), salt, sort(sample_ids[ids])), collapse = "|"))
  h <- 17
  for (ch in chars) h <- (h * 131 + ch) %% 2147483646
  as.integer(h + 1)
}
# Allocate each ranch's triangles as evenly as possible across folds, with
# global balancing for small/singleton ranches in inner training subsets.
make_folds <- function(ids, k, split_seed, strict = FALSE,
                       ranch_ids = ranch, triangle_ids = triangle) {
  stopifnot(length(ids) == length(unique(ids)), length(ids) >= k, k >= 2L)
  set.seed(split_seed)
  assignment <- rep(NA_integer_, length(ids))
  sizes <- integer(k)
  for (r in sort(unique(ranch_ids[ids]))) {
    at <- which(ranch_ids[ids] == r)
    at <- at[sample.int(length(at))]
    start <- order(sizes, stats::runif(k))
    labels <- rep(start, length.out = length(at))
    assignment[at] <- labels
    sizes <- sizes + tabulate(labels, nbins = k)
  }
  if (anyNA(assignment) || any(tabulate(assignment, nbins = k) == 0L)) {
    stop("An empty validation fold was generated.")
  }
  result <- lapply(seq_len(k), function(f) {
    train <- ids[assignment != f]
    test <- ids[assignment == f]
    stopifnot(!length(intersect(triangle_ids[train], triangle_ids[test])))
    if (strict) {
      stopifnot(setequal(unique(ranch_ids[train]), unique(ranch_ids[ids])),
                setequal(unique(ranch_ids[test]), unique(ranch_ids[ids])))
    }
    list(train = train, test = test)
  })
  result
}
prepare_clr <- function(marker, train, test, counts = X[[marker]]) {
  tr <- counts[train, , drop = FALSE]
  te <- counts[test, , drop = FALSE]
  keep <- colMeans(tr > 0) >= prevalence & colSums(tr) >= total_reads
  if (!any(keep)) stop("No ASVs pass training filters for ", marker)
  transform <- function(z) {
    a <- log(z[, keep, drop = FALSE] + pseudocount)
    sweep(a, 1, rowMeans(a), "-")
  }
  tr <- transform(tr)
  te <- transform(te)
  variable <- apply(tr, 2, stats::sd) > 0
  if (!any(variable)) stop("No variable CLR features for ", marker)
  list(train = tr[, variable, drop = FALSE], test = te[, variable, drop = FALSE])
}
predict_spls <- function(fit, data, nc) drop(predict(fit, data)$predict[, 1, nc])
# Match the validation target in both the PLS inner loop and selector tuning.
validation_splits <- function(ids, scheme, salt, ranch_ids = ranch,
                              triangle_ids = triangle, sample_ids = dat$Sample) {
  if (scheme == "Within ranch") {
    return(make_folds(ids, min(desired_inner_folds, length(ids)),
                      subset_seed(ids, salt, sample_ids),
                      ranch_ids = ranch_ids, triangle_ids = triangle_ids))
  }
  sites <- sort(unique(ranch_ids[ids]))
  if (length(sites) < 2L) stop("Too few ranches for ranch-held-out inner validation.")
  lapply(sites, function(site) {
    train <- ids[ranch_ids[ids] != site]
    test <- ids[ranch_ids[ids] == site]
    stopifnot(!length(intersect(ranch_ids[train], ranch_ids[test])),
              !length(intersect(triangle_ids[train], triangle_ids[test])))
    list(train = train, test = test)
  })
}

##### Step 4: Select ASVs Within Each Training Set
# sPLS tuning uses eligible samples from training ranches/triangles only. For
# within-ranch validation, every current outer/inner test triangle is removed
# from the larger pool, including any linked observations sharing that triangle.
selection_cache <- new.env(parent = emptyenv())
select_asvs <- function(marker, train, scheme) {
  pool <- selection_pools[[marker]]
  excluded <- setdiff(seq_len(n), train)
  eligible <- if (scheme == "Ranch held out") {
    pool$ranch %in% unique(ranch[train])
  } else {
    !pool$triangle %in% triangle[excluded]
  }
  pool_train <- which(eligible)
  stopifnot(length(pool_train) >= 3,
            !length(intersect(pool$triangle[pool_train], triangle[excluded])))
  if (scheme == "Ranch held out") {
    stopifnot(!length(intersect(pool$ranch[pool_train], unique(ranch[excluded]))))
  }
  # Any common-cohort samples retained in the larger pool must be training rows.
  matched <- match(pool$sample[pool_train], dat$Sample)
  stopifnot(all(is.na(matched) | matched %in% train))
  key <- paste(scheme, marker, paste(pool$sample[pool_train], collapse = ","), sep = "|")
  if (exists(key, envir = selection_cache, inherits = FALSE)) {
    return(get(key, envir = selection_cache, inherits = FALSE))
  }
  # Tune on the eligible training pool, never on outer/inner PLS test samples.
  folds <- validation_splits(pool_train, scheme, paste0("sPLS_", marker),
                             pool$ranch, pool$triangle, pool$sample)
  prepared <- lapply(folds, function(sp) {
    stopifnot(all(sp$train %in% pool_train), all(sp$test %in% pool_train),
              !length(intersect(pool$triangle[sp$train], triangle[excluded])),
              !length(intersect(pool$triangle[sp$test], triangle[excluded])))
    list(train = sp$train, test = sp$test,
         clr = prepare_clr(marker, sp$train, sp$test, pool$counts))
  })
  grid <- expand.grid(ncomp = spls_components, keepX = keepX_grid)
  grid$RMSE <- NA_real_
  for (g in seq_len(nrow(grid))) {
    observed <- predicted <- numeric()
    for (sp in prepared) {
      nc <- min(grid$ncomp[g], nrow(sp$clr$train) - 1L, ncol(sp$clr$train))
      kx <- min(grid$keepX[g], ncol(sp$clr$train))
      fit <- mixOmics::spls(sp$clr$train, pool$y[sp$train], ncomp = nc,
                            keepX = rep(kx, nc), mode = "regression", scale = TRUE)
      observed <- c(observed, pool$y[sp$test])
      predicted <- c(predicted, predict_spls(fit, sp$clr$test, nc))
    }
    if (any(!is.finite(predicted))) stop("Non-finite sPLS tuning predictions.")
    grid$RMSE[g] <- sqrt(mean((observed - predicted)^2))
  }
  best <- grid[order(grid$RMSE, grid$ncomp, grid$keepX)[1L], , drop = FALSE]
  clr <- prepare_clr(marker, pool_train, pool_train, pool$counts)
  nc <- min(best$ncomp, nrow(clr$train) - 1L, ncol(clr$train))
  kx <- min(best$keepX, ncol(clr$train))
  fit <- mixOmics::spls(clr$train, pool$y[pool_train], ncomp = nc,
                        keepX = rep(kx, nc), mode = "regression", scale = TRUE)
  loadings <- fit$loadings$X[, seq_len(nc), drop = FALSE]
  ids <- rownames(loadings)[rowSums(abs(loadings)) > 0]
  vip <- mixOmics::vip(fit)
  if (is.null(dim(vip))) vip <- matrix(vip, ncol = 1L,
                                       dimnames = list(colnames(clr$train), "comp1"))
  ans <- data.frame(ASV = ids,
                    VIP = apply(vip[ids, seq_len(nc), drop = FALSE], 1, max))
  ans <- head(ans[order(-ans$VIP, method = "radix"), , drop = FALSE], taxa_per_marker)
  if (!nrow(ans)) stop("No sPLS ASVs selected for ", marker)
  ans$Selection_n <- length(pool_train)
  ans$Selection_Samples <- paste(pool$sample[pool_train], collapse = ";")
  assign(key, ans, envir = selection_cache)
  ans
}

make_taxa <- function(train, test, selection) {
  tr <- te <- used <- list()
  for (marker in markers) {
    clr <- prepare_clr(marker, train, test)
    # Match the original PLS rule: rank top five first, then retain those passing
    # the PLS training-only filters. Do not backfill with additional ASVs.
    used[[marker]] <- selection[[marker]][
      selection[[marker]]$ASV %in% colnames(clr$train), , drop = FALSE]
    ids <- used[[marker]]$ASV
    a <- as.data.frame(clr$train[, ids, drop = FALSE])
    b <- as.data.frame(clr$test[, ids, drop = FALSE])
    names(a) <- names(b) <- if (length(ids)) paste0("ASV_", marker, "_", ids) else character()
    if (any(make.names(names(a)) != names(a))) stop("Non-formula-safe ASV names.")
    tr[[marker]] <- a
    te[[marker]] <- b
  }
  list(train = do.call(cbind, unname(tr)), test = do.call(cbind, unname(te)), selected = used)
}
##### Step 5: Define PLS Predictor Sets and Model Fitting
matrix_for <- function(ids, taxa, model) {
  alpha <- as.data.frame(dat[ids, alpha_vars, drop = FALSE])
  baseline <- as.data.frame(dat[ids, baseline_vars, drop = FALSE])
  summaries <- as.data.frame(dat[ids, summary_vars, drop = FALSE])
  ans <- switch(model,
                "Alpha diversity only" = alpha,
                "ASVs only" = taxa,
                "Alpha diversity + ASVs" = cbind(alpha, taxa),
                "Environmental baseline" = baseline,
                "Environmental + ASVs" = cbind(baseline, taxa),
                "Environmental + summaries" = cbind(baseline, summaries),
                "Full integrated model" = cbind(baseline, summaries, taxa))
  if (is.null(ans)) stop("Unknown model: ", model)
  rownames(ans) <- NULL
  ans
}
clean <- function(train, test) {
  keep <- vapply(train, stats::sd, numeric(1), na.rm = TRUE) > 1e-12
  if (!any(keep)) stop("No variable PLS predictors.")
  list(train = train[, keep, drop = FALSE], test = test[, keep, drop = FALSE])
}
fit_pls <- function(X, yy, nc) pls::plsr(
  y ~ ., data = data.frame(y = yy, X), ncomp = nc,
  scale = TRUE, method = "oscorespls", validation = "none")
predict_pls <- function(fit, X, nc) {
  ans <- drop(stats::predict(fit, newdata = X, ncomp = nc))
  if (any(!is.finite(ans))) stop("Non-finite PLS predictions.")
  ans
}
final_vip <- function(fit) {
  SS <- as.numeric(fit$Yloadings)^2 * colSums(fit$scores^2)
  W <- fit$loading.weights
  weighted <- sweep(W^2, 2, SS / colSums(W^2), "*")
  ans <- sqrt(nrow(W) * rowSums(weighted) / sum(SS))
  names(ans) <- rownames(W)
  ans
}

# Inner PLS validation reselects ASVs before this function is called. Component
# candidates share those inner splits; the outer test response is never used.
fit_model <- function(train, test, model, outer_taxa, inner, inner_taxa) {
  outer_X <- clean(matrix_for(train, outer_taxa$train, model),
                   matrix_for(test, outer_taxa$test, model))
  prepared <- lapply(seq_along(inner), function(j) {
    sp <- inner[[j]]
    list(indices = sp, X = clean(matrix_for(sp$train, inner_taxa[[j]]$train, model),
                                 matrix_for(sp$test, inner_taxa[[j]]$test, model)))
  })
  maximum <- min(max_pls_components, ncol(outer_X$train),
                 min(vapply(prepared, function(z) length(z$indices$train) - 1L, integer(1))))
  stopifnot(maximum >= 1)
  tuning <- map_dfr(seq_len(maximum), function(nc) {
    validation <- map_dfr(prepared, function(sp) {
      effective <- min(nc, ncol(sp$X$train), length(sp$indices$train) - 1L)
      fit <- fit_pls(sp$X$train, y[sp$indices$train], effective)
      tibble(observed = y[sp$indices$test],
             predicted = predict_pls(fit, sp$X$test, effective))
    })
    tibble(ncomp = nc, RMSE = sqrt(mean((validation$observed - validation$predicted)^2)))
  }) %>% arrange(RMSE, ncomp)
  best <- tuning$ncomp[1]
  fit <- fit_pls(outer_X$train, y[train], best)
  list(prediction = predict_pls(fit, outer_X$test, best),
       vip = final_vip(fit), ncomp = best)
}


##### Step 6: Run Both Validation Designs
# Predefine all outer splits. Triangle IDs include ranch; each triangle contributes
# one D1 sample. Inner triangle folds balance ranches as sample counts permit.
outer_plan <- list(
  "Within ranch" = lapply(seq_len(repeats), function(rep) {
    make_folds(seq_len(n), outer_k,
               as.integer((as.double(seed) + rep * 1009) %% 2147483646 + 1), strict = TRUE)
  }),
  "Ranch held out" = list(validation_splits(seq_len(n), "Ranch held out", "outer"))
)
predictions <- selections <- full_vips <- list()
for (scheme in schemes) {
  for (rep in seq_along(outer_plan[[scheme]])) {
    folds <- outer_plan[[scheme]][[rep]]
    for (f in seq_along(folds)) {
      tr <- folds[[f]]$train
      te <- folds[[f]]$test
      held_ranch <- if (scheme == "Ranch held out") unique(ranch[te]) else NA_character_
      message(scheme, ": repeat ", rep, ", fold ", f, "/", length(folds))
      stopifnot(!length(intersect(triangle[tr], triangle[te])))
      if (scheme == "Ranch held out") stopifnot(!length(intersect(ranch[tr], ranch[te])))
      
      chosen <- setNames(lapply(markers, function(m) select_asvs(m, tr, scheme)), markers)
      outer_taxa <- make_taxa(tr, te, chosen)
      for (m in markers) {
        used <- outer_taxa$selected[[m]]
        if (!nrow(used)) next
        selections[[paste(scheme, rep, f, m)]] <- tibble(
          Scheme = scheme, Repeat = rep, Fold = f, Held_Out_Ranch = held_ranch,
          Marker = m, ASV = used$ASV, sPLS_VIP = used$VIP,
          Selection_n = used$Selection_n, Selection_Samples = used$Selection_Samples,
          Taxon = unname(raw[[m]]$taxonomy[used$ASV]))
      }
      inner <- validation_splits(tr, scheme, "PLS_inner")
      inner_taxa <- lapply(inner, function(sp) {
        stopifnot(all(sp$train %in% tr), all(sp$test %in% tr),
                  !length(intersect(triangle[sp$train], triangle[te])))
        if (scheme == "Ranch held out") {
          stopifnot(!length(intersect(ranch[sp$train], ranch[te])))
        }
        selected <- setNames(lapply(markers, function(m) select_asvs(m, sp$train, scheme)), markers)
        make_taxa(sp$train, sp$test, selected)
      })
      for (model in models) {
        fitted <- fit_model(tr, te, model, outer_taxa, inner, inner_taxa)
        key <- paste(scheme, rep, f, model)
        predictions[[key]] <- tibble(
          Scheme = scheme, Repeat = rep, Fold = f, Model = model,
          Sample = dat$Sample[te], Ranch = ranch[te], Triangle_Group = triangle[te],
          Held_Out_Ranch = held_ranch, Observed_Log_MAOC = y[te],
          Predicted_Log_MAOC = fitted$prediction, ncomp = fitted$ncomp)
        if (scheme == "Ranch held out" && model == "Full integrated model") {
          full_vips[[key]] <- tibble(Fold = f, Held_Out_Ranch = held_ranch,
                                     Variable = names(fitted$vip), VIP = unname(fitted$vip))
        }
      }
    }
  }
}
prediction_results <- bind_rows(predictions)
asv_selections <- bind_rows(selections)


##### Step 7: Summarize Predictive Performance and Ranch-Held-Out VIPs
# Score one complete set of out-of-fold predictions per repeat, not individual
# small test folds. n is unique samples; repeated predictions do not increase n.
stopifnot(all(prediction_results %>% count(Scheme, Repeat, Model, Sample) %>% pull(n) == 1L))
performance_by_repeat <- prediction_results %>%
  group_by(Scheme, Repeat, Model) %>%
  summarise(n = n(),
            R2_CV = 1 - sum((Observed_Log_MAOC - Predicted_Log_MAOC)^2) /
              sum((Observed_Log_MAOC - mean(Observed_Log_MAOC))^2),
            RMSE = sqrt(mean((Observed_Log_MAOC - Predicted_Log_MAOC)^2)),
            MAE = mean(abs(Observed_Log_MAOC - Predicted_Log_MAOC)), .groups = "drop")
stopifnot(all(performance_by_repeat$n == n))
performance <- performance_by_repeat %>%
  group_by(Scheme, Model) %>%
  summarise(n = first(n), Repeats = n(), SD_R2 = sd(R2_CV),
            Min_R2 = min(R2_CV), Max_R2 = max(R2_CV), R2_CV = mean(R2_CV),
            RMSE = sqrt(mean(RMSE^2)), MAE = mean(MAE), .groups = "drop") %>%
  relocate(SD_R2, Min_R2, Max_R2, .after = MAE) %>%
  group_by(Scheme) %>%
  mutate(Delta_R2_From_Baseline = if_else(Model %in% environmental_models,
                                          R2_CV - R2_CV[Model == "Environmental baseline"], NA_real_),
         Delta_RMSE_From_Baseline = if_else(Model %in% environmental_models,
                                            RMSE - RMSE[Model == "Environmental baseline"], NA_real_)) %>%
  ungroup() %>% arrange(match(Scheme, schemes), match(Model, models))

taxonomy <- map_dfr(markers, function(m) {
  ids <- colnames(X[[m]])
  tibble(Variable = paste0("ASV_", m, "_", ids), Marker = m, ASV = ids,
         Taxon = unname(raw[[m]]$taxonomy[ids]))
})
vip_summary <- bind_rows(full_vips) %>%
  group_by(Variable) %>%
  summarise(Folds_Available = n_distinct(Fold), Median_VIP = median(VIP),
            Lower_VIP = quantile(VIP, 0.25), Upper_VIP = quantile(VIP, 0.75),
            Proportion_VIP_Above_1 = mean(VIP > 1), .groups = "drop") %>%
  left_join(taxonomy, by = "Variable") %>%
  mutate(Shown_In_Figure = !str_starts(Variable, "ASV_") | Folds_Available >= 3,
         Category = case_when(
           Variable %in% c("pom_c", "c_n_ratio", "pH", "percent_clay") ~ "Physicochemical",
           Variable %in% c("mean_annual_temp", "mean_annual_precip") ~ "Climate",
           Variable == "Adaptive_Score" ~ "Management",
           TRUE ~ "Microbiome"))


##### Step 8: Recreate Figure 6 and Save Final Outputs
# Preserve the manuscript's colors, point style, taxonomic labels and VIP = 1
# reference line. Values come from the updated full-model ranch-held-out fits.
variable_labels <- c(
  pom_c = "POC", c_n_ratio = "C:N Ratio", pH = "pH", percent_clay = "%Clay",
  Adaptive_Score = "Adaptive Score", mean_annual_temp = "MAT", mean_annual_precip = "MAP",
  Richness_16s = "Bact/Arch Richness", Evenness_16s = "Bact/Arch Evenness",
  Shannons_16s = "Bact/Arch Diversity", Richness_ITS = "Fungal Richness",
  Evenness_ITS = "Fungal Evenness", Shannons_ITS = "Fungal Diversity",
  Biomass = "Microbial Biomass Proxy")
extract_rank <- function(taxon, rank) {
  str_match(taxon, paste0("(?:^|;)", rank, "_{1,2}([^;]+)"))[, 2]
}
plot_data <- vip_summary %>%
  filter(Shown_In_Figure) %>%
  mutate(Genus = extract_rank(Taxon, "g"), Species = extract_rank(Taxon, "s"),
         Family = extract_rank(Taxon, "f"), Order = extract_rank(Taxon, "o"),
         Class = extract_rank(Taxon, "c"), Phylum = extract_rank(Taxon, "p"),
         Taxon_label = case_when(
           !is.na(Genus) & !is.na(Species) ~ paste0("g_", Genus, ";s_", Species),
           !is.na(Genus) ~ paste0("g_", Genus), !is.na(Family) ~ paste0("f_", Family),
           !is.na(Order) ~ paste0("o_", Order), !is.na(Class) ~ paste0("c_", Class),
           !is.na(Phylum) ~ paste0("p_", Phylum), TRUE ~ "Unclassified")) %>%
  arrange(desc(Median_VIP)) %>%
  group_by(Marker, Taxon_label) %>%
  mutate(Taxon_count = n(), Taxon_number = row_number()) %>%
  ungroup() %>%
  mutate(Figure_Label = case_when(
    str_starts(Variable, "ASV_") & Taxon_count > 1 ~
      paste0(Marker, ": ", Taxon_label, " (ASV ", Taxon_number, ")"),
    str_starts(Variable, "ASV_") ~ paste0(Marker, ": ", Taxon_label),
    TRUE ~ recode(Variable, !!!variable_labels, .default = Variable)),
    Figure_Label = factor(Figure_Label, levels = rev(unique(Figure_Label))))

vip_plot <- ggplot(plot_data, aes(x = Median_VIP, y = Figure_Label, fill = Category)) +
  geom_point(shape = 21, size = 7, color = "black", stroke = 0.7) +
  geom_vline(xintercept = 1, linetype = "dotted", linewidth = 0.6, color = "black") +
  scale_fill_manual(name = "Category", values = c(
    Physicochemical = "#c2c04a", Microbiome = "#74a9b0",
    Climate = "#5b9b62", Management = "#d9962e")) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(x = "Median VIP score", y = NULL) +
  theme_minimal(base_family = "Helvetica") +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title.x = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray80"))

vip_summary <- vip_summary %>%
  left_join(plot_data %>% dplyr::select(Variable, Figure_Label), by = "Variable") %>%
  arrange(desc(Median_VIP))
write_csv(performance, file.path(output_dir, "PLS_performance.csv"))
write_csv(prediction_results, file.path(output_dir, "PLS_predictions.csv"))
write_csv(asv_selections, file.path(output_dir, "PLS_ASV_selections.csv"))
write_csv(vip_summary, file.path(output_dir, "PLS_VIP.csv"))
ggsave(file.path(figure_dir, "Figure_6.pdf"), vip_plot, width = 11, height = 8, bg = "white")
ggsave(file.path(figure_dir, "Figure_6.png"), vip_plot, width = 11, height = 8, dpi = 300, bg = "white")
print(performance)
print(vip_plot)

