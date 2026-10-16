# ============================================================
# Project: Peptide Activity Prediction
# Script:  Analysis Engine (VERSION 1 - modular refactor)
# Based on: peptide_activity_prediction_v7.R (RF + NN + PCA pipeline)
#
# PURPOSE:
#   This is Version 1 of the incremental development plan described in
#   PeptidePredictor_Incremental_Development_Roadmap.txt.
#
#   It converts the original monolithic script into a set of functions
#   culminating in a single entry point, run_analysis(), which can be
#   called from the command line, a test script, or (in Version 2) a
#   Shiny GUI. The scientific/modeling behavior is UNCHANGED from the
#   original script - this is a structural refactor only.
#
# HOW TO RUN (command-line / non-GUI):
#   source("peptide_analysis_engine_v1.R")
#   results <- run_analysis(
#     training_file    = "train_data_capped_300.xlsx",
#     new_file         = "new_peptides.xlsx",
#     target_class     = "AFP",
#     comparison_class = "FUNG",
#     n_runs           = 20,
#     prob_threshold   = 0.90,
#     auc_lower        = 0.60,
#     corr_threshold   = 0.90,
#     output_directory = "output"
#   )
#
# ============================================================

# ============================================================
# PACKAGE CHECKS
# ============================================================

required_packages <- c(
  "readxl", "openxlsx", "caret", "nnet", "ranger", "pROC",
  "recipes", "themis", "dplyr", "plotly", "htmlwidgets"
)

.check_required_packages <- function() {
  missing_packages <- required_packages[
    !sapply(required_packages, requireNamespace, quietly = TRUE)
  ]
  if (length(missing_packages) > 0) {
    stop(
      "ERROR: Missing required R packages: ", paste(missing_packages, collapse = ", "),
      "\nInstall them with: install.packages(c(",
      paste0('"', missing_packages, '"', collapse = ", "), "))"
    )
  }
}

.check_required_packages()

suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(caret)
  library(nnet)
  library(ranger)
  library(pROC)
  library(recipes)
  library(themis)
  library(dplyr)
  library(plotly)
  library(htmlwidgets)
})


# ============================================================
# STAGE 1: LOAD DATA
# ============================================================

#' Force every descriptor column to numeric.
#'
#' Excel exports from descriptor software often store numbers as text
#' (e.g. "935.17") or put error messages in cells. read_excel() guesses a
#' column's type from the first rows only, so such columns would otherwise
#' be read as text and silently dropped later. Here, text is converted to
#' numbers; anything that is not a finite number (error messages, "Inf",
#' blanks) becomes NA and is median-imputed later in the recipes.
#'
#' @param df        data.frame with sanitized column names
#' @param id_cols   columns to leave untouched (ID, Type)
#' @param label     name used in log messages
coerce_descriptors_numeric <- function(df, id_cols, label) {

  desc_cols <- setdiff(colnames(df), id_cols)
  n_converted <- 0L
  n_bad_cells <- 0L

  for (col in desc_cols) {
    x <- df[[col]]
    if (is.numeric(x)) {
      x[!is.finite(x)] <- NA_real_
      df[[col]] <- x
      next
    }
    x_chr <- trimws(as.character(x))
    x_chr[x_chr == ""] <- NA_character_
    x_num <- suppressWarnings(as.numeric(x_chr))
    x_num[!is.finite(x_num)] <- NA_real_
    n_bad_cells <- n_bad_cells + sum(is.na(x_num) & !is.na(x_chr))
    n_converted <- n_converted + 1L
    df[[col]] <- x_num
  }

  cat(label, ": ", n_converted, " descriptor column(s) were stored as text and converted to numeric; ",
      n_bad_cells, " non-numeric cell(s) (e.g. error messages) set to NA.\n", sep = "")
  df
}


#' Load training and new-peptide Excel files and do basic column validation.
#'
#' @param training_file Path to Excel file with columns ID | Type | descriptors
#' @param new_file      Path to Excel file with columns ID | descriptors
#' @return list(train_data, new_data) as data.frames with sanitized column names
load_data <- function(training_file, new_file) {

  if (!file.exists(training_file)) {
    stop("ERROR: Training file not found: ", training_file)
  }
  if (!file.exists(new_file)) {
    stop("ERROR: New-peptide file not found: ", new_file)
  }

  train_data <- tryCatch(
    as.data.frame(readxl::read_excel(training_file, guess_max = Inf)),
    error = function(e) stop("ERROR: Could not read training file '", training_file, "': ", e$message)
  )
  new_data <- tryCatch(
    as.data.frame(readxl::read_excel(new_file, guess_max = Inf)),
    error = function(e) stop("ERROR: Could not read new-peptide file '", new_file, "': ", e$message)
  )

  colnames(train_data) <- make.names(colnames(train_data))
  colnames(new_data)   <- make.names(colnames(new_data))

  # Force descriptor columns to numeric (handles numbers stored as text and
  # error messages inside cells) so no usable descriptor is silently dropped.
  train_id_cols <- intersect(c("ID", "Type"), colnames(train_data))
  new_id_cols   <- intersect("ID", colnames(new_data))
  train_data <- coerce_descriptors_numeric(train_data, train_id_cols, "Training file")
  new_data   <- coerce_descriptors_numeric(new_data,   new_id_cols,   "New-peptide file")

  # Remove descriptors that are entirely empty in the training data
  # (nothing to learn from them) and report it.
  all_na <- vapply(setdiff(colnames(train_data), train_id_cols),
                   function(cn) all(is.na(train_data[[cn]])), logical(1))
  if (any(all_na)) {
    cat("Removed", sum(all_na), "training descriptor(s) that are empty in every row.\n")
    train_data <- train_data[, !(colnames(train_data) %in% names(all_na)[all_na]), drop = FALSE]
  }

  # Tidy labels so stray spaces do not break the class match
  if ("Type" %in% colnames(train_data)) {
    train_data$Type <- trimws(as.character(train_data$Type))
  }

  if (!"Type" %in% colnames(train_data)) {
    stop("ERROR: Training file must contain a 'Type' column.")
  }
  if (!"ID" %in% colnames(train_data)) {
    stop("ERROR: Training file must contain an 'ID' column.")
  }
  if (!"ID" %in% colnames(new_data)) {
    stop("ERROR: New-peptide file must contain an 'ID' column.")
  }
  if (nrow(train_data) == 0) {
    stop("ERROR: Training file contains no rows.")
  }
  if (nrow(new_data) == 0) {
    stop("ERROR: New-peptide file contains no rows.")
  }

  list(train_data = train_data, new_data = new_data)
}


# ============================================================
# STAGE 2: PREPARE DATA (clean, binarize labels, identify
#          descriptors, remove duplicates)
# ============================================================

#' Clean and prepare the raw train/new data for modeling.
#' Mirrors the original script's sections 1-3 exactly.
#'
#' @return list with train_df, new_df, train_data_unique, descriptor_cols,
#'         target_class_r, imbalance_ratio
prepare_data <- function(train_data, new_data, target_class, comparison_class) {

  cat("\n--- DATA LOADING AND CLEANING ---\n")

  if (!is.null(comparison_class)) {
    cat("Filtering to", target_class, "vs", comparison_class, "only\n")
    train_data <- train_data[train_data$Type %in% c(target_class, comparison_class), ]
    cat("Samples after filtering:", nrow(train_data), "\n")
  }

  if (nrow(train_data) == 0) {
    stop("ERROR: No training samples remain after filtering to target_class/comparison_class. ",
         "Check that these class names exist in the 'Type' column.")
  }
  if (sum(train_data$Type == target_class) == 0) {
    stop("ERROR: target_class '", target_class, "' has zero samples in the training data.")
  }

  # Binary labels
  train_data$Type <- ifelse(train_data$Type == target_class, target_class, "other")
  train_data$Type <- factor(train_data$Type, levels = c(target_class, "other"))
  levels(train_data$Type) <- make.names(levels(train_data$Type))
  target_class_r <- make.names(target_class)

  # Empty strings to NA
  train_data[train_data == ""] <- NA
  new_data[new_data == ""]     <- NA

  cat("Class distribution:\n")
  print(table(train_data$Type))

  n_target <- sum(train_data$Type == target_class_r)
  n_other  <- sum(train_data$Type == "other")

  if (n_target == 0 || n_other == 0) {
    stop("ERROR: One of the two classes has zero samples after binarization. ",
         "Check target_class/comparison_class.")
  }

  imbalance_ratio <- max(n_target, n_other) / min(n_target, n_other)
  cat("\nImbalance ratio:", round(imbalance_ratio, 2), "\n")
  if (imbalance_ratio > 3) {
    cat("*** Imbalance > 3 detected - SMOTE will be applied ***\n")
  }

  # Descriptor columns
  descriptor_cols <- setdiff(colnames(train_data), c("ID", "Type"))
  descriptor_cols <- descriptor_cols[
    sapply(train_data[, descriptor_cols, drop = FALSE], is.numeric)
  ]
  cat("\nTotal numeric descriptors:", length(descriptor_cols), "\n")

  if (length(descriptor_cols) == 0) {
    stop("ERROR: No numeric descriptor columns found in training data.")
  }

  missing_new_cols <- setdiff(descriptor_cols, colnames(new_data))
  if (length(missing_new_cols) > 0) {
    stop("ERROR: New-peptide file is missing descriptor column(s) present in training data: ",
         paste(missing_new_cols, collapse = ", "))
  }

  # Remove duplicate samples (within each class separately)
  cat("\n=== Removing Duplicate Samples ===\n")

  train_data_unique <- do.call(rbind, lapply(
    levels(train_data$Type),
    function(cls) {
      rows <- train_data[train_data$Type == cls, ]
      rows[!duplicated(rows[, descriptor_cols]), ]
    }
  ))
  rownames(train_data_unique) <- NULL

  cat("Removed", nrow(train_data) - nrow(train_data_unique), "duplicate rows\n")
  cat("Remaining samples:", nrow(train_data_unique), "\n")
  cat("Class counts after deduplication:\n")
  print(table(train_data_unique$Type))

  train_df <- train_data_unique[, c("Type", descriptor_cols)]
  new_df   <- new_data[, descriptor_cols, drop = FALSE]

  if (nrow(train_df) < 50) {
    stop("ERROR: Too few samples (", nrow(train_df), "). Need at least 50.")
  }

  list(
    train_df          = train_df,
    new_df            = new_df,
    train_data_unique = train_data_unique,
    descriptor_cols   = descriptor_cols,
    target_class_r    = target_class_r,
    imbalance_ratio   = imbalance_ratio
  )
}


# ============================================================
# STAGE 3: TRAIN/TEST SPLIT (before feature selection)
# ============================================================

#' Split train_df into train_inner / test_outer BEFORE feature selection,
#' to avoid feature-selection leakage (matches original script's FIX 1).
split_data <- function(train_df, seed = 123) {

  cat("\n--- Train/Test Split (BEFORE feature selection) ---\n")

  set.seed(seed)
  idx         <- caret::createDataPartition(train_df$Type, p = 0.80, list = FALSE)
  train_inner <- train_df[idx, ]
  test_outer  <- train_df[-idx, ]

  cat("Train:", nrow(train_inner), "| Test:", nrow(test_outer), "\n")
  cat("Train distribution:\n"); print(table(train_inner$Type))
  cat("Test distribution:\n");  print(table(test_outer$Type))

  list(idx = idx, train_inner = train_inner, test_outer = test_outer)
}


# ============================================================
# STAGE 4: FEATURE SELECTION (on train_inner only)
# ============================================================

#' Select descriptors by per-feature AUC on train_inner only, with the
#' original script's fallback-threshold-relax logic (FIX 2).
select_features <- function(train_inner, descriptor_cols, auc_lower) {

  cat("\n--- Feature Selection (train_inner only) ---\n")

  feature_auc <- sapply(descriptor_cols, function(col) {
    tryCatch(
      as.numeric(pROC::auc(pROC::roc(train_inner$Type, train_inner[[col]], quiet = TRUE))),
      error = function(e) NA
    )
  })
  feature_auc <- feature_auc[!is.na(feature_auc)]

  cat("\nAUC distribution:\n")
  cat("  > 0.90  :", sum(feature_auc > 0.90), "\n")
  cat("  0.80-0.90:", sum(feature_auc > 0.80 & feature_auc <= 0.90), "\n")
  cat("  0.70-0.80:", sum(feature_auc > 0.70 & feature_auc <= 0.80), "\n")
  cat("  0.60-0.70:", sum(feature_auc > 0.60 & feature_auc <= 0.70), "\n")
  cat("  < 0.60  :", sum(feature_auc <= 0.60), "\n")

  selected_features <- names(feature_auc[feature_auc >= auc_lower])
  cat("\nFeatures with AUC >=", auc_lower, ":", length(selected_features), "\n")

  if (length(selected_features) < 10) {
    cat("*** Too few features, relaxing threshold to 0.55 ***\n")
    selected_features <- names(feature_auc[feature_auc >= 0.55])
    cat("After relaxing:", length(selected_features), "features\n")
  }

  if (length(selected_features) < 5) {
    stop("ERROR: Only ", length(selected_features), " features passed selection. Cannot proceed.")
  }

  cat("\nTop 10 features by AUC:\n")
  print(round(head(sort(feature_auc[selected_features], decreasing = TRUE), 10), 4))

  list(selected_features = selected_features, feature_auc = feature_auc)
}


# ============================================================
# STAGE 5: DIAGNOSTICS + LEAKAGE CHECK
# ============================================================

#' Apply selected features to all frames, check for train/test leakage,
#' and build the diagnostics table (matches original script's section
#' between feature selection and recipe creation).
apply_feature_selection <- function(train_inner, test_outer, new_df, train_df,
                                     descriptor_cols, target_class_r, imbalance_ratio) {

  train_inner <- train_inner[, c("Type", descriptor_cols)]
  test_outer  <- test_outer[,  c("Type", descriptor_cols)]
  new_df      <- new_df[,      descriptor_cols, drop = FALSE]
  train_df    <- train_df[,    c("Type", descriptor_cols)]

  # Leakage check
  combined  <- rbind(train_inner[, descriptor_cols], test_outer[, descriptor_cols])
  dup_idx   <- duplicated(combined) | duplicated(combined, fromLast = TRUE)
  test_dups <- dup_idx[-(1:nrow(train_inner))]
  cat("\nDuplicate samples between train/test:", sum(test_dups), "\n")
  if (sum(test_dups) > 0) cat("*** WARNING: leakage detected ***\n")

  test_diagnostics <- data.frame(
    Metric = c(
      "Total_Train_Size", "Total_Test_Size",
      paste0("Train_", target_class_r), "Train_other",
      paste0("Test_",  target_class_r), "Test_other",
      paste0("Test_",  target_class_r, "_Proportion"), "Test_other_Proportion",
      "Features_After_AUC_Filter", "Imbalance_Ratio"
    ),
    Value = c(
      nrow(train_inner), nrow(test_outer),
      sum(train_inner$Type == target_class_r),
      sum(train_inner$Type == "other"),
      sum(test_outer$Type  == target_class_r),
      sum(test_outer$Type  == "other"),
      round(sum(test_outer$Type == target_class_r) / nrow(test_outer), 3),
      round(sum(test_outer$Type == "other")         / nrow(test_outer), 3),
      length(descriptor_cols),
      round(imbalance_ratio, 2)
    )
  )
  cat("\n--- DATASET SUMMARY ---\n")
  print(test_diagnostics)

  list(
    train_inner = train_inner, test_outer = test_outer,
    new_df = new_df, train_df = train_df,
    test_diagnostics = test_diagnostics
  )
}


# ============================================================
# STAGE 6: RECIPES (correlation filter + SMOTE if imbalanced)
# ============================================================

#' Build the RF and NN+PCA recipes, with SMOTE added conditionally
#' (matches original script's FIX 3 / FIX 4).
create_recipes <- function(train_inner, corr_threshold, sampling_method) {

  add_sampling <- function(recipe_obj) {
    if (!is.null(sampling_method)) {
      recipe_obj <- recipe_obj %>% themis::step_smote(Type, over_ratio = 1)
    }
    recipe_obj
  }

  rec_rf <- recipe(Type ~ ., data = train_inner) %>%
    step_zv(all_predictors()) %>%
    step_nzv(all_predictors()) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_corr(all_numeric_predictors(), threshold = corr_threshold) %>%
    add_sampling()

  rec_nn_pca <- recipe(Type ~ ., data = train_inner) %>%
    step_zv(all_predictors()) %>%
    step_nzv(all_predictors()) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_corr(all_numeric_predictors(), threshold = corr_threshold) %>%
    step_normalize(all_numeric_predictors()) %>%
    step_pca(all_predictors(), threshold = 0.90) %>%
    add_sampling()

  list(rec_rf = rec_rf, rec_nn_pca = rec_nn_pca, add_sampling = add_sampling)
}


# ============================================================
# STAGE 7: 3D PCA VISUALIZATION
# ============================================================

#' Generate the two 3D PCA plotly widgets and the PCA coordinates table.
#' Saves HTML files into output_directory. Mirrors the original script
#' exactly (same recipe, same num_comp = 3, same plot structure).
generate_pca <- function(train_inner, new_df, train_data_unique, idx, new_data,
                          corr_threshold, target_class_r, output_directory) {

  cat("\n--- GENERATING 3D PCA PLOT ---\n")

  rec_pca_viz <- recipe(Type ~ ., data = train_inner) %>%
    step_zv(all_predictors()) %>%
    step_nzv(all_predictors()) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_corr(all_numeric_predictors(), threshold = corr_threshold) %>%
    step_normalize(all_numeric_predictors()) %>%
    step_pca(all_predictors(), num_comp = 3)

  pca_prepped    <- prep(rec_pca_viz, training = train_inner)
  pca_data_train <- bake(pca_prepped, new_data = train_inner)
  pca_data_new   <- bake(pca_prepped, new_data = new_df)

  n_pcs <- sum(grepl("^PC[0-9]+$", colnames(pca_data_train)))
  cat("Predictors entering PCA recipe:", ncol(train_inner) - 1,
      "| principal components computed:", n_pcs, "\n")
  if (n_pcs < 3) {
    stop("ERROR: Only ", n_pcs, " principal component(s) could be computed for the 3D plot. ",
         "Too few descriptors survive filtering (", ncol(train_inner) - 1,
         " entered the PCA recipe; correlation threshold = ", corr_threshold, "). ",
         "Check that descriptor columns are numeric in the Excel files, raise the ",
         "correlation threshold, or lower the minimum feature AUC.")
  }

  pca_step_idx       <- which(sapply(pca_prepped$steps, function(s) inherits(s, "step_pca")))
  pca_obj            <- pca_prepped$steps[[pca_step_idx]]$res
  variance_explained <- pca_obj$sdev[1:3]^2 / sum(pca_obj$sdev^2) * 100

  cat("Variance explained:\n")
  cat("  PC1:", round(variance_explained[1], 2), "%\n")
  cat("  PC2:", round(variance_explained[2], 2), "%\n")
  cat("  PC3:", round(variance_explained[3], 2), "%\n")
  cat("  Total:", round(sum(variance_explained), 2), "%\n")

  make_axis <- function(label, pct) list(title = paste0(label, " (", round(pct, 1), "%)"))

  fig_train <- plot_ly(
    data = pca_data_train,
    x = ~PC1, y = ~PC2, z = ~PC3,
    color = ~Type,
    colors = c("#1f77b4", "#ff7f0e"),
    type = "scatter3d", mode = "markers",
    marker = list(size = 5, opacity = 0.7),
    text = ~paste("Class:", Type), hoverinfo = "text"
  ) %>% layout(
    title = list(
      text = paste0("3D PCA - Training Data  |  PC1: ",
                    round(variance_explained[1], 1), "%, PC2: ",
                    round(variance_explained[2], 1), "%, PC3: ",
                    round(variance_explained[3], 1), "%"),
      font = list(size = 13)
    ),
    scene = list(
      xaxis = make_axis("PC1", variance_explained[1]),
      yaxis = make_axis("PC2", variance_explained[2]),
      zaxis = make_axis("PC3", variance_explained[3])
    )
  )

  pca_combined <- rbind(
    cbind(pca_data_train[, c("PC1", "PC2", "PC3")],
          Type = as.character(pca_data_train$Type), Dataset = "Training"),
    cbind(pca_data_new[, c("PC1", "PC2", "PC3")],
          Type = "New Sample", Dataset = "New")
  )

  fig_combined <- plot_ly() %>%
    add_trace(
      data = pca_combined[pca_combined$Dataset == "Training" &
                            pca_combined$Type == target_class_r, ],
      x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d", mode = "markers",
      name = paste("Training:", target_class_r),
      marker = list(size = 5, opacity = 0.6, color = "#1f77b4"),
      text = ~paste("Training -", Type), hoverinfo = "text"
    ) %>%
    add_trace(
      data = pca_combined[pca_combined$Dataset == "Training" &
                            pca_combined$Type != target_class_r, ],
      x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d", mode = "markers",
      name = "Training: other",
      marker = list(size = 5, opacity = 0.6, color = "#ff7f0e"),
      text = ~paste("Training - other"), hoverinfo = "text"
    ) %>%
    add_trace(
      data = pca_combined[pca_combined$Dataset == "New", ],
      x = ~PC1, y = ~PC2, z = ~PC3, type = "scatter3d", mode = "markers",
      name = "New Samples",
      marker = list(size = 7, opacity = 0.9, color = "#2ca02c", symbol = "diamond"),
      text = ~paste("New Sample"), hoverinfo = "text"
    ) %>%
    layout(
      title = list(
        text = paste0("3D PCA - Training vs New Samples  |  PC1: ",
                      round(variance_explained[1], 1), "%, PC2: ",
                      round(variance_explained[2], 1), "%, PC3: ",
                      round(variance_explained[3], 1), "%"),
        font = list(size = 13)
      ),
      scene = list(
        xaxis = make_axis("PC1", variance_explained[1]),
        yaxis = make_axis("PC2", variance_explained[2]),
        zaxis = make_axis("PC3", variance_explained[3])
      )
    )

  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  html_train_path    <- file.path(output_directory, paste0(target_class_r, "_PCA_3D_Training.html"))
  html_combined_path <- file.path(output_directory, paste0(target_class_r, "_PCA_3D_Combined.html"))

  saveWidget(fig_train,    html_train_path,    selfcontained = TRUE)
  saveWidget(fig_combined, html_combined_path, selfcontained = TRUE)
  cat("3D PCA plots saved.\n")

  pca_coordinates <- data.frame(
    ID      = c(train_data_unique$ID[idx], new_data$ID),
    PC1     = c(pca_data_train$PC1, pca_data_new$PC1),
    PC2     = c(pca_data_train$PC2, pca_data_new$PC2),
    PC3     = c(pca_data_train$PC3, pca_data_new$PC3),
    Type    = c(as.character(pca_data_train$Type), rep("New Sample", nrow(pca_data_new))),
    Dataset = c(rep("Training", nrow(pca_data_train)), rep("New", nrow(pca_data_new)))
  )

  cat("3D PCA complete.\n")

  list(
    pca_coordinates    = pca_coordinates,
    variance_explained = variance_explained,
    html_train_path    = html_train_path,
    html_combined_path = html_combined_path
  )
}


# ============================================================
# STAGE 8: TRAINING HELPERS (CV control, RF grid, calibration)
# ============================================================

#' Shared cross-validation control used by RF, NN+PCA, and NN+ANOVA.
.make_cv_ctrl <- function() {
  trainControl(
    method          = "repeatedcv",
    number          = 5,
    repeats         = 3,
    classProbs      = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final",
    allowParallel   = FALSE
  )
}

#' Dynamic RF tuning grid (FIX 5 - avoids ranger mtry errors on small
#' post-recipe column counts).
make_rf_grid <- function(n_cols) {
  mtry_vals <- unique(c(
    max(1, floor(sqrt(n_cols))),
    max(1, floor(sqrt(n_cols) * 1.5)),
    max(1, floor(n_cols / 3)),
    max(1, floor(n_cols / 5))
  ))
  mtry_vals <- mtry_vals[mtry_vals < n_cols]
  if (length(mtry_vals) == 0) mtry_vals <- 1
  expand.grid(
    mtry          = mtry_vals,
    splitrule     = c("gini", "extratrees"),
    min.node.size = c(3, 5, 10)
  )
}

#' Platt scaling probability calibration (FIX 7).
calibrate_probs <- function(raw_probs, true_labels, new_raw_probs, positive_class) {
  binary_labels <- as.integer(true_labels == positive_class)
  cal_df        <- data.frame(prob = raw_probs, label = binary_labels)
  platt_model   <- glm(label ~ prob, data = cal_df, family = binomial)
  new_df_cal    <- data.frame(prob = new_raw_probs)
  as.numeric(predict(platt_model, newdata = new_df_cal, type = "response"))
}


# ============================================================
# STAGE 9: SINGLE-RUN MODEL TRAINERS
# ============================================================
# Each of these trains one model for one iteration of the stability
# loop and returns calibrated probabilities + test AUC. Splitting
# these out of the loop makes each model's logic independently
# testable and keeps train_models() readable.

#' Train one Random Forest run.
.train_one_run_rf <- function(run_i, fit_set, cal_set, test_outer, new_df,
                               rec_rf, cv_ctrl, target_class_r) {

  cat("  [RF] Training...\n")

  prepped_rf <- prep(rec_rf, training = fit_set)
  baked_rf   <- bake(prepped_rf, new_data = NULL)
  n_cols_rf  <- ncol(baked_rf) - 1
  rf_grid    <- make_rf_grid(n_cols_rf)
  cat("  [RF] Post-corr columns:", n_cols_rf, "| Grid:", nrow(rf_grid), "combos\n")

  rf_model <- tryCatch(
    train(
      rec_rf, data = fit_set, method = "ranger",
      trControl = cv_ctrl, tuneGrid = rf_grid,
      metric = "ROC", importance = "permutation"
    ),
    error = function(e) stop("ERROR: RF training failed on run ", run_i, ": ", e$message)
  )

  rf_cal_raw <- predict(rf_model, cal_set,    type = "prob")[, target_class_r]
  rf_new_raw <- predict(rf_model, new_df,     type = "prob")[, target_class_r]
  rf_tst_raw <- predict(rf_model, test_outer, type = "prob")[, target_class_r]

  new_prob <- calibrate_probs(rf_cal_raw, cal_set$Type, rf_new_raw, target_class_r)
  tst_cal  <- calibrate_probs(rf_cal_raw, cal_set$Type, rf_tst_raw, target_class_r)
  auc_val  <- as.numeric(auc(roc(test_outer$Type, tst_cal, quiet = TRUE)))

  cat("  [RF] Test AUC:", round(auc_val, 4), "\n")
  cat("  [RF] Confusion matrix:\n")
  print(table(Actual = test_outer$Type, Predicted = tst_cal > 0.5))

  list(new_prob = new_prob, auc = auc_val, model = rf_model)
}

#' Train one NN + PCA run.
.train_one_run_nn_pca <- function(run_i, fit_set, cal_set, test_outer, new_df,
                                   rec_nn_pca, cv_ctrl, target_class_r) {

  cat("\n  [NN-PCA] Training...\n")

  nn_pca_model <- tryCatch(
    train(
      rec_nn_pca, data = fit_set, method = "nnet",
      trControl = cv_ctrl,
      tuneGrid  = expand.grid(size = 2:5, decay = c(0.05, 0.1, 0.3)),
      metric    = "ROC", trace = FALSE, maxit = 1000
    ),
    error = function(e) stop("ERROR: NN+PCA training failed on run ", run_i, ": ", e$message)
  )

  cal_raw <- predict(nn_pca_model, cal_set,    type = "prob")[, target_class_r]
  new_raw <- predict(nn_pca_model, new_df,     type = "prob")[, target_class_r]
  tst_raw <- predict(nn_pca_model, test_outer, type = "prob")[, target_class_r]

  new_prob <- calibrate_probs(cal_raw, cal_set$Type, new_raw, target_class_r)
  tst_cal  <- calibrate_probs(cal_raw, cal_set$Type, tst_raw, target_class_r)
  auc_val  <- as.numeric(auc(roc(test_outer$Type, tst_cal, quiet = TRUE)))

  cat("  [NN-PCA] Test AUC:", round(auc_val, 4), "\n")
  cat("  [NN-PCA] Confusion matrix:\n")
  print(table(Actual = test_outer$Type, Predicted = tst_cal > 0.5))

  list(new_prob = new_prob, auc = auc_val)
}

#' Train one NN + ANOVA-selected-descriptors run.
.train_one_run_nn_anova <- function(run_i, fit_set, cal_set, test_outer, new_df,
                                     descriptor_cols, corr_threshold, add_sampling,
                                     cv_ctrl, target_class_r) {

  cat("\n  [NN-ANOVA] Training...\n")

  anova_pvals <- sapply(
    fit_set[, descriptor_cols, drop = FALSE],
    function(x) tryCatch(
      summary(aov(x ~ fit_set$Type))[[1]][["Pr(>F)"]][1],
      error = function(e) NA
    )
  )

  top_n     <- min(30, sum(!is.na(anova_pvals)))
  top_feats <- names(sort(anova_pvals))[1:top_n]

  train_red <- fit_set[,     c("Type", top_feats)]
  cal_red   <- cal_set[,     c("Type", top_feats)]
  test_red  <- test_outer[,  top_feats, drop = FALSE]
  new_red   <- new_df[,      top_feats, drop = FALSE]

  rec_nn_desc <- recipe(Type ~ ., data = train_red) %>%
    step_zv(all_predictors()) %>%
    step_impute_median(all_predictors()) %>%
    step_corr(all_numeric_predictors(), threshold = corr_threshold) %>%
    step_normalize(all_predictors()) %>%
    add_sampling()

  nn_desc_model <- tryCatch(
    train(
      rec_nn_desc, data = train_red, method = "nnet",
      trControl = cv_ctrl,
      tuneGrid  = expand.grid(size = 2:5, decay = c(0.05, 0.1, 0.3)),
      metric    = "ROC", trace = FALSE, maxit = 1000
    ),
    error = function(e) stop("ERROR: NN+ANOVA training failed on run ", run_i, ": ", e$message)
  )

  cal_raw <- predict(nn_desc_model, cal_red,  type = "prob")[, target_class_r]
  new_raw <- predict(nn_desc_model, new_red,  type = "prob")[, target_class_r]
  tst_raw <- predict(nn_desc_model, test_red, type = "prob")[, target_class_r]

  new_prob <- calibrate_probs(cal_raw, cal_set$Type, new_raw, target_class_r)
  tst_cal  <- calibrate_probs(cal_raw, cal_set$Type, tst_raw, target_class_r)
  auc_val  <- as.numeric(auc(roc(test_outer$Type, tst_cal, quiet = TRUE)))

  cat("  [NN-ANOVA] Test AUC:", round(auc_val, 4), "\n")
  cat("  [NN-ANOVA] Confusion matrix:\n")
  print(table(Actual = test_outer$Type, Predicted = tst_cal > 0.5))

  list(new_prob = new_prob, auc = auc_val)
}


# ============================================================
# STAGE 10: MULTI-RUN TRAINING LOOP + STABILITY SUMMARY
# ============================================================

#' Run the n_runs-iteration stability loop across RF, NN+PCA, and
#' NN+ANOVA, then build the stability_summary table. NOTE: as in the
#' original script, the RF feature-importance model returned here is
#' from the LAST run only (rf_model_last) - this is a known
#' characteristic of the original methodology, preserved deliberately.
train_models <- function(train_inner, test_outer, new_df, new_data,
                          descriptor_cols, recipes, sampling_method,
                          corr_threshold, n_runs, prob_threshold, target_class_r,
                          progress_callback = NULL) {
  # progress_callback, if supplied, is a function(value, detail) called
  # after each of the n_runs completes, with value in [0, 1] and detail
  # a short status string. Optional and NULL by default so command-line
  # use (Version 1) is completely unaffected - this keeps the engine
  # decoupled from any specific GUI (it never references Shiny directly).

  cat("\n--- Multi-Run Training (", n_runs, "runs) ---\n")

  cv_ctrl <- .make_cv_ctrl()

  rf_probs      <- matrix(NA, nrow(new_data), n_runs)
  nn_pca_probs  <- matrix(NA, nrow(new_data), n_runs)
  nn_desc_probs <- matrix(NA, nrow(new_data), n_runs)

  rf_auc      <- numeric(n_runs)
  nn_pca_auc  <- numeric(n_runs)
  nn_desc_auc <- numeric(n_runs)

  rownames(rf_probs) <- rownames(nn_pca_probs) <-
    rownames(nn_desc_probs) <- new_data$ID

  rf_model_last <- NULL

  for (i in seq_len(n_runs)) {

    cat("\n===== RUN", i, "of", n_runs, "=====\n")
    set.seed(1000 + i)

    cal_idx <- createDataPartition(train_inner$Type, p = 0.20, list = FALSE)
    cal_set <- train_inner[ cal_idx, ]
    fit_set <- train_inner[-cal_idx, ]

    rf_result <- .train_one_run_rf(
      i, fit_set, cal_set, test_outer, new_df,
      recipes$rec_rf, cv_ctrl, target_class_r
    )
    rf_probs[, i] <- rf_result$new_prob
    rf_auc[i]     <- rf_result$auc
    rf_model_last <- rf_result$model

    nn_pca_result <- .train_one_run_nn_pca(
      i, fit_set, cal_set, test_outer, new_df,
      recipes$rec_nn_pca, cv_ctrl, target_class_r
    )
    nn_pca_probs[, i] <- nn_pca_result$new_prob
    nn_pca_auc[i]     <- nn_pca_result$auc

    nn_desc_result <- .train_one_run_nn_anova(
      i, fit_set, cal_set, test_outer, new_df,
      descriptor_cols, corr_threshold, recipes$add_sampling,
      cv_ctrl, target_class_r
    )
    nn_desc_probs[, i] <- nn_desc_result$new_prob
    nn_desc_auc[i]     <- nn_desc_result$auc

    if (!is.null(progress_callback)) {
      progress_callback(
        value  = i / n_runs,
        detail = paste0("Completed run ", i, " of ", n_runs)
      )
    }
  }

  cat("\n--- Stability Summary ---\n")

  stability_summary <- data.frame(
    ID = new_data$ID,

    RF_Mean = rowMeans(rf_probs),
    RF_SD   = apply(rf_probs, 1, sd),
    RF_Freq = rowMeans(rf_probs >= prob_threshold),

    NN_PCA_Mean = rowMeans(nn_pca_probs),
    NN_PCA_SD   = apply(nn_pca_probs, 1, sd),
    NN_PCA_Freq = rowMeans(nn_pca_probs >= prob_threshold),

    NN_Desc_Mean = rowMeans(nn_desc_probs),
    NN_Desc_SD   = apply(nn_desc_probs, 1, sd),
    NN_Desc_Freq = rowMeans(nn_desc_probs >= prob_threshold)
  )

  list(
    stability_summary = stability_summary,
    rf_auc = rf_auc, nn_pca_auc = nn_pca_auc, nn_desc_auc = nn_desc_auc,
    rf_model_last = rf_model_last
  )
}


# ============================================================
# STAGE 11: AUC-WEIGHTED ENSEMBLE + MODEL EVALUATION TABLE
# ============================================================

#' Combine per-model mean probabilities into an AUC-weighted ensemble
#' (FIX 6), and build the model evaluation summary table.
calculate_ensemble <- function(stability_summary, rf_auc, nn_pca_auc, nn_desc_auc) {

  cat("\n--- AUC-Weighted Ensemble ---\n")

  mean_rf_auc   <- mean(rf_auc)
  mean_nn_auc   <- mean(nn_pca_auc)
  mean_desc_auc <- mean(nn_desc_auc)
  raw_weights   <- c(mean_rf_auc, mean_nn_auc, mean_desc_auc)
  weights       <- raw_weights / sum(raw_weights)

  cat("Model weights:\n")
  cat("  RF       :", round(weights[1], 4), " (AUC =", round(mean_rf_auc,   4), ")\n")
  cat("  NN_PCA   :", round(weights[2], 4), " (AUC =", round(mean_nn_auc,   4), ")\n")
  cat("  NN_ANOVA :", round(weights[3], 4), " (AUC =", round(mean_desc_auc, 4), ")\n")

  ensemble_prob <- weights[1] * stability_summary$RF_Mean +
                   weights[2] * stability_summary$NN_PCA_Mean +
                   weights[3] * stability_summary$NN_Desc_Mean

  auc_df <- data.frame(
    Model           = c("RF", "NN_PCA", "NN_DESC_ANOVA", "Ensemble_Weight"),
    Mean_AUC        = c(mean_rf_auc, mean_nn_auc, mean_desc_auc, NA),
    SD_AUC          = c(sd(rf_auc),  sd(nn_pca_auc), sd(nn_desc_auc), NA),
    Min_AUC         = c(min(rf_auc), min(nn_pca_auc), min(nn_desc_auc), NA),
    Max_AUC         = c(max(rf_auc), max(nn_pca_auc), max(nn_desc_auc), NA),
    Ensemble_Weight = c(round(weights, 4), NA)
  )

  cat("\n--- MODEL PERFORMANCE ---\n")
  print(auc_df)

  list(
    ensemble_prob = ensemble_prob, weights = weights, auc_df = auc_df,
    mean_rf_auc = mean_rf_auc, mean_nn_auc = mean_nn_auc, mean_desc_auc = mean_desc_auc
  )
}


# ============================================================
# STAGE 12: CONSENSUS + FEATURE IMPORTANCE
# ============================================================

#' Compute consensus confidence tiers and RF feature importance
#' (from the last training run, matching original behavior).
evaluate_models <- function(stability_summary, rf_model_last, train_df, target_class_r) {

  consensus <- stability_summary %>%
    mutate(
      Consensus_Mean  = rowMeans(select(., RF_Freq, NN_PCA_Freq, NN_Desc_Freq)),
      Models_80       = rowSums(select(., RF_Freq, NN_PCA_Freq, NN_Desc_Freq) >= 0.8),
      Confidence_Tier = case_when(
        Consensus_Mean >= 0.9 & Models_80 == 3 ~ "Very High",
        Consensus_Mean >= 0.8 & Models_80 >= 2 ~ "High",
        TRUE ~ "Low"
      )
    ) %>%
    arrange(desc(Confidence_Tier), desc(Consensus_Mean))

  cat("\n--- FEATURE IMPORTANCE (last run RF) ---\n")

  importance_scores <- ranger::importance(rf_model_last$finalModel)
  importance_df     <- data.frame(
    Feature    = names(importance_scores),
    Importance = importance_scores
  ) %>% arrange(desc(Importance))

  cat("\nTop 20 features:\n")
  print(head(importance_df, 20))

  top5_percent <- sum(head(importance_df$Importance, 5)) / sum(importance_df$Importance) * 100
  cat("\nTop 5 features account for", round(top5_percent, 1), "% of total importance\n")
  if (top5_percent > 70) {
    cat("-> Classification relies heavily on a few features\n")
  } else if (top5_percent > 50) {
    cat("-> Moderate feature diversity\n")
  } else {
    cat("-> Information spread across many features\n")
  }

  top5_list <- head(importance_df$Feature, 5)

  feature_stats <- do.call(rbind, lapply(top5_list, function(feat) {
    if (feat %in% colnames(train_df)) {
      tv <- train_df[[feat]][train_df$Type == target_class_r]
      ov <- train_df[[feat]][train_df$Type == "other"]
      data.frame(
        Feature     = feat,
        Target_Mean = mean(tv, na.rm = TRUE),
        Other_Mean  = mean(ov, na.rm = TRUE),
        Difference  = mean(tv, na.rm = TRUE) - mean(ov, na.rm = TRUE),
        Importance  = importance_df$Importance[importance_df$Feature == feat]
      )
    }
  }))

  print(feature_stats)

  list(
    consensus = consensus, importance_df = importance_df,
    feature_stats = feature_stats, top5_percent = top5_percent
  )
}


# ============================================================
# STAGE 13: PREDICTIONS
# ============================================================

#' Build the clear/detailed, simple, and high-confidence prediction
#' tables from the ensemble probability and stability summary.
generate_predictions <- function(new_data, ensemble_prob, stability_summary, target_class_r) {

  cat("\n--- PREDICTIONS ---\n")

  clear_predictions <- data.frame(
    ID                   = new_data$ID,
    Ensemble_Probability = round(ensemble_prob, 4),
    Prediction_50        = ifelse(ensemble_prob > 0.50, target_class_r, "other"),
    Prediction_70        = ifelse(ensemble_prob > 0.70, target_class_r, "other"),
    Prediction_90        = ifelse(ensemble_prob > 0.90, target_class_r, "other"),
    Confidence           = case_when(
      ensemble_prob > 0.9 | ensemble_prob < 0.1 ~ "Very High",
      ensemble_prob > 0.7 | ensemble_prob < 0.3 ~ "High",
      ensemble_prob > 0.6 | ensemble_prob < 0.4 ~ "Medium",
      TRUE ~ "Low"
    ),
    Models_Agree_90 = rowSums(cbind(
      stability_summary$RF_Freq,
      stability_summary$NN_PCA_Freq,
      stability_summary$NN_Desc_Freq
    ) >= 0.9),
    RF_Prob      = round(stability_summary$RF_Mean,      4),
    NN_PCA_Prob  = round(stability_summary$NN_PCA_Mean,  4),
    NN_DESC_Prob = round(stability_summary$NN_Desc_Mean, 4),
    RF_Pred      = ifelse(stability_summary$RF_Mean      > 0.5, target_class_r, "other"),
    NN_PCA_Pred  = ifelse(stability_summary$NN_PCA_Mean  > 0.5, target_class_r, "other"),
    NN_DESC_Pred = ifelse(stability_summary$NN_Desc_Mean > 0.5, target_class_r, "other")
  )

  clear_predictions$Interpretation <- paste0(
    ifelse(clear_predictions$Prediction_50 == target_class_r,
           paste0("POSITIVE for ", target_class_r),
           paste0("NEGATIVE for ", target_class_r)),
    " (", clear_predictions$Confidence, " confidence)"
  )

  cat("50% threshold:", sum(clear_predictions$Prediction_50 == target_class_r),
      target_class_r, "| other:", sum(clear_predictions$Prediction_50 == "other"), "\n")
  cat("70% threshold:", sum(clear_predictions$Prediction_70 == target_class_r),
      target_class_r, "| other:", sum(clear_predictions$Prediction_70 == "other"), "\n")
  cat("90% threshold:", sum(clear_predictions$Prediction_90 == target_class_r),
      target_class_r, "| other:", sum(clear_predictions$Prediction_90 == "other"), "\n")
  cat("Confidence distribution:\n")
  print(table(clear_predictions$Confidence))
  cat("All 3 models agree at 90%:", sum(clear_predictions$Models_Agree_90 == 3), "\n")

  high_confidence_positives <- clear_predictions %>%
    filter(Prediction_70 == target_class_r, Confidence %in% c("High", "Very High")) %>%
    arrange(desc(Ensemble_Probability))

  cat("\n=== HIGH-CONFIDENCE", target_class_r, "CANDIDATES:",
      nrow(high_confidence_positives), "===\n")
  if (nrow(high_confidence_positives) > 0) {
    print(head(high_confidence_positives[,
      c("ID", "Ensemble_Probability", "Confidence", "Models_Agree_90")], 10))
  }

  simple_predictions <- data.frame(
    ID              = new_data$ID,
    Is_Target       = ifelse(ensemble_prob > 0.5, "YES", "NO"),
    Probability_pct = round(ensemble_prob * 100, 1),
    Confidence      = clear_predictions$Confidence
  )
  colnames(simple_predictions)[2] <- paste0("Is_", target_class_r)
  colnames(simple_predictions)[3] <- paste0(target_class_r, "_Probability_%")

  list(
    clear_predictions = clear_predictions,
    simple_predictions = simple_predictions,
    high_confidence_positives = high_confidence_positives
  )
}


# ============================================================
# STAGE 14: EXCEL OUTPUT
# ============================================================

#' Write the 10-sheet Excel workbook to output_directory.
#' Returns the full path to the workbook.
generate_outputs <- function(output_directory, target_class_r, comparison_class,
                              imbalance_ratio, sampling_method, corr_threshold,
                              mean_rf_auc, mean_nn_auc, mean_desc_auc, auc_df,
                              simple_predictions, clear_predictions,
                              high_confidence_positives, consensus,
                              importance_df, feature_stats, test_diagnostics,
                              pca_coordinates) {

  cat("\n--- WRITING EXCEL FILE ---\n")

  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

  wb <- createWorkbook()

  addWorksheet(wb, "1_Summary")
  addWorksheet(wb, "2_Simple_Predictions")
  addWorksheet(wb, "3_Detailed_Predictions")
  addWorksheet(wb, "4_High_Confidence")
  addWorksheet(wb, "5_Consensus")
  addWorksheet(wb, "6_Model_Evaluation")
  addWorksheet(wb, "7_Feature_Importance")
  addWorksheet(wb, "8_Top_Features_Stats")
  addWorksheet(wb, "9_Diagnostics")
  addWorksheet(wb, "10_PCA_Coordinates")

  writeData(wb, "1_Summary", data.frame(
    Section = c(
      "Dataset Info", "", "",
      "Fixes Applied (v2 methodology, v1 engine refactor)", "", "", "", "", "", "",
      "Model Performance", "", "", "",
      "Predictions Summary", "", "", ""
    ),
    Description = c(
      paste("Target Class      :", target_class_r),
      paste("Comparison        :", ifelse(is.null(comparison_class), "All others", comparison_class)),
      paste("Imbalance Ratio   :", round(imbalance_ratio, 2)),
      "FIX 1: Split before feature selection - no leakage",
      "FIX 2: No upper AUC ceiling - all informative features kept",
      paste("FIX 3: Correlation filter (threshold =", corr_threshold, ")"),
      paste("FIX 4: SMOTE -", ifelse(is.null(sampling_method), "skipped (ratio OK)", "applied")),
      "FIX 5: RF mtry grid built dynamically to avoid ranger errors",
      "FIX 6: AUC-weighted ensemble",
      "FIX 7: Platt scaling probability calibration",
      paste("RF Mean AUC       :", round(mean_rf_auc,   4)),
      paste("NN_PCA Mean AUC   :", round(mean_nn_auc,   4)),
      paste("NN_ANOVA Mean AUC :", round(mean_desc_auc, 4)),
      paste("Best Model AUC    :", round(max(auc_df$Mean_AUC, na.rm = TRUE), 4)),
      paste("Total Predictions :", nrow(simple_predictions)),
      paste("Predicted as", target_class_r, "(50%):", sum(simple_predictions[, 2] == "YES")),
      paste("High Conf. Cands  :", nrow(high_confidence_positives)),
      paste("Models Agree 90%  :", sum(clear_predictions$Models_Agree_90 == 3))
    )
  ))

  writeData(wb, "2_Simple_Predictions",   simple_predictions)
  writeData(wb, "3_Detailed_Predictions", clear_predictions)
  writeData(wb, "4_High_Confidence",      high_confidence_positives)
  writeData(wb, "5_Consensus",            consensus)
  writeData(wb, "6_Model_Evaluation",     auc_df)
  writeData(wb, "7_Feature_Importance",   importance_df)
  writeData(wb, "8_Top_Features_Stats",   feature_stats)
  writeData(wb, "9_Diagnostics",          test_diagnostics)
  writeData(wb, "10_PCA_Coordinates",     pca_coordinates)

  output_filename <- file.path(output_directory, paste0(target_class_r, "_Complete_Analysis_v1.xlsx"))
  saveWorkbook(wb, output_filename, overwrite = TRUE)

  cat("Excel workbook saved:", output_filename, "\n")

  output_filename
}


# ============================================================
# STAGE 15: TOP-LEVEL ORCHESTRATOR
# ============================================================

#' Run the complete peptide activity prediction analysis end-to-end.
#'
#' @param training_file    Path to training Excel file (ID | Type | descriptors)
#' @param new_file         Path to new-peptide Excel file (ID | descriptors)
#' @param target_class     Class to predict (e.g. "AFP")
#' @param comparison_class Class to compare against, or NULL for one-vs-all
#' @param n_runs           Number of stability runs (default 20)
#' @param prob_threshold   Threshold for frequency columns (default 0.90)
#' @param auc_lower        Minimum per-feature AUC to keep a descriptor (default 0.60)
#' @param corr_threshold   Correlation filter cutoff (default 0.90)
#' @param output_directory Directory to write Excel/HTML outputs into
#' @param progress_callback Optional function(value, detail) called at
#'        major pipeline stages and after each stability run, with
#'        value in [0, 1]. NULL by default (Version 1 command-line
#'        behavior, unaffected). This is how Version 3's Shiny app
#'        drives a real progress bar without any modeling logic living
#'        in the GUI layer - the engine stays GUI-agnostic since it only
#'        ever calls a generic function pointer, never Shiny itself.
#'
#' @return A named list containing every intermediate and final result
#'         object, plus paths to the generated files. Nothing is left
#'         in the global environment.
run_analysis <- function(training_file, new_file,
                          target_class      = "AFP",
                          comparison_class  = "FUNG",
                          n_runs            = 20,
                          prob_threshold    = 0.90,
                          auc_lower         = 0.60,
                          corr_threshold    = 0.90,
                          output_directory  = "output",
                          progress_callback = NULL) {

  .progress <- function(value, detail) {
    if (!is.null(progress_callback)) progress_callback(value, detail)
  }

  # ---- Parameter validation ----
  if (!is.numeric(n_runs) || n_runs < 1 || n_runs != round(n_runs)) {
    stop("ERROR: n_runs must be a positive integer.")
  }
  if (!is.numeric(prob_threshold) || prob_threshold <= 0 || prob_threshold >= 1) {
    stop("ERROR: prob_threshold must be between 0 and 1 (exclusive).")
  }
  if (!is.numeric(auc_lower) || auc_lower <= 0 || auc_lower >= 1) {
    stop("ERROR: auc_lower must be between 0 and 1 (exclusive).")
  }
  if (!is.numeric(corr_threshold) || corr_threshold <= 0 || corr_threshold >= 1) {
    stop("ERROR: corr_threshold must be between 0 and 1 (exclusive).")
  }
  if (!is.null(comparison_class) && !is.character(comparison_class)) {
    stop("ERROR: comparison_class must be a string or NULL.")
  }
  if (!is.character(target_class) || nchar(target_class) == 0) {
    stop("ERROR: target_class must be a non-empty string.")
  }

  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

  # ---- 1. Load ----
  .progress(0.02, "Loading Excel files...")
  loaded <- load_data(training_file, new_file)

  # ---- 2. Prepare ----
  .progress(0.05, "Cleaning data and identifying descriptors...")
  prepared <- prepare_data(loaded$train_data, loaded$new_data, target_class, comparison_class)

  # ---- 3. Split ----
  .progress(0.08, "Splitting train/test...")
  split <- split_data(prepared$train_df, seed = 123)

  # ---- 4. Feature selection ----
  .progress(0.12, "Selecting features by AUC...")
  feat_sel <- select_features(split$train_inner, prepared$descriptor_cols, auc_lower)

  # ---- 5. Apply selection + diagnostics ----
  .progress(0.15, "Checking for train/test leakage...")
  applied <- apply_feature_selection(
    split$train_inner, split$test_outer, prepared$new_df, prepared$train_df,
    feat_sel$selected_features, prepared$target_class_r, prepared$imbalance_ratio
  )

  sampling_method <- if (prepared$imbalance_ratio > 3) "smote" else NULL

  # ---- 6. Recipes ----
  .progress(0.17, "Building preprocessing recipes...")
  recipes <- create_recipes(applied$train_inner, corr_threshold, sampling_method)

  # ---- 7. PCA visualization ----
  .progress(0.20, "Generating 3D PCA visualization...")
  pca_result <- generate_pca(
    applied$train_inner, applied$new_df, prepared$train_data_unique, split$idx,
    loaded$new_data, corr_threshold, prepared$target_class_r, output_directory
  )

  # ---- 8-10. Multi-run training + stability summary ----
  .progress(0.25, paste0("Starting ", n_runs, " stability run(s)..."))
  training_result <- train_models(
    applied$train_inner, applied$test_outer, applied$new_df, loaded$new_data,
    feat_sel$selected_features, recipes, sampling_method,
    corr_threshold, n_runs, prob_threshold, prepared$target_class_r,
    progress_callback = function(value, detail) .progress(0.25 + 0.60 * value, detail)
  )

  # ---- 11. Ensemble ----
  .progress(0.88, "Computing AUC-weighted ensemble...")
  ensemble_result <- calculate_ensemble(
    training_result$stability_summary,
    training_result$rf_auc, training_result$nn_pca_auc, training_result$nn_desc_auc
  )

  # ---- 12. Consensus + feature importance ----
  .progress(0.92, "Computing consensus and feature importance...")
  evaluation <- evaluate_models(
    training_result$stability_summary, training_result$rf_model_last,
    applied$train_df, prepared$target_class_r
  )

  # ---- 13. Predictions ----
  .progress(0.95, "Building prediction tables...")
  predictions <- generate_predictions(
    loaded$new_data, ensemble_result$ensemble_prob,
    training_result$stability_summary, prepared$target_class_r
  )

  # ---- 14. Excel output ----
  .progress(0.98, "Writing Excel workbook...")
  excel_path <- generate_outputs(
    output_directory, prepared$target_class_r, comparison_class,
    prepared$imbalance_ratio, sampling_method, corr_threshold,
    ensemble_result$mean_rf_auc, ensemble_result$mean_nn_auc, ensemble_result$mean_desc_auc,
    ensemble_result$auc_df, predictions$simple_predictions, predictions$clear_predictions,
    predictions$high_confidence_positives, evaluation$consensus,
    evaluation$importance_df, evaluation$feature_stats, applied$test_diagnostics,
    pca_result$pca_coordinates
  )

  cat("\n--- ANALYSIS COMPLETE ---\n")
  cat("Excel  :", excel_path, "\n")
  cat("PCA    :", pca_result$html_train_path, "\n")
  cat("       ", pca_result$html_combined_path, "\n")

  .progress(1.0, "Analysis complete.")

  list(
    # dataset-level results
    target_class_r    = prepared$target_class_r,
    imbalance_ratio    = prepared$imbalance_ratio,
    descriptor_cols    = feat_sel$selected_features,
    test_diagnostics   = applied$test_diagnostics,

    # model performance
    auc_df             = ensemble_result$auc_df,
    stability_summary  = training_result$stability_summary,
    consensus          = evaluation$consensus,
    importance_df      = evaluation$importance_df,
    feature_stats      = evaluation$feature_stats,

    # predictions
    simple_predictions        = predictions$simple_predictions,
    clear_predictions         = predictions$clear_predictions,
    high_confidence_positives = predictions$high_confidence_positives,

    # PCA
    pca_coordinates    = pca_result$pca_coordinates,
    variance_explained = pca_result$variance_explained,

    # file paths
    excel_path          = excel_path,
    pca_html_train_path = pca_result$html_train_path,
    pca_html_combined_path = pca_result$html_combined_path
  )
}

# ================================ END ==========================================
