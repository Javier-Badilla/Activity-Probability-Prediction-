# ============================================================
# Project: Peptide Activity Prediction
# Script:  Random Forest + NN + PCA Pipeline  (v2 - Clean)
# Author:  Javier Badilla
# Date:    2025 - 2026
#
# WORKFLOW BEFORE RUNNING:
#   1. Cap your training data to max 300 per class using this snippet:
#   2. Load train_data_capped_300.xlsx as train_data
#   3. Load your new peptides as new_data
#   4. Run this script
#
# ---------------------- REQUIREMENTS ------------------------
#
# R version: 4.3+
# Packages:
#   readxl, openxlsx, caret, ranger, nnet, pROC,
#   dplyr, recipes, themis, plotly, htmlwidgets
#
# ---------------------- INPUT FORMAT ------------------------
#
#   train_data : ID | Type | descriptor1 | descriptor2 | ...
#   new_data   : ID | descriptor1 | descriptor2 | ...
#
# ============================================================

# ----------------------- LIBRARIES --------------------------

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


# ----------------------- USER INPUTS ------------------------

target_class     <- "AFP"       # class to predict FUNG CANCER AFP AChE 
comparison_class <- "FUNG"         # NULL = one-vs-all
                                 # e.g. "viral" = AChE vs viral only

n_runs           <- 20           # runs for stability (more = slower)
prob_threshold   <- 0.90         # threshold for frequency columns

auc_lower        <- 0.60         # min AUC to keep a descriptor
corr_threshold   <- 0.90         # correlation filter cutoff


# ============================================================
# 1. CLEANING
# ============================================================

cat("\n--- DATA LOADING AND CLEANING ---\n")

colnames(train_data) <- make.names(colnames(train_data))
colnames(new_data)   <- make.names(colnames(new_data))

# Optional strict binary filter
if (!is.null(comparison_class)) {
  cat("Filtering to", target_class, "vs", comparison_class, "only\n")
  train_data <- train_data[train_data$Type %in% c(target_class, comparison_class), ]
  cat("Samples after filtering:", nrow(train_data), "\n")
}

# Binary labels
train_data$Type <- ifelse(train_data$Type == target_class, target_class, "other")
train_data$Type <- factor(train_data$Type, levels = c(target_class, "other"))
levels(train_data$Type) <- make.names(levels(train_data$Type))
target_class_r  <- make.names(target_class)

# Empty strings to NA
train_data[train_data == ""] <- NA
new_data[new_data == ""]     <- NA

cat("Class distribution:\n")
print(table(train_data$Type))

# Imbalance check
n_target  <- sum(train_data$Type == target_class_r)
n_other   <- sum(train_data$Type == "other")
imbalance_ratio <- max(n_target, n_other) / min(n_target, n_other)
cat("\nImbalance ratio:", round(imbalance_ratio, 2), "\n")
if (imbalance_ratio > 3) {
  cat("*** Imbalance > 3 detected - SMOTE will be applied ***\n")
}


# ============================================================
# 2. DESCRIPTOR COLUMNS
# ============================================================

descriptor_cols <- setdiff(colnames(train_data), c("ID", "Type"))
descriptor_cols <- descriptor_cols[
  sapply(train_data[, descriptor_cols, drop = FALSE], is.numeric)
]
cat("\nTotal numeric descriptors:", length(descriptor_cols), "\n")


# ============================================================
# 3. REMOVE DUPLICATE SAMPLES (within each class separately)
# ============================================================

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


# ============================================================
# FIX 1 - SPLIT FIRST, THEN FEATURE SELECTION
# ============================================================

cat("\n--- Train/Test Split (BEFORE feature selection) ---\n")

set.seed(123)
idx         <- createDataPartition(train_df$Type, p = 0.80, list = FALSE)
train_inner <- train_df[idx, ]
test_outer  <- train_df[-idx, ]

cat("Train:", nrow(train_inner), "| Test:", nrow(test_outer), "\n")
cat("Train distribution:\n"); print(table(train_inner$Type))
cat("Test distribution:\n");  print(table(test_outer$Type))


# ============================================================
# FIX 2 - FEATURE SELECTION ON train_inner ONLY, NO UPPER CEILING
# ============================================================

cat("\n--- Feature Selection (train_inner only) ---\n")

feature_auc <- sapply(descriptor_cols, function(col) {
  tryCatch(
    as.numeric(auc(roc(train_inner$Type, train_inner[[col]], quiet = TRUE))),
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
  stop("ERROR: Only ", length(selected_features), " features. Cannot proceed.")
}

descriptor_cols <- selected_features

cat("\nTop 10 features by AUC:\n")
print(round(head(sort(feature_auc[descriptor_cols], decreasing = TRUE), 10), 4))

# Update all dataframes
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

# Diagnostics table
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


# ============================================================
# FIX 3 - RECIPES WITH CORRELATION FILTER
# FIX 4 - SMOTE IF STILL IMBALANCED
# ============================================================

sampling_method <- if (imbalance_ratio > 3) "smote" else NULL

add_sampling <- function(recipe_obj) {
  if (!is.null(sampling_method)) {
    recipe_obj <- recipe_obj %>% step_smote(Type, over_ratio = 1)
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


# ============================================================
# 3D PCA VISUALIZATION
# ============================================================

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

pca_step_idx       <- which(sapply(pca_prepped$steps, function(s) inherits(s, "step_pca")))
pca_obj            <- pca_prepped$steps[[pca_step_idx]]$res
variance_explained <- pca_obj$sdev[1:3]^2 / sum(pca_obj$sdev^2) * 100

cat("Variance explained:\n")
cat("  PC1:", round(variance_explained[1], 2), "%\n")
cat("  PC2:", round(variance_explained[2], 2), "%\n")
cat("  PC3:", round(variance_explained[3], 2), "%\n")
cat("  Total:", round(sum(variance_explained), 2), "%\n")

make_axis <- function(label, pct) list(title = paste0(label, " (", round(pct,1), "%)"))

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
                  round(variance_explained[1],1), "%, PC2: ",
                  round(variance_explained[2],1), "%, PC3: ",
                  round(variance_explained[3],1), "%"),
    font = list(size = 13)
  ),
  scene = list(
    xaxis = make_axis("PC1", variance_explained[1]),
    yaxis = make_axis("PC2", variance_explained[2]),
    zaxis = make_axis("PC3", variance_explained[3])
  )
)

pca_combined <- rbind(
  cbind(pca_data_train[, c("PC1","PC2","PC3")],
        Type = as.character(pca_data_train$Type), Dataset = "Training"),
  cbind(pca_data_new[,   c("PC1","PC2","PC3")],
        Type = "New Sample",                      Dataset = "New")
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
                    round(variance_explained[1],1), "%, PC2: ",
                    round(variance_explained[2],1), "%, PC3: ",
                    round(variance_explained[3],1), "%"),
      font = list(size = 13)
    ),
    scene = list(
      xaxis = make_axis("PC1", variance_explained[1]),
      yaxis = make_axis("PC2", variance_explained[2]),
      zaxis = make_axis("PC3", variance_explained[3])
    )
  )

saveWidget(fig_train,    paste0(target_class_r, "_PCA_3D_Training.html"),  selfcontained = TRUE)
saveWidget(fig_combined, paste0(target_class_r, "_PCA_3D_Combined.html"),  selfcontained = TRUE)
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


# ============================================================
# CV CONTROL
# ============================================================

cv_ctrl <- trainControl(
  method          = "repeatedcv",
  number          = 5,
  repeats         = 3,
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  allowParallel   = FALSE
)


# ============================================================
# FIX 5 - RF TUNING GRID (dynamic mtry cap to avoid ranger error)
# ============================================================

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


# ============================================================
# FIX 7 - PLATT SCALING CALIBRATION
# ============================================================

calibrate_probs <- function(raw_probs, true_labels, new_raw_probs, positive_class) {
  binary_labels <- as.integer(true_labels == positive_class)
  cal_df        <- data.frame(prob = raw_probs, label = binary_labels)
  platt_model   <- glm(label ~ prob, data = cal_df, family = binomial)
  new_df_cal    <- data.frame(prob = new_raw_probs)
  as.numeric(predict(platt_model, newdata = new_df_cal, type = "response"))
}


# ============================================================
# MULTI-RUN TRAINING LOOP
# ============================================================

cat("\n--- Multi-Run Training (", n_runs, "runs) ---\n")

rf_probs      <- matrix(NA, nrow(new_data), n_runs)
nn_pca_probs  <- matrix(NA, nrow(new_data), n_runs)
nn_desc_probs <- matrix(NA, nrow(new_data), n_runs)

rf_auc      <- numeric(n_runs)
nn_pca_auc  <- numeric(n_runs)
nn_desc_auc <- numeric(n_runs)

rownames(rf_probs) <- rownames(nn_pca_probs) <-
  rownames(nn_desc_probs) <- new_data$ID

for (i in seq_len(n_runs)) {

  cat("\n===== RUN", i, "of", n_runs, "=====\n")
  set.seed(1000 + i)

  # Calibration split: 20% of train_inner held out for Platt scaling
  cal_idx <- createDataPartition(train_inner$Type, p = 0.20, list = FALSE)
  cal_set <- train_inner[ cal_idx, ]
  fit_set <- train_inner[-cal_idx, ]

  # ---- A) Random Forest ----------------------------------------
  cat("  [RF] Training...\n")

  # Build grid dynamically based on actual post-recipe column count
  prepped_rf <- prep(rec_rf, training = fit_set)
  baked_rf   <- bake(prepped_rf, new_data = NULL)
  n_cols_rf  <- ncol(baked_rf) - 1
  rf_grid    <- make_rf_grid(n_cols_rf)
  cat("  [RF] Post-corr columns:", n_cols_rf, "| Grid:", nrow(rf_grid), "combos\n")

  rf_model <- train(
    rec_rf,
    data       = fit_set,
    method     = "ranger",
    trControl  = cv_ctrl,
    tuneGrid   = rf_grid,
    metric     = "ROC",
    importance = "permutation"
  )

  rf_cal_raw <- predict(rf_model, cal_set,    type = "prob")[, target_class_r]
  rf_new_raw <- predict(rf_model, new_df,     type = "prob")[, target_class_r]
  rf_tst_raw <- predict(rf_model, test_outer, type = "prob")[, target_class_r]

  rf_probs[, i] <- calibrate_probs(rf_cal_raw, cal_set$Type, rf_new_raw,  target_class_r)
  rf_tst_cal    <- calibrate_probs(rf_cal_raw, cal_set$Type, rf_tst_raw,  target_class_r)
  rf_auc[i]     <- as.numeric(auc(roc(test_outer$Type, rf_tst_cal, quiet = TRUE)))

  cat("  [RF] Test AUC:", round(rf_auc[i], 4), "\n")
  cat("  [RF] Confusion matrix:\n")
  print(table(Actual = test_outer$Type, Predicted = rf_tst_cal > 0.5))

  # ---- B) Neural Network + PCA ---------------------------------
  cat("\n  [NN-PCA] Training...\n")

  nn_pca_model <- train(
    rec_nn_pca,
    data      = fit_set,
    method    = "nnet",
    trControl = cv_ctrl,
    tuneGrid  = expand.grid(size = 2:5, decay = c(0.05, 0.1, 0.3)),
    metric    = "ROC",
    trace     = FALSE,
    maxit     = 1000
  )

  nn_pca_cal_raw <- predict(nn_pca_model, cal_set,    type = "prob")[, target_class_r]
  nn_pca_new_raw <- predict(nn_pca_model, new_df,     type = "prob")[, target_class_r]
  nn_pca_tst_raw <- predict(nn_pca_model, test_outer, type = "prob")[, target_class_r]

  nn_pca_probs[, i] <- calibrate_probs(nn_pca_cal_raw, cal_set$Type, nn_pca_new_raw, target_class_r)
  nn_pca_tst_cal    <- calibrate_probs(nn_pca_cal_raw, cal_set$Type, nn_pca_tst_raw, target_class_r)
  nn_pca_auc[i]     <- as.numeric(auc(roc(test_outer$Type, nn_pca_tst_cal, quiet = TRUE)))

  cat("  [NN-PCA] Test AUC:", round(nn_pca_auc[i], 4), "\n")
  cat("  [NN-PCA] Confusion matrix:\n")
  print(table(Actual = test_outer$Type, Predicted = nn_pca_tst_cal > 0.5))

  # ---- C) Neural Network + ANOVA descriptors -------------------
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

  train_red <- fit_set[,    c("Type", top_feats)]
  cal_red   <- cal_set[,    c("Type", top_feats)]
  test_red  <- test_outer[,  top_feats, drop = FALSE]
  new_red   <- new_df[,      top_feats, drop = FALSE]

  rec_nn_desc <- recipe(Type ~ ., data = train_red) %>%
    step_zv(all_predictors()) %>%
    step_impute_median(all_predictors()) %>%
    step_corr(all_numeric_predictors(), threshold = corr_threshold) %>%
    step_normalize(all_predictors()) %>%
    add_sampling()

  nn_desc_model <- train(
    rec_nn_desc,
    data      = train_red,
    method    = "nnet",
    trControl = cv_ctrl,
    tuneGrid  = expand.grid(size = 2:5, decay = c(0.05, 0.1, 0.3)),
    metric    = "ROC",
    trace     = FALSE,
    maxit     = 1000
  )

  nn_desc_cal_raw <- predict(nn_desc_model, cal_red,  type = "prob")[, target_class_r]
  nn_desc_new_raw <- predict(nn_desc_model, new_red,  type = "prob")[, target_class_r]
  nn_desc_tst_raw <- predict(nn_desc_model, test_red, type = "prob")[, target_class_r]

  nn_desc_probs[, i] <- calibrate_probs(nn_desc_cal_raw, cal_set$Type, nn_desc_new_raw, target_class_r)
  nn_desc_tst_cal    <- calibrate_probs(nn_desc_cal_raw, cal_set$Type, nn_desc_tst_raw, target_class_r)
  nn_desc_auc[i]     <- as.numeric(auc(roc(test_outer$Type, nn_desc_tst_cal, quiet = TRUE)))

  cat("  [NN-ANOVA] Test AUC:", round(nn_desc_auc[i], 4), "\n")
  cat("  [NN-ANOVA] Confusion matrix:\n")
  print(table(Actual = test_outer$Type, Predicted = nn_desc_tst_cal > 0.5))
}


# ============================================================
# STABILITY SUMMARIES
# ============================================================

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


# ============================================================
# FIX 6 - AUC-WEIGHTED ENSEMBLE
# ============================================================

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


# ============================================================
# MODEL EVALUATION
# ============================================================

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


# ============================================================
# CONSENSUS
# ============================================================

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


# ============================================================
# FEATURE IMPORTANCE
# ============================================================

cat("\n--- FEATURE IMPORTANCE (last run RF) ---\n")

importance_scores <- ranger::importance(rf_model$finalModel)
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


# ============================================================
# PREDICTIONS
# ============================================================

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
    c("ID","Ensemble_Probability","Confidence","Models_Agree_90")], 10))
}

simple_predictions <- data.frame(
  ID              = new_data$ID,
  Is_Target       = ifelse(ensemble_prob > 0.5, "YES", "NO"),
  Probability_pct = round(ensemble_prob * 100, 1),
  Confidence      = clear_predictions$Confidence
)
colnames(simple_predictions)[2] <- paste0("Is_", target_class_r)
colnames(simple_predictions)[3] <- paste0(target_class_r, "_Probability_%")


# ============================================================
# WRITE EXCEL OUTPUT
# ============================================================

cat("\n--- WRITING EXCEL FILE ---\n")

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
    "Fixes Applied (v2)", "", "", "", "", "", "",
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
    paste("Predicted as", target_class_r, "(50%):", sum(simple_predictions[,2] == "YES")),
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

output_filename <- paste0(target_class_r, "_Complete_Analysis_v2_clean.xlsx")
saveWorkbook(wb, output_filename, overwrite = TRUE)

cat("\n--- SCRIPT COMPLETE ---\n")
cat("Excel  :", output_filename, "\n")
cat("PCA    :", paste0(target_class_r, "_PCA_3D_Training.html"), "\n")
cat("       ", paste0(target_class_r, "_PCA_3D_Combined.html"), "\n")

sessionInfo()

# ================================ END ==========================================
