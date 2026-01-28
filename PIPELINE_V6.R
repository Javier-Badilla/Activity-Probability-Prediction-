# ============================================================
# Project: Peptide Activity Prediction
# Script: Random Forest + NN + PCA Pipeline
# Author: Javier BAdilla
# Date: 2025 - 2026
# Description:
#
# This script performs data preprocessing, descriptor selection,
# PCA dimensionality reduction, model training (RF + NN),
# stability analysis and prediction on new samples.
#
# Descripción:
#  
# Este script realiza un preprocesamiento de datos, selección de descriptores.
# Reducción de la dimensionallidad por PCAs, entrenamiento de modelos (RF + NN),
# Análisis de estabilidad y predicción de actividad.
#
# ---------------------- REQUIREMENTS ------------------------
#
# R version: 4.3+
# Packages needed:
# readxl, openxlsx, caret, ranger, nnet, pROC, dplyr, recepies, plotly
#
# ---------------------- INPUT DATA FORMAT -------------------
#
# Two Excel files are needed:
# train_data : labeled data (ID + Type, used for training)
# new_data   : unlabeled data (ID only, used for prediction)
#
# File structure:
# Column 1 (ID): Sample ID
# Column 2 (Type): Class label 
# Columns 3+ : Molecular descriptors (numeric)
# 
# ----------------------- LIBRARIES -------------------------

library(readxl)
library(openxlsx)
library(caret)
library(nnet)
library(ranger)
library(pROC)
library(recipes)
library(dplyr)
library(plotly)

# ----------------------- USER INPUTS/Modifiable Inputs ------------------------

target_class <- "AFP"
comparison_class <- "antiAChE"  # NULL si la idea es uno contra todos

# Run number and probability threshold (aumenta el numero para hacerlo mas robusto, toma mas tiempo)

n_runs <- 20
prob_threshold <- 0.90

# Select moderate features (ajusta para mejor resultado)

auc_lower <- 0.60
auc_upper <- 0.90

# ----------------------- CLEANING --------------------

cat("\n--- DATA LOADING AND CLEANING ---\n")

colnames(train_data) <- make.names(colnames(train_data))
colnames(new_data)   <- make.names(colnames(new_data))

# Filter to specific classes if comparison_class is set

if (!is.null(comparison_class)) {
  cat("Filtering to", target_class, "vs", comparison_class, "only\n")
  train_data <- train_data[train_data$Type %in% c(target_class, comparison_class), ]
  cat("Samples after filtering:", nrow(train_data), "\n")
}

# Create binary classification

train_data$Type <- as.factor(train_data$Type)
train_data$Type <- ifelse(train_data$Type == target_class, target_class, "other")
train_data$Type <- factor(train_data$Type, levels = c(target_class, "other"))
train_data$Type <- droplevels(train_data$Type) 

# Make valid R names

levels(train_data$Type) <- make.names(levels(train_data$Type))
target_class <- make.names(target_class)

# Handle empty strings

train_data[train_data == ""] <- NA
new_data[new_data == ""]     <- NA

cat("Final class distribution:\n")
print(table(train_data$Type))

# ----------------------- DESCRIPTORS -----------------------

descriptor_cols <- setdiff(colnames(train_data), c("ID", "Type"))
descriptor_cols <- descriptor_cols[
  sapply(train_data[, descriptor_cols, drop = FALSE], is.numeric)
]
cat("\nTotal descriptors:", length(descriptor_cols), "\n")

# ----------------------- REMOVE DUPLICATES -----------------

cat("\n=== Removing Duplicate Samples ===\n")
train_data_unique <- train_data[!duplicated(train_data[, descriptor_cols]), ]
cat("Removed", nrow(train_data) - nrow(train_data_unique), "duplicate rows\n")
cat("Remaining samples:", nrow(train_data_unique), "\n")

# Create working dataframes

train_df <- train_data_unique[, c("Type", descriptor_cols)]
new_df   <- new_data[, descriptor_cols, drop = FALSE]

# Check sample number

if (nrow(train_df) < 50) {
  stop("ERROR: Too few samples (", nrow(train_df), "). Need at least 50 samples.")
}

# ----------------------- FEATURE SELECTION -----------------

# Calculate AUC for all descriptors

feature_auc <- sapply(descriptor_cols, function(col) {
  tryCatch(
    auc(roc(train_df$Type, train_df[[col]], quiet = TRUE)),
    error = function(e) NA
  )
})

# Remove features with missing AUC

feature_auc <- feature_auc[!is.na(feature_auc)]

# Show AUC distribution

cat("\nAUC distribution of features:\n")
cat("  AUC > 0.95:", sum(feature_auc > 0.95), "\n")
cat("  AUC 0.85-0.95:", sum(feature_auc > 0.85 & feature_auc <= 0.95), "\n")
cat("  AUC 0.75-0.85:", sum(feature_auc > 0.75 & feature_auc <= 0.85), "\n")
cat("  AUC 0.65-0.75:", sum(feature_auc > 0.65 & feature_auc <= 0.75), "\n")
cat("  AUC 0.55-0.65:", sum(feature_auc > 0.55 & feature_auc <= 0.65), "\n")
cat("  AUC < 0.55:", sum(feature_auc <= 0.55), "\n")

# selection of elements acording the AUC

selected_features <- names(feature_auc[feature_auc >= auc_lower & feature_auc <= auc_upper])

cat("\nFeatures with AUC", auc_lower, "-", auc_upper, ":", length(selected_features), "\n")

# Ensure we have enough features

if (length(selected_features) < 10) {
  cat("*** WARNING: Too few moderate features, expanding range ***\n")
  auc_lower <- 0.50  # this stays as is, unless persistent problems
  auc_upper <- 0.90  # this stays as is, unless persistent problems
  selected_features <- names(feature_auc[feature_auc >= auc_lower & feature_auc <= auc_upper])
  cat("Expanded range", auc_lower, "-", auc_upper, ":", length(selected_features), "features\n")
}

if (length(selected_features) < 5) {
  stop("ERROR: Only ", length(selected_features), " features available. Cannot proceed.")
}

# Update descriptor columns

descriptor_cols <- selected_features

# Show top features by AUC

top_features <- head(sort(feature_auc[descriptor_cols], decreasing = TRUE), 10)
cat("\nTop 10 selected features by AUC:\n")
print(round(top_features, 4))

# Update dataframes with selected features only

train_df <- train_df[, c("Type", descriptor_cols)]
new_df <- new_df[, descriptor_cols, drop = FALSE]

# ----------------------- TRAIN / TEST SPLIT ----------------

cat("\n--- Creating Train/Test Split ---\n")

set.seed(123)
idx <- createDataPartition(train_df$Type, p = 0.80, list = FALSE)
train_inner <- train_df[idx, ]
test_outer  <- train_df[-idx, ]

# Diagnostics

test_diagnostics <- data.frame(
  Metric = c(
    "Total_Train_Size",
    "Total_Test_Size",
    paste0("Train_", target_class),
    "Train_other",
    paste0("Test_", target_class),
    "Test_other",
    paste0("Test_", target_class, "_Proportion"),
    "Test_other_Proportion",
    "Features_Used"
  ),
  Value = c(
    nrow(train_inner),
    nrow(test_outer),
    sum(train_inner$Type == target_class),
    sum(train_inner$Type == "other"),
    sum(test_outer$Type == target_class),
    sum(test_outer$Type == "other"),
    round(sum(test_outer$Type == target_class) / nrow(test_outer), 3),
    round(sum(test_outer$Type == "other") / nrow(test_outer), 3),
    length(descriptor_cols)
  )
)

cat("\n--- DATASET SUMMARY ---\n")
print(test_diagnostics)

# Check for duplicates

train_desc <- train_inner[, descriptor_cols]
test_desc <- test_outer[, descriptor_cols]
combined <- rbind(train_desc, test_desc)
dup_idx <- duplicated(combined) | duplicated(combined, fromLast = TRUE)
test_dups <- dup_idx[-(1:nrow(train_desc))]

cat("\nDuplicate samples between train/test:", sum(test_dups), "\n")
if (sum(test_dups) > 0) {
  cat("*** WARNING: Train/test leakage detected! ***\n")
}

# ----------------------- RECIPES ---------------------------

rec_rf <- recipe(Type ~ ., data = train_inner) %>%
  step_zv(all_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_impute_median(all_numeric_predictors())

rec_nn_pca <- recipe(Type ~ ., data = train_inner) %>%
  step_zv(all_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_pca(all_predictors(), threshold = 0.90)

# -------------------------- 3D PCA VISUALIZATION ----------------------

cat("\n--- GENERATING 3D PCA PLOT ---\n")

# Prepare PCA recipe for visualization (keeping first 3 components explicitly)
rec_pca_viz <- recipe(Type ~ ., data = train_inner) %>%
  step_zv(all_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_impute_median(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_pca(all_predictors(), num_comp = 3)

# Prep and bake the recipe
pca_prepped <- prep(rec_pca_viz, training = train_inner)
pca_data_train <- bake(pca_prepped, new_data = train_inner)
pca_data_new <- bake(pca_prepped, new_data = new_df)

# Get variance explained
pca_obj <- pca_prepped$steps[[5]]$res  # PCA step is the 5th step
variance_explained <- pca_obj$sdev[1:3]^2 / sum(pca_obj$sdev^2) * 100

cat("Variance explained by first 3 PCs:\n")
cat("  PC1:", round(variance_explained[1], 2), "%\n")
cat("  PC2:", round(variance_explained[2], 2), "%\n")
cat("  PC3:", round(variance_explained[3], 2), "%\n")
cat("  Total:", round(sum(variance_explained), 2), "%\n")

# Create 3D plot for training data
fig_train <- plot_ly(
  data = pca_data_train,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~Type,
  colors = c("#1f77b4", "#ff7f0e"),  # Blue for target, orange for other
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 5, opacity = 0.7),
  text = ~paste("Class:", Type),
  hoverinfo = "text"
) %>%
  layout(
    title = list(
      text = paste0("3D PCA Visualization - Training Data\n",
                    "PC1: ", round(variance_explained[1], 1), "%, ",
                    "PC2: ", round(variance_explained[2], 1), "%, ",
                    "PC3: ", round(variance_explained[3], 1), "%"),
      font = list(size = 14)
    ),
    scene = list(
      xaxis = list(title = paste0("PC1 (", round(variance_explained[1], 1), "%)")),
      yaxis = list(title = paste0("PC2 (", round(variance_explained[2], 1), "%)")),
      zaxis = list(title = paste0("PC3 (", round(variance_explained[3], 1), "%)"))
    ),
    legend = list(title = list(text = "Class"))
  )

# Save training data 3D plot as HTML
htmlwidgets::saveWidget(
  fig_train,
  file = paste0(target_class, "_PCA_3D_Training.html"),
  selfcontained = TRUE
)

cat("3D PCA plot for training data saved to:", paste0(target_class, "_PCA_3D_Training.html"), "\n")

# Create 3D plot with both training and new data
pca_combined <- rbind(
  cbind(pca_data_train[, c("PC1", "PC2", "PC3")], 
        Type = as.character(pca_data_train$Type),
        Dataset = "Training"),
  cbind(pca_data_new[, c("PC1", "PC2", "PC3")], 
        Type = "New Sample",
        Dataset = "New")
)

fig_combined <- plot_ly() %>%
  add_trace(
    data = pca_combined[pca_combined$Dataset == "Training" & pca_combined$Type == target_class, ],
    x = ~PC1, y = ~PC2, z = ~PC3,
    type = "scatter3d",
    mode = "markers",
    name = paste("Training:", target_class),
    marker = list(size = 5, opacity = 0.6, color = "#1f77b4"),
    text = ~paste("Training -", Type),
    hoverinfo = "text"
  ) %>%
  add_trace(
    data = pca_combined[pca_combined$Dataset == "Training" & pca_combined$Type != target_class, ],
    x = ~PC1, y = ~PC2, z = ~PC3,
    type = "scatter3d",
    mode = "markers",
    name = "Training: other",
    marker = list(size = 5, opacity = 0.6, color = "#ff7f0e"),
    text = ~paste("Training - other"),
    hoverinfo = "text"
  ) %>%
  add_trace(
    data = pca_combined[pca_combined$Dataset == "New", ],
    x = ~PC1, y = ~PC2, z = ~PC3,
    type = "scatter3d",
    mode = "markers",
    name = "New Samples",
    marker = list(size = 7, opacity = 0.9, color = "#2ca02c", symbol = "diamond"),
    text = ~paste("New Sample"),
    hoverinfo = "text"
  ) %>%
  layout(
    title = list(
      text = paste0("3D PCA - Training vs New Samples\n",
                    "PC1: ", round(variance_explained[1], 1), "%, ",
                    "PC2: ", round(variance_explained[2], 1), "%, ",
                    "PC3: ", round(variance_explained[3], 1), "%"),
      font = list(size = 14)
    ),
    scene = list(
      xaxis = list(title = paste0("PC1 (", round(variance_explained[1], 1), "%)")),
      yaxis = list(title = paste0("PC2 (", round(variance_explained[2], 1), "%)")),
      zaxis = list(title = paste0("PC3 (", round(variance_explained[3], 1), "%)"))
    ),
    legend = list(title = list(text = "Dataset"))
  )

# Save combined 3D plot as HTML
htmlwidgets::saveWidget(
  fig_combined,
  file = paste0(target_class, "_PCA_3D_Combined.html"),
  selfcontained = TRUE
)

cat("3D PCA plot with new samples saved to:", paste0(target_class, "_PCA_3D_Combined.html"), "\n")

# Create a data frame with PCA coordinates for export
pca_coordinates <- data.frame(
  ID = c(train_data_unique$ID[idx], new_data$ID),
  PC1 = c(pca_data_train$PC1, pca_data_new$PC1),
  PC2 = c(pca_data_train$PC2, pca_data_new$PC2),
  PC3 = c(pca_data_train$PC3, pca_data_new$PC3),
  Type = c(as.character(pca_data_train$Type), rep("New Sample", nrow(pca_data_new))),
  Dataset = c(rep("Training", nrow(pca_data_train)), rep("New", nrow(pca_data_new)))
)

cat("\n--- 3D PCA VISUALIZATION COMPLETE ---\n")

# ----------------------- CV CONTROL ------------------------

cv_ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 3,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  allowParallel = FALSE
)

# ---------------- MULTI-RUN ----------------------

cat("\n--- Multi-Run Training ---\n")

rf_probs      <- matrix(NA, nrow(new_data), n_runs)
nn_pca_probs  <- matrix(NA, nrow(new_data), n_runs)
nn_desc_probs <- matrix(NA, nrow(new_data), n_runs)

rf_auc      <- numeric(n_runs)
nn_pca_auc  <- numeric(n_runs)
nn_desc_auc <- numeric(n_runs)

rownames(rf_probs) <- rownames(nn_pca_probs) <-
  rownames(nn_desc_probs) <- new_data$ID

for (i in seq_len(n_runs)) {
  
  cat("\n--- RUN", i, "of", n_runs, "---\n")
  set.seed(1000 + i)
  
  # ------------- Random Forest ------------------
  
  cat("Training Random Forest\n")
  rf_model <- train(
    rec_rf,
    data = train_inner,
    method = "ranger",
    trControl = cv_ctrl,
    tuneGrid = expand.grid(
      mtry = floor(sqrt(length(descriptor_cols))),
      splitrule = "gini",
      min.node.size = 5
    ),
    metric = "ROC",
    importance = "permutation"
  )
  
  rf_probs[, i] <- predict(rf_model, new_df, type = "prob")[, target_class]
  rf_test_pred <- predict(rf_model, test_outer, type = "prob")[, target_class]
  rf_auc[i] <- auc(roc(test_outer$Type, rf_test_pred, quiet = TRUE))
  
  cat("  RF Test AUC:", round(rf_auc[i], 4), "\n")
  cat("  Prediction range:", round(min(rf_test_pred), 4), "to", round(max(rf_test_pred), 4), "\n")
  cat("  Confusion matrix:\n")
  print(table(Actual = test_outer$Type, Predicted = rf_test_pred > 0.5))
  
  # ------------------- NN PCA -------------------------
  
  cat("\nTraining Neural Network (PCA)\n")
  nn_pca_model <- train(
    rec_nn_pca,
    data = train_inner,
    method = "nnet",
    trControl = cv_ctrl,
    tuneGrid = expand.grid(size = 2:4, decay = c(0.1, 0.3)),
    metric = "ROC",
    trace = FALSE,
    maxit = 1000
  )
  
  nn_pca_probs[, i] <- predict(nn_pca_model, new_df, type = "prob")[, target_class]
  nn_pca_test_pred <- predict(nn_pca_model, test_outer, type = "prob")[, target_class]
  nn_pca_auc[i] <- auc(roc(test_outer$Type, nn_pca_test_pred, quiet = TRUE))
  
  cat("  NN_PCA Test AUC:", round(nn_pca_auc[i], 4), "\n")
  cat("  Prediction range:", round(min(nn_pca_test_pred), 4), "to", round(max(nn_pca_test_pred), 4), "\n")
  cat("  Confusion matrix:\n")
  print(table(Actual = test_outer$Type, Predicted = nn_pca_test_pred > 0.5))
  
  # ------------------ ANOVA + NN Descriptor ----------------------------
  
  cat("\nTraining Neural Network (ANOVA features)\n")
  anova_pvals <- sapply(
    train_inner[, descriptor_cols, drop = FALSE],
    function(x) tryCatch(
      summary(aov(x ~ train_inner$Type))[[1]][["Pr(>F)"]][1],
      error = function(e) NA
    )
  )
  
  top_n <- min(30, sum(!is.na(anova_pvals)))
  top_features <- names(sort(anova_pvals))[1:top_n]
  
  train_red <- train_inner[, c("Type", top_features)]
  new_red   <- new_df[, top_features, drop = FALSE]
  test_red  <- test_outer[, top_features, drop = FALSE]
  
  rec_nn_desc <- recipe(Type ~ ., data = train_red) %>%
    step_zv(all_predictors()) %>%
    step_impute_median(all_predictors()) %>%
    step_normalize(all_predictors())
  
  nn_desc_model <- train(
    rec_nn_desc,
    data = train_red,
    method = "nnet",
    trControl = cv_ctrl,
    tuneGrid = expand.grid(size = 2:4, decay = c(0.1, 0.3)),
    metric = "ROC",
    trace = FALSE,
    maxit = 1000
  )
  
  nn_desc_probs[, i] <- predict(nn_desc_model, new_red, type = "prob")[, target_class]
  nn_desc_test_pred <- predict(nn_desc_model, test_red, type = "prob")[, target_class]
  nn_desc_auc[i] <- auc(roc(test_outer$Type, nn_desc_test_pred, quiet = TRUE))
  
  cat("  NN_DESC Test AUC:", round(nn_desc_auc[i], 4), "\n")
  cat("  Prediction range:", round(min(nn_desc_test_pred), 4), "to", round(max(nn_desc_test_pred), 4), "\n")
  cat("  Confusion matrix:\n")
  print(table(Actual = test_outer$Type, Predicted = nn_desc_test_pred > 0.5))
}

# ---------------- Stability SUMMARIES ----------------------

cat("\n--- Stability ---\n")

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

# ---------------- FINAL PREDICTIONS -----------------------

final_predictions <- stability_summary[, c("ID","RF_Mean","NN_PCA_Mean","NN_Desc_Mean")]

# ---------------- AUC (OUTER TEST SET) --------------------

auc_df <- data.frame(
  Model = c("RF", "NN_PCA", "NN_DESC_ANOVA"),
  Mean_AUC = c(mean(rf_auc), mean(nn_pca_auc), mean(nn_desc_auc)),
  SD_AUC = c(sd(rf_auc), sd(nn_pca_auc), sd(nn_desc_auc)),
  Min_AUC = c(min(rf_auc), min(nn_pca_auc), min(nn_desc_auc)),
  Max_AUC = c(max(rf_auc), max(nn_pca_auc), max(nn_desc_auc))
)

# ---------------- CONSENSUS -------------------------------

consensus <- stability_summary %>%
  mutate(
    Consensus_Mean = rowMeans(select(., RF_Freq, NN_PCA_Freq, NN_Desc_Freq)),
    Models_80 = rowSums(select(., RF_Freq, NN_PCA_Freq, NN_Desc_Freq) >= 0.8),
    Confidence_Tier = case_when(
      Consensus_Mean >= 0.9 & Models_80 == 3 ~ "Very High",
      Consensus_Mean >= 0.8 & Models_80 >= 2 ~ "High",
      TRUE ~ "Low"
    )
  ) %>%
  arrange(desc(Confidence_Tier), desc(Consensus_Mean))

# ---------------- FEATURE IMPORTANCE -----------------------

cat("\n--- FEATURE IMPORTANCE ANALYSIS ---\n")

importance_scores <- ranger::importance(rf_model$finalModel)
top_important <- head(sort(importance_scores, decreasing = TRUE), 20)

cat("\nTop 20 most important features for classification:\n")
print(round(top_important, 4))

importance_df <- data.frame(
  Feature = names(importance_scores),
  Importance = importance_scores
) %>% arrange(desc(Importance))

cat("\n--- FEATURE IMPORTANCE ---\n")

top5_importance <- sum(head(importance_df$Importance, 5))
total_importance <- sum(importance_df$Importance)
top5_percent <- (top5_importance / total_importance) * 100

cat("Top 5 features account for", round(top5_percent, 1), "% of total importance\n")

if (top5_percent > 70) {
  cat("→ Classification relies heavily on just a few features\n")
} else if (top5_percent > 50) {
  cat("→ Classification uses moderate feature diversity\n")
} else {
  cat("→ Classification uses information from many features\n")
}

cat("\n=== Top Feature Statistics by Class ===\n")
top_features_list <- head(importance_df$Feature, 5)

for (feat in top_features_list) {
  cat("\n", feat, ":\n")
  
  if (feat %in% colnames(train_df)) {
    target_values <- train_df[[feat]][train_df$Type == target_class]
    other_values <- train_df[[feat]][train_df$Type == "other"]
    
    target_mean <- mean(target_values, na.rm = TRUE)
    other_mean <- mean(other_values, na.rm = TRUE)
    
    cat("  ", target_class, "mean:", round(target_mean, 4), "\n")
    cat("  other mean:", round(other_mean, 4), "\n")
    cat("  Difference:", round(target_mean - other_mean, 4), "\n")
    
    cat("  ", target_class, "range: [", 
        round(min(target_values, na.rm = TRUE), 2), ",", 
        round(max(target_values, na.rm = TRUE), 2), "]\n")
    cat("  other range: [", 
        round(min(other_values, na.rm = TRUE), 2), ",", 
        round(max(other_values, na.rm = TRUE), 2), "]\n")
  } else {
    cat("  Feature not found in training data\n")
  }
}

feature_stats <- data.frame(
  Feature = character(),
  Target_Mean = numeric(),
  Other_Mean = numeric(),
  Difference = numeric(),
  Importance = numeric(),
  stringsAsFactors = FALSE
)

for (feat in top_features_list) {
  if (feat %in% colnames(train_df)) {
    target_values <- train_df[[feat]][train_df$Type == target_class]
    other_values <- train_df[[feat]][train_df$Type == "other"]
    
    feature_stats <- rbind(feature_stats, data.frame(
      Feature = feat,
      Target_Mean = mean(target_values, na.rm = TRUE),
      Other_Mean = mean(other_values, na.rm = TRUE),
      Difference = mean(target_values, na.rm = TRUE) - mean(other_values, na.rm = TRUE),
      Importance = importance_df$Importance[importance_df$Feature == feat]
    ))
  }
}

print(feature_stats)

# ---------------- CLEAR PREDICTIONS ---------------------

cat("\n--- CLEAR PREDICTIONS ---\n")

ensemble_prob <- rowMeans(cbind(
  stability_summary$RF_Mean,
  stability_summary$NN_PCA_Mean,
  stability_summary$NN_Desc_Mean
))

clear_predictions <- data.frame(
  ID = new_data$ID,
  Ensemble_Probability = round(ensemble_prob, 4),
  Prediction_50 = ifelse(ensemble_prob > 0.5, target_class, "other"),
  Prediction_70 = ifelse(ensemble_prob > 0.7, target_class, "other"),
  Prediction_90 = ifelse(ensemble_prob > 0.9, target_class, "other"),
  Confidence = case_when(
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
  RF_Prediction = ifelse(stability_summary$RF_Mean > 0.5, target_class, "other"),
  NN_PCA_Prediction = ifelse(stability_summary$NN_PCA_Mean > 0.5, target_class, "other"),
  NN_DESC_Prediction = ifelse(stability_summary$NN_Desc_Mean > 0.5, target_class, "other")
)

clear_predictions$Interpretation <- paste0(
  ifelse(clear_predictions$Prediction_50 == target_class, 
         paste0("POSITIVE for ", target_class),
         paste0("NEGATIVE for ", target_class)),
  " (", clear_predictions$Confidence, " confidence)"
)

cat("\n--- PREDICTION SUMMARY ---\n")
cat("Using 50% threshold:\n")
cat("  Predicted as", target_class, ":", 
    sum(clear_predictions$Prediction_50 == target_class), "\n")
cat("  Predicted as other:", 
    sum(clear_predictions$Prediction_50 == "other"), "\n")

cat("\nUsing 70% threshold (high confidence):\n")
cat("  Predicted as", target_class, ":", 
    sum(clear_predictions$Prediction_70 == target_class), "\n")
cat("  Predicted as other:", 
    sum(clear_predictions$Prediction_70 == "other"), "\n")

cat("\nUsing 90% threshold (very high confidence):\n")
cat("  Predicted as", target_class, ":", 
    sum(clear_predictions$Prediction_90 == target_class), "\n")
cat("  Predicted as other:", 
    sum(clear_predictions$Prediction_90 == "other"), "\n")

cat("\nConfidence distribution:\n")
print(table(clear_predictions$Confidence))

cat("\nModel agreement (all 3 models agree at 90%):", 
    sum(clear_predictions$Models_Agree_90 == 3), "samples\n")

high_confidence_positives <- clear_predictions %>%
  filter(Prediction_70 == target_class, Confidence %in% c("High", "Very High")) %>%
  arrange(desc(Ensemble_Probability))

cat("\n=== HIGH CONFIDENCE", target_class, "CANDIDATES ===\n")
cat("Found", nrow(high_confidence_positives), "high-confidence", target_class, "candidates\n")

if (nrow(high_confidence_positives) > 0) {
  cat("\nTop 10 candidates:\n")
  print(head(high_confidence_positives[, c("ID", "Ensemble_Probability", "Confidence", "Models_Agree_90")], 10))
}

simple_predictions <- data.frame(
  ID = new_data$ID,
  Is_Target_Class = ifelse(ensemble_prob > 0.5, "YES", "NO"),
  Probability = round(ensemble_prob * 100, 1),
  Confidence = clear_predictions$Confidence
)

colnames(simple_predictions)[2] <- paste0("Is_", target_class)
colnames(simple_predictions)[3] <- paste0(target_class, "_Probability_%")

# ---------------- WRITE EXCEL FILE ----------------

cat("\n--- CREATING FILE ---\n")

# Create workbook

wb <- createWorkbook()

# Add sheets with data

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

# Write data to sheets

writeData(wb, "1_Summary", 
          data.frame(
            Section = c("Dataset Info", "", "Model Performance", "", "Predictions Summary", "", "", ""),
            Description = c(
              paste("Target Class:", target_class),
              paste("Comparison:", ifelse(is.null(comparison_class), "All others", comparison_class)),
              paste("Best Model AUC:", round(max(auc_df$Mean_AUC), 4)),
              paste("Number of runs:", n_runs),
              paste("Total predictions:", nrow(simple_predictions)),
              paste("Predicted as", target_class, "(50% threshold):", sum(simple_predictions[,2] == "YES")),
              paste("High confidence candidates:", nrow(high_confidence_positives)),
              paste("Models agree (90%):", sum(clear_predictions$Models_Agree_90 == 3))
            )
          ))

writeData(wb, "2_Simple_Predictions", simple_predictions)
writeData(wb, "3_Detailed_Predictions", clear_predictions)
writeData(wb, "4_High_Confidence", high_confidence_positives)
writeData(wb, "5_Consensus", consensus)
writeData(wb, "6_Model_Evaluation", auc_df)
writeData(wb, "7_Feature_Importance", importance_df)
writeData(wb, "8_Top_Features_Stats", feature_stats)
writeData(wb, "9_Diagnostics", test_diagnostics)
writeData(wb, "10_PCA_Coordinates", pca_coordinates)

# Save workbook

output_filename <- paste0(target_class, "_Complete_Analysis.xlsx")
saveWorkbook(wb, output_filename, overwrite = TRUE)

cat("\n--- SCRIPT COMPLETE ---\n")
cat("All results saved to:", output_filename, "\n")
cat("\n--- OUTPUT FILES GENERATED ---\n")
cat("1. Excel file:", output_filename, "\n")
cat("2. 3D PCA plot (training):", paste0(target_class, "_PCA_3D_Training.html"), "\n")
cat("3. 3D PCA plot (combined):", paste0(target_class, "_PCA_3D_Combined.html"), "\n")
cat("\n--- PREDICTION SECTION COMPLETE ---\n")


# ---------------------- SESSION INFO ------------------------

sessionInfo()

# ================================ END ======================================