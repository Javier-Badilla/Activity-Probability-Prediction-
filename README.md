# Activity-Probability-Prediction-

# Peptide Activity Prediction - Binary Classification Pipeline

A robust machine learning pipeline for predicting molecular activity using ensemble methods (Random Forest + Neural Networks) with interactive 3D PCA visualization.

## Overview

This repository contains a complete binary classification system for identifying compounds with specific biological activities from unlabeled molecular descriptor datasets. The pipeline uses ensemble machine learning with multiple validation strategies and provides interactive visualizations.

### Key Features

- **Ensemble Learning**: Combines Random Forest, NN-PCA, and NN-ANOVA models
- **Multi-run Stability Analysis**: Configurable runs (2-50) for robust predictions
- **Interactive 3D PCA Visualization**: Explore class separation in principal component space
- **Automated Feature Selection**: AUC-based descriptor filtering
- **Data Leakage Prevention**: Automatic duplicate detection and removal
- **Comprehensive Reports**: Excel outputs with 10 detailed sheets
- **Confidence Scoring**: Multiple threshold levels and model consensus metrics

## Quick Start

The main analysis pipeline is located in:

[Analysis_pipeline](Code/PIPELINE_V6.R)

### Prerequisites

- R version 4.3 or higher
- Required R packages:

```r
install.packages(c(
  "readxl", "openxlsx", "caret", "ranger", 
  "nnet", "pROC", "recipes", "dplyr", "plotly"
))
```

### Basic Usage

1. **Prepare your data**:
   - `train_data`: Labeled samples (ID, Type, descriptors...)
   - `new_data`: Unlabeled samples to classify (ID, descriptors...)

2. **Configure parameters** in the script:
```r
target_class <- "AFP"           # Your target activity
comparison_class <- NULL        # NULL for one-vs-all
n_runs <- 20                    # Number of iterations (2-5 for quick, 20+ for robust)
```

3. **Run the script**:

4. **Review outputs**:
   - Excel file: `(target_class)_Complete_Analysis.xlsx`
   - 3D plots: `(target_class)_PCA_3D_Training.html` and `(target_class)_PCA_3D_Combined.html`
     
## Documentation

### Manuals Available
- Both in Spanish and English (will be added soon)
 
### User Manual Contents
- Complete workflow explanation
- Step-by-step instructions
- Parameter configuration guide
- Troubleshooting section

### Output Guide Contents
- Detailed explanation of all 10 Excel sheets
- Interpretation of 3D PCA visualizations
- How to read confidence scores
- Feature importance analysis
- Model performance metrics

## How It Works

### Pipeline Workflow

```
graph TD
    A[Load Data] --> B[Clean & Filter]
    B --> C[Feature Selection AUC]
    C --> D[Remove Duplicates]
    D --> E[Train/Test Split 80/20]
    E --> F[3D PCA Visualization]
    F --> G[Train 3 Models x N runs]
    G --> H[Stability Analysis]
    H --> I[Consensus Predictions]
    I --> J[Export Results]
```

### Models Used

1. **Random Forest (RF)**
   - Tree-based ensemble
   - Handles non-linear relationships
   - Provides feature importance

2. **Neural Network with PCA (NN_PCA)**
   - Dimensionality reduction (90% variance)
   - Good for high-dimensional data

3. **Neural Network with ANOVA (NN_DESC)**
   - Statistical feature selection
   - Top 30 most discriminative features
   - Most interpretable

## Output Files

### Excel File: `(target_class)_Complete_Analysis.xlsx`

Contains 10 sheets:

1. **Summary** - Executive overview
2. **Simple_Predictions** - YES/NO format with probabilities
3. **Detailed_Predictions** - Multi-threshold predictions with confidence
4. **High_Confidence** - Top candidates for validation
5. **Consensus** - Model agreement analysis
6. **Model_Evaluation** - AUC metrics and stability
7. **Feature_Importance** - Descriptor rankings
8. **Top_Features_Stats** - Statistical analysis of key features
9. **Diagnostics** - Dataset splits and class balance
10. **PCA_Coordinates** - 3D component values

### Interactive HTML Files

- **Training Plot**: Class separation visualization
- **Combined Plot**: Training data + new samples in PCA space

## Data Requirements

### Minimum Requirements
- At least 50 unique samples in training data
- Numeric molecular descriptors
- Class ratio ideally ≤ 3:1

### Data Format

**Training Data:**
| ID | Type | Descriptor1 | Descriptor2 | Etc |
|----|------|-------------|-------------|-----|

**New Data:**
| ID | Descriptor1 | Descriptor2 | Etc |
|----|-------------|-------------|-----|

## Author

**Javier Badilla**

## Acknowledgments

- Built with [caret](https://topepo.github.io/caret/) framework
- Visualizations powered by [plotly](https://plotly.com/r/)
- Random Forest implementation: [ranger](https://github.com/imbs-hl/ranger)

## Documentation

The user manuals are available in Spanish:

- 🇪🇸 [Manual de uso del código](Code_esp/MANUAL_SCRIPT_V6.pdf)
- 🇪🇸 [Manual de interpretación de resultados](Code_esp/MANUAL_RESULTADOS_SCRIPT_V6.pdf)

English versions will be added soon.


## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{badilla2026peptide,
  author = {Badilla, Javier},
  title = {Peptide Activity Prediction - Classification Pipeline},
  year = {2026},
  url = {https://github.com/Javier-Badilla/Activity-Probability-Prediction-}
}
```

**Version:** 2.0  
**Last Updated:** January 2026  
**R Version Required:** 4.3+


