# Activity-Probability-Prediction

# Peptide Activity Prediction — Computational Screening Pipeline

A multi-stage computational pipeline for identifying peptides with specific biological activities. The system combines QSAR-based in silico screening, molecular docking, interaction profiling, and ensemble machine learning to progressively filter candidates before committing to synthesis and in vitro testing.

---

## Overview

This repository contains two independent but complementary pipelines:

**Pipeline A — AFP Candidate Selection (Python)**
Three-stage in silico screening to identify antifreeze peptides (AFPs) from structural models. Each stage acts as a progressively stricter filter before recommending candidates for experimental validation.

**Pipeline B — Activity Prediction (R)**
Binary classification system for predicting whether a peptide belongs to a target activity class (AFP, antifungal, anticancer, AChE inhibitor, etc.) from pre-calculated molecular descriptors, using an ensemble of machine learning models.

---

## Repository Structure

```
├── Codes/
│   ├── afp_sasa_pipeline.py              # AFP Stage 1 — Surface analysis
│   ├── hex_docking_pipeline.py           # AFP Stage 2 — Molecular docking analysis
│   ├── plip_afp_pipeline.py              # AFP Stage 3 — Interaction profiling + Excel output
│   └── peptide_activity_prediction_v7.R  # ML activity prediction pipeline
│
├── Manuales_esp/
│   ├── Manual_Pipeline_AFP_seleccion.pdf  # AFP pipeline user manual (Spanish)
│   └── Manual_Pipeline_Prediccion_Peptidos.pdf  # ML pipeline user manual (Spanish)
```

---

## Pipeline A — AFP Candidate Selection

### How it works

```
Peptide PDB structures
        │
        ▼
[Stage 1] afp_sasa_pipeline.py
  Surface analysis: SASA, hydrophobic patches,
  ice lattice complementarity, B-factor rigidity
  → AFP score 0–10 per peptide
        │
        ▼  (score ≥ 6 advances)
[Stage 2] hex_docking_pipeline.py
  HEX rigid docking (100 orientations vs TIP4P ice slab)
  + PyMOL visual verification of IBS face orientation
  → Top models extracted per peptide
        │
        ▼  (IBS face confirmed toward ice)
[Stage 3] plip_afp_pipeline.py
  PLIP v3.0.0 interaction profiling
  H-bond classification: IBS vs non-IBS residues
  → Excel report with verdict per peptide
        │
        ▼  (IBS ratio ≥ 0.40)
  Candidates recommended for experimental validation
  DLS → IRI assay → Nanoliter osmometer → DSC
```

### Stage 1 — Surface Analysis (`afp_sasa_pipeline.py`)

Analyzes the 3D structure of each peptide and scores it across four criteria:

| Criterion | Points | What it measures |
|-----------|--------|-----------------|
| Hydrophobic ratio | 0–2 | Fraction of surface area from hydrophobic residues |
| IBS candidates | 0–2 | Number of exposed IBS-type residues (THR, ALA, GLY, SER, VAL, LEU, ILE) |
| Ice lattice complementarity | 0–3 | Cβ spacing matching ice Ih lattice (4.52 Å a-axis, 7.36 Å c-axis) |
| Structural rigidity | 0–3 | B-factor analysis of IBS-candidate residues |

**Quick start:**
```bash
# Place PDB files in ./PDB/ then edit PDB_DIR in the script
python afp_sasa_pipeline.py
```

**Key outputs** (`./afp_results/`):
- `summary_ranking.csv` — AFP scores for all peptides
- `summary_ranking.txt` — human-readable report with score breakdown
- `heatmap_scores.png` — visual comparison across all criteria
- `sasa_profiles.png` — per-residue exposure profiles

**Advancement criterion:** AFP score ≥ 6/10. The `ibs_residue_list` column from this output feeds directly into Stage 2.

---

### Stage 2 — Docking Analysis (`hex_docking_pipeline.py`)

Analyzes HEX rigid docking output (100 orientations per peptide against a TIP4P ice slab). HEX must be run manually before this script.

**Folder structure expected:**
```
results_docking/
├── peptide_A/
│   ├── peptide_A-1.pdb
│   ├── peptide_A-2.pdb
│   └── ...  (one PDB per orientation)
└── peptide_B/
    └── ...
```

**Configuration (edit in script):**
```r
HEX_DIR    <- "./results_docking"
OUTPUT_DIR <- "./hex_results"

IBS_RESIDUE_MAP <- list(
  peptide_A = c(2, 3, 5, 7, 8),   # residue numbers from Stage 1 output
  peptide_B = c(1, 4, 6, 9)
)
```

**Key output columns (`summary_all_peptides.csv`):**

| Column | Favorable value |
|--------|----------------|
| `best_hbond_contacts` | ≥ 20 (ideally ≥ 50) |
| `pass_rate_pct` | > 5% acceptable; > 15% good |
| `best_min_dist_A` | < 2.4 Å is normal for rigid docking — not a disqualifying criterion |

**Important:** After running the script, open the top model for each candidate in PyMOL and verify visually that the IBS face is oriented toward the ice slab. A peptide with high H-bond counts but whose IBS face points away from the ice must be discarded regardless of numeric scores.

**Advancement criterion:** IBS face confirmed toward ice in PyMOL + `best_hbond_contacts` ≥ 10 and `pass_rate_pct` > 1%.

---

### Stage 3 — Interaction Profiling (`plip_afp_pipeline.py`)

Runs after PLIP v3.0.0 generates interaction reports for each top model. Classifies every H-bond as IBS or non-IBS and writes a formatted Excel report.

**Run PLIP manually first:**
```bash
plip -f hex_results/peptide_A/top_models/peptide_A-07.pdb \
     -o ./plip_out/peptide_A --txt
```

**Configuration (edit in script):**
```python
PLIP_DIR = "./plip_out"
PEPTIDES = ["peptide_A", "peptide_B", "peptide_C"]

IBS_MAP = {
    "peptide_A": {("LEU", 2), ("ILE", 3), ("THR", 7), ("ALA", 8)},
    "peptide_B": {("GLY", 1), ("ALA", 4), ("VAL", 7)},
}
```

**Run:**
```bash
python plip_afp_pipeline.py
```

**Verdict interpretation:**

| Verdict | IBS ratio | Meaning |
|---------|-----------|---------|
| PASS — IBS dominant | ≥ 0.60 | H-bonds predominantly from IBS residues. Strong candidate. |
| BORDERLINE — advance | 0.40–0.59 | Real IBS signal but competed by non-IBS contacts. Advance with caution. |
| FAIL — non-IBS dominant | < 0.40 | Non-specific contacts dominate. Discard. |

**Excel output** (`./plip_results/PLIP_AFP_Summary.xlsx`):
- Sheet `Summary` — one row per peptide with all metrics and verdict
- Sheet `PLIP Detail` — every individual H-bond with geometry (distances, angles) and IBS classification
- Sheet `IBS Map` — reference table of defined IBS residues per peptide

**Advancement criterion:** Verdict PASS or BORDERLINE, IBS ratio ≥ 0.40, at least 2 IBS residues with confirmed contacts.

---

### AFP Pipeline Dependencies

```bash
pip install freesasa biopython pandas numpy matplotlib seaborn scipy openpyxl
# Also required: HEX 8.x (hex.loria.fr), PyMOL, PLIP v3.0.0 (conda-forge)
```

---

## Pipeline B — Activity Prediction (`peptide_activity_prediction_v7.R`)

### How it works

```
Labeled training data + unlabeled new peptides
(molecular descriptors)
        │
        ▼
Clean: remove duplicates, binary labels, imbalance check
        │
        ▼
Feature selection: AUC ≥ 0.60 per descriptor
+ correlation filter ≤ 0.90
        │
        ▼
Train/test split 80/20 (BEFORE feature selection — no leakage)
        │
        ▼
20 training runs (different seeds each):
  ┌─────────────────────────────┐
  │  A) Random Forest (RF)      │
  │  B) Neural Network + PCA    │  × 20 runs
  │  C) Neural Network + ANOVA  │
  └─────────────────────────────┘
        │
        ▼
Platt scaling calibration per model per run
        │
        ▼
AUC-weighted ensemble
        │
        ▼
Predictions with confidence tiers + stability metrics
        │
        ▼
Excel (10 sheets) + 3D PCA HTML plots
```

### Quick start

```r
# Install dependencies
install.packages(c("readxl","openxlsx","caret","ranger","nnet",
                   "pROC","dplyr","recipes","themis","plotly","htmlwidgets"))

# Load data
library(readxl)
train_data <- read_excel("train_data_capped_300.xlsx")
new_data   <- read_excel("new_peptides.xlsx")

# Configure and run
# Edit USER INPUTS section in the script, then:
source("peptide_activity_prediction_v7.R")
```

### Data format

**Training data** (`train_data`):

| ID | Type | Descriptor1 | Descriptor2 | ... |
|----|------|-------------|-------------|-----|
| PEP_001 | AFP | 1243.5 | 0.43 | ... |
| PEP_002 | FUNG | 987.2 | 0.61 | ... |

**New peptides** (`new_data`):

| ID | Descriptor1 | Descriptor2 | ... |
|----|-------------|-------------|-----|
| NEW_001 | 1105.3 | 0.55 | ... |

### Configuration

```r
target_class     <- "AFP"    # class to predict: AFP, FUNG, CANCER, AChE
comparison_class <- NULL     # NULL = one-vs-all; "FUNG" = binary AFP vs FUNG only

n_runs           <- 20       # training iterations (20 minimum; 50 for publication)
prob_threshold   <- 0.90     # threshold for frequency columns in Consensus sheet
auc_lower        <- 0.60     # minimum AUC to keep a descriptor
corr_threshold   <- 0.90     # correlation filter cutoff
```

### Preparing training data (recommended cap)

If any class has more than 300 samples, cap before running:

```r
library(dplyr); library(openxlsx); set.seed(42)
train_data_capped <- train_data %>%
  group_by(Type) %>%
  slice_sample(n = 300, replace = FALSE) %>%
  ungroup()
write.xlsx(train_data_capped, "train_data_capped_300.xlsx")
```

### Output files

| File | Description |
|------|-------------|
| `AFP_Complete_Analysis_v2_clean.xlsx` | Main Excel with 10 sheets (prefix matches `target_class`) |
| `AFP_PCA_3D_Training.html` | Interactive 3D PCA of training data only |
| `AFP_PCA_3D_Combined.html` | Interactive 3D PCA with new peptides overlaid on training data |

### Excel output — 10 sheets

| Sheet | Primary use |
|-------|-------------|
| `1_Summary` | Executive overview: parameters, model AUCs, prediction counts |
| `2_Simple_Predictions` | **Start here.** YES/NO per peptide with probability % and confidence tier |
| `3_Detailed_Predictions` | Full data: per-model probabilities, three threshold levels, model agreement |
| `4_High_Confidence` | **Priority candidates.** Filtered to probability > 70% and confidence High/Very High |
| `5_Consensus` | Per-run frequency analysis and stability (SD) per model per peptide |
| `6_Model_Evaluation` | Mean AUC, SD, min/max AUC across runs for each model and ensemble weights |
| `7_Feature_Importance` | RF permutation importance — which descriptors drive classification |
| `8_Top_Features_Stats` | Mean values of top 5 descriptors in target class vs others |
| `9_Diagnostics` | Dataset split sizes, class counts, features after filtering |
| `10_PCA_Coordinates` | PC1/PC2/PC3 for all training samples and new peptides |

### Confidence tiers

| Tier | Ensemble probability | Interpretation |
|------|---------------------|----------------|
| Very High | > 90% or < 10% | All three models strongly agree |
| High | > 70% or < 30% | Clear majority agreement |
| Medium | > 60% or < 40% | Tendency present but not strong |
| Low | 40%–60% | Models disagree — interpret with caution |

### Model AUC quality guide

| Mean AUC | Model quality |
|----------|--------------|
| ≥ 0.90 | Excellent |
| 0.80–0.89 | Good |
| 0.70–0.79 | Acceptable |
| 0.60–0.69 | Weak — predictions are indicative only |
| < 0.60 | Near random — review training data quality |

### R dependencies

| Package | Role |
|---------|------|
| caret | Training framework and cross-validation |
| ranger | Random Forest implementation |
| nnet | Neural network (single hidden layer) |
| pROC | AUC calculation and ROC curves |
| recipes + themis | Preprocessing pipeline and SMOTE |
| plotly + htmlwidgets | Interactive 3D PCA plots |
| readxl + openxlsx | Excel I/O |

---

## Documentation

User manuals are available in Spanish:

- 🇪🇸 [Manual del pipeline AFP (selección in silico)](Manuales_esp/Manual_Pipeline_AFP_seleccion.pdf)
- 🇪🇸 [Manual del pipeline de predicción de actividad (R)](Manuales_esp/Manual_Pipeline_Prediccion_Peptidos.pdf)

English versions coming soon.

---

## Requirements Summary

| Tool / Package | Version | Pipeline |
|----------------|---------|----------|
| Python | 3.9+ | A |
| R | 4.3+ | B |
| PLIP | 3.0.0 | A (Stage 3) |
| HEX | 8.x | A (Stage 2) |
| PyMOL | any | A (Stage 2 visual) |
| freesasa | any | A (Stage 1) |
| biopython | any | A (Stage 1) |
| caret | 6.0+ | B |
| ranger | 0.14+ | B |

---

## Author

**Javier Badilla**

## Acknowledgments

- Random Forest: [ranger](https://github.com/imbs-hl/ranger)
- ML framework: [caret](https://topepo.github.io/caret/)
- Visualizations: [plotly](https://plotly.com/r/)
- Interaction profiling: [PLIP](https://github.com/pharmai/plip)
- Surface analysis: [FreeSASA](https://freesasa.github.io/)

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{badilla2026peptide,
  author    = {Badilla, Javier},
  title     = {Peptide Activity Prediction — Computational Screening Pipeline},
  year      = {2026},
  url       = {https://github.com/Javier-Badilla/Activity-Probability-Prediction-}
}
```

**Version:** 2.0 | **Last updated:** March 2026 | **R:** 4.3+ | **Python:** 3.9+
