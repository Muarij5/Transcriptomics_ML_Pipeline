# Transcriptomics ML Pipeline — Epilepsy Multi-Cohort Study

Personal contribution to an ongoing epilepsy transcriptomics research project.
This repository covers end-to-end analysis — data integration, batch correction,
differential expression, and machine learning-based gene ranking, across
7 public RNA-seq cohorts. Part of an unpublished study currently in the
manuscript-writing phase.

## Overview

This pipeline integrates gene expression data from multiple independent
epilepsy studies (mTLE, FCD, TSC vs. Control), corrects for batch effects
across cohorts, identifies differentially expressed genes, and uses a
nested cross-validation machine learning framework to classify samples and
rank the most predictive genes with explainability analysis.

## Datasets

RNA-seq raw counts and metadata integrated from 7 GEO series:
GSE186334, GSE94744, GSE134697, GSE213488, GSE252323 (LEAT), GSE256068, GSE310367

## Pipeline

**1. Metadata curation** (`code/Metada_f2.R`)
Retrieved sample metadata for each GEO series via `GEOquery`, standardized
tissue/condition/platform labels across cohorts, matched sample IDs against
the expression matrix, and merged into one unified metadata table.

**2. Gene annotation & cleaning** (`code/annotate_new.R`)
Mapped Entrez IDs to gene symbols (pinned to Ensembl v112 for reproducibility),
removed pseudogenes, resolved duplicate gene symbols by median expression,
merged all cohorts on their common gene set, and applied QC filtering
(library size and zero-expression thresholds, CPM-based low-expression filter).

**3. Batch correction & differential expression** (`code/Normalization_next.R`)
- Grouped conditions into broad categories (Control, mTLE, FCD, TSC)
- Quantified batch/tissue/platform effects (η²) before and after correction
- Applied ComBat-seq to correct for cross-dataset batch effects on raw counts
- Ran DESeq2 for differential expression (mTLE, FCD, TSC vs. Control)
- Ran a threshold ablation study (log2FC × padj grid) to test result sensitivity
- Generated PCA and volcano plots.
- Identified DEG overlaps across conditions (pairwise + 3-way, with Venn diagram)

**4. ML matrix construction** (`code/ML_matrix_code.R`)
Built a machine learning-ready matrix using the union of DEGs across all
three comparisons, with binary condition labels (Control vs. Disease) for
downstream classification.

**5. ML classification & gene ranking** (`notebook/notebooka35854ddb9.ipynb`)
Built a nested cross-validation pipeline (5 outer folds × 3 inner folds) to
classify Control vs. Disease samples and rank the most predictive genes:
- Two ensemble models: Random Forest and Extra Trees
- Hyperparameter tuning via Optuna (200 trials/fold/model), with an
  overfitting penalty built into the objective function
- SMOTE class balancing and feature filtering applied strictly within inner
  training folds only, to prevent test-fold data leakage
- Gene importance computed via three independent methods — SHAP values,
  permutation importance, and native feature importance — combined into a
  per-model consensus score
- Final gene panel selected via RF ∩ ET intersection (genes ranked highly
  by both models), reducing single-model bias
- Reported accuracy, F1, precision, recall, ROC-AUC, and PR-AUC with
  train/test gap monitoring to flag overfitting per fold

## Repository structure
```
├── data/ raw counts, cleaned expression matrices, sample metadata
├── code/ R scripts — metadata curation, annotation, normalization/DEG analysis, ML matrix construction
├── notebook/ ML classification and gene ranking notebook
├── results/ DEG lists, ablation study, final expression/ML matrices, figures
└── README.md
```
## Status

Full pipeline complete: data integration → batch correction → differential
expression → nested-CV ML classification with consensus gene ranking.
Manuscript currently in preparation.
