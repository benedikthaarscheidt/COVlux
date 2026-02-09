# COVlux-Hybrid: Second Moment Metabolic Analysis Pipeline

This repository contains the code for the Master's Thesis project on **COVLux**. It implements a pipeline integrating single-cell RNA-seq (scRNA-seq) data with Elementary Flux Modes (EFM) to reconstruct metabolic phenotypic plasticity.

## Overview

The pipeline bridges the gap between transcriptomic variability (scRNA-seq) and metabolic flux capabilities (EFMs). It consists of three main stages:

1.  **Clustering (R)**: Groups heterogeneous single-cell data into metabolic phenotypes.
2.  **MRAS Mapping (MATLAB)**: Maps gene expression to reaction activity using GPR rules, computing **Mean** and **Covariance** of reaction activities.
3.  **COVlux Optimization (MATLAB)**: Uses a symmetric, sparsity-inducing optimization to select the minimal set of EFMs that best explain the observed reaction covariances.

## Directory Structure

*   `config/`: Configuration files (controls inputs, parameters, and paths).
*   `data/`: Input storage (Metabolic models, EFM bases, scRNA-seq datasets).
*   `src/`: Core source code.
    *   `clustering/`: R scripts for Seurat-based clustering and data prep.
    *   `mras/`: MATLAB scripts for Gene-to-Reaction mapping.
    *   `optimization/`: COVlux-Hybrid optimization algorithms.
    *   `efm/`: EFM enumeration utilities.
*   `results/`: Pipeline outputs (organized by run name and timestamp).

## Prerequisites

*   **MATLAB** (R2023b or later Recommended)
    *   Optimization Toolbox
*   **R** (4.0 or later)
    *   `Seurat`
    *   `jsonlite`
    *   `Matrix`

## Configuration

The entire pipeline is controlled by `config/config.json`. You do not need to edit the scripts directly.

**Key Fields in `config.json`:**
*   `clustering`: Input file paths (supports `.h5ad`, `.csv`, or lists of files).
*   `model`: Path to the metabolic model (e.g., `iML1515`) and EFM basis files.
*   `params`: Analysis parameters (e.g., `lambda_l21`, `max_iters`, `use_conditions`).
*   `paths`: Directory overrides (optional).

## Usage Guide

### Step 1: Clustering (R)
Processes raw scRNA-seq data and exports expression matrices for each cluster.

*   **Script**: `src/clustering/scRNA_clustering_combined.R`
*   **Input**: Defined in `config.json` (`clustering.input_files`).
*   **Output**: `results/Clustering/<RunName>/data_files/`

### Step 2: MRAS Mapping (MATLAB)
Converts gene expression to reaction activity scores (MRAS) using the metabolic model's GPR rules. It calculates the **Second Moment** (Covariance) of reaction activities for each cluster.

*   **Script**: `src/mras/MRAS_mapping.m`
*   **Input**: Reads automatically from the clustering output defined in `config.json`.
*   **Output**: `MRAS_outputs/*.csv` (Mean and Covariance matrices).

### Step 3: COVlux Optimization (MATLAB)
Reconstructs the metabolic state by selecting EFMs that match the observed reaction covariance.

*   **Script**: `src/optimization/covlux.m`
*   **Input**: Loads `MRAS_outputs` and EFM basis files defined in the config.
*   **Output**: `results/.../COVlux_cov_smallbasis/Run_YYYY-MM-DD.../`
    *   `*_A.csv`: The optimized covariance matrix of EFMs.
    *   `*_metrics.csv`: Quality metrics (Reconstruction error, Sparsity).
    *   `plots/`: Visualization of convergence, EFM length, and sparsity.

## Troubleshooting

*   **"Config file not found"**: Ensure you are running scripts from the project root or that `config.json` is correctly placed in the `config/` folder.
*   **NaN Errors**: The MRAS script includes checks for NaN values. If found, check your scRNA-seq input data for missing values or zero-variance genes.
