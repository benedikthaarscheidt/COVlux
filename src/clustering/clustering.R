# ==============================================================================
# SCRIPT 2: CLUSTERING & ANALYSIS (With QC Plots)
# ==============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(cluster)

# -------------------------------
# 1) SOURCE DATA PREP
# -------------------------------
# This loads 'final_counts', 'final_metadata', and 'config'
source("src/clustering/ETL.R") 

message("--- Starting Clustering Workflow ---")
message("Project: ", config$project_name)

# -------------------------------
# 2) SETUP OUTPUTS
# -------------------------------
run_name <- sprintf("%s_nFeat%d_Dims%d_Res%.1f", config$project_name, config$n_features, max(config$pca_dims), config$cluster_res)
run_dir <- file.path(config$results_base, "Clustering", run_name)
data_dir <- file.path(run_dir, "data_files")
plot_dir <- file.path(run_dir, "plots_and_settings")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

if (config$use_conditions) {
  cond_dir <- file.path(data_dir, "grouped_by_condition")
  dir.create(cond_dir, recursive = TRUE, showWarnings = FALSE)
}

# -------------------------------
# 3) INITIALIZE SEURAT & QC PLOTS
# -------------------------------
message("Initializing Seurat...")

seurat_obj <- CreateSeuratObject(
  counts = final_counts,
  project = config$project_name,
  min.cells = config$min_cells_qc,
  min.features = config$min_features_qc,
  meta.data = final_metadata
)

# --- NEW: QC VIOLIN PLOTS (Before Filtering) ---
message("Generating QC Violin Plots...")
p_vln <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0.1) 
ggsave(file.path(plot_dir, "0_QC_ViolinPlots.png"), p_vln, width = 10, height = 6)

# --- APPLY FILTER ---
raw_counts <- GetAssayData(seurat_obj, layer = "counts")
genes_per_cell <- Matrix::colSums(raw_counts > 0)
cells_to_keep <- names(genes_per_cell)[genes_per_cell >= config$genes_per_cell_qc]
seurat_obj <- subset(seurat_obj, cells = cells_to_keep)

message(sprintf("Cells remaining after QC: %d", ncol(seurat_obj)))

# -------------------------------
# 4) NORMALIZATION & PCA
# -------------------------------
message("Normalizing & PCA...")

# Normalize (Natural Log)
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 1e6, verbose = FALSE)

# Convert to Log2 (Manual Adjustment)
mat <- GetAssayData(seurat_obj, layer = "data") / log(2)
seurat_obj <- SetAssayData(seurat_obj, layer = "data", new.data = mat)

# Find Variable Features & Scale
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = config$n_features)
seurat_obj <- ScaleData(seurat_obj)

# Run PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)

# --- NEW: ELBOW PLOT ---
message("Generating Elbow Plot...")
p_elbow <- ElbowPlot(seurat_obj, ndims = 50)
ggsave(file.path(plot_dir, "1_ElbowPlot.png"), p_elbow, width = 8, height = 6)

# -------------------------------
# 5) CLUSTERING & UMAP
# -------------------------------
message("Clustering...")
seurat_obj <- FindNeighbors(seurat_obj, dims = config$pca_dims)
seurat_obj <- FindClusters(seurat_obj, resolution = config$cluster_res)

# Filter Small Clusters
counts <- table(Idents(seurat_obj))
keep <- names(counts)[counts >= config$min_cluster_size]
if (length(keep) == 0) stop("No clusters met size threshold.")
seurat_obj <- subset(seurat_obj, idents = keep)

# Run UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = config$pca_dims, verbose = FALSE)

# -------------------------------
# 6) PLOTTING RESULTS
# -------------------------------
message("Generating UMAP plots...")

p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle("Clusters")

if (config$use_conditions) {
  p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "condition") + ggtitle("Conditions")
  ggsave(file.path(plot_dir, "2_UMAP_Combined.png"), p1 + p2, width = 14, height = 7)
} else {
  ggsave(file.path(plot_dir, "2_UMAP_Clusters.png"), p1, width = 8, height = 7)
}

# -------------------------------
# 7) EXPORT DATA
# -------------------------------
message("Exporting Data...")
norm_data <- GetAssayData(seurat_obj, layer = "data")
meta_df <- data.frame(
  cell_barcode = colnames(seurat_obj),
  cluster = Idents(seurat_obj),
  condition = seurat_obj$condition, 
  stringsAsFactors = FALSE
)

# A. Generic
fwrite(data.table(gene = rownames(norm_data)), file.path(data_dir, "genes_logCPM.csv"))
fwrite(meta_df, file.path(data_dir, "cell_metadata.csv"))

# B. Per-Cluster
for (cl in unique(meta_df$cluster)) {
  cells <- meta_df$cell_barcode[meta_df$cluster == cl]
  sub_mat <- as.matrix(norm_data[, cells, drop=FALSE])
  
  out_df <- data.table(
    cell_id = cells,
    cluster = cl,
    condition = meta_df$condition[meta_df$cell_barcode %in% cells],
    t(sub_mat)
  )
  fwrite(out_df, file.path(data_dir, sprintf("cluster_%s_logCPM.csv", cl)))
}

# C. Per-Condition (Only if TRUE)
if (config$use_conditions) {
  unique_conds <- unique(na.omit(meta_df$condition))
  for (cond in unique_conds) {
    safe_cond <- gsub("[^A-Za-z0-9]", "_", cond)
    cells <- meta_df$cell_barcode[meta_df$condition == cond]
    
    if(length(cells) > 0) {
      sub_mat <- as.matrix(norm_data[, cells, drop=FALSE])
      out_df <- data.table(
        cell_id = cells,
        cluster = meta_df$cluster[meta_df$cell_barcode %in% cells],
        t(sub_mat)
      )
      fwrite(out_df, file.path(cond_dir, sprintf("Cluster_%s_logCPM.csv", safe_cond)))
    }
  }
}

# -------------------------------
# 8) SAVE SETTINGS
# -------------------------------
settings_text <- paste(
  "--- Settings for run:", run_name, "---",
  "\nProject Name: ", config$project_name,
  "\nDate: ", Sys.time(),
  "\nInput Datasets Processed: ", paste(sapply(config$file_sets, `[[`, "id"), collapse=", "),
  "\nUse Conditions: ", config$use_conditions,
  "\n\n-- PARAMETERS --",
  "\nMin Cells: ", config$min_cells_qc,
  "\nMin Features: ", config$min_features_qc,
  "\nResolution: ", config$cluster_res,
  "\nPCA Dims: ", paste(config$pca_dims, collapse=",")
)

writeLines(settings_text, file.path(plot_dir, "run_settings.txt"))
message("--- Run Complete: ", run_name, " ---")