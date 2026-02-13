# -------------------------------
# 0) LOAD LIBRARIES
# -------------------------------
library(readr)
library(Matrix)
library(Seurat)
library(data.table)
library(ggplot2)
library(cluster)
library(dplyr)
library(SingleCellExperiment)
library(zellkonverter)
library(jsonlite)
# -------------------------------
set.seed(123)
# -------------------------------
# 0.1) LOAD GLOBAL CONFIGURATION & PATHS
# -------------------------------
# Try to determine script location to find Project Root
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg)
  script_dir <- dirname(script_path)
  # Assuming script is in src/clustering/, root is ../../
  project_root <- normalizePath(file.path(script_dir, "../../"))
} else {
  # Interactive mode fallback
  if (dir.exists("config")) {
    project_root <- getwd()
  } else if (dir.exists("../../config")) {
    project_root <- normalizePath("../../")
  } else {
    # Fallback to the hardcoded root from your gold standard if dynamic fails
    project_root <- "~/" 
  }
}

# Load JSON Config
config_file <- file.path(project_root, "config", "config.json")
if (!file.exists(config_file)) {
  # Fallback if running relative to project root manually
  config_file <- "config/config.json"
}
if (!file.exists(config_file)) stop("Config file not found: ", config_file)

global_config <- fromJSON(config_file)

# Resolve directories
datasets_base <- file.path(project_root, global_config$paths$datasets_dir)
results_base <- file.path(project_root, global_config$paths$results_dir)

# -------------------------------
# 0.2) SCRIPT PARAMETERS
# -------------------------------
# Map JSON params to the variables your script expects
param_config <- list(
  project_name = if (!is.null(global_config$params$project_name)) global_config$params$project_name else "Ecoli_eBW4",
  min_cells_qc = if (!is.null(global_config$params$min_cells_qc)) global_config$params$min_cells_qc else 50,
  min_features_qc = if (!is.null(global_config$params$min_features_qc)) global_config$params$min_features_qc else 100,
  genes_per_cell_qc = if (!is.null(global_config$params$genes_per_cell_qc)) global_config$params$genes_per_cell_qc else 10,
  n_features = if (!is.null(global_config$params$n_features)) global_config$params$n_features else 2000,
  pca_dims = if (!is.null(global_config$params$pca_dims)) global_config$params$pca_dims else 1:8,
  cluster_res = if (!is.null(global_config$params$cluster_res)) global_config$params$cluster_res else 0.6,
  min_cluster_size = if (!is.null(global_config$params$min_cluster_size)) global_config$params$min_cluster_size else 10,
  use_conditions = if (!is.null(global_config$params$use_conditions)) global_config$params$use_conditions else FALSE
)

# Parse Input Files from JSON
if (!is.null(global_config$clustering$input_files)) {
  input_raw <- global_config$clustering$input_files
  
  # Helper to resolve path
  resolve_path <- function(p) {
    if (is.null(p) || is.na(p)) return(NULL)
    if (startsWith(p, "/") || startsWith(p, "~")) return(p)
    return(file.path(datasets_base, p))
  }
  
  # Build file sets list
  param_config$file_sets <- lapply(seq_len(nrow(input_raw)), function(i) {
    row <- input_raw[i, ]
    expr_p <- resolve_path(row$expr_path)
    
    list(
      id = row$id,
      expr_path = expr_p,
      gene_path = resolve_path(row$gene_path),
      cell_index_path = resolve_path(row$cell_index_path),
      type = if (endsWith(expr_p, ".h5ad")) "h5ad" else "scifi"
    )
  })
} else {
  stop("No input files defined in config.json.")
}

# --- CREATE UNIQUE RUN NAME & OUTPUT DIRS ---
run_name <- sprintf(
  "%s_nFeat%d_Dims%d_Res%.1f",
  param_config$project_name,
  param_config$n_features,
  max(param_config$pca_dims),
  param_config$cluster_res
)

run_output_dir <- file.path(results_base, "Clustering", run_name)
data_output_dir <- file.path(run_output_dir, "data_files")
plot_output_dir <- file.path(run_output_dir, "plots_and_settings")
condition_output_dir <- file.path(data_output_dir, "grouped_by_condition")

dir.create(data_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)
if (param_config$use_conditions) dir.create(condition_output_dir, recursive = TRUE, showWarnings = FALSE)

message("--- Starting Run: ", run_name, " ---")

# -------------------------------
# 1) Read, Process, and Combine Variable Input Files
# -------------------------------
all_expr_list <- list()
all_gene_indices <- list()
all_metadata_list <- list()

message("Reading and processing input files...")

for (file_set in param_config$file_sets) {
  message("... processing ", file_set$id, " (Type: ", file_set$type, ")")
  
  expr_raw <- NULL
  gene_index <- NULL
  conditions <- NULL
  
  # ==========================
  # H5AD LOGIC (Robust)
  # ==========================
  if (file_set$type == "h5ad") {
    
    # 1. Read in "backed" mode (use_hdf5 = TRUE)
    # This prevents the "string to float" crash by delaying the data read
    message("    ... reading H5AD in backed mode")
    sce <- zellkonverter::readH5AD(file_set$expr_path, use_hdf5 = TRUE)
    
    # 2. Identify the numeric data layer
    assay_names <- assayNames(sce)
    if ("counts" %in% assay_names) {
      target_assay <- "counts"
    } else if ("raw" %in% assay_names) {
      target_assay <- "raw"
    } else {
      target_assay <- 1
    }
    
    # 3. Extract and Force Matrix
    message(sprintf("    ... extracting assay '%s'", target_assay))
    mat_delayed <- assay(sce, target_assay)
    
    # Check first value for string corruption
    first_val <- as.vector(head(mat_delayed, 1)[1])
    
    # Read full matrix into memory
    mat <- as.matrix(mat_delayed)
    
    # 4. Handle "String Matrix" Error (The "ValueError" source)
    if (is.character(first_val) || is.character(mat[1,1])) {
      message("    Warning: Matrix contains strings (likely gene IDs). Forcing to numeric...")
      r_names <- rownames(sce)
      c_names <- colnames(sce)
      
      # Force numeric (non-numbers become NA), then replace NA with 0
      mat <- matrix(suppressWarnings(as.numeric(mat)), 
                    nrow = length(r_names), ncol = length(c_names))
      rownames(mat) <- r_names
      colnames(mat) <- c_names
      mat[is.na(mat)] <- 0
    }
    
    # 5. Orientation Check
    # We want Cells as Rows, Genes as Columns (SCIFI format)
    if (nrow(mat) > ncol(mat)) {
      # This looks like Genes(rows) x Cells(cols) -> Transpose to Cells x Genes
      expr_raw <- as.data.frame(t(mat))
    } else {
      # This looks like Cells(rows) x Genes(cols) -> Keep as is
      expr_raw <- as.data.frame(mat)
    }
    
    # 6. Ensure Gene Names (Columns) and Cell IDs (Rows) are set
    # Logic: If columns match gene count (nrow(sce)), use rownames(sce) for cols
    if (ncol(expr_raw) == nrow(sce)) {
      colnames(expr_raw) <- rownames(sce)  # Genes
      rownames(expr_raw) <- colnames(sce)  # Cells
    } else {
      colnames(expr_raw) <- colnames(sce)
      rownames(expr_raw) <- rownames(sce)
    }
    
    # 7. Extract Metadata (Conditions & Clusters) BEFORE deleting SCE
    # -----------------------------------------------------------
    # Conditions
    if ("Condition" %in% names(colData(sce))) {
      conditions <- as.character(colData(sce)$Condition)
    } else if ("batch" %in% names(colData(sce))) {
      conditions <- as.character(colData(sce)$batch)
    } else {
      conditions <- rep(file_set$id, nrow(expr_raw))
    }
    
    # Clusters (Try common names: 'leiden', 'louvain', 'cluster')
    cd_names <- names(colData(sce))
    cluster_col <- cd_names[grep("leiden|louvain|cluster", cd_names, ignore.case=TRUE)][1]
    
    if (is.na(cluster_col)) {
      message("    Warning: No cluster column found. Saving as single 'Unclustered' file.")
      clusters <- rep("Unclustered", nrow(expr_raw))
    } else {
      clusters <- as.character(colData(sce)[[cluster_col]])
    }
    # -----------------------------------------------------------
    
    # Clean up heavy objects now that we have extracted everything
    rm(sce, mat, mat_delayed); gc()
    
    # 8. Save Cluster Files (Cells x Genes)
    # -----------------------------------------------------------
    message("    ... saving cluster expression files (CSV: Cells x Genes)")
    
    # Create output directory
    out_dir <- file.path(dirname(file_set$expr_path), "cluster_csvs")
    if (!dir.exists(out_dir)) dir.create(out_dir)
    
    unique_clusters <- unique(clusters)
    
    for (cls in unique_clusters) {
      # Subset rows (Cells) belonging to this cluster
      cluster_data <- expr_raw[clusters == cls, , drop = FALSE]
      
      # Clean cluster name for filename
      cls_clean <- gsub("[^A-Za-z0-9]", "_", cls) 
      fname <- file.path(out_dir, paste0(file_set$id, "_cluster_", cls_clean, ".csv"))
      
      # Write CSV: 
      # - quote = FALSE (no quotes around numbers)
      # - row.names = TRUE (preserves Cell IDs as the first column)
      # Result: Rows = Cells, Cols = Genes
      write.csv(cluster_data, file = fname, quote = FALSE, row.names = TRUE)
    }
    # -----------------------------------------------------------
    
    # Dummy gene index (for downstream compatibility with the rest of your script)
    gene_index <- data.frame(Geneid = colnames(expr_raw), Name = colnames(expr_raw))
    
  }else {
    # Explicit readr:: usage to avoid namespace errors
    expr_raw_tmp <- readr::read_csv(file_set$expr_path, col_types = readr::cols())
    expr_raw_tmp <- as.data.frame(expr_raw_tmp)
    
    # Handle index column
    if (is.character(expr_raw_tmp[, 1])) {
      expr_raw <- expr_raw_tmp[, -1, drop = FALSE]
      rownames(expr_raw) <- expr_raw_tmp[, 1]
    } else {
      expr_raw <- expr_raw_tmp
    }
    
    # Force numeric to prevent "string to float" error
    expr_raw <- as.data.frame(lapply(expr_raw, function(x) as.numeric(as.character(x))))
    expr_raw[is.na(expr_raw)] <- 0
    
    # Gene Index
    gene_index <- as.data.frame(readr::read_csv(file_set$gene_path, col_types = readr::cols()))
    
    # Map Columns
    id_col <- grep("id|Geneid", names(gene_index), ignore.case = TRUE, value = TRUE)[1]
    name_col <- grep("name|Symbol", names(gene_index), ignore.case = TRUE, value = TRUE)[1]
    
    if (!is.na(id_col) && !is.na(name_col)) {
      mapping_vector <- setNames(as.character(gene_index[[id_col]]), as.character(gene_index[[name_col]]))
      current_cols <- colnames(expr_raw)
      mapped_ids <- mapping_vector[current_cols]
      valid_map <- !is.na(mapped_ids)
      colnames(expr_raw)[valid_map] <- mapped_ids[valid_map]
      gene_index$Geneid <- gene_index[[id_col]] # Standardize for Section 2
    }
    
    # Metadata / Conditions
    conditions <- rep("Unknown", nrow(expr_raw))
    if (!is.null(file_set$cell_index_path) && file.exists(file_set$cell_index_path)) {
      # explicit readr:: call with col_names = FALSE
      meta_raw <- as.data.frame(readr::read_csv(file_set$cell_index_path, col_names = FALSE, col_types = readr::cols()))
      
      # Find Lib column
      find_treat_col <- function(df) {
        for(i in seq_along(df)) {
          if(any(grepl("^Lib", as.character(df[[i]][1:min(10, nrow(df))])))) return(i)
        }
        return(NULL)
      }
      treat_col_idx <- find_treat_col(meta_raw)
      
      if (!is.null(treat_col_idx)) {
        raw_strings <- as.character(meta_raw[[treat_col_idx]])
        parsed_conditions <- sub(".*_([^_]+)$", "\\1", raw_strings)
        if(length(parsed_conditions) >= nrow(expr_raw)) {
          conditions <- parsed_conditions[1:nrow(expr_raw)]
        } else {
          conditions[1:length(parsed_conditions)] <- parsed_conditions
        }
      }
    }
  }
  
  # Standardize Output for Combine
  cell_ids <- paste0(file_set$id, "_cell_", seq_len(nrow(expr_raw)))
  rownames(expr_raw) <- cell_ids
  
  all_metadata_list[[file_set$id]] <- data.frame(
    cell_barcode = cell_ids,
    dataset_id = file_set$id,
    condition = conditions,
    stringsAsFactors = FALSE
  )
  rownames(all_metadata_list[[file_set$id]]) <- cell_ids
  
  all_expr_list[[file_set$id]] <- expr_raw
  all_gene_indices[[file_set$id]] <- gene_index
}
gc()

# -------------------------------
# 2) Find Common Gene Index & Filter (Safe Version)
# -------------------------------
message("Finding common genes across all datasets...")
all_gene_ids_list <- lapply(all_gene_indices, function(df) as.character(df$Geneid))
common_gene_ids <- Reduce(intersect, all_gene_ids_list)

if (length(common_gene_ids) == 0) stop("No common gene IDs found.")
common_genes_sorted <- sort(unique(common_gene_ids))
message(sprintf("Found %d common genes.", length(common_genes_sorted)))

filtered_expr_list <- lapply(names(all_expr_list), function(id) {
  expr_df <- all_expr_list[[id]]
  # Safe subset: Only select columns that exist
  present_genes <- intersect(common_genes_sorted, colnames(expr_df))
  expr_filtered <- expr_df[, present_genes, drop = FALSE]
  
  # Zero-fill missing common genes if any
  missing_genes <- setdiff(common_genes_sorted, present_genes)
  if (length(missing_genes) > 0) {
    message(sprintf("  ... zero-filling %d genes in %s", length(missing_genes), id))
    for(g in missing_genes) expr_filtered[[g]] <- 0
  }
  
  # Final Sort
  return(expr_filtered[, common_genes_sorted, drop = FALSE])
})
names(filtered_expr_list) <- names(all_expr_list)
rm(all_expr_list, all_gene_indices); gc()

# -------------------------------
# 3) Combine into Final Matrix
# -------------------------------
message("Combining all datasets into one matrix...")
combined_df <- do.call(rbind, filtered_expr_list)
combined_meta <- do.call(rbind, all_metadata_list)
rm(filtered_expr_list, all_metadata_list); gc()

combined_sparse <- Matrix(as.matrix(combined_df), sparse = TRUE)
combined_sparse <- t(combined_sparse) # Genes x Cells
rm(combined_df); gc()

# =========================================================================
# 3b) OPTIONAL: ROBUST GENE MAPPING (BiGG Dictionary)
# =========================================================================
# If you want the advanced mapping from the other script, keep this block.
# Otherwise, skip to Section 4.
message("Fetching iML1515 dictionary from BiGG Models...")
url <- "http://bigg.ucsd.edu/api/v2/models/iML1515/genes"
gene_map_data <- tryCatch({
  json_response <- jsonlite::fromJSON(url)
  if ("results" %in% names(json_response)) json_response$results else NULL
}, error = function(e) return(NULL))

if (!is.null(gene_map_data) && all(c("name", "bigg_id") %in% colnames(gene_map_data))) {
  map_df <- gene_map_data[, c("name", "bigg_id")]
  gene_lookup <- setNames(map_df$bigg_id, tolower(map_df$name))
  
  current_genes <- rownames(combined_sparse)
  new_genes <- current_genes 
  mapped_count <- 0
  
  for (i in seq_along(current_genes)) {
    tokens <- unlist(strsplit(current_genes[i], "[-_ ]"))
    mapped_tokens <- sapply(tokens, function(t) {
      t_low <- tolower(t)
      if (t_low %in% names(gene_lookup)) return(gene_lookup[[t_low]]) else return(t)
    })
    if (any(mapped_tokens != tokens)) {
      new_genes[i] <- paste(mapped_tokens, collapse = "-")
      mapped_count <- mapped_count + 1
    }
  }
  rownames(combined_sparse) <- new_genes
  message(sprintf("Mapped %d genes to iML1515 IDs.", mapped_count))
}

# -------------------------------
# 4) Create Seurat Object & QC
# -------------------------------
message("Creating Seurat object...")
seurat_obj <- CreateSeuratObject(
  counts = combined_sparse,
  project = param_config$project_name,
  min.cells = param_config$min_cells_qc,
  min.features = param_config$min_features_qc,
  meta.data = combined_meta
)

raw_counts <- GetAssayData(seurat_obj, layer = "counts")
genes_per_cell <- Matrix::colSums(raw_counts > 0)
keep_cells <- names(genes_per_cell)[genes_per_cell >= param_config$genes_per_cell_qc]
seurat_obj <- subset(seurat_obj, cells = keep_cells)
message(sprintf("%d cells remain after filtering.", length(keep_cells)))

# -------------------------------
# 5) Normalize, Scale, PCA
# -------------------------------
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 1e6, verbose = FALSE)
mat <- GetAssayData(seurat_obj, layer = "data")
mat <- mat / log(2) # Convert ln to log2
seurat_obj <- SetAssayData(seurat_obj, layer = "data", new.data = mat)
rm(mat); gc()

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = param_config$n_features)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)

# -------------------------------
# 6) Elbow Plot
# -------------------------------
elbow_plot <- ElbowPlot(seurat_obj, ndims = 50)
ggsave(file.path(plot_output_dir, "1_ElbowPlot.png"), elbow_plot, width = 8, height = 6)

# -------------------------------
# 7) Clustering
# -------------------------------
seurat_obj <- FindNeighbors(seurat_obj, dims = param_config$pca_dims)
seurat_obj <- FindClusters(seurat_obj, resolution = param_config$cluster_res)

cluster_counts <- table(Idents(seurat_obj))
clusters_to_keep <- names(cluster_counts)[cluster_counts >= param_config$min_cluster_size]
if (length(clusters_to_keep) == 0) stop("No clusters met size threshold.")
seurat_obj <- subset(seurat_obj, idents = clusters_to_keep)

# -------------------------------
# 8) UMAP & Plots
# -------------------------------
seurat_obj <- RunUMAP(seurat_obj, dims = param_config$pca_dims, verbose = FALSE)

p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle("Cluster")
p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "condition") + ggtitle("Condition")
ggsave(file.path(plot_output_dir, "2_UMAP_Combined.png"), p1 + p2, width = 14, height = 7)

df_meta <- seurat_obj@meta.data
df_summary <- df_meta %>%
  group_by(seurat_clusters, condition) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(seurat_clusters) %>%
  mutate(prop = count / sum(count))

p3 <- ggplot(df_summary, aes(x = seurat_clusters, y = prop, fill = condition)) +
  geom_bar(stat = "identity") + theme_minimal() + ggtitle("Cluster Composition")
ggsave(file.path(plot_output_dir, "3_Cluster_Composition.png"), p3, width = 8, height = 6)

# -------------------------------
# 9) Data Export
# -------------------------------
cell_clusters <- data.frame(
  cell_barcode = names(Idents(seurat_obj)),
  cluster = as.character(Idents(seurat_obj)),
  condition = seurat_obj$condition,
  stringsAsFactors = FALSE
)

norm_data <- GetAssayData(seurat_obj, layer = "data")
norm_data <- norm_data[, cell_clusters$cell_barcode, drop = FALSE]

data.table::fwrite(data.table(gene = rownames(norm_data)), file.path(data_output_dir, "genes_logCPM.csv"))
data.table::fwrite(cell_clusters, file.path(data_output_dir, "cell_metadata.csv"))

for (cl in unique(cell_clusters$cluster)) {
  cols <- cell_clusters$cell_barcode[cell_clusters$cluster == cl]
  sub <- norm_data[, cols, drop = FALSE]
  df <- data.table::data.table(
    cell_id = cols, cluster = cl, 
    condition = cell_clusters$condition[cell_clusters$cluster == cl],
    as.data.frame(t(as.matrix(sub)), check.names = FALSE)
  )
  data.table::fwrite(df, file.path(data_output_dir, sprintf("cluster_%s_logCPM.csv", cl)), sep = ",")
}

# -------------------------------
# 10) CONDITION EXPORT (New Feature)
# -------------------------------
if (param_config$use_conditions) {
  for (cond in unique(cell_clusters$condition)) {
    safe_cond <- gsub("[^A-Za-z0-9]", "_", cond)
    cols <- cell_clusters$cell_barcode[cell_clusters$condition == cond]
    if (length(cols) >= 50) {
      sub <- norm_data[, cols, drop = FALSE]
      df <- data.table::data.table(
        cell_id = cols, cluster = safe_cond, 
        as.data.frame(t(as.matrix(sub)), check.names = FALSE)
      )
      data.table::fwrite(df, file.path(condition_output_dir, sprintf("cluster_%s_logCPM.csv", safe_cond)), sep = ",")
      message("Wrote Condition Group: ", safe_cond)
    }
  }
}

# -------------------------------
# 11) DETAILED QUALITY METRICS (From Gold Standard)
# -------------------------------
message("Computing clustering quality metrics...")

compute_clustering_quality <- function(seurat_obj, pca_dims, cluster_labels, run_name) {
  pca_embedding <- Embeddings(seurat_obj, "pca")[, pca_dims]
  clusters <- as.numeric(as.factor(cluster_labels))
  cluster_ids <- unique(clusters)
  n_clusters <- length(cluster_ids)
  per_cluster_metrics <- data.frame()
  
  message(" Computing per-cluster quality metrics...")
  for (cluster_id in cluster_ids) {
    cluster_points <- pca_embedding[clusters == cluster_id, , drop = FALSE]
    cluster_size <- nrow(cluster_points)
    
    if (cluster_size > 1) {
      centroid <- colMeans(cluster_points)
      dists <- sqrt(rowSums((cluster_points - matrix(centroid, nrow=cluster_size, ncol=length(centroid), byrow=TRUE))^2))
      
      avg_dist <- mean(dists)
      cv_dist <- ifelse(avg_dist > 0, sd(dists) / avg_dist, 0)
      wss <- sum(dists^2)
      
      # Density
      if (cluster_size > 2) {
        # Simple nearest neighbor approximation to avoid huge matrix calc on large data
        if(cluster_size < 2000) {
          cdist <- as.matrix(dist(cluster_points)); diag(cdist) <- Inf
          avg_nn <- mean(apply(cdist, 1, min))
        } else { avg_nn <- NA } 
      } else { avg_nn <- NA }
    } else {
      avg_dist <- 0; cv_dist <- 0; wss <- 0; avg_nn <- NA
    }
    
    per_cluster_metrics <- rbind(per_cluster_metrics, data.frame(
      ClusterID = cluster_id, ClusterSize = cluster_size,
      AvgDistanceToCentroid = avg_dist, DistanceCV = cv_dist, WithinClusterSS = wss,
      AvgNearestNeighborDist = avg_nn, AvgSilhouetteWidth = NA
    ))
  }
  
  # Overall
  wss <- sum(per_cluster_metrics$WithinClusterSS, na.rm = TRUE)
  overall_centroid <- colMeans(pca_embedding)
  bss <- 0
  for (i in 1:nrow(per_cluster_metrics)) {
    cid <- per_cluster_metrics$ClusterID[i]
    pts <- pca_embedding[clusters == cid, , drop = FALSE]
    if (nrow(pts) > 0) bss <- bss + nrow(pts) * sum((colMeans(pts) - overall_centroid)^2)
  }
  
  # Silhouette
  sil_score <- NA
  tryCatch({
    if(nrow(pca_embedding) > 5000) {
      idx <- sample(1:nrow(pca_embedding), 5000)
      sil <- silhouette(clusters[idx], dist(pca_embedding[idx,]))
    } else {
      sil <- silhouette(clusters, dist(pca_embedding))
    }
    sil_score <- summary(sil)$avg.width
  }, error = function(e) message("Sil failed: ", e$message))
  
  return(list(
    overall_metrics = data.frame(
      Config = run_name, Silhouette = sil_score, 
      WSS = wss, BSS = bss, Ratio = bss/wss,
      AvgClusterSize = mean(per_cluster_metrics$ClusterSize)
    ),
    per_cluster_metrics = per_cluster_metrics
  ))
}

quality_results <- compute_clustering_quality(seurat_obj, param_config$pca_dims, Idents(seurat_obj), run_name)

write.csv(quality_results$overall_metrics, file.path(data_output_dir, "clustering_quality_metrics.csv"), row.names = FALSE)
write.csv(quality_results$per_cluster_metrics, file.path(data_output_dir, "per_cluster_quality_metrics.csv"), row.names = FALSE)

# -------------------------------
# 12) Save Settings
# -------------------------------
settings_text <- paste(
  "--- Settings for run:", run_name, "---",
  "\nDate: ", Sys.time(),
  "\nInput Files: ", paste(sapply(param_config$file_sets, `[[`, "id"), collapse=", "),
  "\nMin Cells: ", param_config$min_cells_qc,
  "\nPCA Dims: ", paste(param_config$pca_dims, collapse=","),
  "\nSilhouette: ", quality_results$overall_metrics$Silhouette
)
writeLines(settings_text, file.path(plot_output_dir, "run_settings.txt"))

message("--- Run Complete: ", run_name, " ---")