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
# 0.1) SMART PATH HANDLING & CONFIGURATION
# -------------------------------
# Determine script location for relative paths
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg)
  project_root <- normalizePath(file.path(dirname(script_path), "../../"))
} else {
  # Fallback: check for config directory
  if (dir.exists("config")) {
    project_root <- getwd()
  } else {
    project_root <- "~/"
  }
}

# Try to load config file
config_file <- file.path(project_root, "config", "config.json")
if (!file.exists(config_file)) {
  config_file <- "config/config.json"
}

if (file.exists(config_file)) {
  # Use JSON config if available
  message("Loading configuration from: ", config_file)
  global_config <- fromJSON(config_file)
  
  # Set paths
  datasets_base <- file.path(project_root, global_config$paths$datasets_dir)
  results_base <- file.path(project_root, global_config$paths$results_dir)
  
  # Load parameters
  config <- list(
    project_name = global_config$params$project_name,
    min_cells_qc = global_config$params$min_cells_qc,
    min_features_qc = global_config$params$min_features_qc,
    genes_per_cell_qc = global_config$params$genes_per_cell_qc,
    n_features = global_config$params$n_features,
    pca_dims = global_config$params$pca_dims,
    cluster_res = global_config$params$cluster_res,
    min_cluster_size = global_config$params$min_cluster_size,
    use_conditions = global_config$params$use_conditions
  )
  
  # Resolve input file paths
  resolve_path <- function(p) {
    if (is.null(p) || is.na(p) || p == "") return(NULL)
    if (startsWith(p, "/") || startsWith(p, "~")) return(p)
    return(file.path(datasets_base, p))
  }
  
  # Process input files from config
  if (!is.null(global_config$clustering$input_files)) {
    input_raw <- global_config$clustering$input_files
    config$file_sets <- lapply(seq_len(nrow(input_raw)), function(i) {
      row <- input_raw[i, ]
      
      # Determine file type
      expr_path <- resolve_path(row$expr_path)
      type <- if (endsWith(tolower(expr_path), ".h5ad")) "h5ad" else "scifi"
      
      # Handle condition file - only include if use_conditions is TRUE and file exists
      cell_index_path <- NULL
      if (config$use_conditions && !is.null(row$cell_index_path)) {
        cell_index_path <- resolve_path(row$cell_index_path)
        if (!file.exists(cell_index_path)) {
          warning("Condition file not found: ", cell_index_path, ". Proceeding without conditions.")
          cell_index_path <- NULL
        }
      }
      
      list(
        id = row$id,
        type = type,
        expr_path = expr_path,
        gene_path = resolve_path(row$gene_path),
        cell_index_path = cell_index_path,
        condition_label = if (!is.null(row$condition_label)) row$condition_label else NULL
      )
    })
  } else {
    stop("No input files specified in config.")
  }
  
} else {
  
  message("Config file not found. Using hardcoded configuration.")
  
  # Set default paths relative to project root
  datasets_base <- file.path(project_root, "datasets")
  results_base <- file.path(project_root, "results")
  
  config <- list(
    project_name = "Ecoli_LB_1057",
    min_cells_qc = 50,
    min_features_qc = 100,
    genes_per_cell_qc = 10,
    n_features = 1000,
    pca_dims = 1:8,
    cluster_res = 0.3,
    min_cluster_size = 10,
    use_conditions = TRUE  # Default to using conditions if available
  )
  
  # Define input file sets 
  config$file_sets <- list(
    list(
      id = "LB_1057",
      type = "h5ad",
      expr_path = file.path(datasets_base, "E_coli/ecoli_raw_counts_LB/raw_count_matrix_LB_1057_operons.h5ad"),
      gene_path = NULL,  # Not needed for h5ad
      cell_index_path = NULL,  # Conditions from h5ad metadata
      condition_label = "LB_Growth"
    )
  )
}

# -------------------------------
# 0.2) CREATE OUTPUT DIRECTORIES
# -------------------------------
run_name <- sprintf(
  "%s_nFeat%d_Dims%d_Res%.1f",
  config$project_name, config$n_features, max(config$pca_dims), config$cluster_res
)

base_output_dir <- file.path(results_base, "Clustering")
run_output_dir <- file.path(base_output_dir, run_name)
data_output_dir <- file.path(run_output_dir, "data_files")
plot_output_dir <- file.path(run_output_dir, "plots_and_settings")
condition_output_dir <- file.path(data_output_dir, "grouped_by_condition")

dir.create(data_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)
if (config$use_conditions) {
  dir.create(condition_output_dir, recursive = TRUE, showWarnings = FALSE)
}

message("--- Starting Run: ", run_name, " ---")
message("Project Root: ", project_root)
message("Output Directory: ", run_output_dir)

# -------------------------------
# 1) Read, Process, and Combine Variable Input Files
# -------------------------------
all_expr_list <- list()
all_gene_indices <- list()
all_metadata_list <- list()

message("Reading and processing input files...")
for (file_set in config$file_sets) {
  message("... processing ", file_set$id, " (Type: ", file_set$type, ")")
  cat("\nProcessing:", file_set$id, "\n")
  
  cat("File path:", file_set$expr_path, "\n") 
  cat("File exists?", file.exists(file_set$expr_path), "\n")
  cat("File path class:", class(file_set$expr_path), "\n")
  
  
  tryCatch({
    
    data <- vroom::vroom(file_set$expr_path) 
  }, error = function(e) {
    cat("Error for", file_set$id, ":", e$message, "\n")
   
    cat("Full file path:", normalizePath(file_set$expr_path, mustWork = FALSE), "\n")
  })
  # Initialize variables
  expr_raw <- NULL
  conditions <- NULL
  
  # =========================================================================
  # LOGIC BRANCH: H5AD FILES
  # =========================================================================
  if (file_set$type == "h5ad") {
    
    message("    Reading H5AD file: ", basename(file_set$expr_path))
    sce <- zellkonverter::readH5AD(file_set$expr_path)
    
    # 1. INSPECT LAYERS
    available_assays <- assayNames(sce)
    message("    -> Found layers: ", paste(available_assays, collapse=", "))
    
    # 2. SELECT 'counts' OR 'raw' (Avoid default 'X' if possible)
    if ("counts" %in% available_assays) {
      message("    -> EXPLICITLY LOADING 'counts' layer.")
      mat <- assay(sce, "counts")
    } else if ("raw" %in% available_assays) {
      message("    -> EXPLICITLY LOADING 'raw' layer.")
      mat <- assay(sce, "raw")
    } else {
      message("    -> WARNING: Only found default layer 'X'. Using assay(sce, 1).")
      mat <- assay(sce, 1)
    }
    
    # 3. CHECK FOR "FLAT" CELLS
    # Check if genes within a single cell have identical non-zero values
    mat_sample <- head(as.matrix(mat), 100)
    row_variances <- apply(mat_sample, 1, function(x) var(x[x > 0]))
    
    if (all(row_variances < 1e-6, na.rm=TRUE)) {
      stop("\n!!! CRITICAL ERROR !!!\n",
           "The selected H5AD layer contains FLATTENED data (identical values per cell).\n",
           "This usually means the data was Binarized and then Normalized.\n",
           "You CANNOT cluster this. Please check the python generation of this file\n",
           "or ensure a 'counts' layer exists containing Integers.")
    }
    
    # 4. CHECK FOR DECIMALS 
    sample_vals <- head(as.vector(mat), 100)
    if (any(sample_vals %% 1 != 0, na.rm=TRUE)) {
      message("    -> Detected decimals. Rounding to nearest integer for Seurat...")
      mat <- round(mat)
    }
    
    # 5. ORIENTATION CHECK (Cells x Genes)
    row_names_sample <- rownames(sce)[1:min(50, nrow(sce))]
    # Regex for E. coli: 
    # 1. b-numbers (b followed by 4 digits)
    # 2. Gene symbols (3 lowercase + 1 uppercase, e.g., dnaK)
    # 3. Operons (anything containing a hyphen like ackA-pta)
    rows_are_genes <- any(grepl("^b[0-9]{4}", row_names_sample)) || 
      any(grepl("^[a-z]{3}[A-Z]", row_names_sample)) ||
      any(grepl("-", row_names_sample))
    
    if (rows_are_genes) {
      message("    [Orientation] H5AD Rows look like GENES. Transposing to (Cells x Genes)...")
      # Input is Genes x Cells -> Transpose to Cells x Genes
      expr_raw <- as.data.frame(t(as.matrix(mat)))
      colnames(expr_raw) <- rownames(sce)
      rownames(expr_raw) <- colnames(sce)
    } else {
      message("    [Orientation] H5AD Rows look like CELLS. Keeping as (Cells x Genes).")
      # Input is Cells x Genes -> Keep as is
      expr_raw <- as.data.frame(as.matrix(mat))
      colnames(expr_raw) <- colnames(sce)
      rownames(expr_raw) <- rownames(sce)
    }
    
    # Extract Conditions
    if ("Condition" %in% names(colData(sce))) {
      conditions <- as.character(colData(sce)$Condition)
    } else if (!is.null(file_set$condition_label)) {
      conditions <- rep(file_set$condition_label, nrow(expr_raw))
    } else {
      conditions <- rep(file_set$id, nrow(expr_raw))
    }
    
    # Create dummy gene index
    dummy_gene <- data.frame(Geneid = colnames(expr_raw), Name = colnames(expr_raw))
    all_gene_indices[[file_set$id]] <- dummy_gene
    
    
   
    
    rm(sce,mat); gc()
    
    # =========================================================================
    # LOGIC BRANCH: SCIFI / CSV FILES
    # =========================================================================
  } else if (file_set$type == "scifi") {
    
    message("    Reading SCIFI/CSV file: ", basename(file_set$expr_path))
    
    # 1. Read Expression File
    # Use check.names=FALSE to keep original gene names (e.g., "b0001" not "b0001...1")
    expr_raw <- as.data.frame(read_csv(file_set$expr_path, show_col_types = FALSE))
    
    # --- STEP A: HANDLE ROW NAMES (CELL IDs) ---
    # If the first column is text (e.g. "Cell.Barcode"), move it to rownames
    if (ncol(expr_raw) > 1 && is.character(expr_raw[[1]])) {
      message("    -> Detected cell IDs in first column ('", colnames(expr_raw)[1], "'). Moving to rownames.")
      rownames(expr_raw) <- expr_raw[[1]]
      expr_raw <- expr_raw[, -1, drop = FALSE] # Remove the ID column
    }
    
    # --- STEP B: HANDLE GENE NAMES ---
    # Check if a gene index file was provided in config AND exists
    if (!is.null(file_set$gene_path) && file.exists(file_set$gene_path)) {
      
      message("    -> Gene index file found. Mapping headers...")
      gene_index <- read_csv(file_set$gene_path, show_col_types = FALSE)
      
      # Map column names (Genes) using the index file
      gene_id <- gene_index$Geneid
      gene_name <- gene_index$Name
      
      # Match current columns to the "Name" column in gene_index
      mapped_ids <- gene_id[match(colnames(expr_raw), gene_name)]
      
      # Safety: only rename if mapping worked
      if (any(is.na(mapped_ids))) {
        warning("    Some column names could not be mapped to gene IDs. Keeping originals.")
      } else {
        colnames(expr_raw) <- as.character(mapped_ids)
      }
      
      all_gene_indices[[file_set$id]] <- gene_index
      
    } else {
      
      # --- NO GENE FILE: USE HEADERS AS IS ---
      message("    -> No gene index provided. Using column headers as Gene IDs.")
      
      # Create a dummy gene index from the columns
      current_genes <- colnames(expr_raw)
      all_gene_indices[[file_set$id]] <- data.frame(Geneid = current_genes, Name = current_genes)
    }
    
    # --- STEP C: PROCESS CELL METADATA (Conditions) ---
    conditions <- rep("Unknown", nrow(expr_raw)) 
    
    # Case 1: External Metadata File
    if (!is.null(file_set$cell_index_path) && file.exists(file_set$cell_index_path)) {
      message("    ... parsing SCIFI metadata from index file")
      meta_raw <- read_csv(file_set$cell_index_path, col_names = FALSE, show_col_types = FALSE)
      
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
      
      # Case 2: Config Label
    } else if (!is.null(file_set$condition_label)) {
      conditions <- rep(file_set$condition_label, nrow(expr_raw))
    }
  } else {
    stop("Unknown file type specified in config: ", file_set$type)
  }
  
  # =========================================================================
  # SHARED POST-PROCESSING (Apply Unique IDs & Store)
  # =========================================================================
  
  # Assign unique cell barcodes
  cell_ids <- paste0(file_set$id, "_cell_", seq_len(nrow(expr_raw)))
  rownames(expr_raw) <- cell_ids
  
  # Store Metadata
  all_metadata_list[[file_set$id]] <- data.frame(
    cell_barcode = cell_ids,
    dataset_id = file_set$id,
    condition = conditions,
    stringsAsFactors = FALSE
  )
  rownames(all_metadata_list[[file_set$id]]) <- cell_ids
  
  # Store Expression
  all_expr_list[[file_set$id]] <- expr_raw
}
gc()

# -------------------------------
# 2) Find Common Gene Index & Filter All Datasets
# -------------------------------
message("Finding common genes across all datasets...")
all_gene_ids_list <- lapply(all_gene_indices, function(df) df$Geneid)
common_gene_ids <- Reduce(intersect, all_gene_ids_list)

if (length(common_gene_ids) == 0) stop("No common gene IDs found.")

common_genes_sorted <- sort(common_gene_ids)
message(sprintf("Found %d common genes.", length(common_genes_sorted)))

filtered_expr_list <- lapply(all_expr_list, function(expr_df) {
  expr_filtered <- expr_df[, colnames(expr_df) %in% common_genes_sorted, drop = FALSE]
  expr_filtered <- expr_filtered[, common_genes_sorted, drop = FALSE]
  return(expr_filtered)
})
rm(all_expr_list, all_gene_indices, all_gene_ids_list); gc()

# -------------------------------
# 3) Combine into Final Matrix
# -------------------------------
message("Combining all datasets into one matrix...")
combined_df <- do.call(rbind, filtered_expr_list)
combined_meta <- do.call(rbind, all_metadata_list)
rm(filtered_expr_list, all_metadata_list); gc()

# Convert to sparse matrix (genes x cells)
combined_sparse <- Matrix(as.matrix(combined_df), sparse = TRUE)
combined_sparse <- t(combined_sparse)
rm(combined_df); gc()

# =========================================================================
# 3b) ROBUST MAPPING: Handles Symbols AND Operons (Split -> Map -> Join)
# =========================================================================
message("Fetching iML1515 dictionary from BiGG Models...")

url <- "http://bigg.ucsd.edu/api/v2/models/iML1515/genes"

gene_map_data <- tryCatch({
  json_response <- jsonlite::fromJSON(url)
  
  if ("results" %in% names(json_response)) {
    json_response$results
  } else {
    NULL
  }
}, error = function(e) {
  message("Warning: Could not download gene map. Proceeding with original names.")
  return(NULL)
})

if (!is.null(gene_map_data)) {
  if (is.list(gene_map_data) && !is.data.frame(gene_map_data)) {
    gene_map_data <- dplyr::bind_rows(gene_map_data)
  }
  
  if (all(c("name", "bigg_id") %in% colnames(gene_map_data))) {
    
    map_df <- gene_map_data[, c("name", "bigg_id")]
    colnames(map_df) <- c("Symbol", "LocusTag")
    
    map_df <- map_df[!is.na(map_df$Symbol) & map_df$Symbol != "", ]
    
    gene_lookup <- setNames(map_df$LocusTag, tolower(map_df$Symbol))
    
    current_genes <- rownames(combined_sparse)
    new_genes <- current_genes 
    
    message("Mapping gene symbols/operons to b-numbers...")
    mapped_count <- 0
    
    for (i in seq_along(current_genes)) {
      original_name <- current_genes[i]
      
      tokens <- unlist(strsplit(original_name, "[-_ ]"))
      
      mapped_tokens <- sapply(tokens, function(t) {
        t_low <- tolower(t)
        if (t_low %in% names(gene_lookup)) {
          return(gene_lookup[[t_low]])
        } else {
          return(t)
        }
      })
      
      if (any(mapped_tokens != tokens)) {
        new_name <- paste(mapped_tokens, collapse = "-")
        new_genes[i] <- new_name
        mapped_count <- mapped_count + 1
      }
    }
    
    rownames(combined_sparse) <- new_genes
    
    message(sprintf("Successfully translated %d / %d features to iML1515 IDs.", 
                    mapped_count, length(current_genes)))
    
    rm(map_df, gene_lookup)
  } else {
    message("Error: JSON columns found: ", paste(colnames(gene_map_data), collapse=", "))
    message("Expected 'name' and 'bigg_id'. Skipping mapping.")
  }
  rm(gene_map_data); gc()
}

# -------------------------------
# 4) Create Seurat Object & QC
# -------------------------------
message("Creating single Seurat object...")
seurat_obj <- CreateSeuratObject(
  counts = combined_sparse,
  project = config$project_name,
  min.cells = config$min_cells_qc,
  min.features = config$min_features_qc,
  meta.data = combined_meta
)

# Filter cells by genes per cell
raw_counts <- GetAssayData(seurat_obj, layer = "counts")
genes_per_cell <- Matrix::colSums(raw_counts > 0)
keep_cells <- names(genes_per_cell)[genes_per_cell >= config$genes_per_cell_qc]
seurat_obj <- subset(seurat_obj, cells = keep_cells)

message(sprintf("%d cells remain after initial filtering.", length(keep_cells)))

# -------------------------------
# 5) Normalize, Scale, and Run PCA
# -------------------------------
message("Normalizing and running PCA...")
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 1e6, verbose = FALSE)

# Convert ln to log2
mat <- GetAssayData(seurat_obj, layer = "data")
mat <- mat / log(2)
seurat_obj <- SetAssayData(seurat_obj, layer = "data", new.data = mat)
rm(mat); gc()

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = config$n_features)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)

# -------------------------------
# 6) Save Elbow Plot
# -------------------------------
elbow_plot <- ElbowPlot(seurat_obj, ndims = 50)
ggsave(file.path(plot_output_dir, "1_ElbowPlot.png"), elbow_plot, width = 8, height = 6)

# -------------------------------
# 7) Clustering & Post-Filter
# -------------------------------
message("Clustering...")
seurat_obj <- FindNeighbors(seurat_obj, dims = config$pca_dims)
seurat_obj <- FindClusters(seurat_obj, resolution = config$cluster_res)

# Filter small clusters
cluster_counts <- table(Idents(seurat_obj))
clusters_to_keep <- names(cluster_counts)[cluster_counts >= config$min_cluster_size]
if (length(clusters_to_keep) == 0) stop("No clusters met size threshold.")
seurat_obj <- subset(seurat_obj, idents = clusters_to_keep)
message(sprintf("%d cells remain after cluster filtering.", ncol(seurat_obj)))

# -------------------------------
# 8) Run UMAP and Save PLOTS
# -------------------------------
message("Running UMAP and generating condition plots...")
seurat_obj <- RunUMAP(seurat_obj, dims = config$pca_dims, verbose = FALSE)

# A. Standard Cluster Plot
p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP by Cluster")

# B. Condition Plot (if conditions are available)
if ("condition" %in% colnames(seurat_obj@meta.data) && config$use_conditions) {
  p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "condition", pt.size = 0.5) +
    ggtitle("UMAP by Condition")
  
  # Save Combined UMAPs
  ggsave(file.path(plot_output_dir, "2_UMAP_Clusters_and_Conditions.png"), 
         p1 + p2, width = 14, height = 7)
  
  # C. Cluster Composition Bar Plot
  df_meta <- seurat_obj@meta.data
  df_summary <- df_meta %>%
    group_by(seurat_clusters, condition) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(seurat_clusters) %>%
    mutate(prop = count / sum(count))
  
  p3 <- ggplot(df_summary, aes(x = seurat_clusters, y = count, fill = condition)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Cluster Composition by Condition", x = "Cluster", y = "Cell Count")
  
  p4 <- ggplot(df_summary, aes(x = seurat_clusters, y = prop, fill = condition)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Cluster Composition (Proportion)", x = "Cluster", y = "Proportion")
  
  ggsave(file.path(plot_output_dir, "3_Cluster_Composition_Barplot.png"), 
         p3 + p4, width = 14, height = 7)
} else {
  # Save only cluster plot if conditions not available
  ggsave(file.path(plot_output_dir, "2_UMAP_Clusters.png"), 
         p1, width = 8, height = 7)
}

# -------------------------------
# 9) Extract Final Data
# -------------------------------
# Ensure condition column exists
safe_condition <- if ("condition" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$condition
} else {
  rep("Unspecified", ncol(seurat_obj))
}

cell_clusters <- data.frame(
  cell_barcode = names(Idents(seurat_obj)),
  cluster = as.character(Idents(seurat_obj)),
  condition = safe_condition,
  stringsAsFactors = FALSE
)

norm_data <- GetAssayData(seurat_obj, layer = "data")

# -------------------------------
# 10) Write Data Outputs
# -------------------------------
message("Writing data files...")

norm_data <- norm_data[, cell_clusters$cell_barcode, drop = FALSE]
data.table::fwrite(data.table(gene = rownames(norm_data)), file.path(data_output_dir, "genes_logCPM.csv"))
data.table::fwrite(cell_clusters, file.path(data_output_dir, "cell_metadata.csv"))

clusters <- sort(unique(cell_clusters$cluster))
for (cl in clusters) {
  cols <- cell_clusters$cell_barcode[cell_clusters$cluster == cl]
  sub <- norm_data[, cols, drop = FALSE]
  
  df <- data.table::data.table(
    cell_id = cols,
    cluster = cl,
    condition = cell_clusters$condition[cell_clusters$cluster == cl],
    as.data.frame(t(as.matrix(sub)), check.names = FALSE)
  )
  
  fn <- file.path(data_output_dir, sprintf("cluster_%s_logCPM.csv", cl))
  data.table::fwrite(df, fn, sep = ",")
  message("Wrote: ", fn)
}

# -------------------------------
# 10b) Export Condition-Based Files (Optional)
# -------------------------------
if (config$use_conditions && "condition" %in% colnames(cell_clusters)) {
  message("Exporting data grouped by condition...")
  
  unique_conditions <- unique(cell_clusters$condition)
  
  for (cond in unique_conditions) {
    safe_cond <- gsub("[^A-Za-z0-9]", "_", cond)
    target_cells <- cell_clusters$cell_barcode[cell_clusters$condition == cond]
    
    if (length(target_cells) > 0) {
      sub_data <- norm_data[, target_cells, drop = FALSE]
      
      df_export <- data.table::data.table(
        cell_id = target_cells,
        cluster = cell_clusters$cluster[cell_clusters$condition == cond],
        as.data.frame(t(as.matrix(sub_data)), check.names = FALSE)
      )
      
      out_file <- file.path(condition_output_dir, sprintf("Cluster_%s_logCPM.csv", safe_cond))
      data.table::fwrite(df_export, out_file, sep = ",")
      message(sprintf("   -> Saved: %s (%d cells)", basename(out_file), length(target_cells)))
    }
  }
}

# -------------------------------
# 11) COMPUTE CLUSTERING QUALITY METRICS
# -------------------------------
message("Computing clustering quality metrics...")

# Function to compute comprehensive clustering quality metrics
compute_clustering_quality <- function(seurat_obj, pca_dims, cluster_labels) {
  
  # Get PCA embedding for distance calculations
  pca_embedding <- Embeddings(seurat_obj, "pca")[, pca_dims]
  
  # Get cluster assignments as numeric
  clusters <- as.numeric(as.factor(cluster_labels))
  cluster_ids <- unique(clusters)
  n_clusters <- length(cluster_ids)
  
  # Initialize per-cluster data frame
  per_cluster_metrics <- data.frame()
  
  # 1. Compute per-cluster metrics
  message(" Computing per-cluster quality metrics...")
  for (cluster_id in cluster_ids) {
    cluster_points <- pca_embedding[clusters == cluster_id, , drop = FALSE]
    cluster_size <- nrow(cluster_points)
    
    if (cluster_size > 1) {
      # Cluster centroid
      centroid <- colMeans(cluster_points)
      
      # Distance calculations
      distances_to_centroid <- sqrt(rowSums((cluster_points - matrix(centroid, nrow = cluster_size, ncol = length(centroid), byrow = TRUE))^2))
      avg_distance_to_centroid <- mean(distances_to_centroid)
      max_distance_to_centroid <- max(distances_to_centroid)
      median_distance_to_centroid <- median(distances_to_centroid)
      
      # Compactness metrics
      cv_distance <- ifelse(avg_distance_to_centroid > 0, sd(distances_to_centroid) / avg_distance_to_centroid, 0)
      cluster_wss <- sum(distances_to_centroid^2)
      
      # Density metrics
      if (cluster_size > 2) {
        cluster_dist <- dist(cluster_points)
        cluster_dist_matrix <- as.matrix(cluster_dist)
        diag(cluster_dist_matrix) <- Inf
        avg_nearest_neighbor <- mean(apply(cluster_dist_matrix, 1, min))
        isolation_index <- mean(apply(cluster_dist_matrix, 1, function(x) sort(x)[2])) # Distance to second nearest
      } else {
        avg_nearest_neighbor <- NA
        isolation_index <- NA
      }
      
    } else {
      # For clusters with only 1 point
      avg_distance_to_centroid <- 0
      max_distance_to_centroid <- 0
      median_distance_to_centroid <- 0
      cv_distance <- 0
      cluster_wss <- 0
      avg_nearest_neighbor <- NA
      isolation_index <- NA
    }
    
    # Store per-cluster metrics
    cluster_metrics <- data.frame(
      ClusterID = cluster_id,
      ClusterSize = cluster_size,
      AvgDistanceToCentroid = avg_distance_to_centroid,
      MaxDistanceToCentroid = max_distance_to_centroid,
      MedianDistanceToCentroid = median_distance_to_centroid,
      DistanceCV = cv_distance,
      WithinClusterSS = cluster_wss,
      AvgNearestNeighborDist = avg_nearest_neighbor,
      IsolationIndex = isolation_index,
      AvgSilhouetteWidth = NA, # Will be filled later
      stringsAsFactors = FALSE
    )
    
    per_cluster_metrics <- rbind(per_cluster_metrics, cluster_metrics)
  }
  
  # 2. Compute overall metrics
  message(" Computing overall clustering quality...")
  
  # Overall compactness
  wss <- sum(per_cluster_metrics$WithinClusterSS, na.rm = TRUE)
  
  # Separation (between-cluster)
  overall_centroid <- colMeans(pca_embedding)
  bss <- 0
  for (i in 1:nrow(per_cluster_metrics)) {
    cluster_id <- per_cluster_metrics$ClusterID[i]
    cluster_points <- pca_embedding[clusters == cluster_id, , drop = FALSE]
    if (nrow(cluster_points) > 0) {
      centroid <- colMeans(cluster_points)
      bss <- bss + nrow(cluster_points) * sum((centroid - overall_centroid)^2)
    }
  }
  
  # Silhouette scores
  sil_score <- NA
  message(" Computing silhouette scores...")
  tryCatch({
    max_cells_for_silhouette <- 5000
    if (nrow(pca_embedding) > max_cells_for_silhouette) {
      message("Large dataset detected, sampling for silhouette calculation...")
      set.seed(123)
      sample_indices <- sample(1:nrow(pca_embedding), max_cells_for_silhouette)
      dist_matrix <- dist(pca_embedding[sample_indices, ])
      sil_result <- silhouette(clusters[sample_indices], dist_matrix)
    } else {
      dist_matrix <- dist(pca_embedding)
      sil_result <- silhouette(clusters, dist_matrix)
    }
    sil_score <- summary(sil_result)$avg.width
    
    # Add per-cluster silhouette widths
    sil_summary <- summary(sil_result)
    if (!is.null(sil_summary$clus.avg.widths)) {
      cluster_sil_widths <- sil_summary$clus.avg.widths
      for (i in 1:length(cluster_sil_widths)) {
        cluster_idx <- which(per_cluster_metrics$ClusterID == i)
        if (length(cluster_idx) > 0) {
          per_cluster_metrics$AvgSilhouetteWidth[cluster_idx] <- cluster_sil_widths[i]
        }
      }
    }
  }, error = function(e) {
    message("Silhouette calculation failed: ", e$message)
  })
  
  # Calinski-Harabasz Index
  ch_index <- NA
  n_points <- length(clusters)
  if (n_clusters > 1 && n_points > n_clusters && !is.na(wss) && wss > 0) {
    ch_index <- (bss / (n_clusters - 1)) / (wss / (n_points - n_clusters))
  }
  
  # Davies-Bouldin Index
  db_index <- NA
  if (n_clusters > 1) {
    cluster_centroids <- matrix(0, nrow = n_clusters, ncol = ncol(pca_embedding))
    for (i in 1:n_clusters) {
      cluster_points <- pca_embedding[clusters == i, , drop = FALSE]
      if (nrow(cluster_points) > 0) {
        cluster_centroids[i, ] <- colMeans(cluster_points)
      }
    }
    
    db_sum <- 0
    valid_clusters <- 0
    for (i in 1:n_clusters) {
      cluster_i_metrics <- per_cluster_metrics[per_cluster_metrics$ClusterID == i, ]
      if (nrow(cluster_i_metrics) > 0 && cluster_i_metrics$ClusterSize > 0) {
        max_similarity <- 0
        s_i <- cluster_i_metrics$AvgDistanceToCentroid
        
        for (j in 1:n_clusters) {
          if (i != j) {
            cluster_j_metrics <- per_cluster_metrics[per_cluster_metrics$ClusterID == j, ]
            if (nrow(cluster_j_metrics) > 0 && cluster_j_metrics$ClusterSize > 0) {
              s_j <- cluster_j_metrics$AvgDistanceToCentroid
              d_ij <- sqrt(sum((cluster_centroids[i, ] - cluster_centroids[j, ])^2))
              
              if (d_ij > 0) {
                r_ij <- (s_i + s_j) / d_ij
                max_similarity <- max(max_similarity, r_ij)
              }
            }
          }
        }
        db_sum <- db_sum + max_similarity
        valid_clusters <- valid_clusters + 1
      }
    }
    if (valid_clusters > 0) {
      db_index <- db_sum / valid_clusters
    }
  }
  
  # Cluster size statistics
  cluster_size_stats <- per_cluster_metrics$ClusterSize
  avg_cluster_size <- mean(cluster_size_stats)
  min_cluster_size <- min(cluster_size_stats)
  max_cluster_size <- max(cluster_size_stats)
  size_std_dev <- sd(cluster_size_stats)
  size_cv <- ifelse(avg_cluster_size > 0, size_std_dev / avg_cluster_size, NA)
  
  # Cluster balance
  cluster_proportions <- cluster_size_stats / sum(cluster_size_stats)
  balance_entropy <- -sum(cluster_proportions * log(cluster_proportions)) / log(length(cluster_proportions))
  
  # Compactness statistics
  avg_compactness <- mean(per_cluster_metrics$AvgDistanceToCentroid, na.rm = TRUE)
  compactness_std <- sd(per_cluster_metrics$AvgDistanceToCentroid, na.rm = TRUE)
  
  # Return both overall and per-cluster metrics
  return(list(
    overall_metrics = data.frame(
      ClusterConfiguration = run_name,
      TotalCells = length(clusters),
      NumberOfClusters = n_clusters,
      
      # Quality metrics
      SilhouetteScore = sil_score,
      CalinskiHarabaszIndex = ch_index,
      DaviesBouldinIndex = db_index,
      
      # Variance metrics
      WithinClusterSumSquares = wss,
      BetweenClusterSumSquares = bss,
      VarianceRatio = ifelse(!is.na(wss) && wss > 0, bss / wss, NA),
      
      # Size metrics
      AverageClusterSize = avg_cluster_size,
      MinClusterSize = min_cluster_size,
      MaxClusterSize = max_cluster_size,
      ClusterSizeStdDev = size_std_dev,
      ClusterSizeCV = size_cv,
      ClusterBalanceEntropy = balance_entropy,
      
      # Compactness metrics
      AvgClusterCompactness = avg_compactness,
      CompactnessStdDev = compactness_std,
      
      # Configuration
      PCA_Dimensions_Used = length(pca_dims),
      ClusteringResolution = config$cluster_res,
      VariableFeatures = config$n_features,
      
      stringsAsFactors = FALSE
    ),
    
    per_cluster_metrics = per_cluster_metrics
  ))
}

# Compute quality metrics
cluster_labels <- as.character(Idents(seurat_obj))
quality_results <- compute_clustering_quality(seurat_obj, config$pca_dims, cluster_labels)

# Save quality metrics
quality_file <- file.path(data_output_dir, "clustering_quality_metrics.csv")
write.csv(quality_results$overall_metrics, quality_file, row.names = FALSE)
message("Overall clustering quality metrics saved to: ", quality_file)

per_cluster_file <- file.path(data_output_dir, "per_cluster_quality_metrics.csv")
write.csv(quality_results$per_cluster_metrics, per_cluster_file, row.names = FALSE)
message("Per-cluster quality metrics saved to: ", per_cluster_file)

# Print summary
message("\n--- CLUSTERING QUALITY SUMMARY ---")
message(sprintf("Silhouette Score: %.3f", quality_results$overall_metrics$SilhouetteScore))
message(sprintf("Number of Clusters: %d", quality_results$overall_metrics$NumberOfClusters))
message(sprintf("Avg Cluster Size: %.1f", quality_results$overall_metrics$AverageClusterSize))
message(sprintf("Avg Compactness: %.3f", quality_results$overall_metrics$AvgClusterCompactness))

message("\n--- PER-CLUSTER COMPACTNESS ---")
for (i in 1:nrow(quality_results$per_cluster_metrics)) {
  cm <- quality_results$per_cluster_metrics[i, ]
  message(sprintf("Cluster %d: size=%d, compactness=%.3f, sil=%.3f",
                  cm$ClusterID, cm$ClusterSize, cm$AvgDistanceToCentroid,
                  ifelse(is.na(cm$AvgSilhouetteWidth), NA, cm$AvgSilhouetteWidth)))
}

# Clean up
rm(seurat_obj, combined_sparse, norm_data); gc()

# -------------------------------
# 12) Save the Settings File
# -------------------------------
settings_text <- paste(
  "--- Settings for run:", run_name, "---",
  "\nProject Name: ", config$project_name,
  "\nDate: ", Sys.time(),
  "\nProject Root: ", project_root,
  "\nInput Datasets Processed: ", paste(sapply(config$file_sets, `[[`, "id"), collapse=", "),
  "\nUse Conditions: ", config$use_conditions,
  "\n\n-- QC PARAMETERS --",
  "\nmin.cells (CreateSeurat): ", config$min_cells_qc,
  "\nmin.features (CreateSeurat): ", config$min_features_qc,
  "\nGenes per cell (custom filter): ", config$genes_per_cell_qc,
  "\nMin cluster size: ", config$min_cluster_size,
  "\n\n-- ANALYSIS PARAMETERS --",
  "\nVariable Features (nfeatures): ", config$n_features,
  "\nPCA Dims Used: ", paste(config$pca_dims, collapse=","),
  "\nClustering Resolution: ", config$cluster_res,
  "\n\n-- FINAL OBJECT --",
  "\nCommon genes found: ", length(common_genes_sorted),
  "\nFinal cell count: ", nrow(cell_clusters),
  "\nFinal cluster count: ", length(clusters_to_keep),
  "\nFinal clusters: ", paste(sort(clusters_to_keep), collapse=", ")
)


if (exists("quality_results")) {
  settings_text <- paste(
    settings_text,
    "\n\n-- CLUSTERING QUALITY --",
    sprintf("\nSilhouette Score: %.3f", quality_results$overall_metrics$SilhouetteScore),
    sprintf("\nCalinski-Harabasz: %.1f", quality_results$overall_metrics$CalinskiHarabaszIndex),
    sprintf("\nDavies-Bouldin: %.3f", quality_results$overall_metrics$DaviesBouldinIndex),
    sprintf("\nCluster Balance: %.3f", quality_results$overall_metrics$ClusterBalanceEntropy),
    sprintf("\nAvg Compactness: %.3f", quality_results$overall_metrics$AvgClusterCompactness)
  )
}

writeLines(settings_text, file.path(plot_output_dir, "run_settings.txt"))

message("--- Run Complete: ", run_name, " ---")