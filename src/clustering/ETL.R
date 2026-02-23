# ==============================================================================
# SCRIPT 1: DATA PREPARATION (SOURCED BY CLUSTERING SCRIPT)
# ==============================================================================

# -------------------------------
# 1) SETUP & CONFIG
# -------------------------------
library(readr)
library(Matrix)
library(dplyr)
library(zellkonverter)
library(jsonlite)
library(SingleCellExperiment)

set.seed(123)

# Handle Paths (Robust check for sourcing vs direct run)
if (sys.nframe() == 0) {
  # If run directly
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg)
    project_root <- normalizePath(file.path(dirname(script_path), "../../"))
  } else {
    project_root <- if (dir.exists("config")) getwd() else "~/"
  }
} else {
  # If sourced, assume project root is the working directory
  project_root <- getwd()
}

config_file <- file.path(project_root, "config", "config.json")
if (!file.exists(config_file)) config_file <- "config/config.json"

if (file.exists(config_file)) {
  message("  [Prep] Loading config: ", config_file)
  global_config <- fromJSON(config_file)
  datasets_base <- file.path(project_root, global_config$paths$datasets_dir)
  
  # Load specific prep params
  config <- list(
    project_name = global_config$params$project_name,
    use_conditions = global_config$params$use_conditions,
    
    # Clustering params (needed later in Script 2, but loaded here for consistency)
    min_cells_qc = global_config$params$min_cells_qc,
    min_features_qc = global_config$params$min_features_qc,
    genes_per_cell_qc = global_config$params$genes_per_cell_qc,
    n_features = global_config$params$n_features,
    pca_dims = global_config$params$pca_dims,
    cluster_res = global_config$params$cluster_res,
    min_cluster_size = global_config$params$min_cluster_size,
    
    results_base = file.path(project_root, global_config$paths$results_dir)
  )
  
  resolve_path <- function(p) {
    if (is.null(p) || is.na(p) || p == "") return(NULL)
    if (startsWith(p, "/") || startsWith(p, "~")) return(p)
    return(file.path(datasets_base, p))
  }
  
  input_raw <- global_config$clustering$input_files
  config$file_sets <- lapply(seq_len(nrow(input_raw)), function(i) {
    row <- input_raw[i, ]
    expr_path <- resolve_path(row$expr_path)
    type <- if (endsWith(tolower(expr_path), ".h5ad")) "h5ad" else "csv_scifi"
    
    cell_index_path <- NULL
    if (config$use_conditions && !is.null(row$cell_index_path)) {
      cell_index_path <- resolve_path(row$cell_index_path)
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
  stop("Config file not found.")
}

# -------------------------------
# 2) LOAD & CLEAN DATASETS
# -------------------------------
all_expr_list <- list()
all_gene_indices <- list()
all_metadata_list <- list()

message("  [Prep] Processing ", length(config$file_sets), " datasets...")

for (file_set in config$file_sets) {
  expr_raw <- NULL
  
  # --- H5AD ---
  if (file_set$type == "h5ad") {
    sce <- zellkonverter::readH5AD(file_set$expr_path)
    if ("counts" %in% assayNames(sce)) mat <- assay(sce, "counts") else mat <- assay(sce, 1)
    
    # Orientation Check
    row_sample <- rownames(sce)[1:min(50, nrow(sce))]
    rows_are_genes <- any(grepl("^b[0-9]{4}", row_sample)) || any(grepl("^[a-z]{3}[A-Z]", row_sample))
    
    if (rows_are_genes) {
      expr_raw <- as.data.frame(t(as.matrix(mat)))
      colnames(expr_raw) <- rownames(sce)
      rownames(expr_raw) <- colnames(sce)
    } else {
      expr_raw <- as.data.frame(as.matrix(mat))
      colnames(expr_raw) <- colnames(sce)
      rownames(expr_raw) <- rownames(sce)
    }
    all_gene_indices[[file_set$id]] <- data.frame(Geneid = colnames(expr_raw), Name = colnames(expr_raw))
    rm(sce, mat); gc()
    
    # --- CSV / SCIFI ---
  } else {
    raw_df <- read_csv(file_set$expr_path, show_col_types = FALSE)
    raw_df <- as.data.frame(raw_df)
    
    # Remove Junk Cols
    junk_cols <- grep("^\\.\\.\\.[0-9]+$|^X$|^Unnamed", colnames(raw_df))
    if (length(junk_cols) > 0) raw_df <- raw_df[, -junk_cols, drop = FALSE]
    
    # Handle Cell IDs (First Text Col)
    char_cols <- sapply(raw_df, is.character)
    if (any(char_cols)) {
      idx <- which(char_cols)[1]
      rownames(raw_df) <- raw_df[[idx]]
      expr_raw <- raw_df[, -idx, drop = FALSE]
    } else {
      expr_raw <- raw_df
      rownames(expr_raw) <- paste0(file_set$id, "_cell_", 1:nrow(expr_raw))
    }
    
    # Handle Genes
    if (!is.null(file_set$gene_path) && file.exists(file_set$gene_path)) {
      g_idx <- read_csv(file_set$gene_path, show_col_types = FALSE)
      mapped <- g_idx$Geneid[match(colnames(expr_raw), g_idx$Name)]
      if (!any(is.na(mapped))) colnames(expr_raw) <- as.character(mapped)
      all_gene_indices[[file_set$id]] <- g_idx
    } else {
      all_gene_indices[[file_set$id]] <- data.frame(Geneid = colnames(expr_raw), Name = colnames(expr_raw))
    }
  }
  
  # --- CONDITION LOGIC ---
  if (config$use_conditions) {
    # Default: Use the specific label from config OR the dataset ID
    default_label <- if(!is.null(file_set$condition_label)) file_set$condition_label else file_set$id
    conditions <- rep(default_label, nrow(expr_raw))
    
    # 2. Parse external file if it exists and is valid
    if (!is.null(file_set$cell_index_path) && file.exists(file_set$cell_index_path)) {
      tryCatch({
        # Read the full index file
        meta_raw <- read_csv(file_set$cell_index_path, col_names = FALSE, show_col_types = FALSE)
        
        # Look for the column containing the library/condition identifiers
        for(i in seq_along(meta_raw)) {
          # Check if the column starts with a known prefix like 'Lib'
          if(any(grepl("^Lib", as.character(meta_raw[[i]][1:min(10, nrow(meta_raw))])))) {
            raw_str <- as.character(meta_raw[[i]])
            
            # Extract the unique condition ID (usually the text after the last underscore)
            # Example: "Batch1_WildType" -> "WildType"
            parsed <- sub(".*_([^_]+)$", "\\1", raw_str)
            
            # CRITICAL: Only apply if the length matches the number of cells in the expression matrix
            if(length(parsed) == nrow(expr_raw)) {
              conditions <- parsed
              message(sprintf("    [Prep] Mapped %d granular conditions for %s", length(unique(conditions)), file_set$id))
            }
            break # Exit loop once the correct column is found
          }
        }
      }, error = function(e) warning(sprintf("    [Prep] Condition parsing failed for %s", file_set$id)))
    }
  } else {
    # If conditions are disabled, fill with NA to keep metadata structure consistent
    conditions <- rep(NA, nrow(expr_raw))
  }
  
  # Finalize IDs
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
}

# -------------------------------
# 3) MERGE & MAP
# -------------------------------
message("  [Prep] Merging & Mapping...")
common_genes <- Reduce(intersect, lapply(all_gene_indices, function(x) x$Geneid))
filtered_list <- lapply(all_expr_list, function(df) df[, common_genes, drop=FALSE])
combined_df <- do.call(rbind, filtered_list)
combined_meta <- do.call(rbind, all_metadata_list)

# Create FINAL OUTPUT VARIABLES
final_counts <- t(Matrix(as.matrix(combined_df), sparse = TRUE)) # Genes x Cells
final_metadata <- combined_meta

# BiGG Mapping
tryCatch({
  url <- "http://bigg.ucsd.edu/api/v2/models/iML1515/genes"
  map_data <- jsonlite::fromJSON(url)$results
  if (!is.null(map_data)) {
    map_df <- map_data[!is.na(map_data$name) & map_data$name != "", c("name", "bigg_id")]
    lookup <- setNames(map_df$bigg_id, tolower(map_df$name))
    curr <- rownames(final_counts)
    new_g <- curr
    for (i in seq_along(curr)) {
      toks <- unlist(strsplit(curr[i], "[-_ ]"))
      mapped <- sapply(toks, function(t) if (tolower(t) %in% names(lookup)) lookup[[tolower(t)]] else t)
      if (any(mapped != toks)) new_g[i] <- paste(mapped, collapse = "-")
    }
    rownames(final_counts) <- new_g
  }
}, error = function(e) message("  [Prep] Mapping skipped (API error)."))

message("  [Prep] Done. Variables 'final_counts', 'final_metadata', and 'config' are ready.")