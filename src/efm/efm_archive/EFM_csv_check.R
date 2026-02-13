
# Vector of CSV files to process
files <- c(
  "/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/efm_fluxes_output_gurobi.csv",
  "/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/efm_fluxes_output_builtin.csv",
  "/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/efm_fluxes_output_gurobi_150.csv",
  "/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/efm_fluxes_output_gurobi_500.csv",
  "/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/efm_fluxes_output_gurobi_50.csv",
  "/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/efm_fluxes_output_gurobi_50_2.csv",
  "/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/efm_fluxes_output_gurobi_50_3.csv",
  "/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/efm_fluxes_output_gurobi_50_4.csv",
  "/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/efm_fluxes_output_gurobi_50_5.csv",
  "/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/efm_fluxes_output_gurobi_50_6.csv",
  "/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/efm_fluxes_output_gurobi_50_7.csv",
  "/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/efm_fluxes_output_gurobi_50_8.csv",
  "/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/efm_fluxes_output_gurobi_50_9.csv",
  "/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/efm_fluxes_output_gurobi_100_7.csv",
  "/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/efm_fluxes_output_gurobi_500_7.csv",
  "/Users/benedikthaarscheidt/M.Sc./master_thesis/computed_EFM/efm_fluxes_output_gurobi_1000_7.csv"
  
  # extend as needed
)

# Function to process one file
process_file <- function(filepath) {
  df <- read.csv(filepath, row.names = 1, check.names = FALSE)
  df <- round(df, 3)
  modes <- colnames(df)
  
  # Build nested list: modes[[mode]] = list(indices=…, names=…)
  modes_list <- lapply(modes, function(m) {
    # indices of nonzero entries
    idx <- which(df[[m]] != 0)
    rxn_names <- rownames(df)[idx]
    if (any(duplicated(rxn_names))) {
      warning(sprintf("File '%s', mode '%s': duplicates in support", 
                      basename(filepath), m))
    }
    list(
      indices = sort(idx),
      names   = rxn_names
    )
  })
  names(modes_list) <- modes
  
  # Pairwise difference matrix
  N <- length(modes)
  diff_mat <- matrix(0, nrow = N, ncol = N,
                     dimnames = list(modes, modes))
  for (i in seq_len(N)) {
    for (j in seq_len(N)) {
      a <- modes_list[[i]]$names
      b <- modes_list[[j]]$names
      diff_mat[i, j] <- length(setdiff(a, b))
    }
  }
  
  # Union size
  all_rxns <- unique(unlist(lapply(modes_list, `[[`, "names")))
  union_size <- length(all_rxns)
  
  list(
    file        = basename(filepath),
    modes       = modes_list,
    diff_mat    = diff_mat,
    union_size  = union_size
  )
}

# Process all files
results <- lapply(files, process_file)

# Example: inspect results for the first file
res5 <- results[[6]]
print(res5$file)          # filename
print(res5$union_size)    # total unique reactions
print(res5$modes)         # nested list of indices & names
print(res5$diff_mat)      # difference matrix

# Build a summary table of union sizes
union_summary <- data.frame(
  file       = sapply(results, `[[`, "file"),
  union_size = sapply(results, `[[`, "union_size"),
  stringsAsFactors = FALSE
)
print(union_summary)
