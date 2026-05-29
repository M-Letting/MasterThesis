# Initialize overall start time
overall_start_time <- Sys.time()
# Set to not show warnings
options(warn = -1)

suppressMessages(library(arrow, warn.conflicts = FALSE))
suppressMessages(library(data.table, warn.conflicts = FALSE))
suppressMessages(library(limma, warn.conflicts = FALSE))
# Need discover_complexoforms from ComplexoFinder workflow
setwd(here::here())
source("ComplexoFinder/02-Identification/AssignComplexes.R")

# Set back to Q-Benchmark
setwd(here::here("Q-Benchmark"))

# Globals and Constantss
data_ids <- c("1pep", "2pep", "050pep", "random") # Data IDs to process
score_cutoff <- 0.1
adj_pval_cutoff <- 0.01

# Function to format time difference
format_time_diff <- function(time_diff) {
  if (time_diff < 60) {
    return(paste(round(time_diff, 2), "seconds"))
  } else if (time_diff < 3600) {
    return(paste(round(time_diff / 60, 2), "minutes"))
  } else {
    return(paste(round(time_diff / 3600, 2), "hours"))
  }
}


# Create results directory if it doesn't exist
if (!dir.exists("./01-FindDCF/data/results")) {
  dir.create("./01-FindDCF/data/results", recursive = TRUE)
  cat("Created results directory: ./01-FindDCF/data/results\n")
}

# Initialize results tracking
results_summary <- list()

# Process each data ID
for (i in seq_along(data_ids)) {
  data_id <- data_ids[i]

  cat(sprintf(
    "\n[%d/%d] Processing dataset: %s\n",
    i,
    length(data_ids),
    data_id
  ))
  cat(
    "--------------------------------------------------------------------------------\n"
  )

  # Initialize time for this dataset
  start_time <- Sys.time()

  # Step 1: Data loading
  step_start <- Sys.time()
  cat("Step 1: Loading data...\n")
  # input_data_path <- paste0('./data/processed/bench_', data_id, "_input.feather")
  input_data_path <- paste0(
    './01-FindDCF/data/prepared/bench_',
    data_id,
    "_input.feather"
  )

  # Check if file exists
  if (!file.exists(input_data_path)) {
    cat("  ERROR: File not found:", input_data_path, "\n")
    results_summary[[data_id]] <- list(
      status = "FAILED",
      error = "File not found"
    )
    next
  }

  # Read the input data (feather format)
  input_data <- arrow::read_feather(input_data_path)

  cat("  Data read from:", input_data_path, "\n")
  cat(
    "  Dimensions:",
    nrow(input_data),
    "rows x",
    ncol(input_data),
    "columns\n"
  )
  cat(
    "  Time taken:",
    format_time_diff(as.numeric(difftime(
      Sys.time(),
      step_start,
      units = "secs"
    ))),
    "\n\n"
  )

  # Convert to data.table
  input_data <- as.data.table(input_data)

  # Step 2: Convert from long to wide format
  step_start <- Sys.time()
  cat("Step 2: Converting data from long to wide format...\n")

  # Create a unique sample identifier combining day and filename
  # This will become the column names in wide format
  input_data[, sample_id := paste(day, filename, sep = "_")]

  # Convert from long to wide format using dcast
  # Each row will be a unique protein_id + peptide_id combination
  # Each column will be a sample (day_filename combination)
  wide_data <- dcast(
    input_data,
    protein_id +
      peptide_id +
      n_pep +
      n_perturbed_peptides +
      perturbed_protein +
      perturbed_peptide +
      red_fac ~ sample_id,
    value.var = "intensity",
    fun.aggregate = mean # In case of duplicates, take mean
  )

  # Get intensity column names (all columns that are sample IDs)
  intensity_cols <- grep("^day[0-9]+_", names(wide_data), value = TRUE)

  cat("  Original dimensions (long):", nrow(input_data), "rows\n")
  cat(
    "  Wide format dimensions:",
    nrow(wide_data),
    "rows x",
    ncol(wide_data),
    "columns\n"
  )
  cat("  Intensity columns:", length(intensity_cols), "\n")

  # Step 2b: Log transform and median normalize intensity columns
  cat("  Applying log transformation and median normalization...\n")

  # Log transform intensity columns (log)
  for (col in intensity_cols) {
    wide_data[[col]] <- log(wide_data[[col]])
  }

  # Median normalize: subtract median of each column
  for (col in intensity_cols) {
    col_median <- median(wide_data[[col]], na.rm = TRUE)
    wide_data[[col]] <- wide_data[[col]] - col_median
  }

  cat("  Log transformation and median normalization completed\n")

  # Step 2c: Baseline normalization (flexible - first available condition)
  cat("  Applying flexible baseline normalization...\n")

  # Define the condition regex pattern for extracting conditions
  # Pattern: day1_S1_SW_day1 -> day1
  cond_regex_extract <- "^(day[0-9]+)_.*$"

  # Extract condition from each column name
  cond_from_cols <- sub(cond_regex_extract, "\\1", intensity_cols, perl = TRUE)
  unique_conditions <- unique(cond_from_cols)

  # Order conditions naturally (day1, day3, day5, etc.)
  unique_conditions <- unique_conditions[order(
    as.numeric(gsub("day", "", unique_conditions))
  )]

  cat(
    "    Found",
    length(unique_conditions),
    "conditions:",
    paste(unique_conditions, collapse = ", "),
    "\n"
  )

  # For each peptide, normalize against first available condition
  n_baseline_used <- setNames(
    rep(0, length(unique_conditions)),
    unique_conditions
  )

  for (i in seq_len(nrow(wide_data))) {
    # Get intensity values for this peptide
    intensities <- as.numeric(wide_data[i, intensity_cols, with = FALSE])

    # Find first condition with non-NA values
    baseline_condition <- NULL
    baseline_value <- NA_real_

    for (cond in unique_conditions) {
      # Get columns for this condition
      cond_cols <- intensity_cols[cond_from_cols == cond]
      cond_values <- intensities[match(cond_cols, intensity_cols)]

      # Check if this condition has any non-NA values
      if (any(!is.na(cond_values))) {
        baseline_condition <- cond
        # Use median of available values in baseline condition
        baseline_value <- median(cond_values, na.rm = TRUE)
        n_baseline_used[cond] <- n_baseline_used[cond] + 1
        break
      }
    }

    # Subtract baseline from all samples if baseline exists
    if (!is.na(baseline_value)) {
      for (col in intensity_cols) {
        if (!is.na(wide_data[[col]][i])) {
          wide_data[[col]][i] <- wide_data[[col]][i] - baseline_value
        }
      }
    }
  }

  cat("    Baseline conditions used:\n")
  for (cond in names(n_baseline_used)) {
    cat(sprintf("      %s: %d peptides\n", cond, n_baseline_used[cond]))
  }

  cat(
    "  Time taken:",
    format_time_diff(as.numeric(difftime(
      Sys.time(),
      step_start,
      units = "secs"
    ))),
    "\n\n"
  )

  # Step 3: Run discover_complexoforms  on each protein
  step_start <- Sys.time()
  cat("Step 3: Running discover_complexoforms analysis per protein...\n")

  # Define the condition regex pattern for sample names
  # Pattern matches: day1_S1_SW_day1, day3_S4_SW_day3_rep2, etc.
  # We want to extract the "day" part as the condition
  cond_regex <- "^(day[0-9]+)_.*$"

  # Get unique proteins
  unique_proteins <- unique(wide_data$protein_id)
  n_proteins <- length(unique_proteins)

  cat("  Processing", n_proteins, "proteins...\n")

  # Initialize result storage
  all_results <- list()
  n_processed <- 0
  n_failed <- 0

  # Process each protein separately
  for (protein in unique_proteins) {
    # Extract data for this protein
    protein_data <- wide_data[protein_id == protein]

    # Skip proteins with only 1 peptide (need at least 2 for clustering)
    if (nrow(protein_data) < 2) {
      next
    }

    # Try to run analysis on this protein
    tryCatch(
      {
        result <- discover_complexoforms_new2(
          complex_data = protein_data,
          intensity_cols = intensity_cols,
          group_col = "protein_id",
          peptide_col = "peptide_id",
          alpha = 0.05,
          adjust_method = "BH",
          min_total_non_na_frac = 0.40,
          min_conditions = 2,
          min_reps_per_condition = 2,
          deep_split = 2,
          minClusterSize = 2,
          cond_regex = cond_regex,
          canonical_label = "dCF0",
          singleton_label = "dCF-1",
          return_condition_tests = TRUE,
          verbose = FALSE
        )

        # Store the peptide-level results
        if (
          !is.null(result$peptide_results) && nrow(result$peptide_results) > 0
        ) {
          all_results[[protein]] <- result$peptide_results
          n_processed <- n_processed + 1
        } else {
          n_failed <- n_failed + 1
        }
      },
      error = function(e) {
        # Print first few errors for debugging
        if (n_failed < 5) {
          cat("  Error processing protein", protein, ":", e$message, "\n")
        }
        n_failed <<- n_failed + 1
      }
    )
  }

  # Combine all results into a single data.table
  if (length(all_results) > 0) {
    combined_results <- rbindlist(all_results, fill = TRUE)
  } else {
    stop("No proteins were successfully processed")
  }

  cat("  Processed:", n_processed, "proteins\n")
  cat("  Failed:", n_failed, "proteins\n")
  cat("  Total peptides in results:", nrow(combined_results), "\n")
  cat(
    "  Time taken:",
    format_time_diff(as.numeric(difftime(
      Sys.time(),
      step_start,
      units = "secs"
    ))),
    "\n\n"
  )

  # Step 4: Merge with ground truth and format output
  step_start <- Sys.time()
  cat("Step 4: Formatting output and merging with ground truth...\n")

  # Merge the results with the original metadata to get ground truth labels
  # Keep only the metadata columns we need (matching other methods' output)
  metadata_cols <- wide_data[, .(
    protein_id,
    peptide_id,
    n_pep,
    n_perturbed_peptides,
    perturbed_protein,
    perturbed_peptide,
    red_fac
  )]

  # Merge with results
  final_results <- merge(
    combined_results,
    metadata_cols,
    by = c("protein_id", "peptide_id"),
    all.x = TRUE
  )

  # Ensure column name consistency with other methods
  # Rename key columns to match expected format
  if ("Gene name" %in% names(final_results)) {
    setnames(final_results, "Gene name", "protein_id", skip_absent = TRUE)
  }
  if ("Peptide" %in% names(final_results)) {
    setnames(final_results, "Peptide", "peptide_id", skip_absent = TRUE)
  }

  # Count final statistics
  proteins_processed <- n_processed
  peptides_processed <- nrow(final_results)

  cat(
    "  Processed",
    proteins_processed,
    "proteins with",
    peptides_processed,
    "peptides\n"
  )
  cat(
    "  Final results dimensions:",
    nrow(final_results),
    "rows x",
    ncol(final_results),
    "columns\n"
  )
  cat(
    "  Time taken:",
    format_time_diff(as.numeric(difftime(
      Sys.time(),
      step_start,
      units = "secs"
    ))),
    "\n\n"
  )

  # Step 5: Saving results
  step_start <- Sys.time()
  cat("Step 5: Saving results...\n")

  output_path <- paste0(
    './01-FindDCF/data/results/DCFBaseline_',
    data_id,
    "_result.feather"
  )

  # Save the final results as data.frame (not data.table for consistency)
  arrow::write_feather(as.data.frame(final_results), output_path)

  cat("  Results saved to:", output_path, "\n")
  cat(
    "  Output dimensions:",
    nrow(final_results),
    "rows x",
    ncol(final_results),
    "columns\n"
  )
  cat(
    "  Time taken:",
    format_time_diff(as.numeric(difftime(
      Sys.time(),
      step_start,
      units = "secs"
    ))),
    "\n\n"
  )

  # Calculate total time for this dataset
  end_time <- Sys.time()
  total_time_diff <- as.numeric(difftime(end_time, start_time, units = "secs"))
  total_time_taken <- format_time_diff(total_time_diff)

  # Record success (matching format from other methods)
  results_summary[[data_id]] <- list(
    status = "SUCCESS",
    total_time = total_time_taken,
    total_time_seconds = total_time_diff,
    proteins_processed = proteins_processed,
    peptides_processed = peptides_processed,
    proteins_failed = n_failed,
    input_rows = nrow(input_data),
    output_rows = nrow(final_results),
    output_file = output_path
  )

  cat("Dataset", data_id, "completed successfully in", total_time_taken, "\n")
  cat(
    "================================================================================\n"
  )
}

# Print overall summary
overall_end_time <- Sys.time()
overall_time_diff <- as.numeric(difftime(
  overall_end_time,
  overall_start_time,
  units = "secs"
))
overall_time_taken <- format_time_diff(overall_time_diff)

cat("\n")
cat(
  "================================================================================\n"
)
cat("                    FIND_dCF ANALYSIS COMPLETED\n")
cat(
  "================================================================================\n"
)
cat("\n")

successful <- sum(sapply(results_summary, function(x) x$status == "SUCCESS"))
failed <- length(results_summary) - successful

cat("Summary of processing:\n")
cat("  Total datasets:", length(results_summary), "\n")
cat("  Successful:", successful, "\n")
cat("  Failed:", failed, "\n")
cat("\n")

if (successful > 0) {
  cat("Results saved to: ./01-FindDCF/data/results/\n")
  for (data_id in names(results_summary)) {
    if (results_summary[[data_id]]$status == "SUCCESS") {
      cat(sprintf(
        "  - %s: %d proteins (%d failed), %d peptides, %s\n",
        data_id,
        results_summary[[data_id]]$proteins_processed,
        results_summary[[data_id]]$proteins_failed,
        results_summary[[data_id]]$peptides_processed,
        results_summary[[data_id]]$total_time
      ))
    }
  }
}

cat("\nTotal execution time:", overall_time_taken, "\n")
cat(
  "================================================================================\n"
)
