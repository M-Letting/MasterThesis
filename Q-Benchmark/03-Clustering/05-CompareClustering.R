# Libraries
suppressMessages(library(data.table, warn.conflicts = FALSE))

# Set working directory to Q-Benchmark
setwd(here::here("Q-Benchmark"))

# Set output directory for results
output_dir <- "03-Clustering/00-Data/05-CompareClustering"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat(paste(rep("#", 80), collapse = ""), "\n")
cat("Comparing clustering results across methods...\n")
start_time <- Sys.time()

# Define variants of cannot-link constraint to compare
cannot_link_variants <- c("corr", "rmsd", "ccc")
variant_labels <- c(corr = "Corr", rmsd = "RMSD", ccc = "CCC")

# Define membership threshold mode and value (can be set via options)
membership_threshold_mode <- getOption(
  "qbenchmark.membership_threshold_mode",
  "fixed"
)
membership_threshold_value <- getOption(
  "qbenchmark.membership_threshold_value",
  0.5
)

include_nona <- isTRUE(getOption("qbenchmark.include_nona", TRUE))

filter_nona_datasets <- function(dt) {
  if (isTRUE(include_nona)) {
    return(dt)
  }
  if (!"Dataset" %in% names(dt)) {
    return(dt)
  }
  dt[!grepl("NoNA", Dataset)]
}

get_membership_threshold <- function(n_clusters) {
  if (is.null(n_clusters) || n_clusters <= 0) {
    return(NA_real_)
  }

  mode <- tolower(as.character(membership_threshold_mode))

  if (mode %in% c("one_over_k", "1_over_k", "1/k", "dynamic")) {
    return(1 / n_clusters)
  }

  as.numeric(membership_threshold_value)
}

add_result_if_exists <- function(results_list, method_name, file_path) {
  if (file.exists(file_path)) {
    results_list[[method_name]] <- file_path
  }
  results_list
}

filter_rpc22_result_paths <- function(results_list, context_label) {
  method_names <- names(results_list)
  keep <- !grepl("RPC", method_names)

  filtered <- results_list[keep]

  if (length(filtered) == 0L) {
    stop(
      "No DCF-derived methods available for ",
      context_label,
      "."
    )
  }

  filtered
}

filter_rpc22_metric_rows <- function(data, source_label, require_rows = TRUE) {
  if (!"Method" %in% names(data)) {
    stop("Column 'Method' not found in ", source_label)
  }

  data <- data[!grepl("RPC", Method)]

  if (require_rows && nrow(data) == 0L) {
    stop("No DCF-derived methods found in ", source_label)
  }

  data
}

################################################################################
## Define functions used to evaluate ProteoMaker clustering ####################
################################################################################

# Function to calculate fraction of peptides with Proteoform_ID 1 in top cluster
# Input: PM_result - list of clustering results for a given method
# Output: data.table with fraction for protein for each dataset
PM_evaluate_fraction <- function(PM_result, method_name) {
  # Collect rows in a list and bind once at the end
  result_list <- list()

  # Helper: process one dataset list (accession -> result object)
  process_dataset <- function(dataset_result, dataset_name) {
    if (is.null(dataset_result) || length(dataset_result) == 0L) {
      return(list())
    }

    dataset_rows <- list()
    dataset_idx <- 1L

    # Loop through accessions
    for (accession in names(dataset_result)) {
      accession_result <- dataset_result[[accession]]

      # Blind results: accession_result$ClustOut$Bestcl$membership
      # DCF results: accession_result$vsclust_result$ClustOut$Bestcl$membership
      membership_matrix <- accession_result$ClustOut$Bestcl$membership
      if (is.null(membership_matrix)) {
        membership_matrix <- accession_result$vsclust_result$ClustOut$Bestcl$membership
      }

      if (is.null(membership_matrix)) {
        next
      }

      # Get clustering result for current dataset and accession
      clustering_result <- as.data.table(
        membership_matrix,
        keep.rownames = "Identifier"
      )

      # Split Identifier into components
      clustering_result[,
        c(
          "Accession",
          "Peptide",
          "Position",
          "Proteoform_ID",
          "Peptidoform"
        ) := tstrsplit(Identifier, "_")
      ]

      # Find features where Proteoform_ID contains token "1"
      features <- clustering_result[grepl("(^|\\|)1(\\||$)", Proteoform_ID)]

      membership_columns <- grep(
        "^membership",
        names(clustering_result),
        value = TRUE
      )

      if (nrow(features) == 0L || length(membership_columns) == 0L) {
        next
      }

      threshold <- get_membership_threshold(length(membership_columns))

      fractions <- rbindlist(
        lapply(membership_columns, function(cluster_col) {
          fraction <- features[get(cluster_col) >= threshold, .N] /
            nrow(features)
          data.table(Cluster = cluster_col, Fraction = fraction)
        })
      )

      if (nrow(fractions) == 0L) {
        next
      }

      top_row <- fractions[which.max(Fraction)]

      # Get constraint threshold (ex. 0.5 if dataset name has 05_ prefix)
      constraint_threshold <- NA_real_
      if (grepl("^(\\d+)_", dataset_name)) {
        constraint_threshold <- as.numeric(sub(
          "^(\\d+)_.*",
          "\\1",
          dataset_name
        )) /
          10
      }

      dataset_rows[[dataset_idx]] <- data.table(
        Dataset = dataset_name,
        Accession = accession,
        Method = method_name,
        ConstraintThreshold = constraint_threshold,
        NumClusters = length(membership_columns),
        Top_Cluster = top_row$Cluster,
        Top_Fraction = top_row$Fraction
      )
      dataset_idx <- dataset_idx + 1L
    }

    dataset_rows
  }

  # Loop through top-level elements and support both structures:
  # 1) PM_result[[dataset]][[accession]]
  # 2) PM_result[[constraint_level]][[dataset]][[accession]]
  for (level1_name in names(PM_result)) {
    level1 <- PM_result[[level1_name]]
    if (is.null(level1) || length(level1) == 0L) {
      next
    }

    first_leaf <- level1[[1]]
    first_leaf_names <- names(first_leaf)

    is_direct_dataset <- !is.null(first_leaf_names) &&
      any(c("ClustOut", "vsclust_result") %in% first_leaf_names)

    if (is_direct_dataset) {
      result_list <- c(result_list, process_dataset(level1, level1_name))
    } else {
      for (dataset_name in names(level1)) {
        result_list <- c(
          result_list,
          process_dataset(
            level1[[dataset_name]],
            paste(level1_name, dataset_name, sep = "_")
          )
        )
      }
    }
  }

  if (length(result_list) == 0L) {
    return(data.table(
      Dataset = character(),
      Accession = character(),
      Method = character(),
      ConstraintThreshold = numeric(),
      NumClusters = integer(),
      Top_Cluster = character(),
      Top_Fraction = numeric()
    ))
  }

  rbindlist(result_list, use.names = TRUE, fill = TRUE)
}

# Function to calculate total number of potential proteoforms given a clustering
# Input: PM_result - list of clustering results for a given method
# Output: data.table with total number of potential proteoforms for each dataset
PM_evaluate_total <- function(PM_result, method_name) {
  # Collect rows in a list and bind once at the end
  result_list <- list()

  # Helper: process one dataset list (accession -> result object)
  process_dataset <- function(dataset_result, dataset_name) {
    if (is.null(dataset_result) || length(dataset_result) == 0L) {
      return(list())
    }

    dataset_rows <- list()
    dataset_idx <- 1L

    # Loop through accessions
    for (accession in names(dataset_result)) {
      accession_result <- dataset_result[[accession]]

      # Blind results: accession_result$ClustOut$Bestcl$membership
      # DCF results: accession_result$vsclust_result$ClustOut$Bestcl$membership
      membership_matrix <- accession_result$ClustOut$Bestcl$membership
      if (is.null(membership_matrix)) {
        membership_matrix <- accession_result$vsclust_result$ClustOut$Bestcl$membership
      }

      if (is.null(membership_matrix)) {
        next
      }

      clustering_result <- as.data.table(
        membership_matrix,
        keep.rownames = "Identifier"
      )

      # Split Identifier into components
      clustering_result[,
        c(
          "Accession",
          "Peptide",
          "Position",
          "Proteoform_ID",
          "Peptidoform"
        ) := tstrsplit(Identifier, "_")
      ]

      membership_columns <- grep(
        "^membership",
        names(clustering_result),
        value = TRUE
      )

      if (length(membership_columns) == 0L) {
        next
      }

      # Parse peptide position interval once per accession
      pos_split <- tstrsplit(clustering_result$Position, "-", fixed = TRUE)
      clustering_result[, Start := suppressWarnings(as.numeric(pos_split[[1]]))]
      clustering_result[, Stop := suppressWarnings(as.numeric(pos_split[[2]]))]
      clustering_result <- clustering_result[!is.na(Start) & !is.na(Stop)]

      if (nrow(clustering_result) == 0L) {
        next
      }

      threshold <- get_membership_threshold(length(membership_columns))
      sum_potential_proteoforms <- 0

      # For each cluster, find overlap groups and evaluate potential proteoforms
      for (cluster_col in membership_columns) {
        cluster_dt <- clustering_result[get(cluster_col) >= threshold]
        if (nrow(cluster_dt) == 0L) {
          next
        }

        # Build overlap groups via merged/connected intervals
        setorder(cluster_dt, Start, Stop)
        overlap_group <- integer(nrow(cluster_dt))
        current_group <- 1L
        current_end <- cluster_dt$Stop[1]
        overlap_group[1] <- current_group

        if (nrow(cluster_dt) > 1L) {
          for (i in 2:nrow(cluster_dt)) {
            if (cluster_dt$Start[i] <= current_end) {
              overlap_group[i] <- current_group
              current_end <- max(current_end, cluster_dt$Stop[i])
            } else {
              current_group <- current_group + 1L
              overlap_group[i] <- current_group
              current_end <- cluster_dt$Stop[i]
            }
          }
        }

        cluster_dt[, OverlapGroup := overlap_group]

        # For each overlap group, count unique proteoform IDs across features
        group_pf_counts <- cluster_dt[,
          {
            ids <- unique(unlist(
              strsplit(Proteoform_ID, "\\|", perl = TRUE),
              use.names = FALSE
            ))
            ids <- ids[!is.na(ids) & ids != ""]
            .(n_unique_pf = length(ids))
          },
          by = OverlapGroup
        ]

        if (nrow(group_pf_counts) == 0L) {
          next
        }

        potential_proteoforms_cluster <- prod(group_pf_counts$n_unique_pf)
        sum_potential_proteoforms <- sum_potential_proteoforms +
          potential_proteoforms_cluster
      }

      # Get constraint threshold (ex. 0.5 if dataset name has 05_ prefix)
      constraint_threshold <- NA_real_
      if (grepl("^(\\d+)_", dataset_name)) {
        constraint_threshold <- as.numeric(sub(
          "^(\\d+)_.*",
          "\\1",
          dataset_name
        )) /
          10
      }

      dataset_rows[[dataset_idx]] <- data.table(
        Dataset = dataset_name,
        Accession = accession,
        Method = method_name,
        ConstraintThreshold = constraint_threshold,
        NumClusters = length(membership_columns),
        Sum_Potential_Proteoforms = sum_potential_proteoforms
      )
      dataset_idx <- dataset_idx + 1L
    }

    dataset_rows
  }

  # Loop through top-level elements and support both structures:
  # 1) PM_result[[dataset]][[accession]]
  # 2) PM_result[[constraint_level]][[dataset]][[accession]]
  for (level1_name in names(PM_result)) {
    level1 <- PM_result[[level1_name]]
    if (is.null(level1) || length(level1) == 0L) {
      next
    }

    first_leaf <- level1[[1]]
    first_leaf_names <- names(first_leaf)

    is_direct_dataset <- !is.null(first_leaf_names) &&
      any(c("ClustOut", "vsclust_result") %in% first_leaf_names)

    if (is_direct_dataset) {
      result_list <- c(result_list, process_dataset(level1, level1_name))
    } else {
      for (dataset_name in names(level1)) {
        result_list <- c(
          result_list,
          process_dataset(
            level1[[dataset_name]],
            paste(level1_name, dataset_name, sep = "_")
          )
        )
      }
    }
  }

  if (length(result_list) == 0L) {
    return(data.table(
      Dataset = character(),
      Accession = character(),
      Method = character(),
      ConstraintThreshold = numeric(),
      NumClusters = integer(),
      Sum_Potential_Proteoforms = numeric()
    ))
  }

  rbindlist(result_list, use.names = TRUE, fill = TRUE)
}

# Function to evaluate total number of potential proteoforms given a clustering
# Input: PM_result - list of clustering results for a given method
# Output: data.table with total number of potential proteoforms for each dataset
PM_evaluate_sum <- function(PM_result, method_name) {
  result_list <- list()

  process_dataset <- function(dataset_result, dataset_name) {
    if (is.null(dataset_result) || length(dataset_result) == 0L) {
      return(list())
    }

    dataset_rows <- list()
    dataset_idx <- 1L

    for (accession in names(dataset_result)) {
      accession_result <- dataset_result[[accession]]

      membership_matrix <- accession_result$ClustOut$Bestcl$membership
      if (is.null(membership_matrix)) {
        membership_matrix <- accession_result$vsclust_result$ClustOut$Bestcl$membership
      }
      if (is.null(membership_matrix)) {
        next
      }

      clustering_result <- as.data.table(
        membership_matrix,
        keep.rownames = "Identifier"
      )

      clustering_result[,
        c(
          "Accession",
          "Peptide",
          "Position",
          "Proteoform_ID",
          "Peptidoform"
        ) := tstrsplit(Identifier, "_")
      ]

      membership_columns <- grep(
        "^membership",
        names(clustering_result),
        value = TRUE
      )
      if (length(membership_columns) == 0L) {
        next
      }

      threshold <- get_membership_threshold(length(membership_columns))
      sum_unique_pf <- 0L
      nonempty_clusters <- 0L

      for (cluster_col in membership_columns) {
        cluster_dt <- clustering_result[get(cluster_col) >= threshold]
        if (nrow(cluster_dt) == 0L) {
          next
        }

        pf_tokens <- unique(unlist(
          strsplit(cluster_dt$Proteoform_ID, "\\|", perl = TRUE),
          use.names = FALSE
        ))
        pf_tokens <- pf_tokens[!is.na(pf_tokens) & pf_tokens != ""]

        nonempty_clusters <- nonempty_clusters + 1L
        sum_unique_pf <- sum_unique_pf + length(pf_tokens)
      }

      constraint_threshold <- NA_real_
      if (grepl("^(\\d+)_", dataset_name)) {
        constraint_threshold <- as.numeric(sub(
          "^(\\d+)_.*",
          "\\1",
          dataset_name
        )) /
          10
      }

      dataset_rows[[dataset_idx]] <- data.table(
        Dataset = dataset_name,
        Accession = accession,
        Method = method_name,
        ConstraintThreshold = constraint_threshold,
        NumClusters = length(membership_columns),
        NonEmptyClusters = nonempty_clusters,
        Sum_Unique_Proteoforms = sum_unique_pf
      )
      dataset_idx <- dataset_idx + 1L
    }

    dataset_rows
  }

  for (level1_name in names(PM_result)) {
    level1 <- PM_result[[level1_name]]
    if (is.null(level1) || length(level1) == 0L) {
      next
    }

    first_leaf <- level1[[1]]
    first_leaf_names <- names(first_leaf)

    is_direct_dataset <- !is.null(first_leaf_names) &&
      any(c("ClustOut", "vsclust_result") %in% first_leaf_names)

    if (is_direct_dataset) {
      result_list <- c(result_list, process_dataset(level1, level1_name))
    } else {
      for (dataset_name in names(level1)) {
        result_list <- c(
          result_list,
          process_dataset(
            level1[[dataset_name]],
            paste(level1_name, dataset_name, sep = "_")
          )
        )
      }
    }
  }

  if (length(result_list) == 0L) {
    return(data.table(
      Dataset = character(),
      Accession = character(),
      Method = character(),
      ConstraintThreshold = numeric(),
      NumClusters = integer(),
      NonEmptyClusters = integer(),
      Sum_Unique_Proteoforms = integer()
    ))
  }

  rbindlist(result_list, use.names = TRUE, fill = TRUE)
}

# Function to evaluate cluster purity for each cluster and average across clusters
# Input: PM_result - list of clustering results for a given method
# Output: data.table with mean cluster purity for each dataset and accession
PM_evaluate_purity <- function(PM_result, method_name) {
  result_list <- list()

  process_dataset <- function(dataset_result, dataset_name) {
    if (is.null(dataset_result) || length(dataset_result) == 0L) {
      return(list())
    }

    dataset_rows <- list()
    dataset_idx <- 1L

    for (accession in names(dataset_result)) {
      accession_result <- dataset_result[[accession]]

      membership_matrix <- accession_result$ClustOut$Bestcl$membership
      if (is.null(membership_matrix)) {
        membership_matrix <- accession_result$vsclust_result$ClustOut$Bestcl$membership
      }
      if (is.null(membership_matrix)) {
        next
      }

      clustering_result <- as.data.table(
        membership_matrix,
        keep.rownames = "Identifier"
      )

      clustering_result[,
        c(
          "Accession",
          "Peptide",
          "Position",
          "Proteoform_ID",
          "Peptidoform"
        ) := tstrsplit(Identifier, "_")
      ]

      membership_columns <- grep(
        "^membership",
        names(clustering_result),
        value = TRUE
      )
      if (length(membership_columns) == 0L) {
        next
      }

      threshold <- get_membership_threshold(length(membership_columns))
      cluster_purity <- numeric()
      cluster_sizes <- integer()

      for (cluster_col in membership_columns) {
        cluster_dt <- clustering_result[get(cluster_col) >= threshold]
        if (nrow(cluster_dt) == 0L) {
          next
        }

        pf_tokens <- unlist(
          strsplit(cluster_dt$Proteoform_ID, "\\|", perl = TRUE),
          use.names = FALSE
        )
        pf_tokens <- pf_tokens[!is.na(pf_tokens) & pf_tokens != ""]

        if (length(pf_tokens) == 0L) {
          next
        }

        tab <- table(pf_tokens)
        purity <- as.numeric(max(tab) / sum(tab))

        cluster_purity <- c(cluster_purity, purity)
        cluster_sizes <- c(cluster_sizes, nrow(cluster_dt))
      }

      constraint_threshold <- NA_real_
      if (grepl("^(\\d+)_", dataset_name)) {
        constraint_threshold <- as.numeric(sub(
          "^(\\d+)_.*",
          "\\1",
          dataset_name
        )) /
          10
      }

      dataset_rows[[dataset_idx]] <- data.table(
        Dataset = dataset_name,
        Accession = accession,
        Method = method_name,
        ConstraintThreshold = constraint_threshold,
        NumClusters = length(membership_columns),
        NonEmptyClusters = length(cluster_purity),
        Mean_Cluster_Purity = if (length(cluster_purity) > 0L) {
          mean(cluster_purity)
        } else {
          NA_real_
        },
        Weighted_Mean_Cluster_Purity = if (length(cluster_purity) > 0L) {
          sum(cluster_purity * cluster_sizes) / sum(cluster_sizes)
        } else {
          NA_real_
        }
      )
      dataset_idx <- dataset_idx + 1L
    }

    dataset_rows
  }

  for (level1_name in names(PM_result)) {
    level1 <- PM_result[[level1_name]]
    if (is.null(level1) || length(level1) == 0L) {
      next
    }

    first_leaf <- level1[[1]]
    first_leaf_names <- names(first_leaf)

    is_direct_dataset <- !is.null(first_leaf_names) &&
      any(c("ClustOut", "vsclust_result") %in% first_leaf_names)

    if (is_direct_dataset) {
      result_list <- c(result_list, process_dataset(level1, level1_name))
    } else {
      for (dataset_name in names(level1)) {
        result_list <- c(
          result_list,
          process_dataset(
            level1[[dataset_name]],
            paste(level1_name, dataset_name, sep = "_")
          )
        )
      }
    }
  }

  if (length(result_list) == 0L) {
    return(data.table(
      Dataset = character(),
      Accession = character(),
      Method = character(),
      ConstraintThreshold = numeric(),
      NumClusters = integer(),
      NonEmptyClusters = integer(),
      Mean_Cluster_Purity = numeric(),
      Weighted_Mean_Cluster_Purity = numeric()
    ))
  }

  rbindlist(result_list, use.names = TRUE, fill = TRUE)
}

# Function to evaluate weighted SSE per cluster using condition means
# Input: PM_result - list of clustering results for a given method
# Output: data.table with weighted SSE for each cluster/accession/dataset
PM_evaluate_weighted_sse <- function(PM_result, method_name) {
  result_list <- list()

  process_dataset <- function(dataset_result, dataset_name) {
    if (is.null(dataset_result) || length(dataset_result) == 0L) {
      return(list())
    }

    dataset_rows <- list()
    dataset_idx <- 1L

    for (accession in names(dataset_result)) {
      accession_result <- dataset_result[[accession]]

      # Blind structure
      membership_matrix <- accession_result$ClustOut$Bestcl$membership
      centers_matrix <- accession_result$ClustOut$Bestcl$centers
      original_data <- accession_result$original_data

      # DCF structure
      if (is.null(membership_matrix)) {
        membership_matrix <- accession_result$vsclust_result$ClustOut$Bestcl$membership
      }
      if (is.null(centers_matrix)) {
        centers_matrix <- accession_result$vsclust_result$ClustOut$Bestcl$centers
      }
      if (is.null(original_data)) {
        original_data <- accession_result$vsclust_result$original_data
      }

      if (
        is.null(membership_matrix) ||
          is.null(centers_matrix) ||
          is.null(original_data)
      ) {
        next
      }

      U <- as.data.table(membership_matrix, keep.rownames = "Identifier")
      membership_columns <- grep("^membership", names(U), value = TRUE)
      if (length(membership_columns) == 0L) {
        next
      }

      Xdt <- as.data.table(original_data)
      if (!"Identifier" %in% names(Xdt)) {
        next
      }

      rep_cols <- grep("^C\\d+_R\\d+$", names(Xdt), value = TRUE)
      if (length(rep_cols) == 0L) {
        next
      }

      # Mean of replicates for each condition (C1..Cn)
      condition_map <- sub("_R\\d+$", "", rep_cols)
      condition_names <- unique(condition_map)

      Xcond <- Xdt[, .(Identifier)]
      for (cond in condition_names) {
        cols <- rep_cols[condition_map == cond]
        cond_mean <- Xdt[, rowMeans(.SD, na.rm = TRUE), .SDcols = cols]
        Xcond[, (cond) := cond_mean]
        Xcond[is.nan(get(cond)), (cond) := NA_real_]
      }

      aligned <- merge(U, Xcond, by = "Identifier", all = FALSE)
      if (nrow(aligned) == 0L) {
        next
      }

      X <- as.matrix(aligned[, ..condition_names])
      storage.mode(X) <- "numeric"

      centers <- as.matrix(centers_matrix)
      storage.mode(centers) <- "numeric"

      # Align center columns to condition names from original_data
      if (!is.null(colnames(centers))) {
        common_cols <- intersect(condition_names, colnames(centers))
        if (length(common_cols) == 0L) {
          if (ncol(centers) != length(condition_names)) {
            next
          }
          centers_use <- centers
          center_cond_names <- condition_names
        } else {
          centers_use <- centers[, common_cols, drop = FALSE]
          center_cond_names <- common_cols
          X <- as.matrix(aligned[, ..center_cond_names])
          storage.mode(X) <- "numeric"
        }
      } else {
        if (ncol(centers) != length(condition_names)) {
          next
        }
        centers_use <- centers
        center_cond_names <- condition_names
      }

      constraint_threshold <- NA_real_
      if (grepl("^(\\d+)_", dataset_name)) {
        constraint_threshold <- as.numeric(sub(
          "^(\\d+)_.*",
          "\\1",
          dataset_name
        )) /
          10
      }

      total_weighted_sse <- 0
      accession_row_idx <- integer()

      for (cluster_col in membership_columns) {
        # Extract cluster index from "membership of cluster X"
        cluster_idx <- suppressWarnings(as.integer(sub(
          ".*?(\\d+)\\s*$",
          "\\1",
          cluster_col
        )))
        if (
          is.na(cluster_idx) ||
            cluster_idx < 1L ||
            cluster_idx > nrow(centers_use)
        ) {
          next
        }

        u <- aligned[[cluster_col]]
        if (all(is.na(u)) || sum(u, na.rm = TRUE) == 0) {
          next
        }

        centroid <- centers_use[cluster_idx, ]
        centroid <- as.numeric(centroid)

        diff_mat <- sweep(X, 2, centroid, "-")
        dist2 <- rowSums(diff_mat^2, na.rm = TRUE)
        n_obs <- rowSums(!is.na(diff_mat))
        dist2[n_obs == 0L] <- NA_real_

        weighted_sse <- sum(u * dist2, na.rm = TRUE)
        total_weighted_sse <- total_weighted_sse + weighted_sse

        dataset_rows[[dataset_idx]] <- data.table(
          Dataset = dataset_name,
          Accession = accession,
          Method = method_name,
          ConstraintThreshold = constraint_threshold,
          NumClusters = length(membership_columns),
          Cluster = cluster_col,
          Weighted_SSE = weighted_sse,
          Total_Weighted_SSE = NA_real_
        )
        accession_row_idx <- c(accession_row_idx, dataset_idx)
        dataset_idx <- dataset_idx + 1L
      }

      if (length(accession_row_idx) > 0L) {
        for (ii in accession_row_idx) {
          dataset_rows[[ii]][, Total_Weighted_SSE := total_weighted_sse]
        }
      }
    }

    dataset_rows
  }

  for (level1_name in names(PM_result)) {
    level1 <- PM_result[[level1_name]]
    if (is.null(level1) || length(level1) == 0L) {
      next
    }

    first_leaf <- level1[[1]]
    first_leaf_names <- names(first_leaf)

    is_direct_dataset <- !is.null(first_leaf_names) &&
      any(c("ClustOut", "vsclust_result") %in% first_leaf_names)

    if (is_direct_dataset) {
      result_list <- c(result_list, process_dataset(level1, level1_name))
    } else {
      for (dataset_name in names(level1)) {
        result_list <- c(
          result_list,
          process_dataset(
            level1[[dataset_name]],
            paste(level1_name, dataset_name, sep = "_")
          )
        )
      }
    }
  }

  if (length(result_list) == 0L) {
    return(data.table(
      Dataset = character(),
      Accession = character(),
      Method = character(),
      ConstraintThreshold = numeric(),
      NumClusters = integer(),
      Cluster = character(),
      Weighted_SSE = numeric(),
      Total_Weighted_SSE = numeric()
    ))
  }

  rbindlist(result_list, use.names = TRUE, fill = TRUE)
}

# Function to calculate number of proteoforms based only on position
# Input: PM_result - list of clustering results for a given method
# Output: data.table with total number of potential proteoforms for each dataset
PM_evaluate_position <- function(PM_result, method_name) {
  result_list <- list()

  process_dataset <- function(dataset_result, dataset_name) {
    if (is.null(dataset_result) || length(dataset_result) == 0L) {
      return(list())
    }

    dataset_rows <- list()
    dataset_idx <- 1L

    for (accession in names(dataset_result)) {
      accession_result <- dataset_result[[accession]]

      membership_matrix <- accession_result$ClustOut$Bestcl$membership
      if (is.null(membership_matrix)) {
        membership_matrix <- accession_result$vsclust_result$ClustOut$Bestcl$membership
      }

      if (is.null(membership_matrix)) {
        next
      }

      clustering_result <- as.data.table(
        membership_matrix,
        keep.rownames = "Identifier"
      )

      clustering_result[,
        c(
          "Accession",
          "Peptide",
          "Position",
          "Proteoform_ID",
          "Peptidoform"
        ) := tstrsplit(Identifier, "_")
      ]

      membership_columns <- grep(
        "^membership",
        names(clustering_result),
        value = TRUE
      )

      if (length(membership_columns) == 0L) {
        next
      }

      pos_split <- tstrsplit(clustering_result$Position, "-", fixed = TRUE)
      clustering_result[, Start := suppressWarnings(as.numeric(pos_split[[1]]))]
      clustering_result[, Stop := suppressWarnings(as.numeric(pos_split[[2]]))]
      clustering_result <- clustering_result[!is.na(Start) & !is.na(Stop)]

      if (nrow(clustering_result) == 0L) {
        next
      }

      threshold <- get_membership_threshold(length(membership_columns))
      sum_potential <- 0

      for (cluster_col in membership_columns) {
        cluster_dt <- clustering_result[get(cluster_col) >= threshold]
        if (nrow(cluster_dt) == 0L) {
          next
        }

        setorder(cluster_dt, Start, Stop)
        overlap_group <- integer(nrow(cluster_dt))
        current_group <- 1L
        current_end <- cluster_dt$Stop[1]
        overlap_group[1] <- current_group

        if (nrow(cluster_dt) > 1L) {
          for (i in 2:nrow(cluster_dt)) {
            if (cluster_dt$Start[i] <= current_end) {
              overlap_group[i] <- current_group
              current_end <- max(current_end, cluster_dt$Stop[i])
            } else {
              current_group <- current_group + 1L
              overlap_group[i] <- current_group
              current_end <- cluster_dt$Stop[i]
            }
          }
        }

        cluster_dt[, OverlapGroup := overlap_group]
        group_counts <- cluster_dt[, .N, by = OverlapGroup]

        if (nrow(group_counts) == 0L) {
          next
        }

        sum_potential <- sum_potential + prod(group_counts$N)
      }

      constraint_threshold <- NA_real_
      if (grepl("^(\\d+)_", dataset_name)) {
        constraint_threshold <- as.numeric(sub(
          "^(\\d+)_.*",
          "\\1",
          dataset_name
        )) /
          10
      }

      dataset_rows[[dataset_idx]] <- data.table(
        Dataset = dataset_name,
        Accession = accession,
        Method = method_name,
        ConstraintThreshold = constraint_threshold,
        NumClusters = length(membership_columns),
        Sum_Position_Proteoforms = sum_potential
      )
      dataset_idx <- dataset_idx + 1L
    }

    dataset_rows
  }

  for (level1_name in names(PM_result)) {
    level1 <- PM_result[[level1_name]]
    if (is.null(level1) || length(level1) == 0L) {
      next
    }

    first_leaf <- level1[[1]]
    first_leaf_names <- names(first_leaf)

    is_direct_dataset <- !is.null(first_leaf_names) &&
      any(c("ClustOut", "vsclust_result") %in% first_leaf_names)

    if (is_direct_dataset) {
      result_list <- c(result_list, process_dataset(level1, level1_name))
    } else {
      for (dataset_name in names(level1)) {
        result_list <- c(
          result_list,
          process_dataset(
            level1[[dataset_name]],
            paste(level1_name, dataset_name, sep = "_")
          )
        )
      }
    }
  }

  if (length(result_list) == 0L) {
    return(data.table(
      Dataset = character(),
      Accession = character(),
      Method = character(),
      ConstraintThreshold = numeric(),
      NumClusters = integer(),
      Sum_Position_Proteoforms = numeric()
    ))
  }

  rbindlist(result_list, use.names = TRUE, fill = TRUE)
}

# Function to calculate number of features assigned to each cluster for PM
# Input: PM_result - list of clustering results for a given method
# Output: data.table with number of assigned features for each
#         cluster/accession/dataset
PM_evaluate_assigned <- function(PM_result, method_name) {
  result_list <- list()

  process_dataset <- function(dataset_result, dataset_name) {
    if (is.null(dataset_result) || length(dataset_result) == 0L) {
      return(list())
    }

    dataset_rows <- list()
    dataset_idx <- 1L

    for (accession in names(dataset_result)) {
      accession_result <- dataset_result[[accession]]

      membership_matrix <- accession_result$ClustOut$Bestcl$membership
      if (is.null(membership_matrix)) {
        membership_matrix <- accession_result$vsclust_result$ClustOut$Bestcl$membership
      }
      if (is.null(membership_matrix)) {
        next
      }

      clustering_result <- as.data.table(
        membership_matrix,
        keep.rownames = "Identifier"
      )

      membership_columns <- grep(
        "^membership",
        names(clustering_result),
        value = TRUE
      )

      if (length(membership_columns) == 0L) {
        next
      }

      threshold <- get_membership_threshold(length(membership_columns))
      U <- as.matrix(clustering_result[, ..membership_columns])
      storage.mode(U) <- "numeric"

      assigned <- rowSums(U >= threshold, na.rm = TRUE) > 0
      n_assigned <- sum(assigned, na.rm = TRUE)

      constraint_threshold <- NA_real_
      if (grepl("^(\\d+)_", dataset_name)) {
        constraint_threshold <- as.numeric(sub(
          "^(\\d+)_.*",
          "\\1",
          dataset_name
        )) /
          10
      }

      dataset_rows[[dataset_idx]] <- data.table(
        Dataset = dataset_name,
        Accession = accession,
        Method = method_name,
        ConstraintThreshold = constraint_threshold,
        NumClusters = length(membership_columns),
        Assigned_Features = n_assigned
      )
      dataset_idx <- dataset_idx + 1L
    }

    dataset_rows
  }

  for (level1_name in names(PM_result)) {
    level1 <- PM_result[[level1_name]]
    if (is.null(level1) || length(level1) == 0L) {
      next
    }

    first_leaf <- level1[[1]]
    first_leaf_names <- names(first_leaf)

    is_direct_dataset <- !is.null(first_leaf_names) &&
      any(c("ClustOut", "vsclust_result") %in% first_leaf_names)

    if (is_direct_dataset) {
      result_list <- c(result_list, process_dataset(level1, level1_name))
    } else {
      for (dataset_name in names(level1)) {
        result_list <- c(
          result_list,
          process_dataset(
            level1[[dataset_name]],
            paste(level1_name, dataset_name, sep = "_")
          )
        )
      }
    }
  }

  if (length(result_list) == 0L) {
    return(data.table(
      Dataset = character(),
      Accession = character(),
      Method = character(),
      ConstraintThreshold = numeric(),
      NumClusters = integer(),
      Assigned_Features = integer()
    ))
  }

  rbindlist(result_list, use.names = TRUE, fill = TRUE)
}

# Function to calculate number of features assigned to each cluster for PC
# Input: PC_result - list of clustering results for a given method
# Output: data.table with number of assigned features for each complex/dataset
PC_evaluate_assigned <- function(PC_result, method_name) {
  result_list <- list()

  process_dataset <- function(
    dataset_result,
    dataset_name,
    complex_name = dataset_name
  ) {
    membership_matrix <- dataset_result$ClustOut$Bestcl$membership
    if (is.null(membership_matrix)) {
      membership_matrix <- dataset_result$vsclust_result$ClustOut$Bestcl$membership
    }

    if (is.null(membership_matrix)) {
      return(NULL)
    }

    clustering_result <- as.data.table(
      membership_matrix,
      keep.rownames = "Identifier"
    )

    membership_columns <- grep(
      "^membership",
      names(clustering_result),
      value = TRUE
    )

    if (length(membership_columns) == 0L) {
      return(NULL)
    }

    threshold <- get_membership_threshold(length(membership_columns))
    U <- as.matrix(clustering_result[, ..membership_columns])
    storage.mode(U) <- "numeric"

    assigned <- rowSums(U >= threshold, na.rm = TRUE) > 0
    n_assigned <- sum(assigned, na.rm = TRUE)

    constraint_threshold <- NA_real_
    if (grepl("^(\\d+)_", dataset_name)) {
      constraint_threshold <- as.numeric(sub(
        "^(\\d+)_.*",
        "\\1",
        dataset_name
      )) /
        10
    }

    data.table(
      Dataset = dataset_name,
      Complex = complex_name,
      Method = method_name,
      ConstraintThreshold = constraint_threshold,
      NumClusters = length(membership_columns),
      Assigned_Features = n_assigned
    )
  }

  for (level1_name in names(PC_result)) {
    level1 <- PC_result[[level1_name]]
    if (is.null(level1) || length(level1) == 0L) {
      next
    }

    level1_names <- names(level1)
    is_direct <- !is.null(level1_names) &&
      any(c("ClustOut", "vsclust_result") %in% level1_names)

    if (is_direct) {
      row <- process_dataset(level1, level1_name, level1_name)
      if (!is.null(row)) {
        result_list[[length(result_list) + 1L]] <- row
      }
    } else {
      for (complex_name in names(level1)) {
        row <- process_dataset(
          level1[[complex_name]],
          paste(level1_name, complex_name, sep = "_"),
          complex_name
        )
        if (!is.null(row)) {
          result_list[[length(result_list) + 1L]] <- row
        }
      }
    }
  }

  if (length(result_list) == 0L) {
    return(data.table(
      Dataset = character(),
      Complex = character(),
      Method = character(),
      ConstraintThreshold = numeric(),
      NumClusters = integer(),
      Assigned_Features = integer()
    ))
  }

  rbindlist(result_list, use.names = TRUE, fill = TRUE)
}

# Function to calculate number of proteoforms based only on position for Complex
# Input: PC_result - list of clustering results for a given method
# Output: data.table with total number of potential proteoforms for each dataset
PC_evaluate_position <- function(PC_result, method_name) {
  result_list <- list()

  extract_accession <- function(x) {
    m <- regexpr(
      "([OPQ][0-9][A-Z0-9]{3}[0-9](?:-[0-9]+)?|[A-NR-Z][0-9][A-Z0-9]{3}[0-9](?:-[0-9]+)?)",
      x,
      perl = TRUE
    )
    out <- ifelse(m > 0L, regmatches(x, m), NA_character_)
    out
  }

  process_dataset <- function(
    dataset_result,
    dataset_name,
    complex_name = dataset_name
  ) {
    membership_matrix <- dataset_result$ClustOut$Bestcl$membership
    if (is.null(membership_matrix)) {
      membership_matrix <- dataset_result$vsclust_result$ClustOut$Bestcl$membership
    }

    if (is.null(membership_matrix)) {
      return(NULL)
    }

    clustering_result <- as.data.table(
      membership_matrix,
      keep.rownames = "Identifier"
    )

    membership_columns <- grep(
      "^membership",
      names(clustering_result),
      value = TRUE
    )

    if (length(membership_columns) == 0L) {
      return(NULL)
    }

    clustering_result[, Complex := sub("_.*$", "", Identifier)]
    clustering_result[, PeptideLabel := sub("^[^_]+_", "", Identifier)]
    clustering_result[, Accession := extract_accession(PeptideLabel)]

    pos_match <- regexec(
      "(\\d+)\\s*[-–]\\s*(\\d+)",
      clustering_result$PeptideLabel
    )
    pos_vals <- regmatches(clustering_result$PeptideLabel, pos_match)
    clustering_result[,
      Start := suppressWarnings(as.numeric(vapply(
        pos_vals,
        function(x) if (length(x) >= 3) x[2] else NA_character_,
        character(1)
      )))
    ]
    clustering_result[,
      Stop := suppressWarnings(as.numeric(vapply(
        pos_vals,
        function(x) if (length(x) >= 3) x[3] else NA_character_,
        character(1)
      )))
    ]

    clustering_result <- clustering_result[
      !is.na(Start) & !is.na(Stop) & !is.na(Accession)
    ]

    if (nrow(clustering_result) == 0L) {
      return(NULL)
    }

    threshold <- get_membership_threshold(length(membership_columns))
    sum_potential <- 0

    for (cluster_col in membership_columns) {
      cluster_dt <- clustering_result[get(cluster_col) >= threshold]
      if (nrow(cluster_dt) == 0L) {
        next
      }

      cluster_sum <- cluster_dt[,
        {
          setorder(.SD, Start, Stop)

          overlap_group <- integer(.N)
          current_group <- 1L
          current_end <- Stop[1]
          overlap_group[1] <- current_group

          if (.N > 1L) {
            for (i in 2:.N) {
              if (Start[i] <= current_end) {
                overlap_group[i] <- current_group
                current_end <- max(current_end, Stop[i])
              } else {
                current_group <- current_group + 1L
                overlap_group[i] <- current_group
                current_end <- Stop[i]
              }
            }
          }

          dt_acc <- copy(.SD)
          dt_acc[, OverlapGroup := overlap_group]
          group_counts <- dt_acc[, .N, by = OverlapGroup]
          .(AccessionPotential = prod(group_counts$N))
        },
        by = Accession
      ][, sum(AccessionPotential)]

      sum_potential <- sum_potential + cluster_sum
    }

    constraint_threshold <- NA_real_
    if (grepl("^(\\d+)_", dataset_name)) {
      constraint_threshold <- as.numeric(sub(
        "^(\\d+)_.*",
        "\\1",
        dataset_name
      )) /
        10
    }

    data.table(
      Dataset = dataset_name,
      Complex = complex_name,
      Method = method_name,
      ConstraintThreshold = constraint_threshold,
      NumClusters = length(membership_columns),
      Sum_Position_Proteoforms = sum_potential
    )
  }

  for (level1_name in names(PC_result)) {
    level1 <- PC_result[[level1_name]]
    if (is.null(level1) || length(level1) == 0L) {
      next
    }

    level1_names <- names(level1)
    is_direct <- !is.null(level1_names) &&
      any(c("ClustOut", "vsclust_result") %in% level1_names)

    if (is_direct) {
      row <- process_dataset(level1, level1_name, level1_name)
      if (!is.null(row)) {
        result_list[[length(result_list) + 1L]] <- row
      }
    } else {
      for (complex_name in names(level1)) {
        row <- process_dataset(
          level1[[complex_name]],
          paste(level1_name, complex_name, sep = "_"),
          complex_name
        )
        if (!is.null(row)) {
          result_list[[length(result_list) + 1L]] <- row
        }
      }
    }
  }

  if (length(result_list) == 0L) {
    return(data.table(
      Dataset = character(),
      Complex = character(),
      Method = character(),
      ConstraintThreshold = numeric(),
      NumClusters = integer(),
      Sum_Position_Proteoforms = numeric()
    ))
  }

  rbindlist(result_list, use.names = TRUE, fill = TRUE)
}

################################################################################
## Define functions used to evaluate Protein Complex clustering ################
################################################################################

# Function to calcualte fraction of NMpeptides assigned to cluster with
# most NMpeptides
# Input: PC_result - list of clustering results for a given method
# Output: data.table with fraction for each dataset
PC_evaluate_fraction <- function(PC_result, method_name) {
  result_list <- list()

  process_dataset <- function(
    dataset_result,
    dataset_name,
    complex_name = dataset_name
  ) {
    # Blind: dataset_result contains ClustOut/original_data
    # DCF: dataset_result contains vsclust_result/constraint_matrix
    membership_matrix <- dataset_result$ClustOut$Bestcl$membership
    if (is.null(membership_matrix)) {
      membership_matrix <- dataset_result$vsclust_result$ClustOut$Bestcl$membership
    }

    if (is.null(membership_matrix)) {
      return(NULL)
    }

    clustering_result <- as.data.table(
      membership_matrix,
      keep.rownames = "Identifier"
    )
    membership_columns <- grep(
      "^membership",
      names(clustering_result),
      value = TRUE
    )

    if (length(membership_columns) == 0L) {
      return(NULL)
    }

    # Datatype is the last token after underscore in Identifier
    clustering_result[, Datatype := sub(".*_", "", Identifier)]
    nm_features <- clustering_result[Datatype == "NMpeptide"]

    if (nrow(nm_features) == 0L) {
      return(NULL)
    }

    threshold <- get_membership_threshold(length(membership_columns))
    fractions <- rbindlist(
      lapply(membership_columns, function(cluster_col) {
        frac <- nm_features[get(cluster_col) >= threshold, .N] /
          nrow(nm_features)
        data.table(Cluster = cluster_col, Fraction = frac)
      })
    )

    if (nrow(fractions) == 0L) {
      return(NULL)
    }

    top_row <- fractions[which.max(Fraction)]

    constraint_threshold <- NA_real_
    if (grepl("^(\\d+)_", dataset_name)) {
      constraint_threshold <- as.numeric(sub(
        "^(\\d+)_.*",
        "\\1",
        dataset_name
      )) /
        10
    }

    data.table(
      Dataset = dataset_name,
      Complex = complex_name,
      Method = method_name,
      ConstraintThreshold = constraint_threshold,
      NumClusters = length(membership_columns),
      NumNMpeptides = nrow(nm_features),
      Top_Cluster = top_row$Cluster,
      Top_Fraction = top_row$Fraction
    )
  }

  # Support both structures:
  # 1) PC_result[[complex]]
  # 2) PC_result[[constraint_level]][[complex]]
  for (level1_name in names(PC_result)) {
    level1 <- PC_result[[level1_name]]
    if (is.null(level1) || length(level1) == 0L) {
      next
    }

    level1_names <- names(level1)
    is_direct <- !is.null(level1_names) &&
      any(c("ClustOut", "vsclust_result") %in% level1_names)

    if (is_direct) {
      row <- process_dataset(level1, level1_name, level1_name)
      if (!is.null(row)) {
        result_list[[length(result_list) + 1L]] <- row
      }
    } else {
      for (complex_name in names(level1)) {
        row <- process_dataset(
          level1[[complex_name]],
          paste(level1_name, complex_name, sep = "_"),
          complex_name
        )
        if (!is.null(row)) {
          result_list[[length(result_list) + 1L]] <- row
        }
      }
    }
  }

  if (length(result_list) == 0L) {
    return(data.table(
      Dataset = character(),
      Complex = character(),
      Method = character(),
      ConstraintThreshold = numeric(),
      NumClusters = integer(),
      NumNMpeptides = integer(),
      Top_Cluster = character(),
      Top_Fraction = numeric()
    ))
  }

  rbindlist(result_list, use.names = TRUE, fill = TRUE)
}

# Function to evaluate cluster purity for Protein Complex data using Datatype
# Input: PC_result - list of clustering results for a given method
# Output: data.table with mean cluster purity for each dataset and complex
PC_evaluate_purity <- function(PC_result, method_name) {
  result_list <- list()

  process_dataset <- function(
    dataset_result,
    dataset_name,
    complex_name = dataset_name
  ) {
    membership_matrix <- dataset_result$ClustOut$Bestcl$membership
    if (is.null(membership_matrix)) {
      membership_matrix <- dataset_result$vsclust_result$ClustOut$Bestcl$membership
    }

    if (is.null(membership_matrix)) {
      return(NULL)
    }

    clustering_result <- as.data.table(
      membership_matrix,
      keep.rownames = "Identifier"
    )

    membership_columns <- grep(
      "^membership",
      names(clustering_result),
      value = TRUE
    )
    if (length(membership_columns) == 0L) {
      return(NULL)
    }

    # Datatype is the last token after underscore in Identifier
    clustering_result[, Datatype := sub(".*_", "", Identifier)]

    threshold <- get_membership_threshold(length(membership_columns))
    cluster_purity <- numeric()
    cluster_sizes <- integer()

    for (cluster_col in membership_columns) {
      cluster_dt <- clustering_result[get(cluster_col) >= threshold]
      if (nrow(cluster_dt) == 0L) {
        next
      }

      datatype_tab <- table(cluster_dt$Datatype)
      purity <- as.numeric(max(datatype_tab) / sum(datatype_tab))

      cluster_purity <- c(cluster_purity, purity)
      cluster_sizes <- c(cluster_sizes, nrow(cluster_dt))
    }

    constraint_threshold <- NA_real_
    if (grepl("^(\\d+)_", dataset_name)) {
      constraint_threshold <- as.numeric(sub(
        "^(\\d+)_.*",
        "\\1",
        dataset_name
      )) /
        10
    }

    data.table(
      Dataset = dataset_name,
      Complex = complex_name,
      Method = method_name,
      ConstraintThreshold = constraint_threshold,
      NumClusters = length(membership_columns),
      NonEmptyClusters = length(cluster_purity),
      Mean_Cluster_Purity = if (length(cluster_purity) > 0L) {
        mean(cluster_purity)
      } else {
        NA_real_
      },
      Weighted_Mean_Cluster_Purity = if (length(cluster_purity) > 0L) {
        sum(cluster_purity * cluster_sizes) / sum(cluster_sizes)
      } else {
        NA_real_
      }
    )
  }

  # Support both structures:
  # 1) PC_result[[complex]]
  # 2) PC_result[[constraint_level]][[complex]]
  for (level1_name in names(PC_result)) {
    level1 <- PC_result[[level1_name]]
    if (is.null(level1) || length(level1) == 0L) {
      next
    }

    level1_names <- names(level1)
    is_direct <- !is.null(level1_names) &&
      any(c("ClustOut", "vsclust_result") %in% level1_names)

    if (is_direct) {
      row <- process_dataset(level1, level1_name, level1_name)
      if (!is.null(row)) {
        result_list[[length(result_list) + 1L]] <- row
      }
    } else {
      for (complex_name in names(level1)) {
        row <- process_dataset(
          level1[[complex_name]],
          paste(level1_name, complex_name, sep = "_"),
          complex_name
        )
        if (!is.null(row)) {
          result_list[[length(result_list) + 1L]] <- row
        }
      }
    }
  }

  if (length(result_list) == 0L) {
    return(data.table(
      Dataset = character(),
      Complex = character(),
      Method = character(),
      ConstraintThreshold = numeric(),
      NumClusters = integer(),
      NonEmptyClusters = integer(),
      Mean_Cluster_Purity = numeric(),
      Weighted_Mean_Cluster_Purity = numeric()
    ))
  }

  rbindlist(result_list, use.names = TRUE, fill = TRUE)
}

# Function to calculate weighted SSE per cluster using condition means
# Input: PC_result - list of clustering results for a given method
# Output: data.table with weighted SSE for each cluster/accession/dataset
PC_evaluate_weighted_sse <- function(PC_result, method_name) {
  result_list <- list()

  process_dataset <- function(
    dataset_result,
    dataset_name,
    complex_name = dataset_name
  ) {
    membership_matrix <- dataset_result$ClustOut$Bestcl$membership
    centers_matrix <- dataset_result$ClustOut$Bestcl$centers
    original_data <- dataset_result$original_data

    if (is.null(membership_matrix)) {
      membership_matrix <- dataset_result$vsclust_result$ClustOut$Bestcl$membership
    }
    if (is.null(centers_matrix)) {
      centers_matrix <- dataset_result$vsclust_result$ClustOut$Bestcl$centers
    }
    if (is.null(original_data)) {
      original_data <- dataset_result$vsclust_result$original_data
    }

    if (
      is.null(membership_matrix) ||
        is.null(centers_matrix) ||
        is.null(original_data)
    ) {
      return(NULL)
    }

    U <- as.data.table(membership_matrix, keep.rownames = "Identifier")
    membership_columns <- grep("^membership", names(U), value = TRUE)
    if (length(membership_columns) == 0L) {
      return(NULL)
    }

    Xdt <- as.data.table(original_data)
    if (!"Identifier" %in% names(Xdt)) {
      return(NULL)
    }

    # Build condition means from replicate columns
    rep_cols <- setdiff(names(Xdt), "Identifier")
    rep_cols <- rep_cols[vapply(Xdt[, ..rep_cols], is.numeric, logical(1))]
    if (length(rep_cols) == 0L) {
      return(NULL)
    }

    if (all(grepl("^C\\d+_R\\d+$", rep_cols))) {
      condition_map <- sub("_R\\d+$", "", rep_cols)
    } else if (all(grepl("^B\\d+_D\\d+$", rep_cols))) {
      condition_map <- sub("^B\\d+_", "", rep_cols)
    } else if (all(grepl("^[^_]+_[^_]+$", rep_cols))) {
      condition_map <- sub("^[^_]+_", "", rep_cols)
    } else {
      condition_map <- rep_cols
    }

    condition_names <- unique(condition_map)
    Xcond <- Xdt[, .(Identifier)]
    for (cond in condition_names) {
      cols <- rep_cols[condition_map == cond]
      cond_mean <- Xdt[, rowMeans(.SD, na.rm = TRUE), .SDcols = cols]
      Xcond[, (cond) := cond_mean]
      Xcond[is.nan(get(cond)), (cond) := NA_real_]
    }

    aligned <- merge(U, Xcond, by = "Identifier", all = FALSE)
    if (nrow(aligned) == 0L) {
      return(NULL)
    }

    X <- as.matrix(aligned[, ..condition_names])
    storage.mode(X) <- "numeric"

    centers <- as.matrix(centers_matrix)
    storage.mode(centers) <- "numeric"

    # Align by common names when possible; otherwise rely on column order
    if (!is.null(colnames(centers))) {
      common_cols <- intersect(condition_names, colnames(centers))
      if (length(common_cols) > 0L) {
        centers_use <- centers[, common_cols, drop = FALSE]
        X <- as.matrix(aligned[, ..common_cols])
      } else {
        if (ncol(centers) != ncol(X)) {
          return(NULL)
        }
        centers_use <- centers
      }
    } else {
      if (ncol(centers) != ncol(X)) {
        return(NULL)
      }
      centers_use <- centers
    }

    constraint_threshold <- NA_real_
    if (grepl("^(\\d+)_", dataset_name)) {
      constraint_threshold <- as.numeric(sub(
        "^(\\d+)_.*",
        "\\1",
        dataset_name
      )) /
        10
    }

    rows <- list()
    row_idx <- 1L
    total_weighted_sse <- 0

    for (cluster_col in membership_columns) {
      cluster_idx <- suppressWarnings(as.integer(sub(
        ".*?(\\d+)\\s*$",
        "\\1",
        cluster_col
      )))
      if (
        is.na(cluster_idx) ||
          cluster_idx < 1L ||
          cluster_idx > nrow(centers_use)
      ) {
        next
      }

      u <- aligned[[cluster_col]]
      if (all(is.na(u)) || sum(u, na.rm = TRUE) == 0) {
        next
      }

      centroid <- as.numeric(centers_use[cluster_idx, ])
      diff_mat <- sweep(X, 2, centroid, "-")
      dist2 <- rowSums(diff_mat^2, na.rm = TRUE)
      dist2[rowSums(!is.na(diff_mat)) == 0L] <- NA_real_

      weighted_sse <- sum(u * dist2, na.rm = TRUE)
      total_weighted_sse <- total_weighted_sse + weighted_sse

      rows[[row_idx]] <- data.table(
        Dataset = dataset_name,
        Complex = complex_name,
        Method = method_name,
        ConstraintThreshold = constraint_threshold,
        NumClusters = length(membership_columns),
        Cluster = cluster_col,
        Weighted_SSE = weighted_sse,
        Total_Weighted_SSE = NA_real_
      )
      row_idx <- row_idx + 1L
    }

    if (length(rows) == 0L) {
      return(NULL)
    }

    out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
    out[, Total_Weighted_SSE := total_weighted_sse]
    out
  }

  # Support both structures:
  # 1) PC_result[[complex]]
  # 2) PC_result[[constraint_level]][[complex]]
  for (level1_name in names(PC_result)) {
    level1 <- PC_result[[level1_name]]
    if (is.null(level1) || length(level1) == 0L) {
      next
    }

    level1_names <- names(level1)
    is_direct <- !is.null(level1_names) &&
      any(c("ClustOut", "vsclust_result") %in% level1_names)

    if (is_direct) {
      row <- process_dataset(level1, level1_name, level1_name)
      if (!is.null(row)) {
        result_list[[length(result_list) + 1L]] <- row
      }
    } else {
      for (complex_name in names(level1)) {
        row <- process_dataset(
          level1[[complex_name]],
          paste(level1_name, complex_name, sep = "_"),
          complex_name
        )
        if (!is.null(row)) {
          result_list[[length(result_list) + 1L]] <- row
        }
      }
    }
  }

  if (length(result_list) == 0L) {
    return(data.table(
      Dataset = character(),
      Complex = character(),
      Method = character(),
      ConstraintThreshold = numeric(),
      NumClusters = integer(),
      Cluster = character(),
      Weighted_SSE = numeric(),
      Total_Weighted_SSE = numeric()
    ))
  }

  rbindlist(result_list, use.names = TRUE, fill = TRUE)
}

################################################################################
## Run evaluation and comparison of ProteoMaker clustering results #############
################################################################################

# Define result paths
clustering_results <- list()
clustering_results <- add_result_if_exists(
  clustering_results,
  "BlindDCF",
  "03-Clustering/00-Data/01-BlindVSClust/01-PM_DCF.rds"
)

for (variant in cannot_link_variants) {
  label <- variant_labels[[variant]]

  clustering_results <- add_result_if_exists(
    clustering_results,
    paste0("DCFDCF_", label),
    paste0(
      "03-Clustering/00-Data/03-DCFConstraints/01-PM_DCF_ClusteringResults_",
      variant,
      ".rds"
    )
  )
  clustering_results <- add_result_if_exists(
    clustering_results,
    paste0("DCFRPC2-2_", label),
    paste0(
      "03-Clustering/00-Data/03-DCFConstraints/01-PM_RPC2-2_DCF_ClusteringResults_",
      variant,
      ".rds"
    )
  )
  clustering_results <- add_result_if_exists(
    clustering_results,
    paste0("VSClustDCF_", label),
    paste0(
      "03-Clustering/00-Data/04-VSClustConstraints/01-PM_DCF_VSClust_ClusteringResults_",
      variant,
      ".rds"
    )
  )
  clustering_results <- add_result_if_exists(
    clustering_results,
    paste0("VSClustRPC_", label),
    paste0(
      "03-Clustering/00-Data/04-VSClustConstraints/01-PM_RPC_VSClust_ClusteringResults_",
      variant,
      ".rds"
    )
  )
  clustering_results <- add_result_if_exists(
    clustering_results,
    paste0("VSClustRPC2-2_", label),
    paste0(
      "03-Clustering/00-Data/04-VSClustConstraints/01-PM_RPC2-2_VSClust_ClusteringResults_",
      variant,
      ".rds"
    )
  )
}

clustering_results <- filter_rpc22_result_paths(
  clustering_results,
  "ProteoMaker comparison"
)

#######################################################################
# Evaluate fraction of Proteoform_ID 1 in top cluster for each method #
#######################################################################
fraction_results <- list()
for (method_name in names(clustering_results)) {
  method_path <- clustering_results[[method_name]]
  method_results <- readRDS(method_path)
  fraction_results[[method_name]] <- PM_evaluate_fraction(
    method_results,
    method_name
  )
}

fraction_results_dt <- rbindlist(
  fraction_results,
  use.names = TRUE,
  fill = TRUE
)
fraction_results_dt <- filter_nona_datasets(fraction_results_dt)

# Save results
output_path <- file.path(
  output_dir,
  paste0(
    "ProteoMakerFractionEvaluation_",
    membership_threshold_mode,
    ".csv"
  )
)
fwrite(fraction_results_dt, output_path)

##################################################################
# Evaluate total number of potential proteoforms for each method #
###################################################################
total_results <- list()
for (method_name in names(clustering_results)) {
  method_path <- clustering_results[[method_name]]
  method_results <- readRDS(method_path)
  total_results[[method_name]] <- PM_evaluate_total(
    method_results,
    method_name
  )
}

total_results_dt <- rbindlist(
  total_results,
  use.names = TRUE,
  fill = TRUE
)
total_results_dt <- filter_nona_datasets(total_results_dt)

# Save results
output_path <- file.path(
  output_dir,
  paste0(
    "ProteoMakerTotalPotentialProteoforms_",
    membership_threshold_mode,
    ".csv"
  )
)
fwrite(total_results_dt, output_path)

##############################################################
# Evaluate position-based potential proteoforms for PM data ##
##############################################################
position_results <- list()
for (method_name in names(clustering_results)) {
  method_path <- clustering_results[[method_name]]
  method_results <- readRDS(method_path)
  position_results[[method_name]] <- PM_evaluate_position(
    method_results,
    method_name
  )
}

position_results_dt <- rbindlist(
  position_results,
  use.names = TRUE,
  fill = TRUE
)
position_results_dt <- filter_nona_datasets(position_results_dt)

output_path <- file.path(
  output_dir,
  paste0(
    "ProteoMakerPositionPotentialProteoforms_",
    membership_threshold_mode,
    ".csv"
  )
)
fwrite(position_results_dt, output_path)

######################################################################
# Evaluate sum of unique proteoforms across clusters for each method #
######################################################################
sum_results <- list()
for (method_name in names(clustering_results)) {
  method_path <- clustering_results[[method_name]]
  method_results <- readRDS(method_path)
  sum_results[[method_name]] <- PM_evaluate_sum(
    method_results,
    method_name
  )
}

sum_results_dt <- rbindlist(
  sum_results,
  use.names = TRUE,
  fill = TRUE
)
sum_results_dt <- filter_nona_datasets(sum_results_dt)

# Save results
output_path <- file.path(
  output_dir,
  paste0("ProteoMakerSumUniqueProteoforms_", membership_threshold_mode, ".csv")
)
fwrite(sum_results_dt, output_path)

###########################################
# Evaluate cluster purity for each method #
###########################################
purity_results <- list()
for (method_name in names(clustering_results)) {
  method_path <- clustering_results[[method_name]]
  method_results <- readRDS(method_path)
  purity_results[[method_name]] <- PM_evaluate_purity(
    method_results,
    method_name
  )
}

purity_results_dt <- rbindlist(
  purity_results,
  use.names = TRUE,
  fill = TRUE
)
purity_results_dt <- filter_nona_datasets(purity_results_dt)

# Save results
output_path <- file.path(
  output_dir,
  paste0(
    "ProteoMakerClusterPurityEvaluation_",
    membership_threshold_mode,
    ".csv"
  )
)
fwrite(purity_results_dt, output_path)

#########################################
# Evaluate weighted SSE for each method #
#########################################
sse_results <- list()
for (method_name in names(clustering_results)) {
  method_path <- clustering_results[[method_name]]
  method_results <- readRDS(method_path)
  sse_results[[method_name]] <- PM_evaluate_weighted_sse(
    method_results,
    method_name
  )
}

sse_results_dt <- rbindlist(
  sse_results,
  use.names = TRUE,
  fill = TRUE
)
sse_results_dt <- filter_nona_datasets(sse_results_dt)

# Save results
output_path <- file.path(
  output_dir,
  paste0(
    "ProteoMakerWeightedSSEEvaluation_",
    membership_threshold_mode,
    ".csv"
  )
)
fwrite(sse_results_dt, output_path)

###############################################
# Evaluate assigned features for each method #
###############################################
assigned_results <- list()
for (method_name in names(clustering_results)) {
  method_path <- clustering_results[[method_name]]
  method_results <- readRDS(method_path)
  assigned_results[[method_name]] <- PM_evaluate_assigned(
    method_results,
    method_name
  )
}

assigned_results_dt <- rbindlist(
  assigned_results,
  use.names = TRUE,
  fill = TRUE
)
assigned_results_dt <- filter_nona_datasets(assigned_results_dt)

# Save results
output_path <- file.path(
  output_dir,
  paste0(
    "ProteoMakerAssignedEvaluation_",
    membership_threshold_mode,
    ".csv"
  )
)
fwrite(assigned_results_dt, output_path)

################################################################################
## Run evaluation and comparison of Protein Complex clustering results #########
################################################################################

# Define results path
clustering_results <- list()
clustering_results <- add_result_if_exists(
  clustering_results,
  "BlindDCF",
  "03-Clustering/00-Data/01-BlindVSClust/01-Complex_DCF.rds"
)

for (variant in cannot_link_variants) {
  label <- variant_labels[[variant]]

  clustering_results <- add_result_if_exists(
    clustering_results,
    paste0("DCFDCF_", label),
    paste0(
      "03-Clustering/00-Data/03-DCFConstraints/01-PC_DCF_ClusteringResults_",
      variant,
      ".rds"
    )
  )
  clustering_results <- add_result_if_exists(
    clustering_results,
    paste0("DCFRPC_", label),
    paste0(
      "03-Clustering/00-Data/03-DCFConstraints/01-PC_RPC_DCF_ClusteringResults_",
      variant,
      ".rds"
    )
  )
  clustering_results <- add_result_if_exists(
    clustering_results,
    paste0("DCFRPC2-2_", label),
    paste0(
      "03-Clustering/00-Data/03-DCFConstraints/01-PC_RPC2_2_DCF_ClusteringResults_",
      variant,
      ".rds"
    )
  )
  clustering_results <- add_result_if_exists(
    clustering_results,
    paste0("VSClustDCF_", label),
    paste0(
      "03-Clustering/00-Data/04-VSClustConstraints/01-PC_DCF_VSClust_ClusteringResults_",
      variant,
      ".rds"
    )
  )
  clustering_results <- add_result_if_exists(
    clustering_results,
    paste0("VSClustRPC_", label),
    paste0(
      "03-Clustering/00-Data/04-VSClustConstraints/01-PC_RPC_VSClust_ClusteringResults_",
      variant,
      ".rds"
    )
  )
  clustering_results <- add_result_if_exists(
    clustering_results,
    paste0("VSClustRPC2-2_", label),
    paste0(
      "03-Clustering/00-Data/04-VSClustConstraints/01-PC_RPC2_2_VSClust_ClusteringResults_",
      variant,
      ".rds"
    )
  )
}

clustering_results <- filter_rpc22_result_paths(
  clustering_results,
  "Complex comparison"
)

##################################################################
# Evaluate fraction of NMpeptides in top cluster for each method #
##################################################################
fraction_results <- list()
for (method_name in names(clustering_results)) {
  method_path <- clustering_results[[method_name]]
  method_results <- readRDS(method_path)
  fraction_results[[method_name]] <- PC_evaluate_fraction(
    method_results,
    method_name
  )
}

fraction_results_dt <- rbindlist(
  fraction_results,
  use.names = TRUE,
  fill = TRUE
)
fraction_results_dt <- filter_nona_datasets(fraction_results_dt)

# Save results
output_path <- file.path(
  output_dir,
  paste0("ComplexFractionEvaluation_", membership_threshold_mode, ".csv")
)
fwrite(fraction_results_dt, output_path)

###########################################
# Evaluate cluster purity for each method #
###########################################
purity_results <- list()
for (method_name in names(clustering_results)) {
  method_path <- clustering_results[[method_name]]
  method_results <- readRDS(method_path)
  purity_results[[method_name]] <- PC_evaluate_purity(
    method_results,
    method_name
  )
}

purity_results_dt <- rbindlist(
  purity_results,
  use.names = TRUE,
  fill = TRUE
)
purity_results_dt <- filter_nona_datasets(purity_results_dt)

# Save results
output_path <- file.path(
  output_dir,
  paste0("ComplexClusterPurityEvaluation_", membership_threshold_mode, ".csv")
)
fwrite(purity_results_dt, output_path)

#########################################
# Evaluate weighted SSE for each method #
#########################################
sse_results <- list()
for (method_name in names(clustering_results)) {
  method_path <- clustering_results[[method_name]]
  method_results <- readRDS(method_path)
  sse_results[[method_name]] <- PC_evaluate_weighted_sse(
    method_results,
    method_name
  )
}

sse_results_dt <- rbindlist(
  sse_results,
  use.names = TRUE,
  fill = TRUE
)
sse_results_dt <- filter_nona_datasets(sse_results_dt)

# Save results
output_path <- file.path(
  output_dir,
  paste0("ComplexWeightedSSEEvaluation_", membership_threshold_mode, ".csv")
)
fwrite(sse_results_dt, output_path)

###############################################
# Evaluate assigned features for each method #
###############################################
assigned_results <- list()
for (method_name in names(clustering_results)) {
  method_path <- clustering_results[[method_name]]
  method_results <- readRDS(method_path)
  assigned_results[[method_name]] <- PC_evaluate_assigned(
    method_results,
    method_name
  )
}

assigned_results_dt <- rbindlist(
  assigned_results,
  use.names = TRUE,
  fill = TRUE
)
assigned_results_dt <- filter_nona_datasets(assigned_results_dt)

# Save results
output_path <- file.path(
  output_dir,
  paste0("ComplexAssignedEvaluation_", membership_threshold_mode, ".csv")
)
fwrite(assigned_results_dt, output_path)

###################################################################
# Evaluate position-based potential proteoforms for Complex data ###
###################################################################
position_results <- list()
for (method_name in names(clustering_results)) {
  method_path <- clustering_results[[method_name]]
  method_results <- readRDS(method_path)
  position_results[[method_name]] <- PC_evaluate_position(
    method_results,
    method_name
  )
}

position_results_dt <- rbindlist(
  position_results,
  use.names = TRUE,
  fill = TRUE
)
position_results_dt <- filter_nona_datasets(position_results_dt)

output_path <- file.path(
  output_dir,
  paste0(
    "ComplexPositionPotentialProteoforms_",
    membership_threshold_mode,
    ".csv"
  )
)
fwrite(position_results_dt, output_path)

end_time <- Sys.time()
total_time <- round(end_time - start_time, 2)
cat("Finished comparing clustering methods.\n")
cat("Total time taken:", total_time, "seconds \n")
cat(paste(rep("#", 80), collapse = ""), "\n")
