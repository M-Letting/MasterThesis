library(data.table)

# Lighten a vector of hex colors by mixing toward white.
# factor = 0 leaves the color unchanged; factor = 1 produces pure white.
# Qualitative palette for cluster labels. Uses colorspace::qualitative_hcl()
# which generates perceptually uniform, maximally-separated hues at any k.
cluster_palette <- function(k) {
  colorspace::qualitative_hcl(k, palette = "Dark 3")
}

datatype_palette <- function() {
  c(
    "NM"        = "#00441b",
    "NMpeptide" = "#02d054ff",
    "RmCys"     = "#006d2c",
    "Deglyco"   = "#99d38d",
    "Phospho"   = "#be4e16",
    "FreeCys"   = "#ff8000",
    "LysAc"     = "#fec389"
  )
}

lighten_hex <- function(hex, factor = 0.55) {
  mat <- col2rgb(hex) / 255
  lightened <- mat + (1 - mat) * factor
  apply(lightened, 2, function(x) rgb(x[1], x[2], x[3]))
}

# Extract a "start-end" residue range string from a position annotation such as
# "P54289 [691-704]" or "Q13936 [2191-2201]".
parse_chimerax_range <- function(position_str) {
  m <- regmatches(position_str, regexpr("\\d+-\\d+", position_str))
  if (length(m) == 0L) NA_character_ else m
}

# Parse a range string "start-end" and return the integer midpoint.
range_midpoint <- function(rng) {
  parts <- as.integer(strsplit(rng, "-")[[1]])
  as.integer(round(mean(parts)))
}

# Parse Cα coordinates from an mmCIF structure file.
# Returns a data.table with columns: chain, resno, x, y, z (Cα only).
read_ca_coords_cif <- function(cif_path) {
  lines <- readLines(cif_path, warn = FALSE)

  loop_positions <- which(trimws(lines) == "loop_")
  atom_loop <- NA_integer_
  for (lp in loop_positions) {
    if (
      lp + 1L <= length(lines) &&
        startsWith(trimws(lines[lp + 1L]), "_atom_site.")
    ) {
      atom_loop <- lp
      break
    }
  }
  if (is.na(atom_loop)) {
    stop("No _atom_site loop found in: ", cif_path)
  }

  ci <- atom_loop + 1L
  col_names <- character(0)
  while (ci <= length(lines) && startsWith(trimws(lines[ci]), "_atom_site.")) {
    col_names <- c(col_names, sub("^_atom_site\\.", "", trimws(lines[ci])))
    ci <- ci + 1L
  }

  atom_col <- match("label_atom_id", col_names)
  chain_col <- match("auth_asym_id", col_names)
  resno_col <- match("auth_seq_id", col_names)
  x_col <- match("Cartn_x", col_names)
  y_col <- match("Cartn_y", col_names)
  z_col <- match("Cartn_z", col_names)
  if (is.na(chain_col)) {
    chain_col <- match("label_asym_id", col_names)
  }
  if (is.na(resno_col)) {
    resno_col <- match("label_seq_id", col_names)
  }

  if (anyNA(c(atom_col, chain_col, resno_col, x_col, y_col, z_col))) {
    stop("Required _atom_site columns not found in: ", cif_path)
  }

  raw <- lines[ci:length(lines)]
  raw <- raw[startsWith(raw, "ATOM ") | startsWith(raw, "HETATM ")]
  if (length(raw) == 0L) {
    stop("No ATOM/HETATM records in: ", cif_path)
  }
  raw <- gsub("\\s+", " ", trimws(raw))

  dt <- data.table::fread(
    text = paste(raw, collapse = "\n"),
    header = FALSE,
    sep = " ",
    fill = TRUE,
    quote = ""
  )

  ca <- dt[dt[[atom_col]] == "CA"]
  data.table::data.table(
    chain = as.character(ca[[chain_col]]),
    resno = suppressWarnings(as.integer(as.character(ca[[resno_col]]))),
    x = as.numeric(ca[[x_col]]),
    y = as.numeric(ca[[y_col]]),
    z = as.numeric(ca[[z_col]])
  )
}

# Build a pairwise-distance table from a sites data.table that has columns:
# Gene, Accession, Position, Datatype, mean_res, x, y, z.
# sites_dt columns: Gene, Accession, Position, Datatype, mean_res, x, y, z,
# and optionally coord_source ("reference_frame" | "alphafold").
# When coord_source is present a SameCoordFrame column is added to the output:
#   TRUE  -> both sites share a coordinate frame (distances physically meaningful)
#   FALSE -> sites are from different unaligned frames (distances unreliable)
.ptm_pairwise_distances <- function(sites_dt) {
  n <- nrow(sites_dt)
  has_source <- "coord_source" %in% names(sites_dt)
  rows <- vector("list", n * (n - 1L) / 2L)
  idx <- 0L
  for (i in seq_len(n - 1L)) {
    for (j in seq(i + 1L, n)) {
      d <- sqrt(
        (sites_dt$x[i] - sites_dt$x[j])^2 +
          (sites_dt$y[i] - sites_dt$y[j])^2 +
          (sites_dt$z[i] - sites_dt$z[j])^2
      )
      idx <- idx + 1L
      row <- list(
        Gene1 = sites_dt$Gene[i],
        Accession1 = sites_dt$Accession[i],
        Position1 = sites_dt$Position[i],
        Datatype1 = sites_dt$Datatype[i],
        Residue1 = sites_dt$mean_res[i],
        Gene2 = sites_dt$Gene[j],
        Accession2 = sites_dt$Accession[j],
        Position2 = sites_dt$Position[j],
        Datatype2 = sites_dt$Datatype[j],
        Residue2 = sites_dt$mean_res[j],
        Distance_Angstrom = round(d, 2L)
      )
      if (has_source) {
        s1 <- sites_dt$coord_source[i]
        s2 <- sites_dt$coord_source[j]
        row$SameCoordFrame <-
          (s1 == "reference_frame" & s2 == "reference_frame") |
          (s1 == "alphafold" &
            s2 == "alphafold" &
            sites_dt$Accession[i] == sites_dt$Accession[j])
      }
      rows[[idx]] <- row
    }
  }
  data.table::rbindlist(rows, fill = TRUE)
}

#' Generate ChimeraX .cxc scripts for each complexoform cluster
#'
#' For each cluster in the cluster summary one `.cxc` script is written.
#' The structure chains are colored using a lightened version of
#' `protein_chain_palette()` (defined in ComplexoformStructures.R), and
#' peptides whose membership value meets the threshold are overlaid in their
#' datatype color from `datatype_palette()`.
#'
#' @param complex_obj List. One of the complex definition objects from
#'   `ComplexoformStructures.R`. Must contain:
#'   \describe{
#'     \item{`name`}{Character. Human-readable complex name (used in filenames).}
#'     \item{`cluster_summary`}{Character. Path to the cluster summary CSV.}
#'     \item{`structure_path`}{Character. Path to the structure file (.cif/.pdb).}
#'     \item{`chain_map`}{Named character vector mapping chain IDs (names) to
#'       UniProt accessions (values) as they appear in the Accession column of
#'       the cluster summary.}
#'   }
#' @param output_dir Character. Directory where `.cxc` files are written.
#'   Created if it does not exist.
#' @param image_dir Character. Directory where rendered PNG images are saved.
#'   Defaults to `output_dir` when `NULL`.
#' @param membership_threshold Either `"adaptive"` (threshold = 1/k, where k is
#'   the number of clusters) or a single numeric value (e.g. `0.5`).
#' @param chain_lighten_factor Numeric in \[0, 1\]. How far to mix each chain
#'   color toward white. `0` = original batlow color, `1` = pure white.
#'   Default `0.55` gives a clearly toned-down background.
#' @param save_image Logical. If `TRUE`, each script ends with a `save` command
#'   that renders a PNG to `image_dir`.
#' @param image_width,image_height Integer. PNG dimensions in pixels.
#'
#' @return Invisibly returns a character vector of the written `.cxc` file paths.
#'
#' @export
make_cxc <- function(
  complex_obj,
  output_dir,
  image_dir = NULL,
  membership_threshold = "adaptive",
  chain_lighten_factor = 0.55,
  save_image = TRUE,
  image_width = 2400,
  image_height = 2400
) {
  cs <- data.table::fread(complex_obj$cluster_summary)

  use_cluster_id <- "ClusterID" %in% names(cs)
  if (use_cluster_id) {
    cluster_ids <- sort(unique(cs$ClusterID))
    k <- length(cluster_ids)
  } else {
    mem_cols <- grep("^membership of cluster", names(cs), value = TRUE)
    k <- length(mem_cols)
    if (k == 0L) {
      stop("No membership columns found in cluster summary.")
    }
  }

  threshold <- if (!use_cluster_id) {
    if (identical(membership_threshold, "adaptive")) {
      1 / k
    } else {
      as.numeric(membership_threshold)
    }
  } else {
    NA_real_
  }

  # Chain colors: batlow toned down toward white
  chain_ids <- names(complex_obj$chain_map)
  chain_colors <- lighten_hex(
    protein_chain_palette(length(chain_ids)),
    factor = chain_lighten_factor
  )
  names(chain_colors) <- chain_ids

  dt_pal <- datatype_palette()

  # Reverse lookup: UniProt accession -> chain ID
  acc_to_chain <- setNames(
    names(complex_obj$chain_map),
    complex_obj$chain_map
  )

  struct_path <- normalizePath(complex_obj$structure_path)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  out_dir_abs <- normalizePath(output_dir)

  if (is.null(image_dir)) {
    image_dir <- output_dir
  }
  if (!dir.exists(image_dir)) {
    dir.create(image_dir, recursive = TRUE)
  }
  img_dir_abs <- normalizePath(image_dir)

  # Sanitize complex name for use in filenames
  name_clean <- gsub("[^A-Za-z0-9]+", "_", complex_obj$name)

  written_files <- character(k)

  for (cl in seq_len(k)) {
    cl_id <- if (use_cluster_id) cluster_ids[cl] else cl
    members <- if (use_cluster_id) {
      cs[ClusterID == cl_id & !Datatype %in% c("NM", "NMpeptide")]
    } else {
      mem_col <- paste("membership of cluster", cl)
      cs[get(mem_col) >= threshold & !Datatype %in% c("NM", "NMpeptide")]
    }

    lines <- c(
      paste0("# Complexoform group ", cl_id, " – ", complex_obj$name),
      if (use_cluster_id) {
        "# Cluster assignment: ClusterID"
      } else {
        paste0("# Membership threshold: ", round(threshold, 4))
      },
      "",
      paste0("open ", struct_path),
      ""
    )

    # Cartoon display only
    lines <- c(lines, "show cartoons", "hide atoms", "hide surface", "")

    # Background chain coloring (toned down) + nearly transparent
    for (ch in chain_ids) {
      lines <- c(lines, paste0("color #1/", ch, " ", chain_colors[ch]))
      lines <- c(lines, paste0("transparency #1/", ch, " 75 cartoons"))
    }
    lines <- c(lines, "")

    # Member peptides colored by datatype
    if (nrow(members) > 0L) {
      lines <- c(lines, "# Member peptides")
      for (i in seq_len(nrow(members))) {
        row <- members[i]
        chain <- acc_to_chain[row$Accession]
        rng <- parse_chimerax_range(row$Position)
        dt_color <- dt_pal[row$Datatype]

        if (is.na(chain) || is.na(rng) || is.na(dt_color)) {
          next
        }

        lines <- c(
          lines,
          paste0("color #1/", chain, ":", rng, " ", dt_color)
        )
      }
      lines <- c(lines, "")
    }

    # Visual polish
    lines <- c(
      lines,
      "lighting simple",
      "set bgColor black",
      "graphics silhouettes false",
      ""
    )

    if (save_image) {
      img_path <- file.path(
        img_dir_abs,
        paste0(name_clean, "_cluster_", cl_id, ".png")
      )
      lines <- c(
        lines,
        paste0(
          "save ",
          img_path,
          " width ",
          image_width,
          " height ",
          image_height,
          " supersample 3"
        ),
        ""
      )
    }

    lines <- c(lines)

    cxc_path <- file.path(
      out_dir_abs,
      paste0(name_clean, "_cluster_", cl_id, ".cxc")
    )
    writeLines(lines, cxc_path)
    written_files[cl] <- cxc_path
    message("Written: ", cxc_path)
  }

  invisible(written_files)
}

#' Generate ChimeraX .cxc scripts for complexes without a combined structure
#'
#' Like `make_cxc()`, but for complexes where no multi-chain structure exists.
#' Each subunit is fetched individually from the AlphaFold database. If the
#' `complex_obj` contains a `reference_pdb` entry, that structure is opened
#' first (model #1) and each AlphaFold model is superimposed onto its
#' corresponding chain via `matchmaker` before being colored. Subunits absent
#' from the reference float freely and are left unpositioned.
#'
#' Requires an active internet connection when the generated script is run
#' inside ChimeraX (`alphafold fetch` downloads on demand).
#'
#' @param complex_obj List. Must contain:
#'   \describe{
#'     \item{`name`}{Character. Human-readable complex name.}
#'     \item{`cluster_summary`}{Character. Path to the cluster summary CSV.}
#'     \item{`subunit_af_accessions`}{Named character vector. Names = gene
#'       symbols, values = UniProt accessions. Order determines model numbering
#'       (offset by 1 when a reference is present).}
#'     \item{`reference_pdb`}{Character (optional). PDB ID of a reference
#'       structure used to spatially position the subunits.}
#'     \item{`reference_chain_map`}{Named character vector (required when
#'       `reference_pdb` is set). Names = UniProt accessions from
#'       `subunit_af_accessions`, values = chain IDs in the reference structure.
#'       Only subunits with an entry here are matchmakered; the rest float.}
#'   }
#' @inheritParams make_cxc
#'
#' @return Invisibly returns a character vector of the written `.cxc` file paths.
#'
#' @export
make_cxc_af_subunits <- function(
  complex_obj,
  output_dir,
  image_dir = NULL,
  membership_threshold = "adaptive",
  chain_lighten_factor = 0.55,
  save_image = TRUE,
  image_width = 2400,
  image_height = 2400
) {
  cs <- data.table::fread(complex_obj$cluster_summary)

  use_cluster_id <- "ClusterID" %in% names(cs)
  if (use_cluster_id) {
    cluster_ids <- sort(unique(cs$ClusterID))
    k <- length(cluster_ids)
  } else {
    mem_cols <- grep("^membership of cluster", names(cs), value = TRUE)
    k <- length(mem_cols)
    if (k == 0L) {
      stop("No membership columns found in cluster summary.")
    }
  }

  threshold <- if (!use_cluster_id) {
    if (identical(membership_threshold, "adaptive")) {
      1 / k
    } else {
      as.numeric(membership_threshold)
    }
  } else {
    NA_real_
  }

  subunit_accs <- complex_obj$subunit_af_accessions # named: gene -> accession
  n <- length(subunit_accs)

  has_ref <- !is.null(complex_obj$reference_pdb) &&
    nzchar(complex_obj$reference_pdb)
  # Reference occupies #1 when present; AF models follow sequentially after.
  model_offset <- if (has_ref) 1L else 0L
  acc_to_model <- setNames(seq_len(n) + model_offset, subunit_accs)

  subunit_colors <- lighten_hex(
    protein_chain_palette(n),
    factor = chain_lighten_factor
  )
  names(subunit_colors) <- subunit_accs

  dt_pal <- datatype_palette()

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  out_dir_abs <- normalizePath(output_dir)
  if (is.null(image_dir)) {
    image_dir <- output_dir
  }
  if (!dir.exists(image_dir)) {
    dir.create(image_dir, recursive = TRUE)
  }
  img_dir_abs <- normalizePath(image_dir)

  name_clean <- gsub("[^A-Za-z0-9]+", "_", complex_obj$name)
  written_files <- character(k)

  for (cl in seq_len(k)) {
    cl_id <- if (use_cluster_id) cluster_ids[cl] else cl
    members <- if (use_cluster_id) {
      cs[ClusterID == cl_id & !Datatype %in% c("NM", "NMpeptide")]
    } else {
      mem_col <- paste("membership of cluster", cl)
      cs[get(mem_col) >= threshold & !Datatype %in% c("NM", "NMpeptide")]
    }

    lines <- c(
      paste0("# Complexoform group ", cl_id, " – ", complex_obj$name),
      if (use_cluster_id) {
        "# Cluster assignment: ClusterID"
      } else {
        paste0("# Membership threshold: ", round(threshold, 4))
      },
      ""
    )

    # Open reference scaffold if provided
    if (has_ref) {
      lines <- c(
        lines,
        paste0(
          "# Reference scaffold: PDB ",
          complex_obj$reference_pdb,
          " (model #1)"
        ),
        paste0("open ", complex_obj$reference_pdb),
        ""
      )
    }

    # Fetch AlphaFold structures
    lines <- c(lines, "# Fetch individual AlphaFold subunit structures")
    for (i in seq_len(n)) {
      acc <- subunit_accs[i]
      gene <- names(subunit_accs)[i]
      mnum <- i + model_offset
      lines <- c(lines, paste0("alphafold fetch ", acc))
    }
    lines <- c(lines, "")

    # Cartoon display only
    lines <- c(lines, "show cartoons", "hide atoms", "hide surface", "")

    # Matchmaker each AF model onto its reference chain
    if (has_ref) {
      ref_map <- complex_obj$reference_chain_map # acc -> chain in reference
      lines <- c(lines, "# Superimpose AlphaFold models onto reference chains")
      for (i in seq_len(n)) {
        acc <- subunit_accs[i]
        gene <- names(subunit_accs)[i]
        mnum <- i + model_offset
        ref_chain <- ref_map[acc]
        if (!is.na(ref_chain)) {
          lines <- c(
            lines,
            paste0("matchmaker #", mnum, "/A to #1/", ref_chain)
          )
        } else {
          lines <- c(
            lines,
            paste0("# ", gene, " has no reference chain — left unpositioned")
          )
        }
      }
      lines <- c(lines, "hide #1 models", "")
    }

    # Background coloring: one batlow color per subunit model
    lines <- c(lines, "# Background coloring per subunit")
    for (i in seq_len(n)) {
      acc <- subunit_accs[i]
      mnum <- acc_to_model[acc]
      lines <- c(lines, paste0("color #", mnum, "/A ", subunit_colors[acc]))
      lines <- c(lines, paste0("transparency #", mnum, "/A 75 cartoons"))
    }
    lines <- c(lines, "")

    # Member peptides colored by datatype
    if (nrow(members) > 0L) {
      lines <- c(lines, "# Member peptides")
      for (i in seq_len(nrow(members))) {
        row <- members[i]
        mnum <- acc_to_model[row$Accession]
        rng <- parse_chimerax_range(row$Position)
        dt_color <- dt_pal[row$Datatype]
        if (is.na(mnum) || is.na(rng) || is.na(dt_color)) {
          next
        }
        lines <- c(lines, paste0("color #", mnum, "/A:", rng, " ", dt_color))
        lines <- c(
          lines,
          paste0("transparency #", mnum, "/A:", rng, " 0 cartoons")
        ) # ensure members are fully opaque
      }
      lines <- c(lines, "")
    }

    # Visual polish
    lines <- c(
      lines,
      "lighting simple",
      "set bgColor black",
      "graphics silhouettes false",
      ""
    )

    if (save_image) {
      img_path <- file.path(
        img_dir_abs,
        paste0(name_clean, "_cluster_", cl_id, ".png")
      )
      lines <- c(
        lines,
        paste0(
          "save ",
          img_path,
          " width ",
          image_width,
          " height ",
          image_height,
          " supersample 3"
        ),
        ""
      )
    }

    cxc_path <- file.path(
      out_dir_abs,
      paste0(name_clean, "_cluster_", cl_id, ".cxc")
    )
    writeLines(lines, cxc_path)
    written_files[cl] <- cxc_path
    message("Written: ", cxc_path)
  }

  invisible(written_files)
}

#' Generate a single combined ChimeraX .cxc script showing all clusters at once
#'
#' The structure is colored uniformly in `background_color` (grey by default).
#' Peptides with membership >= `membership_threshold` in cluster N are then
#' overlaid in cluster N's color, so all complexoform groups are visible in one
#' scene. If a peptide qualifies for multiple clusters, the highest cluster
#' index wins (last color command applied).
#'
#' @inheritParams make_cxc
#' @param background_color Character. Hex color for the background chain coloring.
#'
#' @return Invisibly returns the path of the written `.cxc` file.
#'
#' @export
make_cxc_cluster_combined <- function(
  complex_obj,
  output_dir,
  image_dir = NULL,
  membership_threshold = 0.5,
  background_color = "#BBBBBB",
  save_image = TRUE,
  image_width = 2400,
  image_height = 2400
) {
  cs <- data.table::fread(complex_obj$cluster_summary)
  use_cluster_id <- "ClusterID" %in% names(cs)
  if (use_cluster_id) {
    cluster_ids <- sort(unique(cs$ClusterID))
    k <- length(cluster_ids)
  } else {
    mem_cols <- grep("^membership of cluster", names(cs), value = TRUE)
    k <- length(mem_cols)
    if (k == 0L) {
      stop("No membership columns found in cluster summary.")
    }
  }

  threshold <- if (!use_cluster_id) {
    if (identical(membership_threshold, "adaptive")) {
      1 / k
    } else {
      as.numeric(membership_threshold)
    }
  } else {
    NA_real_
  }
  cl_colors <- cluster_palette(k)

  chain_ids <- names(complex_obj$chain_map)
  acc_to_chain <- setNames(names(complex_obj$chain_map), complex_obj$chain_map)

  struct_path <- normalizePath(complex_obj$structure_path)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  out_dir_abs <- normalizePath(output_dir)
  if (is.null(image_dir)) {
    image_dir <- output_dir
  }
  if (!dir.exists(image_dir)) {
    dir.create(image_dir, recursive = TRUE)
  }
  img_dir_abs <- normalizePath(image_dir)

  name_clean <- gsub("[^A-Za-z0-9]+", "_", complex_obj$name)

  lines <- c(
    paste0("# Combined complexoform view – ", complex_obj$name),
    if (use_cluster_id) {
      "# Cluster assignment: ClusterID"
    } else {
      paste0("# Membership threshold: ", round(threshold, 4))
    },
    "",
    paste0("open ", struct_path),
    ""
  )

  # Cartoon display only
  lines <- c(lines, "show cartoons", "hide atoms", "hide surface", "")

  for (ch in chain_ids) {
    lines <- c(lines, paste0("color #1/", ch, " ", background_color))
    lines <- c(lines, paste0("transparency #1/", ch, " 75 cartoons"))
  }
  lines <- c(lines, "")

  for (cl in seq_len(k)) {
    cl_id <- if (use_cluster_id) cluster_ids[cl] else cl
    members <- if (use_cluster_id) {
      cs[ClusterID == cl_id & !Datatype %in% c("NM", "NMpeptide")]
    } else {
      mem_col <- paste("membership of cluster", cl)
      cs[get(mem_col) >= threshold & !Datatype %in% c("NM", "NMpeptide")]
    }
    if (nrow(members) == 0L) {
      next
    }

    lines <- c(lines, paste0("# Cluster ", cl_id, " (", cl_colors[cl], ")"))
    for (i in seq_len(nrow(members))) {
      row <- members[i]
      chain <- acc_to_chain[row$Accession]
      rng <- parse_chimerax_range(row$Position)
      if (is.na(chain) || is.na(rng)) {
        next
      }
      lines <- c(
        lines,
        paste0("color #1/", chain, ":", rng, " ", cl_colors[cl])
      )
      lines <- c(
        lines,
        paste0("transparency #1/", chain, ":", rng, " 0 cartoons") # ensure members are fully opaque
      )
    }
    lines <- c(lines, "")
  }

  lines <- c(
    lines,
    "lighting simple",
    "set bgColor black",
    "graphics silhouettes false",
    ""
  )

  if (save_image) {
    img_path <- file.path(
      img_dir_abs,
      paste0(name_clean, "_cluster_combined.png")
    )
    lines <- c(
      lines,
      paste0(
        "save ",
        img_path,
        " width ",
        image_width,
        " height ",
        image_height,
        " supersample 3"
      ),
      ""
    )
  }

  cxc_path <- file.path(
    out_dir_abs,
    paste0(name_clean, "_cluster_combined.cxc")
  )
  writeLines(lines, cxc_path)
  message("Written: ", cxc_path)
  invisible(cxc_path)
}

#' Generate a combined .cxc script for AF-subunit complexes showing all clusters
#'
#' Like `make_cxc_cluster_combined()` but for complexes without a combined
#' structure: subunits are fetched from AlphaFold and optionally positioned via
#' matchmaker against a reference PDB. The combined coloring (grey background,
#' one color per cluster) is then applied.
#'
#' @inheritParams make_cxc_af_subunits
#' @param background_color Character. Hex color for the background subunit coloring.
#'
#' @return Invisibly returns the path of the written `.cxc` file.
#'
#' @export
make_cxc_af_subunits_cluster_combined <- function(
  complex_obj,
  output_dir,
  image_dir = NULL,
  membership_threshold = 0.5,
  background_color = "#BBBBBB",
  save_image = TRUE,
  image_width = 2400,
  image_height = 2400
) {
  cs <- data.table::fread(complex_obj$cluster_summary)
  use_cluster_id <- "ClusterID" %in% names(cs)
  if (use_cluster_id) {
    cluster_ids <- sort(unique(cs$ClusterID))
    k <- length(cluster_ids)
  } else {
    mem_cols <- grep("^membership of cluster", names(cs), value = TRUE)
    k <- length(mem_cols)
    if (k == 0L) {
      stop("No membership columns found in cluster summary.")
    }
  }

  threshold <- if (!use_cluster_id) {
    if (identical(membership_threshold, "adaptive")) {
      1 / k
    } else {
      as.numeric(membership_threshold)
    }
  } else {
    NA_real_
  }
  cl_colors <- cluster_palette(k)

  subunit_accs <- complex_obj$subunit_af_accessions
  n <- length(subunit_accs)

  has_ref <- !is.null(complex_obj$reference_pdb) &&
    nzchar(complex_obj$reference_pdb)
  model_offset <- if (has_ref) 1L else 0L
  acc_to_model <- setNames(seq_len(n) + model_offset, subunit_accs)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  out_dir_abs <- normalizePath(output_dir)
  if (is.null(image_dir)) {
    image_dir <- output_dir
  }
  if (!dir.exists(image_dir)) {
    dir.create(image_dir, recursive = TRUE)
  }
  img_dir_abs <- normalizePath(image_dir)

  name_clean <- gsub("[^A-Za-z0-9]+", "_", complex_obj$name)

  lines <- c(
    paste0("# Combined complexoform view – ", complex_obj$name),
    if (use_cluster_id) {
      "# Cluster assignment: ClusterID"
    } else {
      paste0("# Membership threshold: ", round(threshold, 4))
    },
    ""
  )

  if (has_ref) {
    lines <- c(
      lines,
      paste0(
        "# Reference scaffold: PDB ",
        complex_obj$reference_pdb,
        " (model #1)"
      ),
      paste0("open ", complex_obj$reference_pdb),
      ""
    )
  }

  lines <- c(lines, "# Fetch individual AlphaFold subunit structures")
  for (i in seq_len(n)) {
    lines <- c(lines, paste0("alphafold fetch ", subunit_accs[i]))
  }
  lines <- c(lines, "")

  # Cartoon display only
  lines <- c(lines, "show cartoons", "hide atoms", "hide surface", "")

  if (has_ref) {
    ref_map <- complex_obj$reference_chain_map
    lines <- c(lines, "# Superimpose AlphaFold models onto reference chains")
    for (i in seq_len(n)) {
      acc <- subunit_accs[i]
      gene <- names(subunit_accs)[i]
      mnum <- i + model_offset
      ref_chain <- ref_map[acc]
      if (!is.na(ref_chain)) {
        lines <- c(lines, paste0("matchmaker #", mnum, "/A to #1/", ref_chain))
      } else {
        lines <- c(
          lines,
          paste0("# ", gene, " has no reference chain — left unpositioned")
        )
      }
    }
    lines <- c(lines, "hide #1 models", "")
  }

  lines <- c(lines, "# Grey background for all subunits")
  for (i in seq_len(n)) {
    lines <- c(
      lines,
      paste0("color #", acc_to_model[subunit_accs[i]], "/A ", background_color)
    )
    lines <- c(
      lines,
      paste0("transparency #", acc_to_model[subunit_accs[i]], "/A 75 cartoons")
    )
  }
  lines <- c(lines, "")

  for (cl in seq_len(k)) {
    cl_id <- if (use_cluster_id) cluster_ids[cl] else cl
    members <- if (use_cluster_id) {
      cs[ClusterID == cl_id & !Datatype %in% c("NM", "NMpeptide")]
    } else {
      mem_col <- paste("membership of cluster", cl)
      cs[get(mem_col) >= threshold & !Datatype %in% c("NM", "NMpeptide")]
    }
    if (nrow(members) == 0L) {
      next
    }

    lines <- c(lines, paste0("# Cluster ", cl_id, " (", cl_colors[cl], ")"))
    for (i in seq_len(nrow(members))) {
      row <- members[i]
      mnum <- acc_to_model[row$Accession]
      rng <- parse_chimerax_range(row$Position)
      if (is.na(mnum) || is.na(rng)) {
        next
      }
      lines <- c(lines, paste0("color #", mnum, "/A:", rng, " ", cl_colors[cl]))
      lines <- c(
        lines,
        paste0("transparency #", mnum, "/A:", rng, " 0 cartoons")
      ) # ensure members are fully opaque
    }
    lines <- c(lines, "")
  }

  lines <- c(
    lines,
    "lighting simple",
    "set bgColor black",
    "graphics silhouettes false",
    ""
  )

  if (save_image) {
    img_path <- file.path(
      img_dir_abs,
      paste0(name_clean, "_cluster_combined.png")
    )
    lines <- c(
      lines,
      paste0(
        "save ",
        img_path,
        " width ",
        image_width,
        " height ",
        image_height,
        " supersample 3"
      ),
      ""
    )
  }

  cxc_path <- file.path(
    out_dir_abs,
    paste0(name_clean, "_cluster_combined.cxc")
  )
  writeLines(lines, cxc_path)
  message("Written: ", cxc_path)
  invisible(cxc_path)
}

#' Compute pairwise Cα distances between PTM sites for each cluster
#'
#' For each cluster, PTM peptides (membership >= threshold, datatype not "NM"
#' or "NMpeptide") are mapped to their approximate modified residue (midpoint
#' of the peptide range). The nearest Cα in the structure is looked up, and
#' pairwise Euclidean distances between all PTM sites in the cluster are saved
#' as a CSV in `output_dir`.
#'
#' @inheritParams make_cxc
#'
#' @return Invisibly NULL.
#'
#' @export
compute_ptm_distances <- function(
  complex_obj,
  output_dir,
  membership_threshold = 0.5
) {
  if (!nzchar(complex_obj$structure_path)) {
    stop(
      "structure_path is empty; use compute_ptm_distances_af_subunits() instead."
    )
  }

  cs <- data.table::fread(complex_obj$cluster_summary)
  use_cluster_id <- "ClusterID" %in% names(cs)
  if (use_cluster_id) {
    cluster_ids <- sort(unique(cs$ClusterID))
    k <- length(cluster_ids)
  } else {
    mem_cols <- grep("^membership of cluster", names(cs), value = TRUE)
    k <- length(mem_cols)
    if (k == 0L) {
      stop("No membership columns found in cluster summary.")
    }
  }

  threshold <- if (!use_cluster_id) {
    as.numeric(membership_threshold)
  } else {
    NA_real_
  }
  nm_types <- c("NM", "NMpeptide")
  acc_to_chain <- setNames(names(complex_obj$chain_map), complex_obj$chain_map)

  message("Reading Cα coordinates from: ", complex_obj$structure_path)
  ca_all <- read_ca_coords_cif(normalizePath(complex_obj$structure_path))

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  out_dir_abs <- normalizePath(output_dir)
  name_clean <- gsub("[^A-Za-z0-9]+", "_", complex_obj$name)

  for (cl in seq_len(k)) {
    cl_id <- if (use_cluster_id) cluster_ids[cl] else cl
    ptm_members <- if (use_cluster_id) {
      cs[ClusterID == cl_id & !Datatype %in% nm_types]
    } else {
      mem_col <- paste("membership of cluster", cl)
      cs[get(mem_col) >= threshold & !Datatype %in% nm_types]
    }

    if (nrow(ptm_members) < 2L) {
      message(
        "Cluster ",
        cl_id,
        ": fewer than 2 PTM members above threshold, skipping."
      )
      next
    }

    site_rows <- list()
    for (i in seq_len(nrow(ptm_members))) {
      row <- ptm_members[i]
      ch_id <- acc_to_chain[row$Accession]
      rng <- parse_chimerax_range(row$Position)
      if (is.na(ch_id) || is.na(rng)) {
        next
      }

      mean_res <- range_midpoint(rng)
      sub_ca <- ca_all[chain == ch_id & !is.na(resno)]
      if (nrow(sub_ca) == 0L) {
        next
      }

      closest <- sub_ca[which.min(abs(sub_ca$resno - mean_res))]
      site_rows[[length(site_rows) + 1L]] <- list(
        Gene = row[["Gene name"]],
        Accession = row$Accession,
        Position = row$Position,
        Datatype = row$Datatype,
        mean_res = mean_res,
        x = closest$x,
        y = closest$y,
        z = closest$z
      )
    }

    sites_dt <- data.table::rbindlist(site_rows, fill = TRUE)
    if (nrow(sites_dt) < 2L) {
      message("Cluster ", cl, ": fewer than 2 mappable PTM sites, skipping.")
      next
    }

    pairs <- .ptm_pairwise_distances(sites_dt)
    csv_path <- file.path(
      out_dir_abs,
      paste0(name_clean, "_cluster_", cl_id, "_ptm_distances.csv")
    )
    data.table::fwrite(pairs, csv_path)
    message("Written: ", csv_path)
  }

  invisible(NULL)
}

#' Compute pairwise PTM-site distances for AF-subunit complexes
#'
#' Like `compute_ptm_distances()` but for complexes without a combined
#' structure. ChimeraX is run headlessly to fetch AlphaFold models, superimpose
#' them onto the reference PDB via matchmaker, and save each subunit's
#' post-matchmaker coordinates to a CIF file. Distances are then computed from
#' those files, so all positioned subunits share the reference coordinate frame
#' and cross-subunit distances are physically meaningful. Floating subunits
#' (no reference chain) are still included but pairs involving them carry
#' `SameCoordFrame = FALSE` in the output.
#'
#' Requires `chimeraX_path` to be defined in the calling environment.
#'
#' @inheritParams make_cxc_af_subunits
#' @param cache_dir Character. Directory where assembled CIF files are cached.
#'
#' @return Invisibly NULL.
#'
#' @export
compute_ptm_distances_af_subunits <- function(
  complex_obj,
  output_dir,
  membership_threshold = 0.5,
  cache_dir = "ComplexoVisualizer/AssembledComplexes",
  af_cache_dir = "ComplexoVisualizer/AfCache"
) {
  cs <- data.table::fread(complex_obj$cluster_summary)
  use_cluster_id <- "ClusterID" %in% names(cs)
  if (use_cluster_id) {
    cluster_ids <- sort(unique(cs$ClusterID))
    k <- length(cluster_ids)
  } else {
    mem_cols <- grep("^membership of cluster", names(cs), value = TRUE)
    k <- length(mem_cols)
    if (k == 0L) {
      stop("No membership columns found in cluster summary.")
    }
  }

  threshold <- if (!use_cluster_id) {
    as.numeric(membership_threshold)
  } else {
    NA_real_
  }
  nm_types <- c("NM", "NMpeptide")
  subunit_accs <- complex_obj$subunit_af_accessions
  n <- length(subunit_accs)

  has_ref <- !is.null(complex_obj$reference_pdb) &&
    nzchar(complex_obj$reference_pdb)
  model_offset <- if (has_ref) 1L else 0L
  ref_map <- if (has_ref) complex_obj$reference_chain_map else character(0)

  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  cache_abs <- normalizePath(cache_dir, mustWork = FALSE)
  if (!dir.exists(af_cache_dir)) {
    dir.create(af_cache_dir, recursive = TRUE)
  }
  af_cache_abs <- normalizePath(af_cache_dir, mustWork = FALSE)
  name_clean <- gsub("[^A-Za-z0-9]+", "_", complex_obj$name)

  # ── Build a headless assembly CXC ─────────────────────────────────────────
  # For each subunit: fetch AF model, matchmaker if positioned, then save the
  # post-matchmaker coordinates to an individual CIF file.  The saved files are
  # all in the reference coordinate frame (for positioned subunits), making
  # cross-subunit distances physically meaningful.

  acc_cif_map <- list() # acc -> list(path, positioned)
  cxc_lines <- character(0)

  if (has_ref) {
    cxc_lines <- c(cxc_lines, paste0("open ", complex_obj$reference_pdb), "")
  }

  cxc_lines <- c(cxc_lines, "# Fetch AlphaFold subunits")
  for (i in seq_len(n)) {
    cxc_lines <- c(cxc_lines, paste0("alphafold fetch ", subunit_accs[i]))
  }
  cxc_lines <- c(cxc_lines, "")

  cxc_lines <- c(cxc_lines, "# Save raw AlphaFold structures")
  for (i in seq_len(n)) {
    acc <- subunit_accs[i]
    mnum <- i + model_offset
    raw_cif <- file.path(af_cache_abs, paste0("AF_", acc, ".cif"))
    cxc_lines <- c(cxc_lines, paste0("save \"", raw_cif, "\" models #", mnum))
  }
  cxc_lines <- c(cxc_lines, "")

  if (has_ref) {
    cxc_lines <- c(cxc_lines, "# Superimpose onto reference chains")
    for (i in seq_len(n)) {
      acc <- subunit_accs[i]
      mnum <- i + model_offset
      ref_chain <- ref_map[acc]
      if (!is.na(ref_chain)) {
        cxc_lines <- c(
          cxc_lines,
          paste0("matchmaker #", mnum, "/A to #1/", ref_chain)
        )
      }
    }
    cxc_lines <- c(cxc_lines, "")
  }

  cxc_lines <- c(
    cxc_lines,
    "# Save each subunit with post-matchmaker coordinates"
  )
  for (i in seq_len(n)) {
    acc <- subunit_accs[i]
    mnum <- i + model_offset
    positioned <- has_ref && !is.na(ref_map[acc])
    suffix <- if (positioned) "assembled" else "unpositioned"
    out_cif <- file.path(cache_abs, paste0("AF_", acc, "_", suffix, ".cif"))
    cxc_lines <- c(cxc_lines, paste0("save \"", out_cif, "\" models #", mnum))
    acc_cif_map[[acc]] <- list(path = out_cif, positioned = positioned)
  }
  cxc_lines <- c(cxc_lines, "", "exit")

  assemble_cxc <- file.path(cache_abs, paste0(name_clean, "_assemble.cxc"))
  writeLines(cxc_lines, assemble_cxc)

  message("Running ChimeraX to assemble: ", complex_obj$name)
  system2(chimeraX_path, args = c("--nogui", normalizePath(assemble_cxc)))

  # ── Read saved coordinates ─────────────────────────────────────────────────
  acc_to_coords <- list()
  for (acc in names(acc_cif_map)) {
    info <- acc_cif_map[[acc]]
    if (!file.exists(info$path)) {
      message("Assembled CIF not found for ", acc, " — skipping.")
      next
    }
    tryCatch(
      acc_to_coords[[acc]] <- list(
        ca = read_ca_coords_cif(info$path),
        positioned = info$positioned
      ),
      error = function(e) {
        message("Could not read coords for ", acc, ": ", e$message)
      }
    )
  }

  # ── Compute per-cluster distances ──────────────────────────────────────────
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  out_dir_abs <- normalizePath(output_dir)

  for (cl in seq_len(k)) {
    cl_id <- if (use_cluster_id) cluster_ids[cl] else cl
    ptm_members <- if (use_cluster_id) {
      cs[ClusterID == cl_id & !Datatype %in% nm_types]
    } else {
      mem_col <- paste("membership of cluster", cl)
      cs[get(mem_col) >= threshold & !Datatype %in% nm_types]
    }

    if (nrow(ptm_members) < 2L) {
      message(
        "Cluster ",
        cl_id,
        ": fewer than 2 PTM members above threshold, skipping."
      )
      next
    }

    site_rows <- list()
    for (i in seq_len(nrow(ptm_members))) {
      row <- ptm_members[i]
      acc <- row$Accession
      info <- acc_to_coords[[acc]]
      if (is.null(info)) {
        next
      }

      rng <- parse_chimerax_range(row$Position)
      if (is.na(rng)) {
        next
      }

      mean_res <- range_midpoint(rng)
      sub_ca <- info$ca[chain == "A" & !is.na(resno)]
      if (nrow(sub_ca) == 0L) {
        next
      }

      closest <- sub_ca[which.min(abs(sub_ca$resno - mean_res))]
      site_rows[[length(site_rows) + 1L]] <- list(
        Gene = row[["Gene name"]],
        Accession = acc,
        Position = row$Position,
        Datatype = row$Datatype,
        mean_res = mean_res,
        coord_source = if (info$positioned) "reference_frame" else "alphafold",
        x = closest$x,
        y = closest$y,
        z = closest$z
      )
    }

    sites_dt <- data.table::rbindlist(site_rows, fill = TRUE)
    if (nrow(sites_dt) < 2L) {
      message("Cluster ", cl, ": fewer than 2 mappable PTM sites, skipping.")
      next
    }

    pairs <- .ptm_pairwise_distances(sites_dt)
    csv_path <- file.path(
      out_dir_abs,
      paste0(name_clean, "_cluster_", cl_id, "_ptm_distances.csv")
    )
    data.table::fwrite(pairs, csv_path)
    message("Written: ", csv_path)
  }

  invisible(NULL)
}

#' Run ChimeraX headlessly on a .cxc script
#'
#' Convenience wrapper around `system2()`. Requires `chimeraX_path` to be
#' defined in the environment (set in ComplexoformStructures.R).
#'
#' @param cxc_path Character. Path to a `.cxc` script file.
#'
#' @return The exit code returned by ChimeraX (0 = success).
#'
#' @export
run_cxc <- function(cxc_path) {
  system2(chimeraX_path, args = normalizePath(cxc_path))
}

#' Compute pairwise Cα distances for all PTM sites across all clusters (combined)
#'
#' Unlike `compute_ptm_distances()`, this function pools every PTM feature in
#' the cluster summary (ignoring cluster membership) and computes a single
#' all-pairs distance file. The membership-based filtering is left to the
#' downstream plotting step, so the output file can be used for both the
#' threshold and highest-membership heatmap modes.
#'
#' @inheritParams compute_ptm_distances
#'
#' @return Invisibly returns the path to the written CSV file.
#'
#' @export
compute_ptm_distances_combined <- function(
  complex_obj,
  output_dir
) {
  if (!nzchar(complex_obj$structure_path)) {
    stop(
      "structure_path is empty; ",
      "use compute_ptm_distances_af_subunits_combined() instead."
    )
  }

  cs <- data.table::fread(complex_obj$cluster_summary)
  nm_types <- c("NM", "NMpeptide")
  ptm_all <- cs[!Datatype %in% nm_types]

  if (nrow(ptm_all) < 2L) {
    message("Fewer than 2 PTM features in summary, skipping.")
    return(invisible(NULL))
  }

  acc_to_chain <- setNames(names(complex_obj$chain_map), complex_obj$chain_map)

  message("Reading Cα coordinates from: ", complex_obj$structure_path)
  ca_all <- read_ca_coords_cif(normalizePath(complex_obj$structure_path))

  site_rows <- list()
  for (i in seq_len(nrow(ptm_all))) {
    row <- ptm_all[i]
    ch_id <- acc_to_chain[row$Accession]
    rng <- parse_chimerax_range(row$Position)
    if (is.na(ch_id) || is.na(rng)) next

    mean_res <- range_midpoint(rng)
    sub_ca <- ca_all[chain == ch_id & !is.na(resno)]
    if (nrow(sub_ca) == 0L) next

    closest <- sub_ca[which.min(abs(sub_ca$resno - mean_res))]
    site_rows[[length(site_rows) + 1L]] <- list(
      Gene = row[["Gene name"]],
      Accession = row$Accession,
      Position = row$Position,
      Datatype = row$Datatype,
      mean_res = mean_res,
      x = closest$x,
      y = closest$y,
      z = closest$z
    )
  }

  sites_dt <- data.table::rbindlist(site_rows, fill = TRUE)
  sites_dt <- unique(sites_dt, by = c("Position", "Datatype"))

  if (nrow(sites_dt) < 2L) {
    message("Fewer than 2 mappable PTM sites, skipping.")
    return(invisible(NULL))
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  out_dir_abs <- normalizePath(output_dir)
  name_clean <- gsub("[^A-Za-z0-9]+", "_", complex_obj$name)

  pairs <- .ptm_pairwise_distances(sites_dt)
  csv_path <- file.path(out_dir_abs, paste0(name_clean, "_all_ptm_distances.csv"))
  data.table::fwrite(pairs, csv_path)
  message("Written: ", csv_path)
  invisible(csv_path)
}

#' Compute pairwise PTM-site distances across all clusters for AF-subunit complexes
#'
#' Like `compute_ptm_distances_combined()` but for complexes without a combined
#' structure. Reuses the assembled CIF files cached by
#' `compute_ptm_distances_af_subunits()`, so that function must have been run
#' first. All PTM features in the summary are pooled regardless of cluster
#' membership; membership-based filtering is left to the plotting step.
#'
#' @inheritParams compute_ptm_distances_af_subunits
#'
#' @return Invisibly returns the path to the written CSV file.
#'
#' @export
compute_ptm_distances_af_subunits_combined <- function(
  complex_obj,
  output_dir,
  cache_dir = "ComplexoVisualizer/AssembledComplexes",
  af_cache_dir = "ComplexoVisualizer/AfCache"
) {
  cs <- data.table::fread(complex_obj$cluster_summary)
  nm_types <- c("NM", "NMpeptide")
  ptm_all <- cs[!Datatype %in% nm_types]

  if (nrow(ptm_all) < 2L) {
    message("Fewer than 2 PTM features in summary, skipping.")
    return(invisible(NULL))
  }

  subunit_accs <- complex_obj$subunit_af_accessions
  n <- length(subunit_accs)
  has_ref <- !is.null(complex_obj$reference_pdb) &&
    nzchar(complex_obj$reference_pdb)
  ref_map <- if (has_ref) complex_obj$reference_chain_map else character(0)

  cache_abs <- normalizePath(cache_dir, mustWork = FALSE)
  name_clean <- gsub("[^A-Za-z0-9]+", "_", complex_obj$name)

  # Reuse CIF files cached by compute_ptm_distances_af_subunits()
  acc_to_coords <- list()
  for (i in seq_len(n)) {
    acc <- subunit_accs[i]
    positioned <- has_ref && !is.na(ref_map[acc])
    suffix <- if (positioned) "assembled" else "unpositioned"
    cif_path <- file.path(cache_abs, paste0("AF_", acc, "_", suffix, ".cif"))
    if (!file.exists(cif_path)) {
      message(
        "Cached CIF not found for ", acc,
        " — run compute_ptm_distances_af_subunits() first."
      )
      next
    }
    tryCatch(
      acc_to_coords[[acc]] <- list(
        ca = read_ca_coords_cif(cif_path),
        positioned = positioned
      ),
      error = function(e) {
        message("Could not read coords for ", acc, ": ", e$message)
      }
    )
  }

  site_rows <- list()
  for (i in seq_len(nrow(ptm_all))) {
    row <- ptm_all[i]
    acc <- row$Accession
    info <- acc_to_coords[[acc]]
    if (is.null(info)) next

    rng <- parse_chimerax_range(row$Position)
    if (is.na(rng)) next

    mean_res <- range_midpoint(rng)
    sub_ca <- info$ca[chain == "A" & !is.na(resno)]
    if (nrow(sub_ca) == 0L) next

    closest <- sub_ca[which.min(abs(sub_ca$resno - mean_res))]
    site_rows[[length(site_rows) + 1L]] <- list(
      Gene = row[["Gene name"]],
      Accession = acc,
      Position = row$Position,
      Datatype = row$Datatype,
      mean_res = mean_res,
      coord_source = if (info$positioned) "reference_frame" else "alphafold",
      x = closest$x,
      y = closest$y,
      z = closest$z
    )
  }

  sites_dt <- data.table::rbindlist(site_rows, fill = TRUE)
  sites_dt <- unique(sites_dt, by = c("Position", "Datatype"))

  if (nrow(sites_dt) < 2L) {
    message("Fewer than 2 mappable PTM sites, skipping.")
    return(invisible(NULL))
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  out_dir_abs <- normalizePath(output_dir)

  pairs <- .ptm_pairwise_distances(sites_dt)
  # Remove cross-frame pairs at source so the CSV only contains valid distances
  if ("SameCoordFrame" %in% names(pairs)) {
    pairs <- pairs[SameCoordFrame == TRUE]
    pairs[, SameCoordFrame := NULL]
  }
  csv_path <- file.path(out_dir_abs, paste0(name_clean, "_all_ptm_distances.csv"))
  data.table::fwrite(pairs, csv_path)
  message("Written: ", csv_path)
  invisible(csv_path)
}
