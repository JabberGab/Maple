
##### Wrappers

# =============================================================================
# S3 generic definitions
# =============================================================================

#' @rdname Maple
#' @export
Maple <- function(x, ...) UseMethod("Maple")

#' @rdname grow_Maple
#' @export
grow_Maple <- function(x, ...) UseMethod("grow_Maple")

#' @rdname cut_Maple
#' @export
cut_Maple <- function(x, ...) UseMethod("cut_Maple")

# =============================================================================
# MapleResult constructor and print method
# =============================================================================

#' Create a MapleResult object
#'
#' Constructs a \code{MapleResult} S3 object from the output of
#' \code{calculate_stability()}. This is the base object used by the Maple
#' pipeline for non-SCE workflows. It is also stored internally in
#' \code{metadata(sce)$Maple} for SCE workflows.
#'
#' @param reference Integer matrix (n_obs x n_grid + 1) of reference
#'   partition labels, with a prepended column of 1s for the virtual root.
#' @param scores Numeric matrix (n_obs x n_grid + 1) of per-cell Jaccard
#'   stability scores, with a prepended column of 1s for the virtual root.
#' @param grid Numeric vector of resolution values.
#'
#' @return A \code{MapleResult} object.
#'
#' @keywords internal
MapleResult <- function(reference, scores, grid) {
  structure(
    list(
      reference      = reference,
      scores         = scores,
      grid           = grid,
      tree           = NULL,
      pruned_tree    = NULL,
      labels         = NULL,
      path_stability = NULL
    ),
    class = "MapleResult"
  )
}

#' Print method for MapleResult
#'
#' @param x A \code{MapleResult} object.
#' @param ... Further arguments (ignored).
#'
#' @export
print.MapleResult <- function(x, ...) {
  cat("MapleResult object\n")
  cat("  Cells:        ", nrow(x$scores), "\n")
  cat("  Grid:         ", paste(round(x$grid, 2), collapse = ", "), "\n")
  cat("  Tree:         ",
      if (!is.null(x$tree)) "computed" else "not yet computed", "\n")
  cat("  Pruned tree:  ",
      if (!is.null(x$pruned_tree)) "computed" else "not yet computed", "\n")
  cat("  Labels:       ",
      if (!is.null(x$labels)) "computed" else "not yet computed", "\n")
  invisible(x)
}

# =============================================================================
# Maple()
# =============================================================================

#' Compute multi-resolution stability scores for single-cell clustering
#'
#' Runs the first step of the Maple pipeline by computing per-cell stability
#' scores across a resolution grid via repeated subsampling and Leiden
#' clustering.
#'
#' Accepts a \code{SingleCellExperiment} object, a numeric matrix, or a
#' data frame of reduced dimensions. Returns the same type as the input —
#' a \code{SingleCellExperiment} with results stored in
#' \code{metadata(sce)$Maple}, or a \code{MapleResult} object.
#'
#' @param x A \code{SingleCellExperiment} object, a numeric matrix, or a
#'   data frame of dimensions (n_obs x dims).
#' @param dimred Character. Name of the reduced dimension to use. Required
#'   when \code{x} is a \code{SingleCellExperiment}.
#' @param dims Integer. Number of dimensions to use. Required when \code{x}
#'   is a \code{SingleCellExperiment}.
#' @param grid Numeric vector of resolution values at which Leiden clustering
#'   is run. Must have at least 2 values. If \code{reference} is provided,
#'   \code{grid} should match the resolutions at which the reference
#'   clustering was obtained.
#' @param nrep Integer. Number of subsampling iterations. Default 25.
#' @param split_prop Numeric in (0, 1). Proportion of cells to subsample.
#'   Default 0.8.
#' @param k Integer. Number of neighbours for SNN graph and majority vote.
#'   Default 25.
#' @param k_max Integer. Number of neighbours in global KNN index. Must
#'   satisfy \code{k_max >= k / split_prop}. Default 50.
#' @param prune_threshold Numeric between 0 and 1. SNN edge weight threshold.
#'   Default \code{1/15}.
#' @param BPPARAM A \code{BiocParallelParam} object. Default
#'   \code{SerialParam()}. Pass \code{MulticoreParam(workers = 4)} for
#'   parallel execution on Linux/macOS, \code{SnowParam(workers = 4)} on
#'   Windows, or
#'   \code{MulticoreParam(workers = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))}
#'   on SLURM.
#' @param reference \code{NULL} (default), a character vector of
#'   \code{colData(sce)} column names (one per resolution in \code{grid},
#'   SCE input only), or a named list of integer vectors (one per resolution
#'   in \code{grid}).
#' @param add_scores Logical. If \code{TRUE} and \code{x} is a
#'   \code{SingleCellExperiment}, per-cell stability scores are written to
#'   \code{colData(sce)} as \code{Maple_stability_R<resolution>}.
#'   Default \code{FALSE}.
#' @param ... Further arguments (ignored).
#'
#' @return A \code{SingleCellExperiment} with results in
#'   \code{metadata(sce)$Maple}, or a \code{MapleResult} object.
#'
#' @seealso \code{\link{grow_Maple}}, \code{\link{cut_Maple}},
#'   \code{\link{calculate_stability}}
#'
#' @name Maple
#' @export
Maple.matrix <- function(x,
                         grid,
                         nrep            = 25,
                         split_prop      = 0.8,
                         k               = 25,
                         k_max           = 50,
                         prune_threshold = 1/15,
                         BPPARAM         = BiocParallel::SerialParam(),
                         reference       = NULL,
                         ...) {

  # --- Input checks ---
  if (is.character(reference)) {
    stop("Character vector 'reference' (colData column names) is only ",
         "supported for SingleCellExperiment input.")
  }
  # all other checks are done in calculate_stability()

  # --- Compute stability ---
  stability_res <- calculate_stability(
    X               = x,
    grid            = grid,
    nrep            = nrep,
    split_prop      = split_prop,
    k               = k,
    k_max           = k_max,
    prune_threshold = prune_threshold,
    BPPARAM         = BPPARAM,
    reference       = reference
  )

  return(MapleResult(
    reference = stability_res$reference,
    scores    = stability_res$scores,
    grid      = stability_res$grid
  ))
}

#' @rdname Maple
#' @export
Maple.data.frame <- function(x, ...) {
  Maple.matrix(as.matrix(x), ...)
}

#' @rdname Maple
#' @export
Maple.default <- function(x, ...) {
  .unsupported_input("Maple", x)
}

#' @rdname Maple
#' @export
Maple.SingleCellExperiment <- function(x,
                                       dimred,
                                       dims,
                                       grid,
                                       nrep            = 25,
                                       split_prop      = 0.8,
                                       k               = 25,
                                       k_max           = 50,
                                       prune_threshold = 1/15,
                                       BPPARAM         = BiocParallel::SerialParam(),
                                       reference       = NULL,
                                       add_scores      = FALSE,
                                       ...) {

  # --- Input checks ---
  if (!is.character(dimred) || length(dimred) != 1L) {
    stop("'dimred' must be a single reducedDim name.")
  }
  if (!dimred %in% SingleCellExperiment::reducedDimNames(x)) {
    stop("'dimred' = '", dimred, "' not found in reducedDimNames(sce). ",
         "Available: ",
         paste(SingleCellExperiment::reducedDimNames(x), collapse = ", "),
         ".")
  }
  if (!.is_count(dims)) {
    stop("'dims' must be a positive integer.")
  }
  if (dims > ncol(SingleCellExperiment::reducedDim(x, dimred))) {
    stop("'dims' = ", dims, " exceeds the number of available dimensions (",
         ncol(SingleCellExperiment::reducedDim(x, dimred)), ").")
  }

  # --- Extract reduced dimensions ---
  X <- SingleCellExperiment::reducedDim(x, dimred)[, seq_len(dims),
                                                   drop = FALSE]

  # --- Convert colData column names to list of integer vectors ---
  if (is.character(reference)) {
    missing_cols <- setdiff(reference,
                            names(SummarizedExperiment::colData(x)))
    if (length(missing_cols) > 0) {
      stop("The following 'reference' column names were not found in ",
           "colData(sce): ",
           paste(missing_cols, collapse = ", "), ".")
    }
    if (length(reference) != length(grid)) {
      stop("'reference' must have the same length as 'grid'.")
    }
    reference <- lapply(reference, function(col) {
      SummarizedExperiment::colData(x)[[col]]
    })
    # coercion to integer ids is done in calculate_stability()
  }

  # --- Run core method on matrix ---
  result <- Maple.matrix(
    x               = X,
    grid            = grid,
    nrep            = nrep,
    split_prop      = split_prop,
    k               = k,
    k_max           = k_max,
    prune_threshold = prune_threshold,
    BPPARAM         = BPPARAM,
    reference       = reference
  )

  # --- Store MapleResult in metadata ---
  S4Vectors::metadata(x)$Maple <- result

  # --- Optionally write per-resolution scores to colData ---
  if (add_scores) {
    for (i in seq_along(result$grid)) {
      col_name      <- paste0("Maple_stability_R", round(result$grid[i], 2))
      x[[col_name]] <- result$scores[, i + 1L]
    }
  }

  return(x)
}

# =============================================================================
# grow_Maple()
# =============================================================================

#' Build the Maple partition tree
#'
#' Runs the second step of the Maple pipeline by constructing a
#' stability-weighted partition tree. Accepts either a
#' \code{SingleCellExperiment} or a \code{MapleResult} object.
#'
#' @param x A \code{SingleCellExperiment} or \code{MapleResult} object on
#'   which \code{Maple()} has been run.
#' @param in_prop_threshold Numeric between 0 and 1. Minimum proportion of a child
#'   cluster's cells that must come from a single parent cluster for the
#'   edge to be retained. Default 0.1.
#' @param qlevel Numeric between 0 and 1. Quantile of per-cell stability scores
#'   used to assign a stability value to each node. Default 0.2.
#' @param verbose Logical. If \code{TRUE}, prints a summary of the tree.
#'   Default \code{TRUE}.
#' @param ... Further arguments (ignored).
#'
#' @return The same type as the input with the tree added.
#'
#' @seealso \code{\link{Maple}}, \code{\link{cut_Maple}},
#'   \code{\link{plotTree}}, \code{\link{plotDiagnostics}}
#'
#' @name grow_Maple
#' @export
grow_Maple.MapleResult <- function(x,
                                   in_prop_threshold = 0.1,
                                   qlevel            = 0.2,
                                   verbose           = TRUE,
                                   ...) {

  # --- Input checks ---
  if (is.null(x$reference)) {
    stop("Maple() must be run before grow_Maple().")
  }
  if (!.is_scalar_num(in_prop_threshold) || in_prop_threshold < 0 ||
      in_prop_threshold > 1) {
    stop("'in_prop_threshold' must be a single value between 0 and 1.")
  }
  if (!.is_scalar_num(qlevel) || qlevel < 0 || qlevel > 1) {
    stop("'qlevel' must be a single value between 0 and 1.")
  }

  # --- Build tree ---
  tree   <- create_tree(x$reference, x$scores, x$grid,
                        in_prop_threshold, qlevel)
  x$tree <- tree

  if (verbose) {
    message("Maple tree built: ", igraph::vcount(tree), " nodes, ",
            igraph::ecount(tree), " edges.")
  }

  return(x)
}

#' @rdname grow_Maple
#' @export
grow_Maple.default <- function(x, ...) .unsupported_input("grow_Maple", x)

#' @rdname grow_Maple
#' @export
grow_Maple.SingleCellExperiment <- function(x,
                                            in_prop_threshold = 0.1,
                                            qlevel            = 0.2,
                                            verbose           = TRUE,
                                            ...) {

  # --- Input checks ---
  if (is.null(S4Vectors::metadata(x)$Maple)) {
    stop("Maple() must be run before grow_Maple(). ",
         "No Maple results found in metadata(sce)$Maple.")
  }

  # --- Run core method on MapleResult ---
  S4Vectors::metadata(x)$Maple <- grow_Maple.MapleResult(
    S4Vectors::metadata(x)$Maple,
    in_prop_threshold = in_prop_threshold,
    qlevel            = qlevel,
    verbose           = verbose
  )

  return(x)
}

# =============================================================================
# cut_Maple()
# =============================================================================

#' Prune the Maple tree to obtain the optimal multi-resolution partition
#'
#' Runs the third step of the Maple pipeline by pruning the
#' stability-weighted partition tree. Accepts either a
#' \code{SingleCellExperiment} or a \code{MapleResult} object.
#'
#' @param x A \code{SingleCellExperiment} or \code{MapleResult} object on
#'   which \code{Maple()} and \code{grow_Maple()} have been run.
#' @param stability_threshold Numeric between 0 and 1 or \code{"adaptive"}.
#'   Minimum stability score for a node to be considered valid. If
#'   \code{"adaptive"}, set to the median stability score across all cells
#'   and resolutions (excluding the virtual root). Use
#'   \code{plotDiagnostics()} to inspect the stability distribution before
#'   choosing a threshold. Default 0.8.
#' @param min_size Integer. Minimum number of cells for a node to be
#'   considered valid. Default 100.
#' @param verbose Logical. If \code{TRUE}, prints a pruning summary.
#'   Default \code{TRUE}.
#' @param ... Further arguments (ignored).
#'
#' @return The same type as the input with labels and path stability added.
#'   For \code{MapleResult}, results are stored in \code{x$labels} and
#'   \code{x$path_stability}. For \code{SingleCellExperiment}, results are
#'   stored in \code{colData(sce)$Maple_labels} and
#'   \code{colData(sce)$Maple_path_stability}.
#'
#' @seealso \code{\link{Maple}}, \code{\link{grow_Maple}},
#'   \code{\link{plotPruning}}, \code{\link{plotDiagnostics}}
#'
#' @name cut_Maple
#' @export
cut_Maple.MapleResult <- function(x,
                                  stability_threshold = 0.8,
                                  min_size            = 100,
                                  verbose             = TRUE,
                                  ...) {

  # --- Input checks ---
  if (is.null(x$tree)) {
    stop("grow_Maple() must be run before cut_Maple().")
  }
  if (!.is_count(min_size)) {
    stop("'min_size' must be a positive integer.")
  }
  if (length(stability_threshold) != 1L) {
    stop("'stability_threshold' must be a single value between 0 and 1 or 'adaptive'.")
  }

  # --- Adaptive threshold ---
  if (is.character(stability_threshold) &&
      stability_threshold == "adaptive") {
    stability_threshold <- median(x$scores[, -1], na.rm = TRUE)
    if (verbose) {
      message("Using adaptive stability threshold: ",
              round(stability_threshold, 3), ".")
    }
  } else if (!.is_scalar_num(stability_threshold) ||
             stability_threshold < 0 ||
             stability_threshold > 1) {
    stop("'stability_threshold' must be a numeric value between 0 and 1 ",
         "or 'adaptive'.")
  }

  # --- Prune tree ---
  final_res <- stability_pruning(
    reference           = x$reference,
    grid                = x$grid,
    scores              = x$scores,
    tree                = x$tree,
    stability_threshold = stability_threshold,
    min_size            = min_size
  )

  # --- Store results ---
  # Pruned tree stored separately to preserve original for plotTree()
  x$pruned_tree    <- final_res$tree
  x$labels         <- final_res$final_labels
  x$path_stability <- final_res$path_stability

  # --- Verbose summary ---
  if (verbose) {
    n_clusters   <- length(levels(final_res$final_labels)) -
      as.integer("unresolved" %in% levels(final_res$final_labels))
    n_unresolved <- sum(final_res$final_labels == "unresolved",
                        na.rm = TRUE)
    message("Maple pruning complete: ", n_clusters, " clusters found ",
            "(stability_threshold = ", round(stability_threshold, 3),
            ", min_size = ", min_size, "). ",
            n_unresolved, " unresolved cells.")
  }

  return(x)
}

#' @rdname cut_Maple
#' @export
cut_Maple.default <- function(x, ...) .unsupported_input("cut_Maple", x)

#' @rdname cut_Maple
#' @export
cut_Maple.SingleCellExperiment <- function(x,
                                           stability_threshold = 0.8,
                                           min_size            = 100,
                                           verbose             = TRUE,
                                           ...) {

  # --- Input checks ---
  if (is.null(S4Vectors::metadata(x)$Maple)) {
    stop("Maple() must be run before cut_Maple(). ",
         "No Maple results found in metadata(sce)$Maple.")
  }
  if (is.null(S4Vectors::metadata(x)$Maple$tree)) {
    stop("grow_Maple() must be run before cut_Maple(). ",
         "No tree found in metadata(sce)$Maple$tree.")
  }

  # --- Run core method on MapleResult ---
  result <- cut_Maple.MapleResult(
    S4Vectors::metadata(x)$Maple,
    stability_threshold = stability_threshold,
    min_size            = min_size,
    verbose             = verbose
  )

  # --- Store results back in SCE ---
  S4Vectors::metadata(x)$Maple     <- result
  x$Maple_labels                    <- result$labels
  x$Maple_path_stability            <- result$path_stability

  return(x)
}
