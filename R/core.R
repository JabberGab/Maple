
##### Core functions

#' Compute per-cell stability scores across a resolution grid
#'
#' Estimates the stability of each cell's cluster assignment by repeatedly
#' subsampling the data, running Leiden clustering on the subsampled graph,
#' projecting labels back to held-out cells via k-NN majority vote, and
#' computing per-cell Jaccard scores against a reference partition. Scores
#' are averaged across \code{nrep} iterations.
#'
#' The stability-via-subsampling framework follows Yu et al. (2019). A global
#' KNN index is built once using \code{k_max} neighbours and filtered to each
#' subsample via \code{.remap_knn()}, avoiding repeated index construction.
#' By a Binomial argument, the recommended rule of thumb is
#' \code{k_max >= k / split_prop}, which ensures the expected number of valid
#' neighbours after filtering equals \code{k}.
#'
#' @param X Numeric matrix of dimensions (n_obs x dims) containing the reduced
#'   dimension representation of the data.
#' @param grid Numeric vector of resolution values at which Leiden clustering
#'   is run on the subsampled graphs. Always required.
#' @param nrep Integer. Number of subsampling iterations. Default 25.
#' @param split_prop Numeric in (0, 1). Proportion of cells to subsample at
#'   each iteration. Default 0.8.
#' @param k Integer. Number of neighbours used to build the SNN graph and
#'   for k-NN majority vote of held-out cells. Default 25.
#' @param k_max Integer. Number of neighbours in the global KNN index.
#'   Must satisfy \code{k_max >= k / split_prop}. Default 50.
#' @param prune_threshold Numeric between 0 and 1. SNN edge weight threshold.
#'   Default \code{1/15}.
#' @param BPPARAM A \code{BiocParallelParam} object controlling parallelisation.
#'   Default is \code{SerialParam()} for serial execution, which works on all
#'   platforms. For parallel execution pass \code{MulticoreParam(workers = 4)}
#'   on Linux/macOS, \code{SnowParam(workers = 4)} on Windows, or
#'   \code{MulticoreParam(workers = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))}
#'   on a SLURM HPC cluster.
#' @param reference A named list of integer vectors of length \code{n_obs},
#'   one per resolution in \code{grid}, providing pre-computed reference
#'   partitions. If NULL, reference partitions are computed via Leiden
#'   clustering on the full graph at each resolution in \code{grid}.
#'   Must have the same length as \code{grid} if provided.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{scores}{Numeric matrix of dimensions (n_obs x n_grid + 1) containing
#'       average per-cell Jaccard stability scores across all iterations,
#'       with a prepended column of 1s for the virtual root at resolution 0.}
#'     \item{reference}{Integer matrix of dimensions (n_obs x n_grid + 1)
#'       containing the reference partition labels at each resolution, with a
#'       prepended column of 1s for the virtual root at resolution 0.}
#'   }
#'
#' @references
#' Yu, H., Chapman, B., Di Florio, A. et al. (2019). Bootstrapping estimates
#' of stability for clusters, observations and model selection.
#' \emph{Computational Statistics}, \strong{34}, 349--372.
#' \doi{10.1007/s00180-018-0830-y}
#'
#' @export
calculate_stability <- function(X,
                                grid,
                                nrep            = 25,
                                split_prop      = 0.8,
                                k               = 25,
                                k_max           = 50,
                                prune_threshold = 1/15,
                                BPPARAM         = BiocParallel::SerialParam(),
                                reference       = NULL) {

  # --- Input checks ---
  checked   <- .check_stability_inputs(X, grid, nrep, split_prop, k, k_max,
                                       prune_threshold, BPPARAM, reference)
  grid      <- checked$grid
  reference <- checked$reference

  # --- Setup ---
  n_obs  <- nrow(X)
  n_grid <- length(grid)

  # Build global KNN index once
  nn_global <- BiocNeighbors::findKNN(X,
                                      k            = k_max,
                                      BNPARAM      = BiocNeighbors::AnnoyParam(),
                                      get.distance = FALSE)$index

  # --- Reference partitions ---
  if (is.null(reference)) {
    message("Computing reference partitions via Leiden...")
    ref_graph  <- .create_graph(nn_global[, 1:k], prune_threshold)
    ref_labels <- BiocParallel::bplapply(grid, function(r) {
      igraph::cluster_leiden(ref_graph,
                             resolution         = r,
                             objective_function = "modularity")$membership
    }, BPPARAM = BPPARAM)
    names(ref_labels) <- as.character(grid)
  } else {
    ref_labels        <- reference
    names(ref_labels) <- as.character(grid)
  }

  # --- Subsampling loop ---
  message("Computing stability scores...")
  results <- BiocParallel::bplapply(seq_len(nrep), function(i) {

    res_matrix <- matrix(NA_real_, nrow = n_obs, ncol = n_grid)
    sub_id     <- sample(seq_len(n_obs), floor(split_prop * n_obs))

    # Remap global KNN to subsample (avoids rebuilding the index)
    remap_res  <- .remap_knn(nn_global, sub_id, n_obs, k)
    nn_sub     <- remap_res$out_nn
    n_fallback <- remap_res$n_fallback

    sub_graph  <- .create_graph(nn_sub, prune_threshold)

    # Query held-out cells against subsample
    nn_map <- BiocNeighbors::queryKNN(X[ sub_id, , drop = FALSE],
                                      X[-sub_id, , drop = FALSE],
                                      k       = k,
                                      BNPARAM = BiocNeighbors::AnnoyParam())

    for (r in seq_along(grid)) {
      star_labels     <- .map_obs(
        n_obs, sub_id,
        igraph::cluster_leiden(
          sub_graph,
          resolution         = grid[r],
          objective_function = "modularity")$membership,
        k, nn_map)
      current_scores  <- .obs_score(star_labels, ref_labels[[r]])
      res_matrix[, r] <- current_scores
    }

    return(list(scores = res_matrix, n_fallback = n_fallback))

  }, BPPARAM = BPPARAM)

  # --- Accumulate fallback warnings ---
  total_fallbacks <- sum(sapply(results, function(x) x$n_fallback))
  if (total_fallbacks > 0) {
    warning("Across all ", nrep, " subsampling iterations, ", total_fallbacks,
            " cell-iteration pairs had fewer than k = ", k,
            " valid neighbours after filtering to the subsample. ",
            "Consider increasing k_max (suggested: k_max >= 2 * k / split_prop).")
  }

  # --- Average scores across iterations ---
  scores <- Reduce("+", lapply(results, function(x) x$scores)) / nrep

  # --- Format reference matrix ---
  ref_matrix           <- do.call(cbind, ref_labels)
  colnames(ref_matrix) <- as.character(round(grid, 2))

  # --- Prepend virtual root (resolution 0): all cells in cluster 1, stability 1 ---
  virtual_ref    <- matrix(1L, nrow = n_obs, ncol = 1,
                           dimnames = list(NULL, "0"))
  virtual_scores <- matrix(1,  nrow = n_obs, ncol = 1,
                           dimnames = list(NULL, "0"))

  ref_matrix <- cbind(virtual_ref,    ref_matrix)
  scores     <- cbind(virtual_scores, scores)

  return(list(scores    = scores,
              reference = ref_matrix,
              grid      = grid))
}

#' Build a stability-weighted partition tree
#'
#' Constructs a directed tree where nodes are (resolution, cluster) pairs and
#' edges connect parent clusters to child clusters across adjacent resolutions,
#' filtered by a minimum cell overlap proportion. Each node is weighted by its
#' stability score, computed as the \code{qlevel}-th quantile of its cells'
#' stability scores. A virtual root node at resolution 0 containing all cells
#' with stability 1 is prepended to the tree, guaranteeing a valid fallback
#' for every cell during pruning.
#'
#' @param reference Integer matrix of dimensions (n_obs x n_grid + 1) as
#'   returned by \code{calculate_stability()}, with the first column being
#'   the virtual root (all 1s) and subsequent columns containing Leiden
#'   partition labels at each resolution in \code{grid}.
#' @param scores Numeric matrix of dimensions (n_obs x n_grid + 1) as
#'   returned by \code{calculate_stability()}, with the first column being
#'   the virtual root stability scores (all 1s) and subsequent columns
#'   containing per-cell Jaccard stability scores at each resolution.
#' @param grid Numeric vector of resolution values, not including resolution 0
#'   (the virtual root). Must match the columns of \code{reference} and
#'   \code{scores} excluding the first column.
#' @param in_prop_threshold Numeric between 0 and 1. Minimum proportion of a child
#'   cluster's cells that must come from a single parent cluster for the
#'   edge to be retained. Higher values produce a cleaner tree with fewer
#'   cross-resolution connections.
#' @param qlevel Numeric between 0 and 1. Quantile of per-cell stability scores
#'   used to assign a stability value to each node. Lower values are more
#'   conservative.
#'
#' @return A directed \code{igraph} object with node attributes \code{name},
#'   \code{resolution}, \code{cluster}, \code{size} and \code{stability},
#'   and edge attributes \code{overlap}, \code{in_prop} and \code{out_prop}.
#'
#' @export
create_tree <- function(reference,
                        scores,
                        grid,
                        in_prop_threshold,
                        qlevel) {

  n_obs  <- nrow(reference)
  n_grid <- length(grid)

  # --- Virtual root node at resolution 0 ---
  virtual_node <- data.frame(
    name       = "0_1",
    resolution = 0,
    cluster    = "1",
    size       = n_obs,
    stability  = 1,
    stringsAsFactors = FALSE
  )

  # --- Nodes for each resolution in grid ---
  # reference and scores columns: 1 = virtual root, 2..(n_grid+1) = grid
  nodes_list <- vector("list", n_grid)
  for (i in seq_len(n_grid)) {
    col      <- i + 1  # offset by 1 for virtual root column
    res_val  <- grid[i]
    res_name <- as.character(round(res_val, 2))
    labels   <- reference[, col]
    clust_sizes <- table(labels)
    clust_stab  <- tapply(scores[, col], labels,
                          function(x) quantile(x, qlevel))
    nodes_list[[i]] <- data.frame(
      name       = paste0(res_name, "_", names(clust_sizes)),
      resolution = res_val,
      cluster    = names(clust_sizes),
      size       = as.numeric(clust_sizes),
      stability  = as.numeric(clust_stab[names(clust_sizes)]),
      stringsAsFactors = FALSE
    )
  }

  nodes_df <- rbind(virtual_node, do.call(rbind, nodes_list))

  # --- Edges from virtual root to grid[1] nodes ---
  first_col   <- reference[, 2]  # grid[1] labels
  first_sizes <- table(first_col)
  virtual_edges <- data.frame(
    from     = "0_1",
    to       = paste0(as.character(round(grid[1], 2)), "_",
                      names(first_sizes)),
    overlap  = as.numeric(first_sizes),
    in_prop  = 1,
    out_prop = as.numeric(first_sizes) / n_obs,
    stringsAsFactors = FALSE
  )

  # --- Edges between adjacent resolutions ---
  edges_list <- vector("list", n_grid - 1)
  for (i in seq_len(n_grid - 1)) {
    col_father    <- i + 1
    col_son       <- i + 2
    res_father    <- as.character(round(grid[i],     2))
    res_son       <- as.character(round(grid[i + 1], 2))
    father_labels <- reference[, col_father]
    son_labels    <- reference[, col_son]
    tt            <- table(father_labels, son_labels)
    overlap_df    <- as.data.frame(as.table(tt), stringsAsFactors = FALSE)
    colnames(overlap_df) <- c("father_clust", "son_clust", "overlap")
    overlap_df    <- overlap_df[overlap_df$overlap > 0, ]
    overlap_df$from     <- paste0(res_father, "_", overlap_df$father_clust)
    overlap_df$to       <- paste0(res_son,    "_", overlap_df$son_clust)
    sizes_from    <- nodes_df$size[match(overlap_df$from, nodes_df$name)]
    sizes_to      <- nodes_df$size[match(overlap_df$to,   nodes_df$name)]
    overlap_df$out_prop <- overlap_df$overlap / sizes_from
    overlap_df$in_prop  <- overlap_df$overlap / sizes_to
    overlap_df    <- overlap_df[overlap_df$in_prop >= in_prop_threshold, ]
    edges_list[[i]] <- overlap_df[, c("from", "to", "overlap",
                                      "in_prop", "out_prop")]
  }

  edges_df <- rbind(virtual_edges, do.call(rbind, edges_list))

  # --- Build igraph tree ---
  tree <- igraph::graph_from_data_frame(
    d        = edges_df,
    vertices = nodes_df,
    directed = TRUE
  )

  return(tree)
}

#' Prune the Maple tree to find the optimal multi-resolution partition
#'
#' Traverses the stability-weighted partition tree via depth-first search to
#' find the optimal frontier of valid nodes — the deepest set of nodes where
#' each cluster meets both the stability and minimum size requirements. The
#' pruning follows a biological principle: a split is accepted if it produces
#' at least 2 valid finer clusters. Cells in failed child nodes are returned
#' to their parent cluster, reflecting the reality that some cells cannot be
#' confidently resolved into finer subtypes at a given resolution. A single
#' split can therefore simultaneously produce finer clusters (for cells that
#' resolve into valid subtypes) and a coarser cluster (for cells whose
#' subtypes are invalid).
#'
#' @param reference Integer matrix of dimensions (n_obs x n_grid + 1) as
#'   returned by \code{calculate_stability()}.
#' @param grid Numeric vector of resolution values, not including resolution 0.
#' @param scores Numeric matrix of dimensions (n_obs x n_grid + 1) as
#'   returned by \code{calculate_stability()}.
#' @param tree A directed \code{igraph} object as returned by
#'   \code{create_tree()}.
#' @param stability_threshold Numeric between 0 and 1. Minimum stability score
#'   for a node to be considered valid.
#' @param min_size Integer. Minimum number of cells for a node to be
#'   considered valid.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{final_labels}{Factor of length n_obs containing the cluster
#'       label for each cell, named as \code{number_Rresolution}. Cells that
#'       could not be resolved into any stable cluster are labelled
#'       \code{"unresolved"}.}
#'     \item{path_stability}{Numeric vector of length n_obs containing the
#'       minimum stability score along the path from the virtual root to the
#'       assigned node for each cell. \code{NA} for unresolved cells.}
#'     \item{tree}{The input \code{igraph} tree with node status annotations
#'       added: \code{"Optimal"} for final frontier nodes, \code{"Path"} for
#'       nodes on the path to the frontier, \code{"Tree"} for all others.}
#'   }
#'
#' @export
stability_pruning <- function(reference,
                              grid,
                              scores,
                              tree,
                              stability_threshold,
                              min_size) {

  # --- 1. Extract data structures from igraph ---
  nodes_df <- igraph::as_data_frame(tree, what = "vertices")
  edges_df <- igraph::as_data_frame(tree, what = "edges")

  # --- 2. DFS setup ---
  isValid <- setNames(
    nodes_df$stability >= stability_threshold & nodes_df$size >= min_size,
    nodes_df$name
  )
  adj <- split(edges_df$to, edges_df$from)

  # Virtual root is always the single root node
  root_node   <- "0_1"
  final_nodes <- character(0)

  # Shared memoisation environment across all DFS calls
  memo_env <- new.env(hash = TRUE, parent = emptyenv())

  # --- 3. DFS pruning from virtual root ---
  res <- .dfs_prune(root_node, adj, isValid, memo_env)

  if (!is.null(res)) {
    final_nodes <- res
  } else {
    final_nodes <- root_node
  }

  # Remove duplicates
  final_nodes <- unique(final_nodes)

  # --- 4. Cluster naming ---
  name_dict <- .assign_maple_names(final_nodes, nodes_df)
  # Always ensure "0_1" maps to "unresolved"
  name_dict["0_1"] <- "unresolved"

  # --- 5. Build assignment matrix ---
  n_cells      <- nrow(reference)
  K            <- length(final_nodes)
  A            <- matrix(FALSE, nrow = n_cells, ncol = K)
  col_idx_node <- integer(K)

  parts      <- strsplit(final_nodes, "_")
  res_vals   <- as.numeric(sapply(parts, `[`, 1))
  clust_vals <- as.numeric(sapply(parts, `[`, 2))

  for (k in seq_len(K)) {
    if (res_vals[k] == 0) {
      # Virtual root: first column of reference
      col_idx <- 1L
    } else {
      col_idx <- which(abs(grid - res_vals[k]) < 1e-9) + 1L
    }
    col_idx_node[k] <- col_idx
    A[, k]          <- reference[, col_idx] == clust_vals[k]
  }

  # --- 6. Conflict resolution and assignment ---
  n_total          <- rowSums(A)
  final_labels_raw <- character(n_cells)

  # Unresolved cells (no matching final node)
  final_labels_raw[n_total == 0] <- "0_1"

  # Unique assignment
  idx_1 <- which(n_total == 1)
  if (length(idx_1) > 0) {
    final_labels_raw[idx_1] <- final_nodes[max.col(A[idx_1, , drop = FALSE])]
  }

  # Conflict resolution via stability scores
  idx_m <- which(n_total > 1)
  if (length(idx_m) > 0) {
    S_nodes  <- matrix(-1, nrow = length(idx_m), ncol = K)
    for (k in seq_len(K)) {
      S_nodes[, k] <- scores[idx_m, col_idx_node[k]]
    }
    S_masked <- ifelse(A[idx_m, , drop = FALSE], S_nodes, -1)
    col_win  <- max.col(S_masked, ties.method = "first")
    final_labels_raw[idx_m] <- final_nodes[col_win]
  }

  # --- 7. Path stability ---
  # NA for unresolved cells (virtual root or n_total == 0)
  path_stability <- rep(NA_real_, n_cells)

  for (k in seq_len(K)) {
    if (final_nodes[k] == "0_1") next
    idx_cells <- which(final_labels_raw == final_nodes[k])
    if (length(idx_cells) > 0) {
      res_max <- col_idx_node[k]
      path_stability[idx_cells] <- matrixStats::rowMins(
        scores[idx_cells, 1:res_max, drop = FALSE],
        na.rm = TRUE
      )
    }
  }

  # --- 8. Translate raw node names to human-readable labels ---
  final_labels <- unname(name_dict[final_labels_raw])

  # --- 9. Tree annotation ---
  if (length(final_nodes) > 0) {
    path_nodes <- final_nodes
    curr       <- final_nodes

    while (length(curr) > 0) {
      parents    <- edges_df$from[edges_df$to %in% curr]
      curr       <- unique(parents)
      path_nodes <- unique(c(path_nodes, curr))
    }

    node_names  <- igraph::V(tree)$name
    node_status <- rep("Tree", length(node_names))
    node_status[node_names %in% path_nodes]  <- "Path"
    node_status[node_names %in% final_nodes] <- "Optimal"

    igraph::V(tree)$node_status <- node_status

    maple_label <- rep(NA_character_, length(node_names))
    is_final    <- node_names %in% final_nodes
    maple_label[is_final] <- name_dict[node_names[is_final]]
    igraph::V(tree)$maple_label <- maple_label
    igraph::E(tree)$is_path     <- edges_df$to %in% path_nodes
  }

  return(list(
    final_labels   = factor(final_labels,
                            levels = intersect(c(name_dict[name_dict != "unresolved"],
                                                 "unresolved"),
                                               final_labels)),
    path_stability = path_stability,
    tree           = tree
  ))
}
