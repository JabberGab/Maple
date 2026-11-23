
##### Internal helpers

#' Classify held-out cells via k-NN majority vote
#'
#' Assigns cluster labels to held-out cells by majority vote over their
#' k nearest neighbours within the subsampled partition.
#'
#' @param n_obs Integer. Total number of cells in the full dataset.
#' @param id Integer vector. Indices of the subsampled cells.
#' @param labels Integer vector. Cluster labels for the subsampled cells,
#'   as returned by Leiden clustering on the subsampled graph.
#' @param k Integer. Number of nearest neighbours to use for majority vote.
#' @param nns List output of \code{BiocNeighbors::queryKNN()}, containing
#'   an \code{index} matrix of dimensions (n_held_out x k) with local
#'   indices into the subsampled cells for each held-out cell.
#'
#' @return Integer vector of length \code{n_obs} containing cluster labels
#'   for all cells — subsampled cells retain their Leiden label, held-out
#'   cells receive the majority vote label from their k nearest neighbours
#'   in the subsample.
#'
#' @keywords internal
.map_obs <- function(n_obs, id, labels, k, nns) {
  new_labels <- numeric(n_obs)
  new_labels[id] <- labels
  new_labels[-id] <- collapse::fmode(
    matrix(labels[nns$index], nrow = k, byrow = TRUE)
  )
  return(new_labels)
}

#' Compute per-cell Jaccard stability scores
#'
#' For each cell, computes the Jaccard similarity between its cluster in the
#' subsampled partition and its cluster in the reference partition. This score
#' reflects how consistently the cell is assigned across subsampling iterations.
#' The per-cell Jaccard scoring approach is based on Yu et al. (2019).
#'
#' @param star_labels Integer vector. Cluster labels for all cells from the
#'   subsampled partition (after majority vote for held-out cells).
#' @param ref_labels Integer vector. Cluster labels for all cells from the
#'   reference partition at the same resolution.
#'
#' @return Numeric vector of length \code{n_obs} containing per-cell Jaccard
#'   scores between 0 and 1, where 1 indicates perfect agreement between the
#'   subsampled and reference partitions for that cell's clusters.
#'
#' @references
#' Yu, H., Chapman, B., Di Florio, A. et al. (2019). Bootstrapping estimates
#' of stability for clusters, observations and model selection.
#' \emph{Computational Statistics}, \strong{34}, 349--372.
#' \doi{10.1007/s00180-018-0830-y}
#'
#' @keywords internal
.obs_score <- function(star_labels, ref_labels) {
  tab            <- collapse::qtab(star_labels, ref_labels)
  size_star      <- rowSums(tab)
  size_ref       <- colSums(tab)
  union_matrix   <- outer(size_star, size_ref, "+") - tab
  jaccard_matrix <- tab / union_matrix
  indices        <- cbind(as.character(star_labels), as.character(ref_labels))
  scores         <- jaccard_matrix[indices]
  return(scores)
}

#' Build a pruned shared nearest neighbour graph
#'
#' Constructs a Shared Nearest Neighbour (SNN) graph from a KNN index matrix,
#' weights edges by the Jaccard similarity between neighbourhood sets, and
#' removes edges whose weight falls below a pruning threshold. Called on both
#' the full KNN index (to build the reference graph) and the filtered subsample
#' KNN index (to build each subsampled graph). Note that the subsampled graph is
#' not a subgraph of the full graph — both edge presence and weights change when
#' the neighbourhood structure is recomputed within the subsample.
#'
#'
#' With larger k values the SNN graph becomes denser compared to typical
#' scran/bluster usage (k = 10). Pruning via \code{prune_threshold} is
#' therefore important to remove spurious weak connections. The larger k
#' is motivated by the subsampling design: each cell needs enough neighbours
#' to retain at least k valid connections after filtering to the subsample
#' (see \code{Maple()} for the binomial argument underlying this choice).
#'
#' @param nn Integer matrix of dimensions (n_cells x k) containing the indices
#'   of the k nearest neighbours for each cell, as returned by
#'   \code{BiocNeighbors::findKNN()} or \code{.remap_knn()}.
#' @param prune_threshold Numeric between 0 and 1. Edges with Jaccard weight strictly
#'   below this value are removed from the SNN graph. Higher values produce
#'   sparser graphs and coarser clusters. Default in \code{Maple()} is
#'   \code{1/15}, matching Seurat's convention, empirically chosen for k = 25.
#'
#' @return An undirected weighted \code{igraph} object with SNN Jaccard weights
#'   on remaining edges.
#'
#' @keywords internal
.create_graph <- function(nn, prune_threshold) {
  graph <- bluster::neighborsToSNNGraph(nn, type = "jaccard")
  adj   <- igraph::as_adjacency_matrix(graph, attr = "weight")
  adj@x[adj@x < prune_threshold] <- 0
  adj   <- Matrix::drop0(adj)
  graph <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected",
                                               weighted = TRUE)
  return(graph)
}

#' Filter and remap a global KNN index to local subsample indices
#'
#' Takes the global KNN index matrix and filters it to retain only neighbours
#' that belong to the subsample, remapping global cell indices to local
#' subsample indices (1..n_sub). This avoids rebuilding the KNN index at
#' each subsampling iteration.
#'
#' The critical assumption underlying this approach is that the nearest
#' neighbours of a subsampled cell within the subsample are exactly the
#' subsampled cells among its global nearest neighbours — i.e. the neighbours
#' of the subsample are the subsample of the neighbours. This holds when
#' k_max is large enough that the k nearest subsampled neighbours of every
#' cell all appear within its k_max global neighbours. By a Binomial argument,
#' the recommended rule of thumb is k_max >= k / split_prop, which ensures
#' the expected number of valid neighbours after filtering equals k.
#'
#' Cells with fewer than k valid neighbours after filtering fall back to
#' self-assignment (the cell is its own neighbour for the missing slots).
#' The number of such cells is returned so \code{calculate_stability()} can
#' report a summary warning to the user.
#'
#' @param nn_global Integer matrix of dimensions (n_obs x k_max) containing
#'   the global KNN index for all cells, as returned by
#'   \code{BiocNeighbors::findKNN()}.
#' @param sub_id Integer vector. Global indices of the subsampled cells.
#' @param n_obs Integer. Total number of cells in the full dataset.
#' @param k Integer. Number of neighbours to retain after filtering.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{out_nn}{Integer matrix of dimensions (n_sub x k) containing
#'       the remapped KNN index for subsampled cells, with local indices
#'       into the subsample (1..n_sub).}
#'     \item{n_fallback}{Integer. Number of subsampled cells that had fewer
#'       than k valid neighbours after filtering and fell back to
#'       self-assignment.}
#'   }
#'
#' @keywords internal
.remap_knn <- function(nn_global, sub_id, n_obs, k) {
  n_sub       <- length(sub_id)
  mapping     <- rep(NA_integer_, n_obs)
  mapping[sub_id] <- seq_len(n_sub)
  nn_raw      <- nn_global[sub_id, , drop = FALSE]
  nn_mapped   <- matrix(mapping[nn_raw], nrow = n_sub, ncol = ncol(nn_raw))
  valid_mask  <- !is.na(nn_mapped)
  valid_counts <- rowSums(valid_mask)
  valid_vals  <- t(nn_mapped)[t(valid_mask)]
  col_indices <- sequence(valid_counts)
  keep        <- col_indices <= k
  final_vals  <- valid_vals[keep]
  final_cols  <- col_indices[keep]
  final_rows  <- rep(seq_len(n_sub), times = valid_counts)[keep]
  out_nn      <- matrix(seq_len(n_sub), nrow = n_sub, ncol = k)
  out_nn[cbind(final_rows, final_cols)] <- final_vals
  n_fallback  <- sum(valid_counts < k)
  return(list(out_nn = out_nn, n_fallback = n_fallback))
}

#' Depth-first search pruning of the Maple tree
#'
#' Traverses the partition tree post-order to find the optimal frontier of
#' valid nodes — the deepest set of nodes where each cluster meets both the
#' stability and minimum size requirements. The pruning is guided by a
#' biological principle: a split is accepted if it produces at least 2 valid
#' subtypes, meaning genuine new substructure is discovered. Cells in failed
#' child nodes are not lost — they are returned to the parent cluster,
#' reflecting the biological reality that some cells cannot be confidently
#' resolved into finer subtypes at a given resolution.
#'
#' A key consequence of this design is that a single split can simultaneously
#' produce finer clusters (for cells that resolve into valid subtypes) and a
#' coarser cluster (for cells whose subtypes are invalid). Maple therefore
#' finds the locally appropriate resolution for each cell population rather
#' than imposing a globally uniform granularity.
#'
#' Memoisation via a shared environment ensures each node is visited at most
#' once across all root nodes, which is important since multiple coarser
#' clusters can share finer descendants in the tree.
#'
#' @param node Character. Name of the current node being evaluated.
#' @param adj Named list. Adjacency list of the tree, mapping each node name
#'   to a character vector of its children, as produced by
#'   \code{split(edges_df$to, edges_df$from)}.
#' @param isValid Named logical vector. Whether each node meets the stability
#'   and minimum size requirements, as defined in \code{stability_pruning()}.
#' @param memo An \code{environment} used as a hash table for memoisation.
#'   Shared across all root node calls to avoid redundant computation.
#'
#' @return A character vector of node names forming the valid frontier of the
#'   subtree rooted at \code{node}, or \code{NULL} if no valid node exists
#'   in the subtree.
#'
#' @keywords internal
.dfs_prune <- function(node, adj, isValid, memo) {

  # Check memoisation cache first
  node_str <- as.character(node)
  if (exists(node_str, envir = memo)) {
    return(memo[[node_str]])
  }

  children <- adj[[node]]

  # Condition 1 — Leaf node: return itself if valid, NULL otherwise
  if (is.null(children) || length(children) == 0) {
    res <- if (isTRUE(isValid[node])) node else NULL
    memo[[node_str]] <- res
    return(res)
  }

  # Recurse post-order: evaluate all children before the parent
  frontier         <- character(0)
  n_children_total <- length(children)
  n_children_valid <- 0

  for (child in children) {
    child_frontier <- .dfs_prune(child, adj, isValid, memo)
    if (!is.null(child_frontier)) {
      frontier         <- c(frontier, child_frontier)
      n_children_valid <- n_children_valid + 1
    }
  }

  # Condition 2 — Accept split: at least 2 valid child subtrees (genuine new
  # substructure), or all children valid (handles single-child nodes correctly)
  if (length(frontier) > 0) {
    is_valid_split   <- n_children_valid >= 2
    is_perfect_split <- n_children_valid == n_children_total
    if (is_valid_split || is_perfect_split) {
      memo[[node_str]] <- frontier
      return(frontier)
    }
  }

  # Condition 3 — Parsimony fallback: split failed downstream, return this
  # node if valid (cells remain at the coarser level), NULL otherwise
  res <- if (isTRUE(isValid[node])) node else NULL
  memo[[node_str]] <- res
  return(res)
}

#' Assign human-readable names to Maple final cluster nodes
#'
#' Generates unique cluster identifiers for the final nodes returned by the
#' DFS pruning. Each cluster is named as \code{number_Rresolution}, where
#' the number reflects the depth order (coarser clusters get lower numbers,
#' finer clusters get higher numbers) and the resolution is the Leiden
#' resolution at which the cluster was found, rounded to 2 decimal places.
#' Within the same resolution, clusters are ordered by Leiden cluster index.
#' The virtual root cluster at resolution 0 is named \code{"unresolved"}.
#'
#' @param final_nodes Character vector. Names of the final nodes returned by
#'   the DFS pruning, in the format \code{"resolution_clusterid"}.
#' @param nodes_df Data frame. Node attributes from the Maple tree, containing
#'   at least \code{name} and \code{resolution} columns.
#'
#' @return A named character vector mapping each final node name to its
#'   human-readable cluster label.
#'
#' @keywords internal
.assign_maple_names <- function(final_nodes, nodes_df) {

  # Match final nodes to their resolution
  idx      <- match(final_nodes, nodes_df$name)
  res_vals <- nodes_df$resolution[idx]

  # Virtual root (resolution 0) is always named "unresolved"
  is_root   <- res_vals == 0
  name_dict <- setNames(rep(NA_character_, length(final_nodes)), final_nodes)
  name_dict[is_root] <- "unresolved"

  # For remaining nodes order by resolution ascending, then by Leiden index
  regular_nodes <- final_nodes[!is_root]
  if (length(regular_nodes) == 0) return(name_dict)

  regular_res   <- res_vals[!is_root]
  ord           <- order(regular_res)
  regular_nodes <- regular_nodes[ord]
  regular_res   <- regular_res[ord]

  # Assign numbers 1..n in depth order
  # Use as.character(round(..., 2)) to avoid trailing zeros (e.g. "0.1" not "0.10")
  res_labels <- paste0(seq_along(regular_nodes), "_R",
                       as.character(round(regular_res, 2)))

  name_dict[regular_nodes] <- res_labels

  return(name_dict)
}
##### Input validation

.is_scalar_num <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
.is_count      <- function(x) .is_scalar_num(x) && x >= 1 && x == round(x)

#' Validate and normalise inputs of calculate_stability()
#'
#' @return A list with \code{grid} sorted in increasing order and
#'   \code{reference} normalised to integer vectors and reordered
#'   consistently with \code{grid} (or NULL).
#' @keywords internal
.check_stability_inputs <- function(X, grid, nrep, split_prop, k, k_max,
                                    prune_threshold, BPPARAM, reference) {

  # X
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("'X' must be a numeric matrix (cells x dimensions).")
  }
  if (anyNA(X) || any(!is.finite(X))) {
    stop("'X' contains NA or non-finite values.")
  }
  n_obs <- nrow(X)

  # grid
  if (!is.numeric(grid) || length(grid) < 1L || anyNA(grid) ||
      any(!is.finite(grid))) {
    stop("'grid' must be a numeric vector of finite resolution values.")
  }
  if (any(grid <= 0)) {
    stop("'grid' values must be strictly positive (0 is reserved for the ",
         "virtual root).")
  }
  if (anyDuplicated(grid)) {
    stop("'grid' contains duplicated values.")
  }
  if (any(abs(grid - round(grid, 2)) > 1e-9)) {
    stop("'grid' values must have at most 2 decimal places.")
  }

  # scalar parameters
  if (!.is_count(nrep))  stop("'nrep' must be a positive integer.")
  if (!.is_count(k))     stop("'k' must be a positive integer.")
  if (!.is_count(k_max)) stop("'k_max' must be a positive integer.")
  if (!.is_scalar_num(split_prop) || split_prop <= 0 || split_prop >= 1) {
    stop("'split_prop' must be a single value in (0, 1).")
  }
  if (!.is_scalar_num(prune_threshold) || prune_threshold < 0 ||
      prune_threshold > 1) {
    stop("'prune_threshold' must be a single value between 0 and 1.")
  }
  if (k_max < k) {
    stop("'k_max' (", k_max, ") must be >= 'k' (", k, ").")
  }
  if (k_max >= n_obs) {
    stop("'k_max' (", k_max, ") must be smaller than the number of cells (",
         n_obs, ").")
  }
  if (floor(split_prop * n_obs) <= k || n_obs - floor(split_prop * n_obs) < 1) {
    stop("'split_prop' leaves too few subsampled or held-out cells for k = ",
         k, ".")
  }
  if (k_max < k / split_prop) {
    warning("k_max < k / split_prop: some subsampled cells may have fewer ",
            "than k valid neighbours after filtering. ",
            "Consider increasing k_max (suggested: k_max >= 2 * k / split_prop).")
  }
  if (!is(BPPARAM, "BiocParallelParam")) {
    stop("'BPPARAM' must be a BiocParallelParam object.")
  }

  # reference
  if (!is.null(reference)) {
    if (is.data.frame(reference) || is.matrix(reference)) {
      reference <- lapply(seq_len(ncol(reference)), function(j) reference[, j])
    }
    if (!is.list(reference)) {
      stop("'reference' must be NULL, a list, a matrix or a data frame ",
           "with one partition per resolution in 'grid'.")
    }
    if (length(reference) != length(grid)) {
      stop("'reference' must contain one partition per resolution in 'grid' ",
           "(got ", length(reference), ", expected ", length(grid), ").")
    }
    reference <- lapply(seq_along(reference), function(j) {
      lab <- reference[[j]]
      if (length(lab) != n_obs) {
        stop("Reference partition ", j, " has length ", length(lab),
             ", expected ", n_obs, " (number of cells).")
      }
      if (anyNA(lab)) {
        stop("Reference partition ", j, " contains NA labels.")
      }
      # keep integer ids when possible, recode anything else
      if (is.numeric(lab) && all(lab == round(lab))) {
        as.integer(lab)
      } else {
        as.integer(factor(lab))
      }
    })
  }

  # sort grid; reference is paired positionally with grid, so reorder it too
  ord <- order(grid)
  if (is.unsorted(grid)) {
    message("'grid' is not in increasing order: sorting it",
            if (!is.null(reference)) " (and 'reference' accordingly)" else "",
            ".")
  }
  grid <- grid[ord]
  if (!is.null(reference)) reference <- reference[ord]

  return(list(grid = grid, reference = reference))
}

#' Error out on unsupported input classes
#' @keywords internal
.unsupported_input <- function(fun, x) {
  supported <- if (fun == "Maple") {
    "SingleCellExperiment, matrix, data.frame"
  } else {
    "SingleCellExperiment, MapleResult"
  }
  stop(fun, "() does not support objects of class '",
       paste(class(x), collapse = "/"), "'. Supported: ", supported, ".",
       call. = FALSE)
}
