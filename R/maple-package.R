# =============================================================================
# MAPLE — Multi-resolution Automated Partitioning via Leiden Ensemble
# =============================================================================
#
# ALGORITHM OVERVIEW
# ------------------
# Maple builds a quasi-hierarchical clustering of single cells by running Leiden
# clustering across a resolution grid, measuring the stability of each
# partition via subsampling, and pruning the resulting tree to find the
# optimal multi-resolution partition.
#
# The algorithm proceeds in three steps, mirroring the three user-facing
# wrapper functions:
#
# STEP 1 — Maple()  [wrappers.R]
#   Computes a stability score matrix (cells x resolutions) by:
#   (a) Building a global KNN graph on the reduced dimensions using k_max
#       neighbours — deliberately larger than k to ensure each subsampled
#       cell retains at least k valid neighbours after filtering (see
#       Performance Design below).
#   (b) Running Leiden on the full graph at each resolution in `grid` to
#       obtain the reference partitions.
#   (c) Repeating `nrep` times:
#       - Subsampling `split_prop` of cells.
#       - Filtering the global KNN index to the subsample via .remap_knn()
#         to obtain the subsampled graph without rebuilding the index.
#       - Running Leiden on the subsampled graph at each resolution.
#       - Projecting labels back to held-out cells via k-NN majority vote,
#         using a fresh queryKNN() call against the subsample (see
#         Performance Design for why this is not replaced by filtering).
#       - Computing per-cell Jaccard-based stability scores against the
#         reference partition at each resolution via .obs_score().
#   (d) Averaging scores across repetitions.
#   Calls: calculate_stability()  [core.R]
#     Calls: .create_graph()      [internals.R]
#             .remap_knn()        [internals.R]
#             .map_obs()          [internals.R]
#             .obs_score()        [internals.R]
#
# STEP 2 — grow_Maple()  [wrappers.R]
#   Builds a directed tree where nodes are (resolution, cluster) pairs and
#   edges connect parent clusters to child clusters across adjacent
#   resolutions, filtered by a minimum cell overlap proportion
#   (in_prop_threshold). Each node is coloured by its stability score,
#   computed as the qlevel-th quantile of its cells' stability scores.
#   Calls: create_tree()          [core.R]
#
# STEP 3 — cut_Maple()  [wrappers.R]
#   Prunes the tree via a depth-first search to find the optimal frontier:
#   the set of nodes where stability >= stability_threshold and
#   size >= min_size, maximising resolution depth. A shared memoisation
#   environment ensures each node is visited at most once across all
#   root nodes. Returns cell labels and a per-cell path stability score
#   (minimum stability along the path from root to the assigned node),
#   which reflects the weakest link in each cell's assignment chain.
#   Calls: stability_pruning()    [core.R]
#     Calls: .dfs_prune()         [internals.R]
#             .assign_maple_names() [internals.R]
#
# PERFORMANCE DESIGN
# ------------------
# The dominant computational cost is KNN index construction. Maple builds
# this index ONCE on the full dataset using k_max neighbours (default 50),
# which is deliberately larger than k (default 25). This guarantees that
# after subsampling a proportion `split_prop` of cells, each subsampled
# cell retains at least k valid neighbours within the subsample — because
# the global neighbours are filtered to subsample members via .remap_knn(),
# which is a cheap O(n) operation requiring no distance computations.
#
# Rule of thumb:  k_max >= k / split_prop
# This follows from a Binomial argument: each of the k_max global neighbours
# of a cell is independently retained in the subsample with probability
# split_prop, so the number of valid neighbours is Bin(k_max, split_prop)
# with expectation k_max * split_prop. Setting k_max = k / split_prop
# ensures the expected number of valid neighbours equals k. For stronger
# guarantees, k_max can be chosen via the Binomial CDF to bound the
# per-cell failure probability at a desired level epsilon:
#   P(Bin(k_max, split_prop) < k) <= epsilon
# Default: k_max=50, k=25, split_prop=0.8 -> E[valid] = 40 >= 25 ✓
#
# This avoids rebuilding the KNN index at each of the `nrep` subsampling
# iterations for the subsampled cells, multiplying the index cost by nrep.
#
# NOTE on held-out cells: an analogous filtering approach was considered
# to also replace the queryKNN() call used to project labels onto held-out
# cells. This is theoretically valid under the same k_max assumption, but
# was deliberately not implemented: when k_max is insufficient, filtering
# produces silent label corruption in the majority vote with no recoverable
# fallback, whereas queryKNN() is exact by construction regardless of k_max.
# The queryKNN() call is therefore retained for correctness and safety.
#
# NOTE ON REFERENCE CONSISTENCY
# ------------------------------
# When a pre-computed reference is provided, the stability scores measure
# how consistently cells are assigned relative to that reference across
# subsampled Leiden runs. For stability scores to be meaningful, the
# reference should ideally have been obtained using a graph built with
# the same parameters (k, prune_threshold) as passed to Maple(). A future
# version will support passing a pre-built reference graph directly.
#
# FILE STRUCTURE
# --------------
#   internals.R  — internal helper functions (not exported):
#                  .map_obs, .obs_score, .create_graph,
#                  .remap_knn, .dfs_prune, .assign_maple_names
#   core.R       — core algorithm functions (exported):
#                  calculate_stability, create_tree, stability_pruning
#   wrappers.R   — user-facing SCE wrappers (exported):
#                  Maple, grow_Maple, cut_Maple
#   plots.R       — visualisation functions (exported):
#                  plotTree, plotPruning, plotDiagnostics
#
# DEPENDENCIES
# ------------
#   BiocNeighbors  — KNN index construction and query
#                    (calculate_stability, .create_graph)
#   BiocParallel   — parallel subsampling loops (calculate_stability)
#   bluster        — SNN graph construction (.create_graph)
#   igraph         — graph operations (create_tree, stability_pruning,
#                    .dfs_prune)
#   collapse       — fast mode computation (.map_obs)
#   Matrix         — sparse matrix operations (.create_graph)
#   ggraph         — tree visualisation (plotTree, plotPruning)
#   ggplot2        — diagnostic plots (plotDiagnostics)
#   SingleCellExperiment — SCE class and accessors (wrappers.R)
#
# =============================================================================

#' Maple: Multi-resolution Automated Partitioning via Leiden Ensemble
#'
#' @description
#' Maple identifies a multi-resolution clustering of single-cell data by
#' building a stability-weighted partition tree across a Leiden resolution
#' grid and pruning it to find the optimal simultaneous partition at
#' multiple scales.
#'
#' @author Gabriele Dolfi
#'
#' @importFrom methods is
#' @importFrom stats aggregate median quantile setNames
#'
#' @keywords internal
"_PACKAGE"

# Suppress R CMD check notes for ggplot2 aes() variables
utils::globalVariables(c(
  "Metric", "Value", "resolution", "stability",
  "in_prop", "size", "cluster",
  "is_path", "node_status"
))
