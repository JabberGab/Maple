##### Plot functions

# =============================================================================
# Helper to extract Maple objects from SCE or MapleResult
# =============================================================================

.get_maple_result <- function(x) {
  if (is(x, "SingleCellExperiment")) {
    result <- S4Vectors::metadata(x)$Maple
    if (is.null(result)) {
      stop("Maple() must be run first. ",
           "No Maple results found in metadata(sce)$Maple.")
    }
    return(result)
  } else if (inherits(x, "MapleResult")) {
    return(x)
  } else {
    stop("'x' must be a SingleCellExperiment or a MapleResult object.")
  }
}

.get_tree <- function(x) {
  result <- .get_maple_result(x)
  if (is.null(result$tree)) {
    stop("grow_Maple() must be run before plotting the tree.")
  }
  return(result$tree)
}

.get_pruned_tree <- function(x) {
  result <- .get_maple_result(x)
  if (is.null(result$pruned_tree)) {
    stop("cut_Maple() must be run before plotting the pruning.")
  }
  return(result$pruned_tree)
}

# =============================================================================
# plotTree()
# =============================================================================

#' Plot the Maple partition tree
#'
#' Visualises the stability-weighted partition tree built by
#' \code{grow_Maple()}, with nodes sized by cluster size and coloured by
#' stability score. The virtual root node at resolution 0 is shown in grey
#' to indicate it is an internal anchor node.
#'
#' @param x A \code{SingleCellExperiment} or \code{MapleResult} object on
#'   which \code{Maple()} and \code{grow_Maple()} have been run.
#' @param colours Character vector of colours for the stability gradient,
#'   from low to high stability. Default is a red-to-green palette.
#'
#' @return A \code{ggplot2} object.
#'
#' @seealso \code{\link{grow_Maple}}, \code{\link{plotPruning}},
#'   \code{\link{plotDiagnostics}}
#'
#' @export
plotTree <- function(x,
                     colours = c("#7A111E", "#D7263D", "#E76F51",
                                 "#E9C46A", "#A3B14A", "#4C9A2A")) {

  tree     <- .get_tree(x)
  igraph::V(tree)["0_1"]$stability <- NA

  p <- ggraph::ggraph(tree, layout = "sugiyama") +
    ggraph::geom_edge_link(ggplot2::aes(alpha = in_prop),
                           colour      = "grey50",
                           show.legend = TRUE) +
    ggraph::geom_node_point(ggplot2::aes(size  = size,
                                         color = stability)) +
    ggraph::geom_node_text(ggplot2::aes(label = cluster,
                                        size  = size * 0.6),
                           vjust    = 0.4,
                           color    = "black",
                           fontface = "bold") +
    ggplot2::scale_color_gradientn(colours  = colours,
                                   limits   = c(0, 1),
                                   na.value = "grey80",
                                   name     = "Stability") +
    ggraph::theme_graph()

  return(p)
}

# =============================================================================
# plotPruning()
# =============================================================================

#' Plot the Maple pruning selection
#'
#' Visualises the result of \code{cut_Maple()} by highlighting the optimal
#' frontier nodes, the path leading to them, and the remaining tree nodes
#' in distinct colours.
#'
#' @param x A \code{SingleCellExperiment} or \code{MapleResult} object on
#'   which \code{Maple()}, \code{grow_Maple()} and \code{cut_Maple()} have
#'   been run.
#' @param path_color Character. Colour for path nodes and edges.
#'   Default \code{"#D7263D"}.
#' @param tree_color Character. Colour for non-path tree nodes and edges.
#'   Default \code{"grey80"}.
#' @param optimal_color Character. Colour for optimal frontier nodes.
#'   Default \code{"#7A111E"}.
#'
#' @return A \code{ggplot2} object.
#'
#' @seealso \code{\link{cut_Maple}}, \code{\link{plotTree}}
#'
#' @export
plotPruning <- function(x,
                        path_color    = "#D7263D",
                        tree_color    = "grey80",
                        optimal_color = "#7A111E") {

  tree <- .get_pruned_tree(x)

  p <- ggraph::ggraph(tree, layout = "sugiyama") +
    ggraph::geom_edge_link(
      ggplot2::aes(edge_colour = is_path),
      show.legend = FALSE) +
    ggraph::geom_node_point(
      ggplot2::aes(size  = size,
                   color = node_status),
      shape       = 16,
      show.legend = FALSE) +
    ggplot2::scale_color_manual(
      values = c("Tree"    = tree_color,
                 "Path"    = path_color,
                 "Optimal" = optimal_color)) +
    ggraph::scale_edge_colour_manual(
      values = c("TRUE"  = path_color,
                 "FALSE" = tree_color)) +
    ggraph::theme_graph()

  return(p)
}

# =============================================================================
# plotDiagnostics()
# =============================================================================

#' Plot Maple diagnostic distributions
#'
#' Visualises key Maple metrics to inform parameter choices before running
#' \code{grow_Maple()} and \code{cut_Maple()}. Returns a two-panel figure:
#'
#' \itemize{
#'   \item Left panel: violin plots with overlaid boxplots of the
#'     in-proportion distribution (to inform \code{in_prop_threshold}) and
#'     the per-cell observational stability distribution (to inform
#'     \code{stability_threshold}).
#'   \item Right panel: minimum node stability per resolution across the
#'     grid, with a dashed horizontal line at the adaptive threshold
#'     (\code{median(scores)}). Helps identify which resolutions produce
#'     reliable clusters and whether the adaptive threshold is appropriate.
#' }
#'
#' The virtual root node and its edges (resolution 0) are excluded from all
#' distributions.
#'
#' @param x A \code{SingleCellExperiment} or \code{MapleResult} object on
#'   which \code{Maple()} has been run.
#' @param qlevel Numeric between 0 and 1. Quantile used to compute node stability
#'   scores for the diagnostic tree. Default 0.2.
#' @param colours Character vector of 2 colours for the left panel
#'   distributions. Default is an orange-blue palette.
#'
#' @return A combined \code{patchwork} figure.
#'
#' @seealso \code{\link{Maple}}, \code{\link{grow_Maple}},
#'   \code{\link{cut_Maple}}
#'
#' @export
plotDiagnostics <- function(x,
                            qlevel  = 0.2,
                            colours = c("#E76F51", "#457B9D")) {

  Maple_res <- .get_maple_result(x)

  # Build temporary raw tree (in_prop_threshold = 0) for diagnostics only
  raw_tree <- create_tree(
    reference         = Maple_res$reference,
    scores            = Maple_res$scores,
    grid              = Maple_res$grid,
    in_prop_threshold = 0,
    qlevel            = qlevel
  )

  edges_df <- igraph::as_data_frame(raw_tree, what = "edges")
  nodes_df <- igraph::as_data_frame(raw_tree, what = "vertices")

  # Exclude virtual root edges and nodes
  edges_df <- edges_df[edges_df$from != "0_1", ]
  nodes_df <- nodes_df[nodes_df$resolution != 0, ]

  # Exclude virtual root column from observational stability
  obs_scores <- as.numeric(Maple_res$scores[, -1])

  # Adaptive threshold
  adaptive_threshold <- median(obs_scores, na.rm = TRUE)

  # --- Panel 1: violin/boxplot of In Proportion and Obs Stability ---
  plot_df <- do.call(rbind, list(
    data.frame(Metric = "In Proportion", Value = edges_df$in_prop),
    data.frame(Metric = "Obs Stability", Value = obs_scores)
  ))
  plot_df        <- plot_df[!is.na(plot_df$Value), ]
  plot_df$Metric <- factor(plot_df$Metric,
                           levels = c("In Proportion", "Obs Stability"))

  p1 <- ggplot2::ggplot(plot_df,
                        ggplot2::aes(x = Metric, y = Value,
                                     fill = Metric)) +
    ggplot2::geom_violin(trim         = TRUE,
                         alpha        = 0.7,
                         scale        = "width",
                         color        = NA) +
    ggplot2::geom_boxplot(width         = 0.1,
                          fill          = "white",
                          color         = "black",
                          outlier.shape = NA,
                          alpha         = 0.8) +
    ggplot2::scale_fill_manual(values = colours) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                breaks = seq(0, 1, 0.2)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x    = ggplot2::element_blank(),
      axis.text.x     = ggplot2::element_text(size = 11, face = "bold"),
      axis.title.y    = ggplot2::element_text(size = 11, face = "bold"),
      plot.title      = ggplot2::element_text(face = "bold", hjust = 0.5)
    ) +
    ggplot2::labs(y     = "Score / Proportion",
                  title = "Distributions")

  # --- Panel 2: minimum node stability per resolution ---
  min_stab_df <- aggregate(stability ~ resolution,
                           data = nodes_df,
                           FUN  = min)
  min_stab_df <- min_stab_df[order(min_stab_df$resolution), ]

  p2 <- ggplot2::ggplot(min_stab_df,
                        ggplot2::aes(x = resolution, y = stability)) +
    ggplot2::geom_line(color     = "#1D3557",
                       linewidth = 0.8) +
    ggplot2::geom_point(color = "#1D3557",
                        size  = 2) +
    ggplot2::geom_hline(yintercept = adaptive_threshold,
                        linetype   = "dashed",
                        color      = "#D7263D",
                        linewidth  = 0.7) +
    ggplot2::annotate("text",
                      x     = max(min_stab_df$resolution),
                      y     = adaptive_threshold,
                      vjust = -0.5,
                      label = paste0("adaptive = ",
                                     round(adaptive_threshold, 3)),
                      color = "#D7263D",
                      hjust = 1,
                      size  = 3.5) +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                breaks = seq(0, 1, 0.2)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 11, face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    ) +
    ggplot2::labs(x     = "Resolution",
                  y     = "Min Node Stability",
                  title = "Stability across Resolution Grid")

  # --- Combine panels ---
  patchwork::wrap_plots(p1, p2)
}

