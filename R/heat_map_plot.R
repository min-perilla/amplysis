#' @title Heatmap Visualization
#'
#' @description
#' After performing heatmap analysis using the `heat_map()` function,
#' you will obtain a dataframe. This function allows you to visualize the
#' dataframe, thereby generating a high-quality heatmap.
#' If you want to customize the x-axis order of the heatmap, please input a
#' vector into the parameter cluster_cols. For example,
#' cluster_cols = c(4, 5, 6, 1, 2, 3), indicates sorting according to column
#' numbers 4 5 6 1 2 3.
#'
#' @param mat (data.frame) Plotting Data
#' @param scale (character) scale is used to set normalization. "row" represents
#' row-wise normalization, "column" represents column-wise normalization, and
#' "none" represents no normalization.
#' @param cellwidth (numeric) Translation: Represents the width of a single
#' cell, default is "NA".
#' @param cellheight (numeric) Translation: Represents the height of a single
#' cell, default is "NA".
#' @param color (character) The heatmap cell colors are generated automatically
#' with a gradient.
#' For example: c("#2196f3", "#a8d1f2", "#f4faff", "#ec9fa2", "#ec1c24")
#' @param cluster_cols (logical or integer vector) Whether to enable column
#' clustering (when clustering is enabled,
#' custom sorting is not possible). If you want to customize the x-axis order of
#' the heatmap,
#' please input a vector into the parameter cluster_cols.
#' For example, cluster_cols = c(4, 5, 6, 1, 2, 3), indicates sorting according
#' to column numbers 4 5 6 1 2 3.
#' @param clustering_method (character) "Represents clustering methods,
#' including: 'ward.D', ward.D2', 'single', 'complete', 'average', 'mcquitty',
#' 'median', 'centroid'"
#' @param treeheight_row (numeric) Row clustering tree height
#' @param treeheight_col (numeric) Col clustering tree height
#'
#' @param filename (character) File name for saving
#' @param file_width (numeric) Image width
#' @param file_height (numeric) Image height
#'
#' @return pheatmap
#' @export
#'
#' @examples
#' \dontrun{heat_map_plot(mat = heatmap1, scale = "row", cellwidth = NA,
#' cellheight = NA, color =  c("#2196f3", "#a8d1f2", "#f4faff", "#ec9fa2",
#' "#ec1c24"), cluster_cols = F, clustering_method = "ward.D",
#' treeheight_row = 100, treeheight_col = 50, filename = "heatmap.png")}
#'
#' @importFrom stats dist hclust
#'
# Please use tools:: showNonASCIIfile(file.R) to check for the presence of non ASCII characters.
# tools::showNonASCIIfile(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "heat_map_plot.R"))
heat_map_plot <- function(mat, scale = "row", cellwidth = NA, cellheight = NA,
         color =  c("#2196f3", "#a8d1f2", "#f4faff", "#ec9fa2", "#ec1c24"),
         cluster_cols = F, clustering_method = "ward.D", treeheight_row = 50,
         treeheight_col = 50, filename = "heatmap.png", file_width = 12,
         file_height = 12)
{
  valid_methods <- c("ward.D", "ward.D2", "single", "complete", "average",
                     "mcquitty", "median", "centroid")

  if (!(clustering_method %in% valid_methods)) {
    stop("Invalid clustering method. Please choose from \n",
         "'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty',
         'median', or 'centroid'.")
  } else {}

  color =  grDevices::colorRampPalette(color)(100)

  if(isTRUE(cluster_cols)) {
    cat("Column clustering: Enabled\n")
    cat("\033[32m", "If you want to customize the x-axis order of the heatmap,
        please input a vector into the parameter cluster_cols.\n", "For example,
        cluster_cols = c(4, 5, 6, 1, 2, 3), indicates sorting according to
        column numbers 4 5 6 1 2 3.", "\033[0m", sep = "")

  } else if(isFALSE(cluster_cols)) {
    cat("Column clustering: Disabled\n")
    cat("\033[32m", "If you want to customize the x-axis order of the heatmap,
        please input a vector into the parameter cluster_cols.\n", "For example,
        cluster_cols = c(4, 5, 6, 1, 2, 3), indicates sorting according to
        column numbers 4 5 6 1 2 3.", "\033[0m", sep = "")

  } else if(is.vector(cluster_cols)) {
    cat("Column clustering: Disabled\n", "custom x-axis order: ", sep = "")
    cat(cluster_cols, "\n")

    matrix2 <- round(t(apply(mat, MARGIN = 1, base::scale)), 2)
    colnames(matrix2) <- colnames(mat)
    exprTable <- as.data.frame(t(matrix2))
    row_dist = stats::dist(exprTable,method = "euclidean")
    hclust_1 <- stats::hclust(row_dist)
    hclust_1[[3]] <- cluster_cols
    cluster_cols <- hclust_1
  }



  ###
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("The 'pheatmap' package is required but not installed. Please install
         it using 'install.packages(\"pheatmap\")'.")
  }


  ##
  p <- pheatmap::pheatmap(
    mat = mat,
    scale = scale,
    border_color = "white",
    cellwidth = cellwidth, cellheight = cellheight,
    color = color,

    cluster_cols = cluster_cols,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = clustering_method,
    cutree_rows = NA, cutree_cols = 4,
    treeheight_row = treeheight_row,
    treeheight_col = treeheight_col,
    gaps_row = NULL,

    labels_row = NULL,

    angle_col='0',

    fontsize = 16,
    fontsize_row = 16,
    fontsize_col = 16,
    legend_labels = NA,

    filename = filename, width = file_width, height = file_height
  )

  cat("\033[32mheatmap: success!\033[0m\n")
  cat("\033[0;32m", "The file \"", filename, "\" has been saved to \n",
      getwd(), "/", filename, "\033[0m\n", sep = "")

  return(p)
}
