#' @title R6 Class: class_heat_map
#'
#' @description
#' The class_heat_map class encapsulates the heatmap analysis function
#' `heat_map()` and its accompanying visualization function `heat_map_plot()`,
#' enabling convenient analysis.
#'
#' @return class_heat_map: otu, tax, metadata, data, chart
#' @export
#'
#' @importFrom R6 R6Class
class_heat_map = R6::R6Class(
  "class_heat_map",

  inherit = otu_tax_metadata_data_chart,


  ##
  public = list(

    #' @description
    #' The `plot()` function integrates the `heat_map()` and `heat_map_plot()`
    #' functions, allowing for one-click plotting of species stacked plots.
    #'
    #' @param id_col (integer) The column number of the OTU ID column in
    #' the OTU table, by default, is 1.
    #' @param tax_cla (character) Taxonomic level. Only column names from the
    #' tax table can be entered, for example tax_cla = 'genus'.
    #' @param group1 (character) Group 1, please enter the column name or column
    #' number of the grouping information in the metadata table.
    #' @param group2 (character) Group 2 for facetting plots, please enter the
    #' column name or column number of the grouping information in the metadata
    #' table.
    #' @param parallel_method (character) Parallel sample processing method,
    #' defaulting to mean. Options: mean (average), sum (summation),
    #' median (median).
    #' @param row_n (integer) Preserve the top N taxa (including the Nth) based
    #' on abundance.
    #'
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
    #' @param filename (character) File name for saving
    #' @param file_width (numeric) Image width
    #' @param file_height (numeric) Image height
    #'
    #' @return data and chart.
    #' @export
    #'
    #' @examples
    #' \dontrun{
    #' a <- class_heat_map$new(otu, tax, metadata)
    #' a$plot(group1 = "group2")
    #' }
    #'
    plot = function(id_col = 1,
                    tax_cla = "genus",
                    group1 = "group",
                    group2 = NULL,
                    parallel_method = "mean",
                    row_n = 35,

                    scale = "row",
                    cellwidth = NA,
                    cellheight = NA,
                    color = c("#2196f3", "#a8d1f2", "#f4faff", "#ec9fa2", "#ec1c24"),
                    cluster_cols = F,
                    clustering_method = "ward.D",

                    treeheight_row = 50,
                    treeheight_col = 50,

                    filename = "heatmap.png",
                    file_width = 12,
                    file_height = 12
    ) {

      # amplysis::heat_map
      heat_data <- amplysis::heat_map(
        otu = private$.otu,
        tax = private$.tax,
        metadata = private$.metadata,
        id_col = id_col,
        tax_cla = tax_cla,
        group1 = group1,
        group2 = group2,
        parallel_method = parallel_method,
        row_n = row_n
      )


      # amplysis::heat_map_plot
      p <- amplysis::heat_map_plot(
        mat = heat_data,
        scale = scale,
        cellwidth = cellwidth,
        cellheight = cellheight,
        color = color,
        cluster_cols = cluster_cols,
        clustering_method = clustering_method,
        treeheight_row = treeheight_row,
        treeheight_col = treeheight_col,
        filename = filename,
        file_width = file_width,
        file_height = file_height
      )


      #
      print(p)

      #
      private$.data <- heat_data
      private$.chart <- p

      return("plot success!")
    }
  ),


  ##
  private = list()
)
