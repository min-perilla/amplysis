#' @title R6 Class: class_taxa_bar
#'
#' @description
#' The class_taxa_bar class encapsulates the species composition analysis
#' function `taxa_bar()` and its accompanying visualization function
#' `taxa_bar_plot()`, enabling convenient analysis.
#'
#' @return class_taxa_bar: otu, tax, metadata, data, chart
#' @export
#'
#' @importFrom R6 R6Class
class_taxa_bar <- R6::R6Class(
  "class_taxa_bar",

  #
  inherit = otu_tax_metadata_data_chart,


  ##
  public = list(

    #' @description
    #' The `plot()` function integrates the `taxa_bar()` and `taxa_bar_plot()`
    #' functions, allowing for one-click plotting of species stacked plots.
    #'
    #' @param id_col (integer)The column number of the OTU ID column in the
    #' OTU table, by default, is 1.
    #' @param tax_cla (character) Taxonomic level. Only column names from the tax
    #' table can be entered, for example tax_cla = 'genus'.
    #' @param group1 (Required, character) Group 1, please enter the column name of
    #' the grouping information in the metadata table.
    #' @param group2 (Optional, character) Group 2 for facetting plots, please enter
    #' the column name of the grouping information in the metadata table.
    #' @param parallel_method (character) Parallel sample processing method,
    #' defaulting to mean. Options: mean (average), sum (summation), median (median).
    #' @param row_n (integer) Preserve the top N taxa (including the Nth) based on
    #' abundance, while merging taxa with lower abundance into "others".
    #'
    #'
    #' @param color_scheme (character) Color scheme. Please enter a hexadecimal
    #' format color scheme character vector, such as color_scheme = c("#7FC97F",
    #' "#BEAED4", "#FDC086", "#FFFF99", "#386CB0"). If NULL, the default color
    #' scheme will be used.
    #'
    #' @param y_abun (character) Abundance information. Please enter the column name
    #' for "Abundance Information" in metadata.
    #' @param custom_order (character) Customize the order of horizontal axis grouping.
    #' @param custom_order_F (character) (Facet plot) Customize the order of
    #' horizontal axis grouping.
    #'
    #' @param bar_type (character) Type of species stacked plot. Options: "fill"
    #' (relative abundance) or "stack" (absolute abundance).
    #' @param bar_width (numeric) Width of each column in the stacked plot.
    #' @param grid_line (logical) Grid lines on the background of the stacked plot.
    #' @param size_point_legend (numeric) Size of legend points.
    #' @param spacing_legend_point (numeric) The spacing inside the legend.
    #' @param spacing_legend_title (numeric) The spacing between the legend title
    #' and the main text.
    #' @param legend_ncol (integer) Number of columns in the legend.
    #'
    #' @param title (character) Main title.
    #' @param title_x (character) Horizontal axis title.
    #' @param title_y (character) Vertical axis title.
    #' @param title_legend (character) Legend title.
    #'
    #' @param size_title (numeric) Font size of the main title.
    #' @param size_title_x (numeric) Font size of the horizontal axis title.
    #' @param size_title_y (numeric) Font size of the vertical axis title.
    #' @param size_title_legned (numeric) Font size of legend title.
    #' @param size_title_facet (numeric) Font size of the facet plot titles.
    #'
    #' @param size_x (numeric) Font size of the horizontal axis tick labels.
    #' @param size_y (numeric) Font size of the vertical axis tick labels.
    #' @param size_legned (numeric) Font size of the legend.
    #'
    #' @param filename (character) File name for saving.
    #' @param file_width (numeric) Width of the image.
    #' @param file_height (numeric) Height of the image.
    #'
    #' @return data and chart.
    #' @export
    #'
    #' @examples
    #' \dontrun{
    #' a <- class_taxa_bar$new(otu, tax, metadata)
    #' a$plot()
    #' }
    #'
    plot = function(id_col = 1,
                    tax_cla = "phylum",
                    group1 = "group",
                    group2 = NULL,
                    parallel_method = "mean",
                    row_n = 8,

                    color_scheme = NULL,
                    y_abun = "abun",
                    custom_order = NULL,
                    custom_order_F = NULL,
                    bar_type = "fill",
                    bar_width = 0.7,
                    grid_line = F,
                    size_point_legend = 0.75,
                    spacing_legend_point = 0.75,
                    spacing_legend_title = 0.5,
                    legend_ncol = 1,

                    title = NULL,
                    title_x = "Groups",
                    title_y = NULL,
                    title_legend = NULL,

                    size_title = 28,
                    size_title_x = 20,
                    size_title_y = 20,
                    size_title_legned = 24,
                    size_title_facet = 32,

                    size_x = 14,
                    size_y = 18,
                    size_legned = 16,

                    filename = NULL,
                    file_width = 16,
                    file_height = 9
    )
    {
      # amplysis::taxa_bar
      taxa <- amplysis::taxa_bar(
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


      #
      if(is.null(title_legend)){
        title_legend <- tax_cla
        title_legend <- tolower(title_legend)
        title_legend <- stringr::str_to_title(title_legend)
      }


      # amplysis::taxa_bar_plot
      p <- amplysis::taxa_bar_plot(
        data = taxa,
        color_scheme = color_scheme,
        tax_cla = tax_cla,
        x_group = group1,
        facet_group = group2,
        y_abun = y_abun,
        custom_order = custom_order,
        custom_order_F = custom_order_F,

        bar_type = bar_type,
        bar_width = bar_width,
        grid_line = grid_line,
        size_point_legend = spacing_legend_title,
        spacing_legend_point = spacing_legend_title,
        spacing_legend_title = spacing_legend_title,
        legend_ncol = legend_ncol,

        title = title,
        title_x = title_x,
        title_y = title_y,
        title_legend = title_legend,

        size_title = size_title,
        size_title_x = size_title_x,
        size_title_y = size_title_y,
        size_title_legned = size_title_legned,
        size_title_facet = size_title_facet,

        size_x = size_x,
        size_y = size_y,
        size_legned = size_legned,

        filename = filename,
        file_width = file_width,
        file_height = file_height
      )

      #
      print(p)

      #
      private$.data <- taxa
      private$.chart <- p

      return("plot success!")
    }
  ),


  ##
  private = list(),
)
