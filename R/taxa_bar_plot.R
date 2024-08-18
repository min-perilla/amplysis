#' @title Visualization of Species Stacked Plot
#'
#' @description
#' After processing with the `taxa_bar()` function, you can use this function
#' for visualization. Note that the parameters `group1` and `group2` in
#' `taxa_bar()` correspond to the parameters `x_group` and `facet_group` in
#' `taxa_bar_plot()`, respectively, just with different names.
#'
#' @param data Plot data. Please use the function `taxa_bar()` to generate plot
#' data.
#' @param color_scheme (character) Color scheme. Please enter a hexadecimal
#' format color scheme character vector, such as color_scheme = c("#7FC97F",
#' "#BEAED4", "#FDC086", "#FFFF99", "#386CB0"). If NULL, the default color
#' scheme will be used.
#' @param tax_cla (character) Taxonomic level. Only column names from the tax
#' table can be entered, for example tax_cla = 'genus'.
#' @param x_group (Required, character) Group 1, please enter the column name of
#' the grouping information in the metadata table.
#' @param facet_group (Optional, character) Group 2 for facetting plots, please
#' enter the column name of the grouping information in the metadata table.
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
#' @return A graph of ggplot2 class
#' @export
#'
#' @examples
#' \dontrun{taxa_bar_plot(data = taxa, color_scheme = NULL, tax_cla = "genus",
#' x_group = "group", facet_group = "group2", y_abun = "abun",
#' custom_order = c("AA", "AB", "AR", "AD", "BA", "BB", "BR", "BD", "Soil"),
#' custom_order_F = c("A", "B", "S"),
#' bar_type = "fill", bar_width = 0.7, grid_line = F, legend_ncol = 1,
#' title = "Species Stacked Plot", title_x = "Groups", title_y = NULL,
#' title_legend = "Enrichment Culture \nTop 20 Genera",
#' size_title = 28, size_title_x = 20, size_title_y = 20,
#' size_title_legned = 24, size_title_facet = 32, size_x = 14, size_y = 18,
#' size_legned = 16,
#' filename = NULL, file_width = 16, file_height = 9)
#' }
#'
#' @importFrom dplyr arrange mutate
#' @importFrom forcats fct_relevel
#' @importFrom ggplot2 aes element_blank element_line element_text facet_grid
#' geom_bar ggplot ggsave guides guide_legend labs margin scale_fill_manual
#' scale_y_continuous theme_bw unit
#' @importFrom rlang quo_text sym :=
#' @importFrom stringr str_to_title
#'
# Please use tools:: showNonASCIIfile(file.R) to check for the presence of non ASCII characters.
# tools::showNonASCIIfile(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "taxa_bar_plot.R"))
taxa_bar_plot <- function(data, color_scheme = NULL, tax_cla, x_group = "group",
         facet_group = NULL, y_abun = "abun", custom_order = NULL,
         custom_order_F = NULL,

         bar_type = "fill", bar_width = 0.7, grid_line = F,
         size_point_legend = 0.75, spacing_legend_point = 0.75,
         spacing_legend_title = 0.5, legend_ncol = 1,

         title = NULL, title_x = "Groups", title_y = NULL, title_legend = NULL,

         size_title = 28, size_title_x = 20, size_title_y = 20,
         size_title_legned = 24, size_title_facet = 32, size_x = 14,
         size_y = 18, size_legned = 16,
         filename = NULL, file_width = 16, file_height = 9)
{
  ##
  if(!is.null(custom_order)) {
    x_group_sym <- rlang::sym(x_group)

    data <- data %>%
      dplyr::mutate(!!x_group_sym := factor(!!x_group_sym,
                                            levels = custom_order)) %>%
      dplyr::arrange(!!x_group_sym)

    cat("Custom horizontal axis order: ", custom_order, "\n")
  }


  ##
  if(is.null(facet_group)){
    facet_grid2 = ". ~ ."


  } else {
    facet_grid2 = paste0(". ~ ", facet_group)

    if(!is.null(custom_order_F)) {
      facet_group_sym <- rlang::sym(facet_group)

      data <- data %>%
        dplyr::mutate(!!facet_group_sym := factor(!!facet_group_sym,
                                                  levels = custom_order_F)) %>%
        dplyr::arrange(!!facet_group_sym)

      cat("Custom faceted plot horizontal axis order: ", custom_order_F, "\n")
    }
  }


  ##
  if(bar_type == "fill") {
    label_y = ggplot2::scale_y_continuous(
      name = if(is.null(title_y)){ "Relative abundance (%)" } else { title_y },
      limits = c(0,1),
      breaks = seq(0, 1, 0.25),
      labels = paste(seq(0, 100, 25), "%"))


  } else if(bar_type == "stack") {
    label_y = ggplot2::scale_y_continuous(
      name = if(is.null(title_y)){ "Relative abundance (%)" } else { title_y })


  } else {
    stop("The formal parameter 'bar_type' is entered incorrectly. \n",
         "Please enter the correct parameter: \"fill\" or \"stack\"!")
  }


  ##
  if(isFALSE(grid_line)){
    panel_grid = ggplot2::element_blank()


  } else {
    panel_grid = ggplot2::element_line(color = "gray", linewidth = 0.2,
                                       linetype = "dashed")
  }


  ##
  if(is.null(title_legend)) {
    title_legend = deparse(substitute(tax_cla))
    title_legend <- tolower(title_legend)
    title_legend <- stringr::str_to_title(title_legend)

  } else { }


  ##
  tax_cla <- rlang::sym(tax_cla)
  x_group <- rlang::sym(x_group)
  y_abun <- rlang::sym(y_abun)


  ##
  if(is.null(color_scheme)) {
    color_scheme = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0",
                     "#F0027F", "#1ba7e5", "#66C2A5", "#FC8D62", "#8DA0CB",
                     "#E78AC3", "#FFD92F", "#E5C494", "#95beae", "#FF7F00",
                     "#E31A1C", "#FB9A99", "#33A02C", "#A6CEE3", "#95be3e",
                     "#7F217F", "#2EAED4", "#ADC086", "#BFFF99", "#C86CB0",
                     "#DE023F", "#E8a745", "#F6C2A5", "#AC2D62", "#ADA0CB",
                     "#B78AC3", "#CFD92F", "#D5C494", "#E5beae", "#AF7F00",
                     "#B31A1C", "#CB9A99", "#D3A02C", "#EE6CE3", "#F12e3e",
                     "#413dAC", "#2EAEA4", "#B43086", "#BFDF99", "#C84A30"
    )

    color_n <- length((unique(data[[tax_cla]])))
    color_scheme2 <- color_scheme[1:color_n]
    color_scheme2[color_n] <- "#7f7f7f"
  }


  ##
  p <- ggplot2::ggplot() +


    ##
    ggplot2::geom_bar(data = data,
                      ggplot2::aes(
                        x = !!x_group,
                        weight = !!y_abun,
                        fill = forcats::fct_relevel(!!tax_cla, after = Inf,
                                                    "others")),
                      position = bar_type,
                      width = bar_width) +


    ##
    ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(values = color_scheme) +
    ggplot2::theme(panel.grid = panel_grid) +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 30, r = 40, b = 30, l = 40,
                                                 unit = "pt")) +


    ##
    ggplot2::labs(title = title) +
    ggplot2::theme(plot.title = ggplot2::element_text(
      face = "bold", size = size_title, hjust = 0.5)) +
    ggplot2::theme(plot.title = ggplot2::element_text(
      margin = ggplot2::margin(t = 0, r = 0, b = 15, l = 0, unit = "pt"))) +


    ##
    ggplot2::labs(x = title_x) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(
      face = "bold", size = size_title_x, color = "black")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      face = "bold", size = size_x, color="black")) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(
      margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0, unit = "pt"))) +


    ##
    label_y +
    ggplot2::theme(axis.title.y = ggplot2::element_text(
      face = "bold", size = size_title_y, color = "black")) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(
      face = "bold", size = size_y, color="black")) +
    ggplot2::theme(axis.title.y = ggplot2::element_text(
      margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0, unit = "pt"))) +


    ##
    ggplot2::guides(fill = ggplot2::guide_legend(
      title = title_legend, ncol = legend_ncol)) +
    ggplot2::theme(legend.position = "right") +
    ggplot2::theme(legend.title = ggplot2::element_text(
      face = "bold", size = size_title_legned, color="black")) +
    ggplot2::theme(legend.text = ggplot2::element_text(
      face = "bold.italic", size = size_legned, color="black")) +
    ggplot2::theme(legend.key.size = unit(size_point_legend, "cm")) +
    ggplot2::theme(legend.key.spacing.y = ggplot2::unit(
      spacing_legend_point, "pt")) +
    ggplot2::theme(legend.title = element_text(margin = ggplot2::margin(
      b = spacing_legend_title, unit = 'cm'))) +
    ggplot2::theme(legend.margin = ggplot2::margin(
      t = 0, r = 0, b = 0, l = 20, unit = "pt")) +


    ##
    ggplot2::facet_grid(facet_grid2, scales = "free") +
    ggplot2::theme(strip.text.x = ggplot2::element_text(
      face = "bold", size = size_title_facet, color="black"))


  ##
  if(is.null(filename) || filename == "") {
    filename <- as.character(rlang::quo_text(tax_cla))

  } else {  }
  ggplot2::ggsave(filename = paste0(filename, ".png"), plot = p,
                  width = file_width, height = file_height)
  ggplot2::ggsave(filename = paste0(filename, ".pdf"), plot = p,
                  width = file_width, height = file_height)

  ##
  cat("\033[32mtaxa_bar: success!\033[0m\n")
  cat("\033[0;32m", "The file has been saved to \n",
      getwd(), "\033[0m\n", sep = "")
  return(p)
}
