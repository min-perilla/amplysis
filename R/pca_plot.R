#' @title Visualization of PCA
#'
#' @description
#' After processing with the `pca()` function, you can use this function
#' for visualization.
#'
#' @param data Plotting data.
#' @param color_scheme (character) Color scheme.
#' @param custom_order (character) Custom legend order.
#' @param seed Seed
#'
#' @param size_point (numeric) The size of points.
#' @param size_point_legend (numeric) The size of legend points.
#' @param spacing_legend_point (numeric) The internal padding of the legend.
#' @param spacing_legend_title (numeric) The spacing between the legend title and the body
#' @param legend_ncol (integer) Number of columns in the legend.
#' @param label_is (logical) Whether to display data labels.
#' @param size_label (numeric) Font size of the data label.
#' @param label_font_color (character) The label font color should default to the group color.
#' @param ellipse_type (character / numeric) The method for calculating the confidence ellipse.
#' Three options are available (can also be selected using numeric values):
#'
#' 1 - "t" (default): Computes the ellipse based on the t-distribution,
#'     suitable for small samples, using robust estimation (MASS::cov.trob()).
#'
#' 2 - "norm": Computes the ellipse based on the normal distribution,
#'     using the covariance matrix estimation without robust estimation.
#'
#' 3 - "euclid": Draws a fixed-radius circle based on Euclidean distance,
#'     which is related to the data scale.
#'
#' @param title (character) Main title.
#' @param title_sub (character) Subtitle.
#' @param title_legend (character) Legend title.
#'
#' @param size_title (numeric) Font size of the main title.
#' @param size_title_sub (numeric) Font size of the subtitle.
#' @param size_title_x (numeric) Font size of the horizontal axis title.
#' @param size_title_y (numeric) Font size of the vertical axis title.
#' @param size_title_legend (numeric) Font size of legend title.
#'
#' @param size_x (numeric) Font size of the horizontal axis tick labels.
#' @param size_y (numeric) Font size of the vertical axis tick labels.
#' @param size_legend (numeric) Font size of the legend.
#'
#' @param filename (character) File name for saving.
#' @param file_width (numeric) Width of the image.
#' @param file_height (numeric) Height of the image.
#'
#' @return A graph of ggplot2 class
#' @export
#'
#' @examples
#' \dontrun{
#' pca1 <- pca(otu = otu, metadata = metadata, id_col = 1,group = "group",
#'             parallel_method = "none")
#'
#' pca_plot(data = pca1,
#'          color_scheme = c("#00b0f6", "#FFC24B", "#f8766d",
#'                           "#ae876d", "#AFC24B"),
#'          seed = 123, custom_order = c("A", "B", "R", "D", "S"),
#'
#'          size_point = 8, size_point_legend = 8, spacing_legend_point = 1.2,
#'          spacing_legend_title = 0.5, legend_ncol = 1, label_is = T,
#'          size_label = 5, label_font_color = NULL, ellipse_type = "t",
#'
#'          title = "RDA", title_sub = NULL, title_legend = "Group",
#'
#'          size_title = 28, size_title_sub = 16, size_title_x = 20,
#'          size_title_y = 20, size_title_legend = 24,
#'
#'          size_x = 16, size_y = 16, size_legend = 16,
#'
#'          filename = "RDA", file_width = 12, file_height = 9)
#' }
#'
#' @importFrom dplyr arrange
#' @importFrom ggplot2 aes element_blank element_text ggplot ggsave guides
#' guide_legend labs margin scale_fill_manual stat_ellipse theme_bw unit
#' @importFrom rlang sym
#'
# Please use tools:: showNonASCIIfile(file.R) to check for the presence of non ASCII characters.
# tools::showNonASCIIfile(file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "pca_plot.R"))
pca_plot <- function(
    data,
    color_scheme = NULL,
    custom_order = NULL,
    seed = 123,

    size_point = 4.5,
    size_point_legend = 8,
    spacing_legend_point = 1.2,
    spacing_legend_title = 0.5,
    legend_ncol = 1,
    label_is = T,
    size_label = 5,
    label_font_color = NULL,
    ellipse_type = "t",

    title = "PCA",
    title_sub = NULL,
    title_legend = "Group",

    size_title = 28,
    size_title_sub = 16,
    size_title_x = 20,
    size_title_y = 20,
    size_title_legend = 24,

    size_x = 18,
    size_y = 18,
    size_legend = 16,

    filename = "PCA",
    file_width = 12,
    file_height = 9)
{
  #
  set.seed(seed = seed)
  group = "group"

  #
  if(isTRUE(label_is)) {
    labelNum = 200
  }else{
    labelNum = 0
  }

  #
  if(!is.null(custom_order)){
    data[["PCA"]][["group"]] <- factor(data[["PCA"]][["group"]],
                                       levels = custom_order)
    data[["PCA"]] <- dplyr::arrange(data[["PCA"]], group)

    cat("Custom legend order: ", custom_order, sep = "")
  }


  # ----------------------------------------------------------------------------
  ## Confidence ellipse calculation method
  # Three methods for calculating confidence ellipses
  ellipse_methods <- c("t", "norm", "euclid")

  # Descriptions of the calculation methods
  ellipse_descriptions <- list(
    "t" = "Ellipse based on t-distribution, suitable for small samples, using robust estimation",
    "norm" = "Ellipse based on normal distribution, using covariance matrix estimation, without robust estimation",
    "euclid" = "Circle with a fixed radius based on Euclidean distance, related to data scale"
  )

  # If the input is numeric, determine the corresponding method using modulo operation
  if (is.numeric(ellipse_type)) {
    index <- (ellipse_type - 1) %% length(ellipse_methods) + 1
    ellipse_type <- ellipse_methods[index]
  }

  # If the input is not a valid method, default to "t"
  if (!(ellipse_type %in% ellipse_methods)) {
    ellipse_type <- "t"
  }

  # Output the selected confidence ellipse calculation method
  cat("\nConfidence ellipse calculation method: ellipse_type = \"", ellipse_type, "\" ",
      "(", ellipse_descriptions[[ellipse_type]], ")\n", sep = "")

  # Output other available methods (excluding the selected one)
  remaining_methods <- setdiff(ellipse_methods, ellipse_type)
  if (length(remaining_methods) > 0) {
    cat("Other available methods (can also be selected by entering the corresponding number):\n")
    for (method in ellipse_methods) {
      if (method != ellipse_type) {
        index <- which(ellipse_methods == method)
        cat(index, ":", method, "-", ellipse_descriptions[[method]], "\n")
      }
    }
  }
  cat("\n")
  # ----------------------------------------------------------------------------


  pc = NULL

  #
  pc[1] <- round(data[["PoA"]][1], 2)
  pc[2] <- round(data[["PoA"]][2], 2)

  #
  xName = paste0("PC1 (", pc[1], "%)")
  yName = paste0("PC2 (", pc[2], "%)")


  #
  group_sym <- rlang::sym(group)


  ##
  #
  p1 <- ggplot2::ggplot(
    #
    data = data[["PCA"]],
    ggplot2::aes(x = data[["PCA"]][["PC1"]],
                 y = data[["PCA"]][["PC2"]],
                 color = !!group_sym, shape = !!group_sym)) +

    #
    ggplot2::theme_bw() +

    #
    ggplot2::geom_vline(xintercept = 0, lty = "dashed", alpha = 0.2) +
    ggplot2::geom_hline(yintercept = 0, lty = "dashed", alpha = 0.2) +
    ggplot2::geom_point(size = size_point) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())


  if(is.null(label_font_color)) {
    p1 <- p1 +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = sample), size = size_label,
        box.padding = ggplot2::unit(0.6, "lines"),
        point.padding = ggplot2::unit(0.5, "lines"),
        max.overlaps = labelNum, alpha = 0.8, show.legend = F, seed = seed)
  } else {
    p1 <- p1 +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = sample),

        color = label_font_color,  #

        size = size_label,
        box.padding = ggplot2::unit(0.6, "lines"),
        point.padding = ggplot2::unit(0.5, "lines"),
        max.overlaps = labelNum, alpha = 0.8, show.legend = F, seed = seed)
  }



  #
  p1 <- p1 +
    ggplot2::labs(x = xName,
                  y = yName) +

    #
    ggplot2::stat_ellipse(data = data[["PCA"]],
                          geom = "polygon",
                          level = 0.95,
                          linetype = 2,
                          linewidth = 0.4,
                          ggplot2::aes(fill = group),
                          alpha = 0.15,
                          show.legend = F,   # Do not display legend
                          # type = "t"       # Default: Ellipse based on t-distribution (for small samples, using robust estimation MASS::cov.trob())
                          # type = "norm"    # Ellipse based on normal distribution (calculated using covariance matrix, without robust estimation)
                          # type = "euclid"  # Circle with a fixed radius drawn using Euclidean distance (depends on data scale)
                          type = ellipse_type  # Calculation method for the confidence ellipse

    ) +

    #
    ggplot2::labs(title = title) +
    ggplot2::theme(plot.title = ggplot2::element_text(
      face = "bold", size = size_title, hjust = 0.5)) +

    #
    ggplot2::labs(subtitle = title_sub) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(
      face = "bold", size = size_title_sub, hjust = 0)) +

    #
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = size_title_x),
                   axis.title.y = ggplot2::element_text(size = size_title_y, angle = 90),
                   axis.text.x = ggplot2::element_text(size = size_x),
                   axis.text.y = ggplot2::element_text(size = size_y)) +

    #
    ggplot2::guides(
      shape = "none",
      color = ggplot2::guide_legend(
        title = title_legend,
        ncol = legend_ncol,
        override.aes = list(size = size_point_legend))) +
    ggplot2::theme(legend.title = ggplot2::element_text(
      face = "bold", size = size_title_legend, color = "black")) +
    ggplot2::theme(legend.text = ggplot2::element_text(
      face = "bold", size = size_legend, color = "black")) +

    #
    ggplot2::theme(legend.text = ggplot2::element_text(
      margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"))) +
    ggplot2::theme(legend.title = ggplot2::element_text(hjust = 0.5)) +

    #
    ggplot2::theme(legend.key.height = ggplot2::unit(
      spacing_legend_point, "cm")) +
    #
    ggplot2::theme(legend.title = element_text(
      margin = ggplot2::margin(b = spacing_legend_title, unit = 'cm'))) +


    #
    ggplot2::theme(plot.title = ggplot2::element_text(
      margin = ggplot2::margin(t = 0, r = 0, b = 15, l = 0, unit = "pt"))) +
    #
    ggplot2::theme(axis.title.x = ggplot2::element_text(
      margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0, unit = "pt"))) +
    #
    ggplot2::theme(axis.title.y = ggplot2::element_text(
      margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0, unit = "pt"))) +
    #
    ggplot2::theme(legend.margin = ggplot2::margin(
      t = 0, r = 0, b = 0, l = 20, unit = "pt")) +

    #
    ggplot2::theme(plot.margin = ggplot2::margin(
      t = 20, r = 30, b = 20, l = 30, unit = "pt"))


  ##
  if(!is.null(color_scheme)) {
    #
    color_scheme_point <- color_scheme
    color_scheme_ellipse <- color_scheme_point

    p1 <- p1 +
      ggplot2::scale_color_manual(values = color_scheme_point) +
      ggplot2::scale_fill_manual(values = color_scheme_ellipse)
  }


  ##
  ggplot2::ggsave(filename = paste0(filename, ".png"), plot = p1,
                  width = file_width, height = file_height)
  ggplot2::ggsave(filename = paste0(filename, ".pdf"), plot = p1,
                  width = file_width, height = file_height)

  ##
  cat("\033[32mtaxa_bar: success!\033[0m\n")
  cat("\033[0;32m", "The file has been saved to \n",
      getwd(), "\033[0m\n", sep = "")

  return(p1)
}
