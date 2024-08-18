#' @title Visualization of CCA
#'
#' @description
#' After processing with the `CCA()` function, you can use this function
#' for visualization.
#'
#' @param data Plotting data.
#' @param color_scheme (character) Color scheme.
#' @param group (character) Grouping information.
#' @param custom_order Custom legend order.
#' @param seed Seed
#'
#' @param size_point (numeric) The size of points.
#' @param size_point_legend (numeric) The size of legend points.
#' @param spacing_legend_point (numeric) The internal padding of the legend.
#' @param spacing_legend_title (numeric) The spacing between the legend title and the body
#' @param legend_ncol (integer) Number of columns in the legend.
#'
#' @param title (character) Main title.
#' @param title_legend (character) Legend title.
#'
#' @param size_title (numeric) Font size of the main title.
#' @param size_title_x (numeric) Font size of the horizontal axis title.
#' @param size_title_y (numeric) Font size of the vertical axis title.
#' @param size_title_legned (numeric) Font size of legend title.
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
#' \dontrun{
#' cca1 <- CCA(otu, env, metadata, 1, "group")
#' CCA_plot(data = cca1, group = "group", custom_order = NULL,
#'          color_scheme = c("#00b0f6", "#FFC24B", "#f8766d"), seed = 123,
#'
#'          size_point = 8, size_point_legend = 8, spacing_legend_point = 1.2,
#'          spacing_legend_title = 0.5, legend_ncol = 1,
#'
#'          title = "CCA", title_legend = "Group",
#'
#'          size_title = 28, size_title_x = 20, size_title_y = 20,
#'          size_title_legned = 24,
#'
#'          size_x = 16, size_y = 16, size_legned = 16,
#'
#'          filename = "CCA", file_width = 16, file_height = 9)
#' }
#'
#' @importFrom dplyr arrange
#' @importFrom ggplot2 aes element_blank element_text ggplot ggsave guides
#' guide_legend labs margin scale_fill_manual stat_ellipse theme_bw unit
#' @importFrom rlang sym
#'
CCA_plot = function(data, color_scheme = NULL, group = "group", custom_order = NULL, seed = 123,

         size_point = 4.5, size_point_legend = 8, spacing_legend_point = 1.2,
         spacing_legend_title = 0.5, legend_ncol = 1,

         title = "RDA", title_legend = "Group",

         size_title = 28, size_title_x = 20, size_title_y = 20,
         size_title_legned = 24,

         size_x = 18, size_y = 18, size_legned = 16,
         filename = "RDA", file_width = 16, file_height = 9)
{
  set.seed(seed = seed)

  x = "CCA1"
  y = "CCA2"

  x_contrib <- data[[x]]
  y_contrib <- data[[y]]

  group_sym <- rlang::sym(group)
  x_sym <- rlang::sym(x)
  y_sym <- rlang::sym(y)

  p1 <- ggplot2::ggplot(
    data = data[["data"]],
    ggplot2::aes(
      x = !!x_sym,   # x
      y = !!y_sym,   # y
      color = !!group_sym)) +

    ##
    ggplot2::theme_bw() +
    ggplot2::geom_vline(xintercept = 0, lty = "dashed") +
    ggplot2::geom_hline(yintercept = 0, lty = "dashed") +
    ggplot2::geom_point(size = size_point, shape = 16, alpha = 0.7) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = sample), size = 5,
      box.padding = ggplot2::unit(0.6, "lines"),
      point.padding = ggplot2::unit(0.5, "lines"), max.overlaps = 100,
      alpha = 0.8, show.legend = F) +
    ggplot2::labs(x = paste0("RDA1 (", x_contrib, "%)"),
                  y = paste0("RDA2 (", y_contrib, "%)")) +

    #
    ggplot2::stat_ellipse(
      data = data[["data"]],
      geom = "polygon",
      level = 0.95,
      linetype = 2,
      linewidth = 0.4,
      ggplot2::aes(fill = group),
      alpha = 0.1,
      show.legend = F
    ) +

    #
    ggplot2::labs(title = title) +
    ggplot2::theme(plot.title = ggplot2::element_text(
      face = "bold", size = size_title, hjust = 0.5)) +

    #
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = size_title_x),
                   axis.title.y = ggplot2::element_text(size = size_title_y, angle = 90),
                   axis.text.x = ggplot2::element_text(size = size_x),
                   axis.text.y = ggplot2::element_text(size = size_y),
                   panel.grid = ggplot2::element_blank()) +

    ggplot2::guides(
      shape = "none",
      color = ggplot2::guide_legend(
        title = title_legend,
        ncol = legend_ncol,
        override.aes = list(size = size_point_legend))) +
    ggplot2::theme(legend.title = ggplot2::element_text(
      face = "bold", size = size_title_legned, color = "black")) +
    ggplot2::theme(legend.text = ggplot2::element_text(
      face = "bold", size = size_legned, color = "black")) +


    ggplot2::theme(legend.text = ggplot2::element_text(
      margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"))) +
    ggplot2::theme(legend.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(legend.title = ggplot2::element_text(
      face = "bold", size = 24, color="black")) +
    ggplot2::theme(legend.text = ggplot2::element_text(
      face = "bold", size = 24, color="black")) +
    ggplot2::theme(legend.key.height = ggplot2::unit(
      spacing_legend_point, "cm")) +
    ggplot2::theme(legend.title = element_text(
      margin = ggplot2::margin(b = spacing_legend_title, unit = 'cm'))) +
    ggplot2::theme(plot.title = ggplot2::element_text(
      margin = ggplot2::margin(t = 0, r = 0, b = 15, l = 0, unit = "pt"))) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(
      margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0, unit = "pt"))) +
    ggplot2::theme(axis.title.y = ggplot2::element_text(
      margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0, unit = "pt"))) +
    ggplot2::theme(legend.margin = ggplot2::margin(
      t = 0, r = 0, b = 0, l = 20, unit = "pt")) +
    ggplot2::theme(plot.margin = ggplot2::margin(
      t = 20, r = 30, b = 20, l = 30, unit = "pt"))


  ##
  p1 <- p1 +
    ggplot2::geom_segment(
      data = data[["env"]],
      ggplot2::aes(x = 0,           # X
                   y = 0,           # Y
                   xend = data[["env"]][,1],
                   yend = data[["env"]][,2]),
      color = "#585858",
      linewidth = 0.8,
      alpha = 0.6,
      arrow = ggplot2::arrow(angle = 35, length = ggplot2::unit(0.3, "cm")))


  ##
  p1 <- p1 +
    ggrepel::geom_text_repel(
      data = data[["env"]],
      ggplot2::aes(
        x = data[["env"]][, 1],
        y = data[["env"]][, 2],
        label = rownames(data[["env"]])
      ),
      size = 5,
      color = "#000000",
      box.padding = ggplot2::unit(0.45, "lines"),
      alpha = 0.75
    )


  ##
  if(!is.null(color_scheme)) {
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
