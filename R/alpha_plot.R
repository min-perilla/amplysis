#' @title Visualization of Alpha Diversity Analysis
#'
#' @description
#' After processing with the `alpha()` function, you can use this function
#' for visualization.
#'
#' @param data Plotting data.
#' @param color_scheme (character) Color scheme.
#' @param custom_order (character) Custom legend order.
#'
#' @param size_point (numeric) The size of points.
#' @param size_differ (numeric) The size of significance marker letters.
#' @param errorbar_width (numeric) The width of the horizontal lines on the error bars.
#' @param errorbar_linewidth (numeric) The width of the vertical lines on the error bars.
#'
#' @param title (character) The title of the chart.
#' @param title_x (character) The title of the X-axis.
#' @param title_y (character) The title of the Y-axis.
#'
#' @param size_title (numeric) Font size of the main title.
#' @param size_title_x (numeric) Font size of the horizontal axis title.
#' @param size_title_y (numeric) Font size of the vertical axis title.
#'
#' @param size_x (numeric) Font size of the horizontal axis tick labels.
#' @param size_y (numeric) Font size of the vertical axis tick labels.
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
#' alpha_plot(data,
#' color_scheme = c('#aaf200','#0082ff',"#d23aa4","#c777ff", "#79ff79"),
#' custom_order = c("A", "B", "R", "D", "S"))
#' }
#'
#' @importFrom dplyr arrange
alpha_plot <- function(data, color_scheme = NULL, custom_order = NULL,
                       size_point = 5, size_differ = 14, errorbar_width = 0.15,
                       errorbar_linewidth = 0.8,
                       title = NULL, title_x = NULL, title_y = NULL,
                       size_title = 40, size_title_x = 28, size_title_y = 28,
                       size_x = 28, size_y = 28,
                       filename = "alpha", file_width = 12, file_height = 9)
{
  customOrder <- function(data, col, custom_order){
    data[, col] <- factor(data[, col], levels = custom_order)
    data <- dplyr::arrange(data, col)
    return(data)
  }

  if(!is.null(custom_order)) {
    group = "group"
    for (i in 1:length(data)) {
      data[[i]] = customOrder(data[[i]], group, custom_order)
    }
    cat("Custom legend order: ", custom_order, "\n", sep = " ")
  }


  result <- list()

  if(is.null(title)) {
    title = names(data)
  }

  result_name = names(data)

  #
  for(i in 1:length(data)) {
    title1 = title[i]    # 绘图标题
    x = data[[i]][, "group"]
    y = data[[i]][, 1]
    filename1 = NULL
    filename1 = paste0(filename, "_", names(data)[i])

    #
    p1 <- box_plot(
      data = data[[i]], index_type = 1, x = x, y = y,
      color_scheme = color_scheme,

      size_point = size_point, size_differ = size_differ, errorbar_width = errorbar_width,
      errorbar_linewidth = errorbar_linewidth,

      title = title1, title_x = title_x, title_y = title_y,

      size_title = size_title, size_title_x = size_title_x, size_title_y = size_title_y,
      size_x = size_x, size_y = size_y,

      filename = filename1, file_width = file_width, file_height = file_height)

    result[[result_name[i]]] <- p1

    print(p1)
  }


  ##
  cat("\033[32mtaxa_bar: success!\033[0m\n")
  cat("\033[0;32m", "The file has been saved to \n",
      getwd(), "/result\033[0m\n", sep = "")

  ##
  return(result)
}


################################################################################
#' @title Box Plot
#' @description
#' This function is specifically designed to provide box plot drawing
#' capabilities for the alpha_plot() function.
#'
#' @param data Plotting data.
#' @param index_type (character) Type of Alpha diversity indicex
#' @param x The data on the X-axis of the box plot
#' @param y The data on the X-axis of the box plot
#' @param color_scheme (character) Color scheme.
#'
#' @param size_point (numeric) The size of points.
#' @param size_differ (numeric) The size of significance marker letters.
#' @param errorbar_width (numeric) The width of the horizontal lines on the error bars.
#' @param errorbar_linewidth (numeric) The width of the vertical lines on the error bars.
#'
#' @param title (character) Main title.
#' @param title_x (character) The title of the X-axis.
#' @param title_y (character) The title of the Y-axis.
#'
#' @param size_title (numeric) Font size of the main title.
#' @param size_title_x (numeric) Font size of the horizontal axis title.
#' @param size_title_y (numeric) Font size of the vertical axis title.
#'
#' @param size_x (numeric) Font size of the horizontal axis tick labels.
#' @param size_y (numeric) Font size of the vertical axis tick labels.
#'
#' @param filename (character) File name for saving.
#' @param file_width (numeric) Width of the image.
#' @param file_height (numeric) Height of the image.
#'
#' @return A graph of ggplot2 class
#'
#' @examples
#' \dontrun{
#' box_plot(data = data, index_type = "Shannon", x = data[[index_type]][, "group"],
#' y = data[[index_type]][, 1], color_scheme = NULL,
#' size_point = 5, size_differ = 14, errorbar_width = 0.15,
#' errorbar_linewidth = 0.8, title_x = NULL, title_y = NULL,
#' size_title = 40, size_title_x = 28, size_title_y = 28, size_x = 28,
#' size_y = 28, filename = "alpha", file_width = 12, file_height = 9)
#' }
#'
#' @importFrom ggplot2 aes element_text geom_boxplot geom_jitter geom_text
#' ggsave labs margin theme theme_bw scale_color_manual scale_fill_manual
#' stat_boxplot
box_plot <- function(data, index_type, x, y, color_scheme,
         size_point, size_differ, errorbar_width, errorbar_linewidth, title,
         title_x, title_y, size_title, size_title_x, size_title_y, size_x,
         size_y, filename, file_width, file_height) {

  group <- rlang::sym("group")
  differ_y = rlang::sym("differ_y")
  differ = rlang::sym("differ")

  p1 <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = x, y = y)) +
    ggplot2::theme_bw() +

    ggplot2::stat_boxplot(
      geom = "errorbar",
      width = errorbar_width, linewidth = errorbar_linewidth,
      ggplot2::aes(color = factor(x))) +

    ggplot2::geom_boxplot(
      ggplot2::aes(color = x),
      outlier.colour = "red",
      outlier.size = size_point,
      outlier.alpha = 0,
      size = 1.1,
      fill = "transparent",
      alpha = 1
    ) +

    ggplot2::geom_jitter(width = 0.2,
                         size = size_point,
                         ggplot2::aes(color = factor(x)),
                         alpha = 0.5) +

    # 显著性标记
    ggplot2::geom_text(data = data,
                       ggplot2::aes(x = !!group, y = !!differ_y, color = !!group,
                                    label = !!differ),
                       size = size_differ) +

    ggplot2::labs(title = title) +
    ggplot2::theme(plot.title = ggplot2::element_text(
      face = "bold", size = size_title, hjust = 0.5)) +
    ggplot2::labs(x = title_x, y = title_y) +
    ggplot2::theme(
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold")
    ) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = size_title_x),
                   axis.title.y = ggplot2::element_text(size = size_title_y, angle = 90),
                   axis.text.x = ggplot2::element_text(size = size_x),
                   axis.text.y = ggplot2::element_text(size = size_y)) +
    theme(legend.position = "none") +
    ggplot2::theme(plot.title = ggplot2::element_text(
      margin = ggplot2::margin(t = 0, r = 0, b = 15, l = 0, unit = "pt"))) +
    ggplot2::theme(axis.title.x = ggplot2::element_text(
      margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0, unit = "pt"))) +
    ggplot2::theme(axis.title.y = ggplot2::element_text(
      margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0, unit = "pt"))) +
    ggplot2::theme(legend.margin = ggplot2::margin(
      t = 0, r = 0, b = 0, l = 20, unit = "pt")) +
    ggplot2::theme(plot.margin = ggplot2::margin(
      t = 40, r = 50, b = 40, l = 50, unit = "pt"))

  if(!is.null(color_scheme)) {
    p1 <- p1 +
      ggplot2::scale_color_manual(values = color_scheme) +
      ggplot2::scale_fill_manual(values = color_scheme)
  }


  ##
  folder_path <- "relust"
  if (!dir.exists(folder_path)) {
    dir.create(folder_path)
  }
  png_folder_path <- file.path(folder_path, "PNG")
  if (!dir.exists(png_folder_path)) {
    dir.create(png_folder_path)
  }

  #
  png_file_path <- file.path(getwd(), png_folder_path, paste0(filename, ".png"))
  ggplot2::ggsave(filename = png_file_path, plot = p1,
                  width = file_width, height = file_height)


  ##
  pdf_folder_path <- file.path(folder_path, "PDF")
  if (!dir.exists(pdf_folder_path)) {
    dir.create(pdf_folder_path)
  }

  #
  pdf_file_path <- file.path(getwd(), pdf_folder_path, paste0(filename, ".pdf"))
  ggplot2::ggsave(filename = pdf_file_path, plot = p1,
                  width = file_width, height = file_height)


  ##
  return(p1)
}
