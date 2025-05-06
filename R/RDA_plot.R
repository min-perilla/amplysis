#' @title Visualization of RDA
#'
#' @description
#' After processing with the `RDA()` function, you can use this function
#' for visualization.
#'
#' @param data Plotting data.
#' @param color_scheme (character) Color scheme.
#' @param custom_order (character) Custom legend order.
#' @param seed (numeric) Seed
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
#' rda1 <- RDA(otu = otu, env = env, metadata = metadata, id_col = 1,
#' group = "group", replicate_method = "none")
#'
#' RDA_plot(data = rda1, custom_order = NULL,
#'          color_scheme = c("#00b0f6", "#FFC24B", "#f8766d"), seed = 123,
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
RDA_plot = function(data, color_scheme = NULL, custom_order = NULL, seed = 123,

                    size_point = 4.5, size_point_legend = 8, spacing_legend_point = 1.2,
                    spacing_legend_title = 0.5, legend_ncol = 1, label_is = T,
                    size_label = 5, label_font_color = NULL, ellipse_type = "t",

                    title = "RDA", title_sub = NULL, title_legend = "Group",

                    size_title = 28, size_title_sub = 16, size_title_x = 20,
                    size_title_y = 20, size_title_legend = 24,

                    size_x = 18, size_y = 18, size_legend = 16,

                    filename = "RDA", file_width = 12, file_height = 9)
{
  # Seed setting
  set.seed(seed = seed)

  group = "group"    # Group information

  # Check if data labels should be displayed
  if(isTRUE(label_is)) {
    labelNum = 200
  } else {
    labelNum = 0
  }

  # Custom legend order
  if(!is.null(custom_order)){
    data[["data"]][[group]] <- factor(data[["data"]][[group]], levels = custom_order)
    data[["data"]] <- dplyr::arrange(data[["data"]], group)

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


  # x and y axes
  x = "RDA1"
  y = "RDA2"


  # Get the contribution percentages for RDA1 and RDA2
  x_contrib <- data[[x]]
  y_contrib <- data[[y]]


  # Convert to symbol objects
  group_sym <- rlang::sym(group)
  x_sym <- rlang::sym(x)
  y_sym <- rlang::sym(y)


  ##
  # Plotting
  p1 <- ggplot2::ggplot(
    data = data[["data"]],     # Plot data
    ggplot2::aes(
      x = !!x_sym,   # X-axis
      y = !!y_sym,   # Y-axis
      color = !!group_sym)) +

    ##
    # Theme color settings
    ggplot2::theme_bw() +  # White background and black lines

    # Add dashed lines through the origin
    ggplot2::geom_vline(xintercept = 0, lty = "dashed") +
    ggplot2::geom_hline(yintercept = 0, lty = "dashed") +

    # Plot points with size setting
    ggplot2::geom_point(size = size_point, shape = 16, alpha = 0.7) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())  # Remove grid lines


  ## Add adaptive labels
  if (is.null(label_font_color)) {  # Default label font color
    p1 <- p1 +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = sample), size = size_label,
        box.padding = ggplot2::unit(0.6, "lines"),
        point.padding = ggplot2::unit(0.5, "lines"),
        max.overlaps = labelNum, alpha = 0.8, show.legend = F, seed = seed)
  } else {
    # Custom label color
    p1 <- p1 +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = sample),

        color = label_font_color,  # Label font color

        size = size_label,
        box.padding = ggplot2::unit(0.6, "lines"),
        point.padding = ggplot2::unit(0.5, "lines"),
        max.overlaps = labelNum, alpha = 0.8, show.legend = F, seed = seed)

  }


  p1 <- p1 +

    # Change X and Y axis titles to contribution percentages
    ggplot2::labs(x = paste0("RDA1 (", x_contrib, "%)"),
                  y = paste0("RDA2 (", y_contrib, "%)")) +

    # Add confidence ellipses
    ggplot2::stat_ellipse(
      data = data[["data"]],
      geom = "polygon",
      level = 0.95,      # level: confidence level
      linetype = 2,      # Line style
      linewidth = 0.4,
      ggplot2::aes(fill = group),
      alpha = 0.15,
      show.legend = F,   # Do not display legend
      # type = "t"       # Default: Ellipse based on t-distribution (for small samples, using robust estimation MASS::cov.trob())
      # type = "norm"    # Ellipse based on normal distribution (calculated using covariance matrix, without robust estimation)
      # type = "euclid"  # Circle with a fixed radius drawn using Euclidean distance (depends on data scale)
      type = ellipse_type  # Calculation method for the confidence ellipse
    ) +

    # Set main title
    ggplot2::labs(title = title) +
    # Main title font size
    ggplot2::theme(plot.title = ggplot2::element_text(
      face = "bold", size = size_title, hjust = 0.5)) +   # hjust parameter (range 0 to 1) controls horizontal alignment of the text

    # Set subtitle
    ggplot2::labs(subtitle = title_sub) +
    # Subtitle font size
    ggplot2::theme(plot.subtitle = ggplot2::element_text(
      face = "bold", size = size_title_sub, hjust = 0)) +  # hjust parameter (range 0 to 1) controls horizontal alignment of the text


    # Tick label font size settings
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = size_title_x),     # Modify X-axis title text
                   axis.title.y = ggplot2::element_text(size = size_title_y, angle = 90),  # Modify Y-axis title text
                   axis.text.x = ggplot2::element_text(size = size_x),                     # Modify X-axis tick label text
                   axis.text.y = ggplot2::element_text(size = size_y)) +

    # Set legend
    ggplot2::guides(
      shape = "none",
      color = ggplot2::guide_legend(
        title = title_legend,              # Set legend title
        ncol = legend_ncol,                # Number of columns in legend
        override.aes = list(size = size_point_legend))) +

    # Legend font size
    ggplot2::theme(legend.title = ggplot2::element_text(
      face = "bold", size = size_title_legend, color = "black")) + # Legend title size
    ggplot2::theme(legend.text = ggplot2::element_text(
      face = "bold", size = size_legend, color = "black")) +       # Legend text size


    # Set legend text margin
    ggplot2::theme(legend.text = ggplot2::element_text(
      margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"))) +
    ggplot2::theme(legend.title = ggplot2::element_text(hjust = 0.5)) +                                 # Center legend title


    # Legend internal spacing
    ggplot2::theme(legend.key.height = ggplot2::unit(
      spacing_legend_point, "cm")) +
    # Spacing between legend title and body
    ggplot2::theme(legend.title = element_text(
      margin = ggplot2::margin(b = spacing_legend_title, unit = 'cm'))) +


    # Margin settings: t for top, b for bottom, r for right, l for left
    # Main title margin
    ggplot2::theme(plot.title = ggplot2::element_text(
      margin = ggplot2::margin(t = 0, r = 0, b = 15, l = 0, unit = "pt"))) +
    # X-axis margin
    ggplot2::theme(axis.title.x = ggplot2::element_text(
      margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0, unit = "pt"))) +
    # Y-axis margin
    ggplot2::theme(axis.title.y = ggplot2::element_text(
      margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0, unit = "pt")))


  ## Saving as an image
  ggplot2::ggsave(
    filename = paste0(filename, ".png"),
    plot = p1,
    dpi = 300, width = file_width, height = file_height)

  return(p1)  # Return the plot object
}
