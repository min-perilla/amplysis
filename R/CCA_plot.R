#' @title Visualization of CCA
#'
#' @description
#' After processing with the `CCA()` function, you can use this function
#' for visualization.
#'
#' @param data Plotting data.
#' @param color_scheme (character) Color scheme.
#' @param custom_order Custom legend order.
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
#' cca1 <- CCA(otu = otu, env = env, metadata = metadata, id_col = 1,
#' group = "group", replicate_method = "none")
#'
#' CCA_plot(data = cca1, custom_order = NULL,
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
CCA_plot = function(data, color_scheme = NULL, custom_order = NULL, seed = 123,

                    size_point = 4.5, size_point_legend = 8, spacing_legend_point = 1.2,
                    spacing_legend_title = 0.5, legend_ncol = 1, label_is = T,
                    size_label = 5, label_font_color = NULL, ellipse_type = "t",

                    title = "RDA", title_sub = NULL, title_legend = "Group",

                    size_title = 28, size_title_sub = 16, size_title_x = 20,
                    size_title_y = 20, size_title_legend = 24,

                    size_x = 18, size_y = 18, size_legend = 16,

                    filename = "CCA", file_width = 12, file_height = 9)
{
  # Set seed
  set.seed(seed = seed)

  group = "group"    # Grouping information

  # Determine whether to display data labels
  if(isTRUE(label_is)) {
    labelNum = 200
  } else {
    labelNum = 0
  }

  # Custom legend order
  if(!is.null(custom_order)){
    data[["data"]][[group]] <- factor(data[["data"]][[group]], levels = custom_order)
    data[["data"]] <- dplyr::arrange(data[["data"]], group)

    cat("Custom legend order: ", custom_order, "\n", sep = "")
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

  # X and Y axes
  x = "CCA1"
  y = "CCA2"

  # Get the percentage contribution of CCA1 and CCA2
  x_contrib <- data[[x]]
  y_contrib <- data[[y]]

  # Convert to symbol objects
  group_sym <- rlang::sym(group)
  x_sym <- rlang::sym(x)
  y_sym <- rlang::sym(y)

  ## Plotting
  p1 <- ggplot2::ggplot(
    data = data[["data"]],     # Plot data
    ggplot2::aes(
      x = !!x_sym,   # X-axis
      y = !!y_sym,   # Y-axis
      color = !!group_sym)) +

    ##
    # Theme color settings
    ggplot2::theme_bw() +  # White background with black lines

    # Add dashed lines through the origin
    ggplot2::geom_vline(xintercept = 0, lty = "dashed") +
    ggplot2::geom_hline(yintercept = 0, lty = "dashed") +

    # Draw scatter plot with size settings
    ggplot2::geom_point(size = size_point, shape = 16, alpha = 0.7) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())  # Remove grid lines


  ## Add adaptive labels
  if (is.null(label_font_color)) {  # Default font color for labels
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

    # Change x and y axis titles to contribution values
    ggplot2::labs(x = x_contrib,
                  y = y_contrib) +

    # Add confidence ellipse
    ggplot2::stat_ellipse(
      data = data[["data"]],
      geom = "polygon",
      level = 0.95,      # Confidence level
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
    # The hjust parameter (range 0 to 1) controls the horizontal alignment of the text
    ggplot2::theme(plot.title = ggplot2::element_text(
      face = "bold", size = size_title, hjust = 0.5)) +

    # Set subtitle
    ggplot2::labs(subtitle = title_sub) +
    # Subtitle font size
    ggplot2::theme(plot.subtitle = ggplot2::element_text(
      face = "bold", size = size_title_sub, hjust = 0)) +  # The hjust parameter (range 0 to 1) controls the horizontal alignment of the text

    # Axis label font size settings
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = size_title_x),     # Modify X-axis title text
                   axis.title.y = ggplot2::element_text(size = size_title_y, angle = 90),  # Modify Y-axis title text
                   axis.text.x = ggplot2::element_text(size = size_x),                     # Modify X-axis tick labels
                   axis.text.y = ggplot2::element_text(size = size_y)                      # Modify Y-axis tick labels
    ) +

    # Set legend
    ggplot2::guides(
      shape = "none",
      color = ggplot2::guide_legend(
        title = title_legend,              # Set legend title
        ncol = legend_ncol,                # Number of legend columns
        override.aes = list(size = size_point_legend))) +
    # Legend title font size
    ggplot2::theme(legend.title = ggplot2::element_text(
      face = "bold", size = size_title_legend, color = "black")) + # Title size
    ggplot2::theme(legend.text = ggplot2::element_text(
      face = "bold", size = size_legend, color = "black")) +       # Font and size


    # Set legend text margins
    ggplot2::theme(legend.text = ggplot2::element_text(
      margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"))) +
    ggplot2::theme(legend.title = ggplot2::element_text(hjust = 0.5)) +  # Center legend title


    # Internal spacing of the legend
    ggplot2::theme(legend.key.height = ggplot2::unit(
      spacing_legend_point, "cm")) +
    # Spacing between legend title and content
    ggplot2::theme(legend.title = element_text(
      margin = ggplot2::margin(b = spacing_legend_title, unit = 'cm'))) +


    # Margin settings: t = top, b = bottom, r = right, l = left
    # Main title margin
    ggplot2::theme(plot.title = ggplot2::element_text(
      margin = ggplot2::margin(t = 0, r = 0, b = 15, l = 0, unit = "pt"))) +
    # X-axis margin
    ggplot2::theme(axis.title.x = ggplot2::element_text(
      margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0, unit = "pt"))) +
    # Y-axis margin
    ggplot2::theme(axis.title.y = ggplot2::element_text(
      margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0, unit = "pt"))) +
    # Legend margin
    ggplot2::theme(legend.margin = ggplot2::margin(
      t = 0, r = 0, b = 0, l = 20, unit = "pt")) +

    # Adjust the margins of this plot
    ggplot2::theme(plot.margin = ggplot2::margin(
      t = 20, r = 30, b = 20, l = 30, unit = "pt"))


  ##
  # Add environmental factor data
  p1 <- p1 +
    # Add environmental factor arrows
    ggplot2::geom_segment(
      data = data[["env"]],           # Data for plotting
      ggplot2::aes(x = 0,           # X-axis
                   y = 0,           # Y-axis
                   xend = data[["env"]][,1],   # X-axis endpoint
                   yend = data[["env"]][,2]),  # Y-axis endpoint
      color = "#585858",
      linewidth = 0.8,
      alpha = 0.6,
      arrow = ggplot2::arrow(angle = 35, length = ggplot2::unit(0.3, "cm")))


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

    # Change x and y axis titles to contribution percentages
    ggplot2::labs(x = x_contrib,
                  y = y_contrib) +

    # Add confidence ellipse
    ggplot2::stat_ellipse(
      data = data[["data"]],
      geom = "polygon",
      level = 0.95,      # level: confidence level
      linetype = 2,      # Line type style
      linewidth = 0.4,
      ggplot2::aes(fill = group),
      alpha = 0.15,
      show.legend = F,   # Hide legend
      # type = "t"       # Default, ellipse based on t-distribution (suitable for small samples, using robust estimation MASS::cov.trob())
      # type = "norm"    # Ellipse based on normal distribution (computed using covariance matrix, without robust estimation)
      # type = "euclid"  # Draws a fixed-radius circle based on Euclidean distance (scale-dependent)
      type = ellipse_type  # Method for calculating the confidence ellipse
    ) +


    # Set main title
    ggplot2::labs(title = title) +
    # The `hjust` parameter (range: 0 to 1) controls the horizontal alignment of the text
    ggplot2::theme(plot.title = ggplot2::element_text(
      face = "bold", size = size_title, hjust = 0.5)) +

    # Set subtitle
    ggplot2::labs(subtitle = title_sub) +
    # Subtitle font size
    ggplot2::theme(plot.subtitle = ggplot2::element_text(
      face = "bold", size = size_title_sub, hjust = 0)) +  # `hjust` parameter (range: 0 to 1) controls horizontal alignment

    # Set font size for axis labels and tick marks
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = size_title_x),     # Modify X-axis title text
                   axis.title.y = ggplot2::element_text(size = size_title_y, angle = 90),  # Modify Y-axis title text
                   axis.text.x = ggplot2::element_text(size = size_x),                     # Modify X-axis tick labels
                   axis.text.y = ggplot2::element_text(size = size_y)                      # Modify Y-axis tick labels
    ) +

    # Set legend
    ggplot2::guides(
      shape = "none",
      color = ggplot2::guide_legend(
        title = title_legend,              # Set legend title
        ncol = legend_ncol,                # Number of legend columns
        override.aes = list(size = size_point_legend))) +
    # Legend font size
    ggplot2::theme(legend.title = ggplot2::element_text(
      face = "bold", size = size_title_legend, color = "black")) + # Title size
    ggplot2::theme(legend.text = ggplot2::element_text(
      face = "bold", size = size_legend, color = "black")) +       # Font style and size


    # Set legend text margin
    ggplot2::theme(legend.text = ggplot2::element_text(
      margin = ggplot2::margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"))) +
    ggplot2::theme(legend.title = ggplot2::element_text(hjust = 0.5)) +                                 # Center-align legend title


    # Legend internal spacing
    ggplot2::theme(legend.key.height = ggplot2::unit(
      spacing_legend_point, "cm")) +
    # Spacing between legend title and content
    ggplot2::theme(legend.title = element_text(
      margin = ggplot2::margin(b = spacing_legend_title, unit = 'cm'))) +


    # Margin settings: `t` for top margin, `b` for bottom margin, `r` for right margin, `l` for left margin
    # Main title margin
    ggplot2::theme(plot.title = ggplot2::element_text(
      margin = ggplot2::margin(t = 0, r = 0, b = 15, l = 0, unit = "pt"))) +
    # X-axis margin
    ggplot2::theme(axis.title.x = ggplot2::element_text(
      margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0, unit = "pt"))) +
    # Y-axis margin
    ggplot2::theme(axis.title.y = ggplot2::element_text(
      margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0, unit = "pt"))) +
    # Legend margin
    ggplot2::theme(legend.margin = ggplot2::margin(
      t = 0, r = 0, b = 0, l = 20, unit = "pt")) +

    # Adjust the overall plot margin
    ggplot2::theme(plot.margin = ggplot2::margin(
      t = 20, r = 30, b = 20, l = 30, unit = "pt"))


  ##
  # Add environmental factor data
  p1 <- p1 +
    # Add environmental factor arrows
    ggplot2::geom_segment(
      data = data[["env"]],           # Plot data
      ggplot2::aes(x = 0,           # X-axis
                   y = 0,           # Y-axis
                   xend = data[["env"]][,1],   # End position on X-axis
                   yend = data[["env"]][,2]),  # End position on Y-axis
      color = "#585858",
      linewidth = 0.8,
      alpha = 0.6,
      arrow = ggplot2::arrow(angle = 35, length = ggplot2::unit(0.3, "cm")))


  ##
  # Add labels to the arrows
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
      box.padding = ggplot2::unit(0.45, "lines"),  # Control the space around labels
      alpha = 0.75
    )


  ## Color scheme
  if(!is.null(color_scheme)) {

    # Custom color settings for points
    color_scheme_point <- color_scheme          # Point fill color
    color_scheme_ellipse <- color_scheme_point  # Confidence ellipse fill color

    p1 <- p1 +
      # Point color
      ggplot2::scale_color_manual(values = color_scheme_point) +

      # Outline color
      ggplot2::scale_fill_manual(values = color_scheme_ellipse)
  }



  ##
  # Save file
  #
  ggplot2::ggsave(filename = paste0(filename, ".png"), plot = p1, width = file_width, height = file_height)  # Save as PNG file
  ggplot2::ggsave(filename = paste0(filename, ".pdf"), plot = p1, width = file_width, height = file_height)  # Save as PDF file

  ##
  cat("\033[32mtaxa_bar: success!\033[0m\n")
  cat("\033[0;32m", "The file has been saved to \n",
      getwd(), "\033[0m\n", sep = "")

  return(p1)
}
