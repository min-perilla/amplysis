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
#' custom_order = c("A", "B", "R", "D", "S"))}
#'
#' @importFrom dplyr arrange
alpha_plot <- function(data, color_scheme = NULL, custom_order = NULL,
                       size_point = 5, size_differ = 18, errorbar_width = 0.15,
                       errorbar_linewidth = 0.8, title_x = NULL, title_y = NULL,
                       size_title = 48, size_title_x = 32, size_title_y = 32,
                       size_x = 32, size_y = 32,
                       filename = "alpha", file_width = 12, file_height = 9)
{
  ## Custom order function
  customOrder <- function(data,  # Input data
                          col,   # Column name to be converted to a factor
                          custom_order  # Custom order
  ) {
    # Convert the specified column to a factor and order by custom levels
    data[, col] <- factor(data[, col], levels = custom_order)
    data <- dplyr::arrange(data, col)
    return(data)
  }

  # Apply custom order to the x-axis
  if (!is.null(custom_order)) {
    # Set group information
    group = "group"

    # Apply custom order
    for (i in 1:length(data)) {
      data[[i]] = customOrder(data[[i]], group, custom_order)
    }
    cat("Custom legend order: ", custom_order, "\n", sep = " ")
  }

  # List to store results
  result <- list()

  title = NULL
  # Set title
  if (is.null(title)) {
    title = tools::toTitleCase(names(data))
  }
  # Convert "pd" to uppercase "PD"
  title <- gsub("pd", "PD", title, ignore.case = TRUE)

  # Set result names
  result_name = names(data)

  # Generate plots
  for (i in 1:length(data)) {
    # Parameter settings
    title1 = title[i]    # Plot title
    x = data[[i]][, "group"]  # x-axis
    y = data[[i]][, 1]        # y-axis
    filename1 = paste0(filename, "_", names(data)[i])  # File name

    # Create plot
    p1 <- box_plot(data = data[[i]], index_type = 1, x = x, y = y,
                   color_scheme = color_scheme,

                   size_point = size_point, size_differ = size_differ, errorbar_width = errorbar_width,
                   errorbar_linewidth = errorbar_linewidth,

                   title = title1, title_x = title_x, title_y = title_y,

                   size_title = size_title, size_title_x = size_title_x, size_title_y = size_title_y,
                   size_x = size_x, size_y = size_y,

                   filename = filename1, file_width = file_width, file_height = file_height)

    # Save plot to result list
    result[[result_name[i]]] <- p1

    print(p1) # Preview result
  }

  # Print success message
  cat("\033[32mtaxa_bar: success!\033[0m\n")
  cat("\033[0;32m", "The file has been saved to \n",
      getwd(), "/result\033[0m\n", sep = "")

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
#' size_title = 48, size_title_x = 32, size_title_y = 32, size_x = 32,
#' size_y = 32, filename = "alpha", file_width = 12, file_height = 9)
#' }
#'
#' @importFrom ggplot2 aes element_text geom_boxplot geom_jitter geom_text
#' ggsave labs margin theme theme_bw scale_color_manual scale_fill_manual
#' stat_boxplot
#'
box_plot <- function(data, index_type, x, y, color_scheme,
         size_point, size_differ, errorbar_width, errorbar_linewidth, title,
         title_x, title_y, size_title, size_title_x, size_title_y, size_x,
         size_y, filename, file_width, file_height)
{
  # Convert to symbols
  group <- rlang::sym("group")
  differ_y = rlang::sym("differ_y")
  differ = rlang::sym("differ")

  p1 <-
    ggplot2::ggplot(
      data,
      ggplot2::aes(x = x, y = y)) +

    # Theme settings
    ggplot2::theme_bw() +

    # Plot settings
    # Add error bars
    ggplot2::stat_boxplot(
      geom = "errorbar",
      width = errorbar_width, linewidth = errorbar_linewidth,  # Width and size
      ggplot2::aes(color = factor(x))) + # Color

    # Box plot
    ggplot2::geom_boxplot(
      ggplot2::aes(color = x),    # Set border color and fill color based on group
      outlier.colour = "red",     # Specify outlier color for visibility
      outlier.size = size_point,  # Outlier point size
      outlier.alpha = 0,
      size = 1.1,                 # Line thickness of the box plot
      fill = "transparent",       # Transparent fill color
      alpha = 1                   # Border transparency
    ) +

    # Jitter points
    ggplot2::geom_jitter(width = 0.2,             # Jitter range
                         size = size_point,       # Jitter point size
                         ggplot2::aes(color = factor(x)),  # Jitter point color
                         alpha = 0.5) +           # Jitter point transparency

    # Significance labels
    ggplot2::geom_text(data = data,   # data = data[[index_type]],
                       ggplot2::aes(x = !!group, y = !!differ_y, color = !!group, label = !!differ),
                       size = size_differ) +

    # Set main title
    ggplot2::labs(title = title) +
    # Set X and Y axis titles
    ggplot2::labs(x = title_x, y = title_y) +

    # Unified theme() to avoid overriding
    ggplot2::theme(
      # Main title font size
      plot.title = ggplot2::element_text(
        face = "bold",
        size = size_title, hjust = 0.5,
        margin = ggplot2::margin(t = 0, r = 0, b = 30, l = 0, unit = "pt")
      ),

      # X-axis title
      axis.title.x = ggplot2::element_text(
        size = size_title_x,
        margin = ggplot2::margin(t = 20, r = 0, b = 0, l = 0, unit = "pt")
      ),

      # Y-axis title
      axis.title.y = ggplot2::element_text(
        size = size_title_y, angle = 90,
        margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0, unit = "pt")
      ),

      # Axis tick label font size
      # Modify X-axis tick labels
      axis.text.x = ggplot2::element_text(
        face = "bold",
        size = size_x,
        margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")
      ),
      # Modify Y-axis tick labels
      axis.text.y = ggplot2::element_text(
        size = size_y,
        margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0, unit = "pt")
      ),

      # Remove legend
      legend.position = "none",

      # Plot margins
      plot.margin = ggplot2::margin(
        t = 50, r = 120, b = 50, l = 50, unit = "pt"
      )
    )

  ## Apply color scheme
  if (!is.null(color_scheme)) {
    p1 <- p1 +
      ggplot2::scale_color_manual(values = color_scheme) +
      ggplot2::scale_fill_manual(values = color_scheme)
  }

  ## Save the file
  # Define folder path
  folder_path <- "result"  # Corrected folder name

  # Check if the folder exists; if not, create it
  if (!dir.exists(folder_path)) {
    dir.create(folder_path)
  }

  # Save as PNG
  png_file_path <- file.path(getwd(), folder_path, paste0(filename, ".png"))
  ggplot2::ggsave(filename = png_file_path, plot = p1,
                  width = file_width, height = file_height, dpi = 300)

  # Check if the PDF folder exists; if not, create it
  pdf_folder_path <- file.path(folder_path, "PDF")
  if (!dir.exists(pdf_folder_path)) {
    dir.create(pdf_folder_path)
  }

  # Save as PDF
  pdf_file_path <- file.path(getwd(), pdf_folder_path, paste0(filename, ".pdf"))
  ggplot2::ggsave(filename = pdf_file_path, plot = p1,
                  width = file_width, height = file_height, dpi = 300)

  return(p1)
}

