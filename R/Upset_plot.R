#' @title Visualization of Upset
#'
#' @description
#' After processing with the `Upset()` function, you can use this function
#' for visualization.
#'
#'
#' @param data Plotting data.
#' @param color_matrix (character) Color scheme of the bar chart in the matrix
#' plot.
#' @param n (integer) Number of bars displayed in the bar chart.
#' @param custom_order (character) Order of the x-axis in the bar chart.
#' @param order_by (character) Sorting order of the matrix plot, options are
#' "freq" or "degree".
#'
#' @param mb.ratio (character) The proportion of the heights between the bar
#' plot and the matrix plot.
#' @param size_point (numeric) The size of the points in the matrix plot.
#'
#' @param color_point (character) The color of the points in the matrix plot.
#' @param color_matrix_shade (character) The color of the shaded areas in the
#' matrix plot.
#' @param color_bar (character) The color of the bars in the y-axis bar plot.
#'
#' @param title_matrix_x (character) The labels on the x-axis of the matrix plot.
#' @param title_bar_y (character) The title of the Y-axis in the bar plot.
#'
#' @param size_title_matrix_x (numeric) The font size of the x-axis title in the
#' matrix plot.
#' @param size_title_matrix_y (numeric) The font size of the y-axis title in the
#' matrix plot.
#' @param size_title_bar_y (numeric) The font size of the vertical axis title in
#' the bar plot.
#'
#' @param size_matrix_x (numeric) The size of the tick labels in the matrix plot.
#' @param size_bar_y (numeric) The font size of the characters on the vertical
#' axis scale of the bar plot.
#' @param size_bar_label (numeric) The size of the numbers on the bars in the
#' bar plot.
#'
#' @param queries (a list) Highlight specific elements in the matrix plot or the
#' bar plot.
#'
#' @param filename (character) File name for saving.
#' @param file_width (numeric) Width of the image.
#' @param file_height (numeric) Height of the image.
#'
#' @return upset plot
#' @export
#'
#' @examples
#' \dontrun{
#' upset1 <- Upset(otu, metadata, 1, "group", "mean")
#' Upset_plot(upset1,
#'            queries = list(
#'              list(query = intersects,
#'                   params = list("AA", "AB", "AD"), color="#f06676", active = T),
#'              list(query = intersects,
#'                   params = list("AA", "AB"), color="#f06676", active = T)))
#' }
#' \dontrun{
#' upset1 <- Upset(otu, metadata, 1, "group", "mean")
#' Upset_plot(
#'   upset1, color_matrix = NULL, n = 40, custom_order = NULL, order_by = "freq",
#'   mb.ratio = c(0.6, 0.4),
#'   size_point = 2, color_point = "#505050", color_matrix_shade = "#f89c9f",
#'   color_bar = '#505050',
#'   title_matrix_x = "Set Size", title_bar_y = "Intersection Size",
#'   size_title_matrix_x = 2, size_title_matrix_y = 1.7, size_title_bar_y = 2,
#'   size_matrix_x = 2, size_bar_y = 2, size_bar_label = 1.75)
#' }
#'
#'
#@importFrom UpSetR fromList upset
#' @importFrom grDevices dev.off pdf
#'
Upset_plot <- function(
    data,
    color_matrix = NULL,

    n = 40,
    custom_order = NULL,
    order_by = "freq",

    mb.ratio = c(0.6, 0.4),
    size_point = 2,

    color_point = "#505050",
    color_matrix_shade = "#f89c9f",
    color_bar = '#505050',

    title_matrix_x = "Set Size",
    title_bar_y = "Intersection Size",

    size_title_matrix_x = 2,
    size_title_matrix_y = 1.7,
    size_title_bar_y = 2,

    size_matrix_x = 2,
    size_bar_y = 2,
    size_bar_label = 1.75,

    queries = NULL,

    filename = "Upset",
    file_width = 16,
    file_height = 9
) {
  order_by <- tolower(order_by)

  if(!order_by %in% c("freq", "degree")) {
    cat("Please enter 'freq' or 'degree'!\n")
    order_by <- "freq"
    cat("The formal argument `order_by` has been reset to \"freq\".\n")
  }

  custom_order = rev(custom_order)


  ##
  df <- list()


  ##
  for (i in 1:length(colnames(data))){
    group <- colnames(data)[i]
    df[[group]] <- rownames(data)[which(data[,i] != 0)]
  }


  ##
  if(is.null(color_matrix)) {
    color_matrix = c("#3d97cb", "#57b9e5", "#8ae0ff", "#b0eaff",
                     "#ff4f24", "#ff7752", "#ffab7c", "#ffd9be",
                     "#038d6f", "#2ebe8c", "#48e092", "#70ffa9",
                     "#8f38ff", "#ae75ed", "#c886f2", "#e7a3c6",
                     "#ffca18", "#ffdb63", "#ffe899", "#fff3cc",
                     "#00fff6", "#4aebe5", "#8aefeb", "#d3f5f4",
                     "#ff2aba", "#fe70d0", "#ecaad7", "#ffe6f2")

    num_df <- length(df)
    num_color <- length(color_matrix)

    if(num_color >= num_df) {
      color_matrix2 <- color_matrix[1: num_df]

    } else {
      color_matrix <- c(color_matrix, color_matrix, color_matrix, color_matrix)
      color_matrix2 <- color_matrix[1: num_df]
    }

  } else {
    num_df <- length(df)
    num_color <- length(color_matrix)
    color_matrix2 <- color_matrix[1: num_df]
  }

  ###
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    stop("The 'UpSetR' package is required but not installed. Please install
         it using 'install.packages(\"UpSetR\")'.")
  }

  ###
  p1 <- UpSetR::upset(
    data = UpSetR::fromList(df),
    nsets = length(df),
    nintersects = n,

    #
    keep.order = TRUE,
    sets = custom_order,

    #
    number.angles = 0,
    point.size = size_point,
    line.size = 1,
    mb.ratio = mb.ratio,

    #
    text.scale = c(size_title_bar_y,
                   size_bar_y,
                   size_title_matrix_x,
                   size_matrix_x,
                   size_title_matrix_y,
                   size_bar_label),

    shade.color = color_matrix_shade,

    # y
    mainbar.y.label = title_bar_y,
    main.bar.color = color_bar,
    order.by = order_by,
    decreasing = c(T, F),

    # x
    sets.x.label = title_matrix_x,
    matrix.color = color_point,

    #
    sets.bar.color = color_matrix2,

    #
    queries = queries
  )
  print(p1)
  ###


  ##
  grDevices::pdf(file = paste0(filename, ".pdf"), width = file_width,
                 height = file_height)
  print(p1)
  grDevices::dev.off()


  ##
  cat("\033[32mtaxa_bar: success!\033[0m\n")
  cat("\033[0;32m", "The file has been saved to \n",
      getwd(), "\033[0m\n", sep = "")


  return(p1)
}
